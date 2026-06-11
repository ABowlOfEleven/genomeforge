//! The eframe application: layout, worker wiring, and all the Phase 1 panels
//! (genome browser, variant detail, variant table, search, health highlights).

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use eframe::egui;
use egui_extras::{Column, TableBuilder};

use gx_core::{Assembly, GenomicRange, Strand};

use crate::ancestry::AncestryState;
use crate::browser::{self, BrowserData};
use crate::crispr::CrisprState;
use crate::document::{Document, RefData};
use crate::health;
use crate::pgx::PgxState;
use crate::phenotype::{PhenoRequest, PhenotypeState};
use crate::plasmid::PlasmidState;
use crate::settings::Settings;
use crate::theme::{self, palette};
use crate::tutorial::{self, Section};
use crate::ux::{self, Tier};
use crate::view::GenomeView;
use crate::worker::{Request, Response, Worker};

const FEATURE_TILE: u64 = 1_000_000;

/// Which tool/workspace is active.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    Genome,
    Plasmid,
    Crispr,
    Phenotype,
    Ancestry,
    Pharma,
}

impl Tool {
    fn section(self) -> Section {
        match self {
            Tool::Genome => Section::Genome,
            Tool::Plasmid => Section::Plasmid,
            Tool::Crispr => Section::Crispr,
            Tool::Phenotype => Section::Phenotype,
            Tool::Ancestry => Section::Ancestry,
            Tool::Pharma => Section::Pharma,
        }
    }
}

pub struct GenomeForgeApp {
    view: GenomeView,
    settings: Settings,
    last_online: bool,
    show_settings: bool,
    status: String,

    worker: Worker,
    document: Option<Document>,
    refdata: RefData,
    plasmid: Option<PlasmidState>,
    crispr: CrisprState,
    phenotype: PhenotypeState,
    ancestry: AncestryState,
    pgx: PgxState,
    tool: Tool,
    show_tutorial: bool,
    tutorial_section: Section,
    tutorial_step: usize,
    /// Number of fetch requests awaiting a reply (drives the busy spinner).
    inflight: u32,

    selection: Option<(String, u64)>,
    search: String,
    show_table: bool,
    table_filter: String,
    table_only_notable: bool,
    requested_annot: HashSet<String>,
    /// PubMed articles fetched per rsID (None entry = fetch in flight).
    literature: HashMap<String, Option<Vec<gx_annotate::Article>>>,
    /// Last liftover result for the selected variant: (target build, mapped coord).
    lift: Option<(Assembly, Option<(String, u64)>)>,
    /// Newer release found by the update check, if any.
    update: Option<gx_annotate::UpdateInfo>,
    show_about: bool,
    show_shortcuts: bool,
    /// Local annotation-cache entry counts (refreshed when Settings opens).
    cache_stats: Option<gx_annotate::CacheStats>,
    /// On-disk path of the annotation cache (for size/age display).
    cache_file: PathBuf,
}

/// The GitHub repo the update check queries.
const REPO: &str = "ABowlOfEleven/genomeforge";

impl GenomeForgeApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        crate::theme::apply(&cc.egui_ctx);
        let settings: Settings = cc
            .storage
            .and_then(|s| eframe::get_value(s, "settings"))
            .unwrap_or_default();
        let cache_file = cache_path();
        let worker = Worker::spawn(
            cc.egui_ctx.clone(),
            cache_file.clone(),
            settings.online,
            settings.api_key_opt(),
        );
        // CLI: first non-flag arg = file to import; `--pgs <ID>` precomputes a score.
        let args: Vec<String> = std::env::args().skip(1).collect();
        let initial = args.iter().find(|a| !a.starts_with("--")).map(PathBuf::from);
        let mut startup_requests = 0u32; // pre-count so the in-flight tally is exact
        let status = match &initial {
            Some(p) => {
                worker.send(Request::Import(p.clone()));
                startup_requests += 1;
                format!("Importing {}…", p.display())
            }
            None => "Ready. Open a genome file to begin (File ▸ Open / Ctrl+O).".to_string(),
        };
        if let Some(i) = args.iter().position(|a| a == "--pgs")
            && let Some(id) = args.get(i + 1)
        {
            worker.send(Request::FetchPgs(id.clone()));
            startup_requests += 1;
        }
        // Quietly check for a newer release on startup (online only).
        if settings.online {
            worker.send(Request::CheckUpdate(REPO.to_string()));
            startup_requests += 1;
        }
        Self {
            view: GenomeView::default(),
            last_online: settings.online,
            settings,
            show_settings: false,
            status,
            worker,
            document: None,
            refdata: RefData::default(),
            plasmid: None,
            crispr: CrisprState::default(),
            phenotype: PhenotypeState::default(),
            ancestry: AncestryState::default(),
            pgx: PgxState::default(),
            tool: Tool::Genome,
            show_tutorial: false,
            tutorial_section: Section::Genome,
            tutorial_step: 0,
            inflight: startup_requests,
            selection: None,
            search: String::new(),
            show_table: true,
            table_filter: String::new(),
            table_only_notable: false,
            requested_annot: HashSet::new(),
            literature: HashMap::new(),
            lift: None,
            update: None,
            show_about: false,
            show_shortcuts: false,
            cache_stats: None,
            cache_file,
        }
    }

    /// Open the Settings window and refresh the cache statistics shown in it.
    fn open_settings(&mut self) {
        self.show_settings = true;
        self.send(Request::CacheStats);
    }

    fn current_assembly(&self) -> Assembly {
        self.document
            .as_ref()
            .and_then(|d| d.variant_doc())
            .map(|d| d.assembly)
            .unwrap_or_else(|| self.settings.assembly_override.unwrap_or_default())
    }

    // ---- worker plumbing ---------------------------------------------------

    /// Send a request, counting the ones that produce exactly one reply so the
    /// status bar can show a busy spinner. (Set*/Shutdown are fire-and-forget.)
    fn send(&mut self, req: Request) {
        if !matches!(
            req,
            Request::SetOnline(_) | Request::SetAssembly(_) | Request::Shutdown
        ) {
            self.inflight += 1;
        }
        self.worker.send(req);
    }

    fn drain_worker(&mut self) {
        while let Ok(resp) = self.worker.rx.try_recv() {
            self.inflight = self.inflight.saturating_sub(1);
            match resp {
                Response::Imported(Ok(imported)) => self.on_imported(imported),
                Response::Imported(Err(e)) => self.status = format!("⚠ Import failed: {e}"),
                Response::Annotations(map) => {
                    if let Some(doc) = self.document.as_mut().and_then(|d| d.variant_doc_mut()) {
                        for (k, v) in map {
                            doc.annotations.insert(k, v);
                        }
                    }
                }
                Response::Features { features, .. } => self.refdata.add_features(features),
                Response::Sequence { range, seq } => {
                    self.refdata.current_seq = Some((range.contig, range.start, seq));
                }
                Response::Gene { symbol, location } => match location {
                    Some(g) => {
                        self.view.go_to(g.contig.clone(), g.start as f64, g.end as f64, 1200.0);
                        self.status = format!("Jumped to {symbol} · {}", g.range().display_region());
                    }
                    None => self.status = format!("Gene '{symbol}' not found."),
                },
                Response::Pgs(Ok(sf)) => {
                    if let Some(doc) = self.document.as_ref().and_then(|d| d.variant_doc()) {
                        let result = gx_pgs::apply(&sf, &doc.store);
                        log::info!(
                            "PRS {} applied: {}/{} matched, percentile={:?}",
                            sf.trait_name,
                            result.matched,
                            result.total,
                            result.percentile
                        );
                        self.phenotype.status = format!(
                            "{}: {}/{} variants matched",
                            sf.trait_name, result.matched, result.total
                        );
                        self.phenotype.results.push(result);
                        self.phenotype.selected = Some(self.phenotype.results.len() - 1);
                    } else {
                        self.phenotype.status = "Load a genome before computing scores.".into();
                    }
                }
                Response::Pgs(Err(e)) => self.phenotype.status = format!("PGS error: {e}"),
                Response::Literature { rsid, articles } => {
                    self.literature.insert(rsid, Some(articles));
                }
                Response::Liftover { to, mapped } => {
                    self.status = match &mapped {
                        Some((c, p)) => format!("Lifted to {}: {c}:{}", to.label(), p + 1),
                        None => format!("Position does not map to {}.", to.label()),
                    };
                    self.lift = Some((to, mapped));
                }
                Response::Update(info) => {
                    if let Some(u) = &info
                        && gx_annotate::is_newer(&u.tag, env!("CARGO_PKG_VERSION"))
                    {
                        self.status = format!("Update available: {} (Help ▸ About).", u.tag);
                        self.update = info;
                    }
                }
                Response::CacheStats(s) => self.cache_stats = Some(s),
                Response::Notice(s) => self.status = s,
                Response::Failed(e) => {
                    self.status = format!("⚠ {e}");
                    // A transient failure shouldn't permanently suppress retries:
                    // re-enable annotation fetches for the variants still in view.
                    self.requested_annot.clear();
                }
                Response::Idle => {}
            }
        }
    }

    fn on_imported(&mut self, imported: gx_io::Imported) {
        self.refdata.clear();
        self.requested_annot.clear();
        self.selection = None;
        // PRS results belong to the previous sample — never carry them over.
        self.phenotype = PhenotypeState::default();
        // Haplogroups are sample-specific too; recompute lazily for the new one.
        self.ancestry = AncestryState::default();
        // Same for the pharmacogenomics report.
        self.pgx = PgxState::default();
        let doc = Document::from_imported(imported, self.settings.assembly_override);
        match &doc {
            Document::Variants(v) => {
                self.status = format!(
                    "Loaded {} variants · {} · {}",
                    v.store.len(),
                    v.format_label,
                    v.assembly.label()
                );
                if let Some(c) = v.store.contigs().next() {
                    let len = gx_core::chrom_length(v.assembly, c).unwrap_or(1);
                    self.view.contig = c.to_string();
                    self.view.center = (len / 2) as f64;
                    self.view.bp_per_px = (len as f64 / 1200.0).max(1.0);
                }
                self.send(Request::SetAssembly(v.assembly));
            }
            Document::Sequences(s) => {
                self.status = format!(
                    "Loaded {} sequence record(s). Plasmid designer ready.",
                    s.len()
                );
            }
        }
        self.document = Some(doc);
        self.plasmid = match &self.document {
            Some(Document::Sequences(recs)) if !recs.is_empty() => {
                Some(PlasmidState::new(0, &recs[0]))
            }
            _ => None,
        };
        // Switch to the matching workspace, but DON'T auto-open the tutorial over
        // the freshly loaded data — first-run tours fire only on an explicit tab click.
        self.tool = if matches!(self.document, Some(Document::Sequences(_))) {
            Tool::Plasmid
        } else {
            Tool::Genome
        };
    }

    /// Fetch the gene models / sequence / annotations needed for the current
    /// view, de-duplicated so we never re-request the same thing.
    fn request_visible(&mut self, width: f32) {
        if !self.settings.online {
            return;
        }
        let (lo, hi) = self.view.visible_range(width);
        let span = hi - lo;
        let contig = self.view.contig.clone();

        // Gene-model tiles, only when zoomed in enough to be meaningful.
        if self.document.is_some() && span <= 5_000_000.0 {
            let mut tile = (lo.max(0.0) as u64) / FEATURE_TILE * FEATURE_TILE;
            let end = hi.max(0.0) as u64;
            let mut sent = 0;
            while tile <= end && sent < 6 {
                let range = GenomicRange::new(contig.clone(), tile, tile + FEATURE_TILE, Strand::Unknown);
                if self.refdata.fetched_tiles.insert(range.ensembl_region()) {
                    self.send(Request::Features(range));
                    sent += 1;
                }
                tile += FEATURE_TILE;
            }
        }

        // Reference sequence, only at base resolution.
        if self.view.shows_bases() {
            let start = lo.max(0.0) as u64;
            let end = hi.max(0.0) as u64 + 1;
            let range = GenomicRange::new(contig.clone(), start, end, Strand::Unknown);
            let covered = matches!(&self.refdata.current_seq,
                Some((c, st, sq)) if *c == contig && *st <= start && st + sq.len() as u64 >= end);
            if !covered && self.refdata.requested_seq.insert(range.ensembl_region()) {
                self.send(Request::Sequence(range));
            }
        }

        // Annotate the variants currently in view (bounded, rate-limit friendly).
        if let Some(doc) = self.document.as_ref().and_then(|d| d.variant_doc()) {
            let visible = doc
                .store
                .variants_in(&contig, lo.max(0.0) as u64, hi.max(0.0) as u64);
            if visible.len() <= 500 {
                let mut batch = Vec::new();
                for v in visible {
                    if let Some(id) = &v.rsid
                        && !doc.annotations.contains_key(id) && self.requested_annot.insert(id.clone()) {
                            batch.push(id.clone());
                            if batch.len() >= 300 {
                                break;
                            }
                        }
                }
                if !batch.is_empty() {
                    self.send(Request::Annotate(batch));
                }
            }
        }
    }

    fn select_variant(&mut self, contig: String, pos: u64, rsid: Option<String>) {
        self.view.contig = contig.clone();
        self.view.center = pos as f64;
        self.selection = Some((contig, pos));
        self.lift = None; // a liftover result belongs to the previously selected variant
        if let Some(id) = rsid {
            let needed = self
                .document
                .as_ref()
                .and_then(|d| d.variant_doc())
                .map(|d| !d.annotations.contains_key(&id))
                .unwrap_or(false);
            if needed && self.settings.online && self.requested_annot.insert(id.clone()) {
                self.send(Request::Annotate(vec![id]));
            }
        }
    }

    fn select_by_rsid(&mut self, rsid: String) {
        if let Some(v) = self
            .document
            .as_ref()
            .and_then(|d| d.variant_doc())
            .and_then(|d| d.store.find_by_rsid(&rsid))
        {
            let (c, p) = (v.contig.clone(), v.pos);
            self.view.bp_per_px = self.view.bp_per_px.min(2.0);
            self.select_variant(c, p, Some(rsid));
        } else {
            self.status = format!("rsID {rsid} not present in this dataset.");
        }
    }

    fn run_search(&mut self) {
        let q = self.search.trim().to_string();
        if q.is_empty() {
            return;
        }
        if q.starts_with("rs") || q.starts_with('i') {
            self.select_by_rsid(q);
        } else if let Some((contig, pos)) = parse_locus(&q) {
            self.view.contig = contig;
            self.view.center = pos as f64;
            self.view.bp_per_px = self.view.bp_per_px.min(500.0);
        } else {
            self.send(Request::LookupGene(q.clone()));
            self.status = format!("Looking up gene {q}…");
        }
    }

    fn open_file_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter(
                "Genome / sequence",
                &["txt", "vcf", "gz", "fa", "fasta", "fna", "gb", "gbk", "gbff"],
            )
            .add_filter("All files", &["*"])
            .pick_file()
        {
            self.import_file(path);
        }
    }

    /// Import a file (from the dialog, the recent list, or a drag-and-drop),
    /// recording it in the recent-files list.
    fn import_file(&mut self, path: PathBuf) {
        self.settings.push_recent(&path);
        self.status = format!("Importing {}…", path.display());
        self.send(Request::Import(path));
    }

    // ---- panels ------------------------------------------------------------

    fn menu_bar(&mut self, ui: &mut egui::Ui) {
        use theme::tokens;
        ui.horizontal(|ui| {
            // Wordmark — drawn double-helix mark + two-tone lockup. The mark is
            // painted (not an emoji) so it always renders and matches the app icon.
            ui.add_space(3.0);
            draw_brand_mark(ui);
            ui.add_space(2.0);
            let mut job = egui::text::LayoutJob::default();
            job.append(
                "Genome",
                0.0,
                egui::TextFormat::simple(egui::FontId::proportional(15.0), tokens::TEXT),
            );
            job.append(
                "Forge",
                0.0,
                egui::TextFormat::simple(egui::FontId::proportional(15.0), tokens::ACCENT),
            );
            ui.label(job).on_hover_text("GenomeForge: local genome studio");
            ui.separator();

            ui.menu_button("File", |ui| {
                if ui.button("Open genome…  (Ctrl+O)").clicked() {
                    self.open_file_dialog();
                    ui.close();
                }
                let mut reopen: Option<PathBuf> = None;
                ui.menu_button("Open Recent", |ui| {
                    if self.settings.recent_files.is_empty() {
                        ui.label(egui::RichText::new("(nothing yet)").weak());
                    }
                    for f in &self.settings.recent_files {
                        if ui.button(f).clicked() {
                            reopen = Some(PathBuf::from(f));
                            ui.close();
                        }
                    }
                    if !self.settings.recent_files.is_empty() {
                        ui.separator();
                        if ui.button("Clear recent").clicked() {
                            self.settings.recent_files.clear();
                            ui.close();
                        }
                    }
                });
                if let Some(p) = reopen {
                    self.import_file(p);
                }
            });
            ui.menu_button("View", |ui| {
                ui.checkbox(&mut self.show_table, "Variant table");
                ui.separator();
                ui.label("Jump to chromosome:");
                egui::ScrollArea::vertical().max_height(240.0).show(ui, |ui| {
                    let assembly = self.current_assembly();
                    for c in gx_core::chrom_order() {
                        if ui.button(format!("chr{c}")).clicked() {
                            let len = gx_core::chrom_length(assembly, c).unwrap_or(1);
                            self.view.contig = c.to_string();
                            self.view.center = (len / 2) as f64;
                            self.view.bp_per_px = (len as f64 / 1200.0).max(1.0);
                        }
                    }
                });
            });
            if ui.button("⚙ Settings").clicked() {
                self.open_settings();
            }

            ui.separator();
            // Tool tabs.
            let has_variants = matches!(self.document, Some(Document::Variants(_)));
            let has_seq = matches!(self.document, Some(Document::Sequences(_)));
            let genome_ok = has_variants || self.document.is_none();
            if ui
                .add_enabled(
                    genome_ok,
                    egui::Button::selectable(self.tool == Tool::Genome, "Genome"),
                )
                .on_hover_text("Browse your genome and explore variants")
                .on_disabled_hover_text("A sequence file is loaded; open variant data to use this")
                .clicked()
            {
                self.set_tool(Tool::Genome);
            }
            if ui
                .add_enabled(
                    has_seq,
                    egui::Button::selectable(self.tool == Tool::Plasmid, "Plasmid"),
                )
                .on_hover_text("Design and visualise plasmids / sequences")
                .on_disabled_hover_text("Open a FASTA or GenBank file to enable")
                .clicked()
            {
                self.set_tool(Tool::Plasmid);
            }
            if ui
                .add(egui::Button::selectable(self.tool == Tool::Crispr, "CRISPR"))
                .on_hover_text("Design CRISPR guides and simulate edits")
                .clicked()
            {
                self.set_tool(Tool::Crispr);
            }
            if ui
                .add_enabled(
                    has_variants,
                    egui::Button::selectable(self.tool == Tool::Phenotype, "Phenotype"),
                )
                .on_hover_text("Polygenic risk scores and Mendelian findings")
                .on_disabled_hover_text("Load a genome (variant data) to enable")
                .clicked()
            {
                self.set_tool(Tool::Phenotype);
            }
            if ui
                .add_enabled(
                    has_variants,
                    egui::Button::selectable(self.tool == Tool::Ancestry, "Ancestry"),
                )
                .on_hover_text("Predicted maternal & paternal haplogroups")
                .on_disabled_hover_text("Load a genome (variant data) to enable")
                .clicked()
            {
                self.set_tool(Tool::Ancestry);
            }
            if ui
                .add_enabled(
                    has_variants,
                    egui::Button::selectable(self.tool == Tool::Pharma, "Pharma"),
                )
                .on_hover_text("Pharmacogenomics: how your variants relate to drug response")
                .on_disabled_hover_text("Load a genome (variant data) to enable")
                .clicked()
            {
                self.set_tool(Tool::Pharma);
            }

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ux::tier_picker(ui, &mut self.settings.tier);
                ui.menu_button("Help", |ui| {
                    if ui.button("Tutorial for this workspace").clicked() {
                        let section = self.tool.section();
                        self.open_tutorial(section);
                        ui.close();
                    }
                    if ui.button("Keyboard shortcuts").clicked() {
                        self.show_shortcuts = true;
                        ui.close();
                    }
                    if ui.button("Check for updates").clicked() {
                        self.send(Request::CheckUpdate(REPO.to_string()));
                        self.status = "Checking for updates…".into();
                        ui.close();
                    }
                    ui.separator();
                    if ui.button("About GenomeForge").clicked() {
                        self.show_about = true;
                        ui.close();
                    }
                });
                ui.separator();
                let (mode, mode_color) = if self.settings.online {
                    ("online", palette::BENIGN)
                } else {
                    ("offline", palette::VUS)
                };
                // Drawn status dot (a painted circle, not a glyph, so it always renders).
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 5.0;
                    let (dot, _) = ui.allocate_exact_size(egui::vec2(8.0, 8.0), egui::Sense::hover());
                    ui.painter().circle_filled(dot.center(), 4.0, mode_color);
                    ui.label(egui::RichText::new(mode).color(mode_color));
                })
                .response
                .on_hover_text("Toggle network access in Settings");
                if has_variants {
                    ui.separator();
                    ui.label(
                        egui::RichText::new(self.current_assembly().label())
                            .color(palette::RULER_TEXT),
                    );
                }
            });
        });
    }

    fn sidebar(&mut self, ui: &mut egui::Ui) {
        let online = self.settings.online;
        let tier = self.settings.tier;
        ui.add_space(4.0);
        ui.heading("Genome browser");
        ui.label(
            egui::RichText::new("Variant explorer")
                .size(11.0)
                .color(palette::RULER_TEXT),
        );
        ui.separator();

        let resp = ui.add(
            egui::TextEdit::singleline(&mut self.search)
                .hint_text("Search: gene / rsID / chr:pos")
                .desired_width(f32::INFINITY),
        );
        if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
            self.run_search();
        }
        ui.separator();

        if !online {
            ui.label(
                egui::RichText::new(
                    "Offline: gene models and new annotations won't load. Enable Online in Settings.",
                )
                .size(11.0)
                .color(palette::VUS),
            );
            ui.separator();
        }

        match self.document.as_ref() {
            None => {
                ui.label(
                    egui::RichText::new(
                        "No genome loaded.\nFile ▸ Open to import 23andMe/AncestryDNA, VCF, FASTA, or GenBank.",
                    )
                    .size(11.0)
                    .color(palette::RULER_TEXT),
                );
            }
            Some(Document::Sequences(s)) => {
                ui.label(format!("{} sequence record(s) loaded.", s.len()));
                for rec in s.iter().take(20) {
                    ui.label(format!(
                        "• {} ({} bp{})",
                        rec.name,
                        rec.len(),
                        if rec.circular { ", circular" } else { "" }
                    ));
                }
            }
            Some(Document::Variants(doc)) => {
                ui.label(format!("{} variants · {}", doc.store.len(), doc.format_label));
                ui.label(
                    egui::RichText::new(format!(
                        "{} annotated · {} notable",
                        doc.annotations.len(),
                        doc.notable_count()
                    ))
                    .size(11.0)
                    .color(palette::RULER_TEXT),
                );
                for w in &doc.warnings {
                    ui.label(egui::RichText::new(format!("⚠ {w}")).size(11.0).color(palette::VUS));
                }
                ui.add_space(6.0);
                ux::explain_beginner(
                    ui,
                    tier,
                    "The colours below show how serious each variant is, based on the ClinVar \
                     database. Click a lollipop in the browser to learn more.",
                );
                variant_legend(ui);
                ui.add_space(6.0);
                ui.label(egui::RichText::new("Health highlights").strong());

                let notable: Vec<(String, String, String)> = doc
                    .annotations
                    .values()
                    .filter(|a| a.is_clinically_notable())
                    .take(40)
                    .map(|a| {
                        (
                            a.rsid.clone(),
                            a.gene.clone().unwrap_or_default(),
                            a.clinical_significance.clone().unwrap_or_default(),
                        )
                    })
                    .collect();

                if notable.is_empty() {
                    ui.label(
                        egui::RichText::new("Browse / search to annotate variants. Notable ClinVar hits appear here.")
                            .size(11.0)
                            .color(palette::RULER_TEXT),
                    );
                } else {
                    let mut jump: Option<String> = None;
                    for (rsid, gene, sig) in &notable {
                        if ui
                            .button(format!("{rsid}  {gene}\n{sig}"))
                            .on_hover_text("Jump to this variant")
                            .clicked()
                        {
                            jump = Some(rsid.clone());
                        }
                    }
                    if let Some(rsid) = jump {
                        self.select_by_rsid(rsid);
                    }
                }

                // Curated quick-jumps to famous health/trait SNPs (the genome
                // analogue of the Phenotype "Example scores"). Jumps to the SNP
                // in the loaded data and pulls its annotation; if a chip didn't
                // genotype it, the status bar says so.
                ui.add_space(6.0);
                let picked = egui::CollapsingHeader::new("Common health SNPs")
                    .default_open(false)
                    .show(ui, |ui| {
                        ux::explain_beginner(
                            ui,
                            tier,
                            "Well-known SNPs people often look up. Click one to jump to it in your \
                             genome. If your data doesn't include it, you'll see a note. These are \
                             starting points for learning, not a medical screen.",
                        );
                        let mut pick = None;
                        for s in health::HEALTH_SNPS {
                            let label = if tier.is_expert() {
                                format!("{}  {}\n{}", s.rsid, s.gene, s.trait_)
                            } else {
                                format!("{}\n{}", s.gene, s.trait_)
                            };
                            if ui
                                .button(label)
                                .on_hover_text(format!("{} · {}: jump to this SNP", s.rsid, s.gene))
                                .clicked()
                            {
                                pick = Some(s.rsid.to_string());
                            }
                        }
                        pick
                    })
                    .body_returned
                    .flatten();
                if let Some(rsid) = picked {
                    self.select_by_rsid(rsid);
                }
            }
        }
    }

    fn detail_panel(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);
        ui.heading("Variant details");
        ui.separator();
        let tier = self.settings.tier;
        let online = self.settings.online;
        ux::explain_beginner(
            ui,
            tier,
            "Each variant is one spot where your DNA was read. Click a lollipop in the browser \
             to see what is known about it.",
        );

        let Some(doc) = self.document.as_ref().and_then(|d| d.variant_doc()) else {
            ui.label(egui::RichText::new("Load a genome, then click a variant.").color(palette::RULER_TEXT));
            return;
        };
        let Some((contig, pos)) = self.selection.clone() else {
            ui.label(egui::RichText::new("Click a variant in the browser to inspect it.").color(palette::RULER_TEXT));
            return;
        };
        let Some(v) = doc.store.variants_in(&contig, pos, pos + 1).first() else {
            ui.label("Variant no longer in view.");
            return;
        };

        egui::Grid::new("variant_detail").num_columns(2).striped(true).show(ui, |ui| {
            row(ui, "ID", &v.display_id());
            row(ui, "Location", &format!("{}:{}", contig, pos + 1));
            if let Some(gt) = &v.genotype {
                row(ui, "Your genotype", &format!("{}  ({:?})", gt.display(), gt.zygosity()));
            }
            if let Some(r) = &v.ref_allele {
                row(ui, "Ref / Alt", &format!("{} / {}", r, v.alt_alleles.join(",")));
            }

            match doc.annotation_for(v) {
                None => {
                    ui.label("Annotation");
                    let msg = if v.rsid.is_none() {
                        "no rsID, can't annotate"
                    } else if online {
                        "fetching…"
                    } else {
                        "offline; enable Online in Settings"
                    };
                    ui.label(egui::RichText::new(msg).color(palette::RULER_TEXT));
                    ui.end_row();
                }
                Some(ann) => {
                    if let Some(g) = &ann.gene {
                        row(ui, "Gene", g);
                    }
                    if let Some(c) = &ann.consequence {
                        let k = if tier.is_beginner() { "Effect on protein" } else { "Consequence" };
                        row(ui, k, c);
                    }
                    if let Some(sig) = &ann.clinical_significance {
                        ui.label(if tier.is_beginner() { "Clinical significance" } else { "ClinVar" });
                        ui.label(egui::RichText::new(sig).color(theme::sig_color(sig)).strong());
                        ui.end_row();
                    }
                    if !ann.conditions.is_empty() {
                        row(ui, "Condition(s)", &ann.conditions.join("; "));
                    }
                    if let Some(af) = ann.gnomad_af {
                        let k = if tier.is_beginner() { "How common (population)" } else { "gnomAD freq" };
                        row(ui, k, &format_af(af));
                    }
                    if tier.at_least(Tier::Intermediate)
                        && let Some(cadd) = ann.cadd_phred {
                            row(ui, "CADD (phred)", &format!("{cadd:.1}"));
                        }
                    if !ann.gwas_traits.is_empty() {
                        row(ui, "GWAS traits", &ann.gwas_traits.join("; "));
                    }
                }
            }
        });
        if doc.annotation_for(v).is_some() {
            ui.label(
                egui::RichText::new("via MyVariant.info, cached locally")
                    .size(10.0)
                    .color(palette::RULER_TEXT),
            );
        }

        // Pull owned identifiers, then the `doc`/`v` borrow ends so the
        // &mut self drill-downs (literature, liftover) are legal.
        let rsid = v.rsid.clone();
        let gene = doc.annotation_for(v).and_then(|a| a.gene.clone());
        let gnomad_id = v
            .ref_allele
            .clone()
            .zip(v.alt_alleles.first().cloned())
            .map(|(r, a)| format!("{}-{}-{}-{}", contig.trim_start_matches("chr"), pos + 1, r, a));

        resource_links(ui, rsid.as_deref(), gene.as_deref(), gnomad_id.as_deref(), tier);
        if let Some(rsid) = &rsid {
            self.literature_section(ui, rsid, tier, online);
        }
        self.liftover_section(ui, &contig, pos, tier);
        ui.add_space(8.0);
        ux::disclaimer(ui);
    }

    /// "Related publications" drill-down: an on-demand PubMed search for the
    /// selected variant, rendered as clickable references.
    fn literature_section(&mut self, ui: &mut egui::Ui, rsid: &str, tier: Tier, online: bool) {
        ui.add_space(10.0);
        ui.label(egui::RichText::new("Related publications").strong());
        ux::explain_beginner(
            ui,
            tier,
            "Research papers from NIH's PubMed that mention this variant. A starting point for \
             reading the primary literature. Quality and relevance vary.",
        );
        // Clone the small entry so the match doesn't hold a borrow of
        // `self.literature` while the button arm mutates it / calls `self.send`.
        match self.literature.get(rsid).cloned() {
            Some(None) => {
                ui.horizontal(|ui| {
                    ui.add(egui::Spinner::new().size(13.0));
                    ui.label(
                        egui::RichText::new("Searching PubMed…")
                            .size(11.0)
                            .color(palette::RULER_TEXT),
                    );
                });
            }
            Some(Some(arts)) if arts.is_empty() => {
                ui.label(
                    egui::RichText::new("No PubMed articles mention this rsID.")
                        .size(11.0)
                        .color(palette::RULER_TEXT),
                );
            }
            Some(Some(arts)) => {
                for a in &arts {
                    let title = ui
                        .add(
                            egui::Label::new(egui::RichText::new(&a.title).color(palette::ACCENT))
                                .sense(egui::Sense::click()),
                        )
                        .on_hover_text(format!("Open PubMed {} in your browser", a.pmid));
                    if title.clicked() {
                        ui.ctx().open_url(egui::OpenUrl::new_tab(a.url()));
                    }
                    let meta = [a.journal.as_str(), a.year.as_str()]
                        .iter()
                        .filter(|s| !s.is_empty())
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(" · ");
                    if !meta.is_empty() {
                        ui.label(egui::RichText::new(meta).size(10.0).color(palette::RULER_TEXT));
                    }
                    ui.add_space(4.0);
                }
            }
            None => {
                if ui
                    .add_enabled(online, egui::Button::new("Find related papers").small())
                    .on_disabled_hover_text("Enable Online in Settings to search PubMed")
                    .clicked()
                {
                    self.literature.insert(rsid.to_string(), None);
                    self.send(Request::Literature(rsid.to_string()));
                }
            }
        }
    }

    /// Liftover control: convert this position to the other human build.
    fn liftover_section(&mut self, ui: &mut egui::Ui, contig: &str, pos: u64, tier: Tier) {
        let from = self.current_assembly();
        let to = other_assembly(from);
        ui.add_space(10.0);
        ui.label(egui::RichText::new("Liftover").strong());
        ux::explain_beginner(
            ui,
            tier,
            "Genome builds number the same base differently. This converts the position to \
             the other build (GRCh37 ↔ GRCh38) using Ensembl.",
        );
        if ui
            .add_enabled(
                self.settings.online,
                egui::Button::new(format!("Lift {contig}:{} to {}", pos + 1, to.label())).small(),
            )
            .on_disabled_hover_text("Enable Online in Settings to use liftover")
            .clicked()
        {
            self.lift = None;
            self.send(Request::Liftover { from, to, contig: contig.to_string(), pos });
            self.status = "Lifting coordinate…".into();
        }
        if let Some((t, mapped)) = &self.lift {
            match mapped {
                Some((c, p)) => {
                    ui.label(
                        egui::RichText::new(format!("→ {}  {c}:{}", t.label(), p + 1))
                            .color(palette::BENIGN),
                    );
                }
                None => {
                    ui.label(
                        egui::RichText::new(format!("Does not map to {}.", t.label()))
                            .size(11.0)
                            .color(palette::VUS),
                    );
                }
            }
        }
    }

    /// Returns an rsID if a table row was clicked.
    fn table_panel(&mut self, ui: &mut egui::Ui) -> Option<String> {
        let GenomeForgeApp {
            document,
            table_filter,
            table_only_notable,
            ..
        } = self;
        let Some(doc) = document.as_ref().and_then(|d| d.variant_doc()) else {
            ui.label("No variants loaded.");
            return None;
        };

        ui.horizontal(|ui| {
            ui.label("Filter:");
            ui.add(egui::TextEdit::singleline(table_filter).hint_text("gene / rsID / condition").desired_width(180.0));
            ui.checkbox(table_only_notable, "Notable only");
            ui.label(
                egui::RichText::new(format!("{} annotated", doc.annotations.len()))
                    .size(11.0)
                    .color(palette::RULER_TEXT),
            );
        });

        let needle = table_filter.to_ascii_lowercase();
        let mut rows: Vec<&gx_annotate::VariantAnnotation> = doc
            .annotations
            .values()
            .filter(|a| !*table_only_notable || a.is_clinically_notable())
            .filter(|a| {
                if needle.is_empty() {
                    return true;
                }
                a.rsid.to_ascii_lowercase().contains(&needle)
                    || a.gene.as_deref().unwrap_or("").to_ascii_lowercase().contains(&needle)
                    || a.conditions.join(";").to_ascii_lowercase().contains(&needle)
            })
            .collect();
        rows.sort_by(|a, b| {
            b.is_clinically_notable()
                .cmp(&a.is_clinically_notable())
                .then(a.gene.cmp(&b.gene))
        });

        let mut clicked = None;
        TableBuilder::new(ui)
            .striped(true)
            .column(Column::auto().at_least(90.0))
            .column(Column::auto().at_least(70.0))
            .column(Column::auto().at_least(120.0))
            .column(Column::remainder())
            .header(18.0, |mut h| {
                h.col(|ui| { ui.strong("rsID"); });
                h.col(|ui| { ui.strong("Gene"); });
                h.col(|ui| { ui.strong("ClinVar"); });
                h.col(|ui| { ui.strong("Condition / trait"); });
            })
            .body(|mut body| {
                for a in &rows {
                    body.row(18.0, |mut r| {
                        r.col(|ui| {
                            if ui.link(&a.rsid).clicked() {
                                clicked = Some(a.rsid.clone());
                            }
                        });
                        r.col(|ui| { ui.label(a.gene.as_deref().unwrap_or("-")); });
                        r.col(|ui| {
                            match &a.clinical_significance {
                                Some(s) => { ui.label(egui::RichText::new(s).color(theme::sig_color(s))); }
                                None => { ui.label("-"); }
                            }
                        });
                        r.col(|ui| {
                            let txt = if !a.conditions.is_empty() {
                                a.conditions.join("; ")
                            } else {
                                a.gwas_traits.join("; ")
                            };
                            // Truncate long condition/trait lists to the column,
                            // full text on hover, so the panel can't blow out.
                            let resp = ui.add(egui::Label::new(&txt).truncate());
                            if !txt.is_empty() {
                                resp.on_hover_text(txt);
                            }
                        });
                    });
                }
            });

        ui.add_space(4.0);
        ui.horizontal(|ui| {
            if ui
                .button("Export table (CSV)")
                .on_hover_text("Save the currently filtered variants to a CSV file")
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .set_file_name("genomeforge-variants.csv")
                    .add_filter("CSV", &["csv"])
                    .save_file()
                && let Err(e) = std::fs::write(&path, variants_to_csv(&rows))
            {
                log::warn!("CSV export failed: {e}");
            }
            ui.label(
                egui::RichText::new(format!("{} rows", rows.len()))
                    .size(11.0)
                    .color(palette::RULER_TEXT),
            );
        });
        clicked
    }

    fn about_window(&mut self, ctx: &egui::Context) {
        let mut open = self.show_about;
        egui::Window::new("About GenomeForge")
            .open(&mut open)
            .resizable(false)
            .collapsible(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    draw_brand_mark(ui);
                    ui.add_space(2.0);
                    ui.heading("GenomeForge");
                });
                ui.label(
                    egui::RichText::new(concat!("Version ", env!("CARGO_PKG_VERSION")))
                        .color(palette::RULER_TEXT),
                );
                ui.add_space(4.0);
                ui.label("A local, native genome browser, variant explorer, plasmid designer,");
                ui.label("CRISPR studio, and polygenic-risk tool. Your data stays on your machine.");
                ui.add_space(6.0);
                if let Some(u) = self.update.clone() {
                    egui::Frame::new()
                        .fill(theme::tokens::ACCENT.gamma_multiply(0.12))
                        .stroke(egui::Stroke::new(1.0, theme::tokens::ACCENT))
                        .corner_radius(egui::CornerRadius::same(6))
                        .inner_margin(egui::Margin::same(8))
                        .show(ui, |ui| {
                            ui.label(
                                egui::RichText::new(format!("Update available: {}", u.tag)).strong(),
                            );
                            ui.hyperlink_to("Download the latest release", u.url);
                        });
                    ui.add_space(6.0);
                }
                ui.horizontal_wrapped(|ui| {
                    ui.spacing_mut().item_spacing.x = 10.0;
                    ui.hyperlink_to("Project", format!("https://github.com/{REPO}"));
                    ui.hyperlink_to("Releases", format!("https://github.com/{REPO}/releases"));
                    ui.hyperlink_to("Report an issue", format!("https://github.com/{REPO}/issues"));
                });
                ui.add_space(4.0);
                ui.label(egui::RichText::new("MIT licensed.").size(11.0).color(palette::RULER_TEXT));
                ux::disclaimer(ui);
            });
        self.show_about = open;
    }

    fn settings_window(&mut self, ctx: &egui::Context) {
        let mut open = self.show_settings;
        egui::Window::new("Settings")
            .open(&mut open)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.checkbox(&mut self.settings.online, "Online: query live databases");
                ui.horizontal(|ui| {
                    ui.label("MyVariant API key:");
                    ui.text_edit_singleline(&mut self.settings.api_key);
                });
                ui.add_space(6.0);
                ui.label("Assembly build (applied on next import):");
                ui.horizontal(|ui| {
                    let mut choice = self.settings.assembly_override;
                    if ui.selectable_label(choice.is_none(), "Auto-detect").clicked() {
                        choice = None;
                    }
                    if ui.selectable_label(choice == Some(Assembly::Grch38), "GRCh38").clicked() {
                        choice = Some(Assembly::Grch38);
                    }
                    if ui.selectable_label(choice == Some(Assembly::Grch37), "GRCh37").clicked() {
                        choice = Some(Assembly::Grch37);
                    }
                    self.settings.assembly_override = choice;
                });

                ui.separator();
                ui.label(egui::RichText::new("Local cache").strong());
                let meta = std::fs::metadata(&self.cache_file).ok();
                let size = meta.as_ref().map(|m| m.len()).unwrap_or(0);
                let age = meta
                    .as_ref()
                    .and_then(|m| m.modified().ok())
                    .and_then(|t| t.elapsed().ok());
                let stats = self.cache_stats;
                egui::Grid::new("cache_stats").num_columns(2).spacing([12.0, 2.0]).show(ui, |ui| {
                    ui.label("On disk:");
                    ui.label(human_bytes(size));
                    ui.end_row();
                    if let Some(age) = age {
                        ui.label("Last updated:");
                        ui.label(humanize_ago(age));
                        ui.end_row();
                    }
                    if let Some(s) = stats {
                        ui.label("Annotations:");
                        ui.label(s.annotations.to_string());
                        ui.end_row();
                        ui.label("Gene / sequence regions:");
                        ui.label((s.feature_regions + s.sequence_regions).to_string());
                        ui.end_row();
                    }
                });
                if ui
                    .button("Clear cache")
                    .on_hover_text(
                        "Delete all cached annotations, gene models, and sequence. They \
                         re-download on demand when Online.",
                    )
                    .clicked()
                {
                    self.send(Request::ClearCache);
                    self.cache_stats = None;
                    self.requested_annot.clear();
                    self.refdata.fetched_tiles.clear();
                    self.refdata.requested_seq.clear();
                    self.status = "Cache cleared.".into();
                }

                ui.separator();
                ui.label(egui::RichText::new("Data sources").strong());
                ui.label(
                    egui::RichText::new(
                        "Queried live when Online, then cached locally. Versions track each \
                         upstream service. Your raw genotypes are never sent; only rsIDs and \
                         genomic regions are.",
                    )
                    .size(11.0)
                    .color(palette::RULER_TEXT),
                );
                for line in [
                    "Variant annotations: MyVariant.info (ClinVar, dbSNP, gnomAD, CADD, GWAS Catalog)",
                    "Gene models & sequence: Ensembl REST",
                    "Polygenic scores: PGS Catalog",
                    "Publications: NCBI PubMed",
                    "Liftover: Ensembl assembly map",
                ] {
                    ui.label(egui::RichText::new(format!("• {line}")).size(11.0));
                }

                ui.separator();
                ux::disclaimer(ui);
            });
        self.show_settings = open;
    }

    fn draw_central(&mut self, ui: &mut egui::Ui) {
        let width = ui.available_size().x.max(1.0);
        self.request_visible(width);

        let assembly = self.current_assembly();
        let (lo, hi) = self.view.visible_range(width);
        let doc = self.document.as_ref().and_then(|d| d.variant_doc());
        let features =
            self.refdata
                .features_in(&self.view.contig, lo.max(0.0) as u64, hi.max(0.0) as u64);
        let sequence = match &self.refdata.current_seq {
            Some((c, start, seq))
                if *c == self.view.contig
                    && (*start as f64) <= lo
                    && (start + seq.len() as u64) as f64 >= hi =>
            {
                Some((*start, seq.as_str()))
            }
            _ => None,
        };
        let data = BrowserData {
            doc,
            features,
            sequence,
            selected: self.selection.clone(),
        };
        let clicked = browser::show(ui, &mut self.view, assembly, &data);
        drop(data);
        if let Some(cv) = clicked {
            self.select_variant(cv.contig, cv.pos, cv.rsid);
        }
    }

    // ---- workspace routing & tools -----------------------------------------

    fn set_tool(&mut self, tool: Tool) {
        self.tool = tool;
        let already = match tool {
            Tool::Genome => self.settings.seen_genome,
            Tool::Plasmid => self.settings.seen_plasmid,
            Tool::Crispr => self.settings.seen_crispr,
            Tool::Phenotype => self.settings.seen_phenotype,
            Tool::Ancestry => self.settings.seen_ancestry,
            Tool::Pharma => self.settings.seen_pharma,
        };
        if !already {
            match tool {
                Tool::Genome => self.settings.seen_genome = true,
                Tool::Plasmid => self.settings.seen_plasmid = true,
                Tool::Crispr => self.settings.seen_crispr = true,
                Tool::Phenotype => self.settings.seen_phenotype = true,
                Tool::Ancestry => self.settings.seen_ancestry = true,
                Tool::Pharma => self.settings.seen_pharma = true,
            }
            self.open_tutorial(tool.section());
        }
    }

    fn open_tutorial(&mut self, section: Section) {
        self.tutorial_section = section;
        self.tutorial_step = 0;
        self.show_tutorial = true;
    }

    fn plasmid_central(&mut self, ui: &mut egui::Ui) {
        let tier = self.settings.tier;
        if let (Some(Document::Sequences(recs)), Some(state)) =
            (self.document.as_ref(), self.plasmid.as_mut())
        {
            let idx = state.record_index.min(recs.len().saturating_sub(1));
            crate::plasmid::central(ui, &recs[idx], state, tier);
        }
    }

    fn plasmid_sidebar(&mut self, ui: &mut egui::Ui) {
        let tier = self.settings.tier;
        if let (Some(Document::Sequences(recs)), Some(state)) =
            (self.document.as_ref(), self.plasmid.as_mut())
            && let Some(idx) = crate::plasmid::sidebar(ui, recs, state, tier) {
                state.record_index = idx;
                state.recompute(&recs[idx]);
            }
    }

    fn plasmid_detail(&mut self, ui: &mut egui::Ui) {
        let tier = self.settings.tier;
        if let (Some(Document::Sequences(recs)), Some(state)) =
            (self.document.as_ref(), self.plasmid.as_ref())
        {
            let idx = state.record_index.min(recs.len().saturating_sub(1));
            crate::plasmid::detail(ui, &recs[idx], state, tier);
        }
    }

    fn crispr_sidebar(&mut self, ui: &mut egui::Ui) {
        let tier = self.settings.tier;
        let action = {
            let loaded: Option<(&str, &[u8])> = match &self.document {
                Some(Document::Sequences(recs)) if !recs.is_empty() => {
                    let idx = self
                        .plasmid
                        .as_ref()
                        .map(|p| p.record_index)
                        .unwrap_or(0)
                        .min(recs.len() - 1);
                    Some((recs[idx].name.as_str(), recs[idx].seq.as_slice()))
                }
                _ => None,
            };
            crate::crispr::sidebar(ui, &mut self.crispr, tier, loaded)
        };
        if let Some(crate::crispr::CrisprAction::LoadReference) = action {
            self.load_crispr_reference();
        }
    }

    fn load_crispr_reference(&mut self) {
        let Some(path) = rfd::FileDialog::new()
            .add_filter("FASTA reference", &["fa", "fasta", "fna", "gz", "txt"])
            .add_filter("All files", &["*"])
            .pick_file()
        else {
            return;
        };
        match gx_io::load_fasta(&path) {
            Ok(records) => {
                let name = records
                    .first()
                    .map(|r| r.name.clone())
                    .unwrap_or_else(|| "reference".to_string());
                // Concatenate records with N-spacers so off-targets can't span
                // a chromosome junction.
                let mut seq = Vec::new();
                for (i, r) in records.iter().enumerate() {
                    if i > 0 {
                        seq.extend_from_slice(&[b'N'; 25]);
                    }
                    seq.extend_from_slice(&r.seq);
                }
                let total = seq.len();
                let label = if records.len() > 1 {
                    format!("{name} +{} more", records.len() - 1)
                } else {
                    name
                };
                self.crispr.set_reference(label, seq);
                self.status = format!("Loaded off-target reference: {total} bp.");
            }
            Err(e) => self.status = format!("⚠ Reference load failed: {e}"),
        }
    }

    fn crispr_central(&mut self, ui: &mut egui::Ui) {
        crate::crispr::central(ui, &mut self.crispr, self.settings.tier);
    }

    fn crispr_detail(&mut self, ui: &mut egui::Ui) {
        crate::crispr::detail(ui, &mut self.crispr, self.settings.tier);
    }

    fn phenotype_sidebar(&mut self, ui: &mut egui::Ui) {
        let tier = self.settings.tier;
        match crate::phenotype::sidebar(ui, &mut self.phenotype, tier) {
            Some(PhenoRequest::Fetch(id)) => {
                if self.document.as_ref().and_then(|d| d.variant_doc()).is_none() {
                    self.phenotype.status = "Load a genome first (File ▸ Open).".into();
                } else if !self.settings.online {
                    self.phenotype.status = "Enable Online in Settings to download scores.".into();
                } else {
                    self.phenotype.status = format!("Fetching {id}…");
                    self.send(Request::FetchPgs(id));
                }
            }
            Some(PhenoRequest::LocalFile) => self.load_local_pgs(),
            None => {}
        }
    }

    fn phenotype_central(&mut self, ui: &mut egui::Ui) {
        crate::phenotype::central(ui, &self.phenotype, self.settings.tier);
    }

    fn phenotype_detail(&mut self, ui: &mut egui::Ui) {
        let tier = self.settings.tier;
        let doc = self.document.as_ref().and_then(|d| d.variant_doc());
        crate::phenotype::detail(ui, doc, tier);
    }

    /// Compute haplogroups and ancestry composition once for the current
    /// document. Cheap (a few hundred lookups), so we run it lazily and cache.
    fn ensure_ancestry(&mut self) {
        if self.ancestry.done {
            return;
        }
        if let Some(doc) = self.document.as_ref().and_then(|d| d.variant_doc()) {
            self.ancestry.maternal = gx_haplo::classify_maternal(&doc.store);
            self.ancestry.paternal = gx_haplo::classify_paternal(&doc.store);
            self.ancestry.composition = gx_ancestry::estimate(&doc.store);
            self.ancestry.done = true;
        }
    }

    fn ancestry_sidebar(&mut self, ui: &mut egui::Ui) {
        crate::ancestry::sidebar(ui, self.settings.tier);
    }

    fn ancestry_central(&mut self, ui: &mut egui::Ui) {
        self.ensure_ancestry();
        let has_doc = self.document.as_ref().and_then(|d| d.variant_doc()).is_some();
        crate::ancestry::central(ui, &self.ancestry, has_doc, self.settings.tier);
    }

    fn ancestry_detail(&mut self, ui: &mut egui::Ui) {
        self.ensure_ancestry();
        crate::ancestry::detail(ui, &self.ancestry, self.settings.tier);
    }

    /// Build the pharmacogenomics report once for the current document, lazily.
    fn ensure_pgx(&mut self) {
        if self.pgx.done {
            return;
        }
        if let Some(doc) = self.document.as_ref().and_then(|d| d.variant_doc()) {
            self.pgx.results = gx_pgx::report(&doc.store);
            self.pgx.selected = (!self.pgx.results.is_empty()).then_some(0);
            self.pgx.done = true;
        }
    }

    fn pgx_sidebar(&mut self, ui: &mut egui::Ui) {
        self.ensure_pgx();
        let has_doc = self.document.as_ref().and_then(|d| d.variant_doc()).is_some();
        crate::pgx::sidebar(ui, &mut self.pgx, has_doc, self.settings.tier);
    }

    fn pgx_central(&mut self, ui: &mut egui::Ui) {
        self.ensure_pgx();
        let has_doc = self.document.as_ref().and_then(|d| d.variant_doc()).is_some();
        crate::pgx::central(ui, &self.pgx, has_doc, self.settings.tier);
    }

    fn pgx_detail(&mut self, ui: &mut egui::Ui) {
        crate::pgx::detail(ui, self.settings.tier);
    }

    fn load_local_pgs(&mut self) {
        let Some(path) = rfd::FileDialog::new()
            .add_filter("PGS scoring file", &["txt", "gz", "tsv"])
            .add_filter("All files", &["*"])
            .pick_file()
        else {
            return;
        };
        let text = match read_pgs_text(&path) {
            Ok(t) => t,
            Err(e) => {
                self.phenotype.status = format!("Read error: {e}");
                return;
            }
        };
        match gx_pgs::parse(&text) {
            Ok(sf) => {
                if let Some(doc) = self.document.as_ref().and_then(|d| d.variant_doc()) {
                    let result = gx_pgs::apply(&sf, &doc.store);
                    self.phenotype.status =
                        format!("{}: {}/{} matched", sf.trait_name, result.matched, result.total);
                    self.phenotype.results.push(result);
                    self.phenotype.selected = Some(self.phenotype.results.len() - 1);
                } else {
                    self.phenotype.status = "Load a genome first.".into();
                }
            }
            Err(e) => self.phenotype.status = format!("Parse error: {e}"),
        }
    }
}

/// External "Learn more" links for a variant — authoritative NCBI/NIH/EBI
/// resources, keyed off the rsID (build-independent), the variant locus (gnomAD),
/// and the gene where known. These open in the user's browser; nothing here
/// leaves the machine until clicked.
fn resource_links(
    ui: &mut egui::Ui,
    rsid: Option<&str>,
    gene: Option<&str>,
    gnomad_id: Option<&str>,
    tier: ux::Tier,
) {
    ui.add_space(8.0);
    ui.label(egui::RichText::new("Learn more").strong());
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing.x = 10.0;
        if let Some(rsid) = rsid {
            ui.hyperlink_to("dbSNP", format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}"));
            ui.hyperlink_to("ClinVar", format!("https://www.ncbi.nlm.nih.gov/clinvar/?term={rsid}"));
            ui.hyperlink_to("PubMed", format!("https://pubmed.ncbi.nlm.nih.gov/?term={rsid}"));
            ui.hyperlink_to("GWAS Catalog", format!("https://www.ebi.ac.uk/gwas/variants/{rsid}"));
            ui.hyperlink_to("SNPedia", format!("https://www.snpedia.com/index.php/{rsid}"));
            if tier.at_least(ux::Tier::Intermediate) {
                ui.hyperlink_to(
                    "Ensembl",
                    format!("https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={rsid}"),
                );
            }
        }
        if let Some(id) = gnomad_id {
            ui.hyperlink_to(
                "gnomAD",
                format!("https://gnomad.broadinstitute.org/variant/{id}?dataset=gnomad_r4"),
            );
        }
    });

    // Gene-level resources on their own line.
    if let Some(g) = gene.filter(|g| !g.is_empty()) {
        ui.horizontal_wrapped(|ui| {
            ui.spacing_mut().item_spacing.x = 10.0;
            ui.label(egui::RichText::new(format!("{g}:")).color(palette::RULER_TEXT));
            ui.hyperlink_to(
                "MedlinePlus",
                format!("https://medlineplus.gov/genetics/gene/{}/", g.to_lowercase()),
            );
            ui.hyperlink_to(
                "NCBI Gene",
                format!("https://www.ncbi.nlm.nih.gov/gene/?term={g}%5Bsym%5D+AND+human%5Borgn%5D"),
            );
            ui.hyperlink_to(
                "ClinVar",
                format!("https://www.ncbi.nlm.nih.gov/clinvar/?term={g}%5Bgene%5D"),
            );
            ui.hyperlink_to("PubMed", format!("https://pubmed.ncbi.nlm.nih.gov/?term={g}"));
        });
    }
}

/// Serialise the (already-filtered) variant rows to CSV.
fn variants_to_csv(rows: &[&gx_annotate::VariantAnnotation]) -> String {
    fn field(s: &str) -> String {
        // Neutralise spreadsheet formula injection: annotation text comes from
        // remote databases, and a leading =, +, -, @ (or tab/CR) is executed as
        // a formula by Excel / Sheets / LibreOffice. Prefix such values with '.
        let s = if s.starts_with(['=', '+', '-', '@', '\t', '\r']) {
            format!("'{s}")
        } else {
            s.to_string()
        };
        if s.contains([',', '"', '\n']) {
            format!("\"{}\"", s.replace('"', "\"\""))
        } else {
            s
        }
    }
    let mut out =
        String::from("rsid,gene,consequence,clinvar_significance,conditions,gnomad_af,gwas_traits\n");
    for a in rows {
        let cols = [
            a.rsid.clone(),
            a.gene.clone().unwrap_or_default(),
            a.consequence.clone().unwrap_or_default(),
            a.clinical_significance.clone().unwrap_or_default(),
            a.conditions.join("; "),
            a.gnomad_af.map(|x| format!("{x}")).unwrap_or_default(),
            a.gwas_traits.join("; "),
        ];
        out.push_str(&cols.iter().map(|c| field(c)).collect::<Vec<_>>().join(","));
        out.push('\n');
    }
    out
}

/// The other human build (GRCh37 ⇄ GRCh38).
fn other_assembly(a: Assembly) -> Assembly {
    match a {
        Assembly::Grch38 => Assembly::Grch37,
        Assembly::Grch37 => Assembly::Grch38,
    }
}

/// A small reference of keyboard shortcuts and interactions.
fn shortcuts_window(ctx: &egui::Context, open: &mut bool) {
    let mut keep = true;
    egui::Window::new("Keyboard shortcuts")
        .open(&mut keep)
        .resizable(false)
        .collapsible(false)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            egui::Grid::new("shortcuts").num_columns(2).striped(true).spacing([24.0, 6.0]).show(
                ui,
                |ui| {
                    for (k, what) in [
                        ("Ctrl/Cmd + O", "Open a genome / sequence file"),
                        ("Drag & drop", "Open a file by dropping it on the window"),
                        ("Drag in browser", "Pan the genome / plasmid view"),
                        ("Scroll", "Zoom the view in and out"),
                        ("Click a lollipop", "Inspect that variant"),
                        ("Esc", "Close this window / Settings / a tutorial"),
                    ] {
                        ui.label(egui::RichText::new(k).strong().monospace());
                        ui.label(what);
                        ui.end_row();
                    }
                },
            );
        });
    if !keep {
        *open = false;
    }
}

/// Paint the GenomeForge brand mark: a small two-strand DNA helix with rungs,
/// in the brand accent. Drawn (not an emoji) so it always renders and stays
/// crisp at any DPI, echoing the app icon.
fn draw_brand_mark(ui: &mut egui::Ui) {
    use theme::tokens;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(15.0, 19.0), egui::Sense::hover());
    let p = ui.painter();
    let cx = rect.center().x;
    let (top, bot) = (rect.top() + 2.5, rect.bottom() - 2.5);
    let amp = 5.0;
    let turns = std::f32::consts::TAU * 1.25;
    let n = 18;
    let strand = |sign: f32| -> Vec<egui::Pos2> {
        (0..=n)
            .map(|i| {
                let t = i as f32 / n as f32;
                let ph = t * turns;
                egui::pos2(cx + sign * amp * ph.sin(), top + t * (bot - top))
            })
            .collect()
    };
    let a = strand(1.0);
    let b = strand(-1.0);
    // base-pair rungs first (behind the strands)
    for i in (2..n).step_by(4) {
        p.line_segment([a[i], b[i]], egui::Stroke::new(1.0, tokens::TEXT_MUTED));
    }
    for w in a.windows(2) {
        p.line_segment([w[0], w[1]], egui::Stroke::new(1.8, tokens::ACCENT));
    }
    for w in b.windows(2) {
        p.line_segment([w[0], w[1]], egui::Stroke::new(1.8, tokens::ACCENT_HOVER));
    }
}

impl eframe::App for GenomeForgeApp {
    // egui 0.34 deprecated the `SidePanel`/`TopBottomPanel` aliases and
    // `default_width`/`default_height` in favour of `Panel::*`/`default_size`.
    // We pin egui 0.34, so the deprecated-but-present forms are fine until a
    // deliberate egui bump.
    #[allow(deprecated)]
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.drain_worker();

        if self.settings.online != self.last_online {
            self.send(Request::SetOnline(self.settings.online));
            if self.settings.online {
                // Reconnected: let regions/annotations skipped while offline retry.
                self.refdata.fetched_tiles.clear();
                self.refdata.requested_seq.clear();
                self.requested_annot.clear();
            }
            self.last_online = self.settings.online;
        }

        // Ctrl/Cmd+O opens a file from any workspace.
        if ui.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::O)) {
            self.open_file_dialog();
        }
        // Drag-and-drop a genome/sequence file onto the window to open it.
        if let Some(path) = ui
            .input(|i| i.raw.dropped_files.iter().find_map(|f| f.path.clone()))
        {
            self.import_file(path);
        }
        // Esc dismisses the transient windows.
        if (self.show_settings || self.show_tutorial || self.show_about || self.show_shortcuts)
            && ui.input(|i| i.key_pressed(egui::Key::Escape))
        {
            self.show_settings = false;
            self.show_tutorial = false;
            self.show_about = false;
            self.show_shortcuts = false;
        }

        egui::TopBottomPanel::top("menubar").show_inside(ui, |ui| self.menu_bar(ui));
        let busy = self.inflight > 0;
        egui::TopBottomPanel::bottom("statusbar").show_inside(ui, |ui| {
            ui.horizontal(|ui| {
                if busy {
                    ui.add(egui::Spinner::new().size(13.0));
                }
                ui.label(
                    egui::RichText::new(&self.status)
                        .size(11.0)
                        .color(palette::RULER_TEXT),
                );
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        egui::RichText::new(concat!("GenomeForge v", env!("CARGO_PKG_VERSION")))
                            .size(11.0)
                            .color(theme::tokens::TEXT_FAINT),
                    );
                });
            });
        });
        if busy {
            ui.ctx().request_repaint(); // animate the spinner between replies
        }
        let tool = self.tool;
        egui::SidePanel::left("sidebar")
            .resizable(true)
            .default_width(250.0)
            .show_inside(ui, |ui| {
                // Scroll when a workspace's controls run taller than the window.
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| match tool {
                        Tool::Genome => self.sidebar(ui),
                        Tool::Plasmid => self.plasmid_sidebar(ui),
                        Tool::Crispr => self.crispr_sidebar(ui),
                        Tool::Phenotype => self.phenotype_sidebar(ui),
                        Tool::Ancestry => self.ancestry_sidebar(ui),
                        Tool::Pharma => self.pgx_sidebar(ui),
                    });
            });
        egui::SidePanel::right("detail")
            .resizable(true)
            .default_width(320.0)
            .show_inside(ui, |ui| {
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| match tool {
                        Tool::Genome => self.detail_panel(ui),
                        Tool::Plasmid => self.plasmid_detail(ui),
                        Tool::Crispr => self.crispr_detail(ui),
                        Tool::Phenotype => self.phenotype_detail(ui),
                        Tool::Ancestry => self.ancestry_detail(ui),
                        Tool::Pharma => self.pgx_detail(ui),
                    });
            });

        if tool == Tool::Genome
            && self.show_table
            && self.document.as_ref().and_then(|d| d.variant_doc()).is_some()
        {
            let clicked = egui::TopBottomPanel::bottom("table")
                .resizable(true)
                .default_height(200.0)
                .show_inside(ui, |ui| self.table_panel(ui))
                .inner;
            if let Some(rsid) = clicked {
                self.select_by_rsid(rsid);
            }
        }

        egui::CentralPanel::default().show_inside(ui, |ui| match tool {
            Tool::Genome => self.draw_central(ui),
            Tool::Plasmid => self.plasmid_central(ui),
            Tool::Crispr => self.crispr_central(ui),
            Tool::Phenotype => self.phenotype_central(ui),
            Tool::Ancestry => self.ancestry_central(ui),
            Tool::Pharma => self.pgx_central(ui),
        });

        if self.show_settings {
            let ctx = ui.ctx().clone();
            self.settings_window(&ctx);
        }

        if self.show_tutorial {
            let ctx = ui.ctx().clone();
            tutorial::show(
                &ctx,
                self.tutorial_section,
                &mut self.tutorial_step,
                &mut self.show_tutorial,
            );
        }

        if self.show_about {
            let ctx = ui.ctx().clone();
            self.about_window(&ctx);
        }
        if self.show_shortcuts {
            let ctx = ui.ctx().clone();
            shortcuts_window(&ctx, &mut self.show_shortcuts);
        }
    }

    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        eframe::set_value(storage, "settings", &self.settings);
    }
}

// ---- helpers ---------------------------------------------------------------

fn row(ui: &mut egui::Ui, label: &str, value: &str) {
    ui.label(label);
    ui.label(value);
    ui.end_row();
}

/// Significant-figure-aware allele-frequency formatting (no trailing junk).
fn format_af(af: f64) -> String {
    let pct = af * 100.0;
    if pct <= 0.0 {
        "0%".to_string()
    } else if pct >= 1.0 {
        format!("{pct:.2}%")
    } else if pct >= 0.01 {
        format!("{pct:.3}%")
    } else {
        format!("{pct:.4}%")
    }
}

/// A compact colour key for the genome-browser variant track.
fn variant_legend(ui: &mut egui::Ui) {
    for (label, color) in [
        ("Pathogenic", palette::PATHOGENIC),
        ("Likely pathogenic", palette::LIKELY_PATH),
        ("Uncertain / conflicting", palette::VUS),
        ("Benign", palette::BENIGN),
        ("Not yet annotated", palette::VARIANT),
    ] {
        ui.horizontal(|ui| {
            let (rect, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
            ui.painter().rect_filled(rect, 2.0, color);
            ui.label(egui::RichText::new(label).size(11.0));
        });
    }
}


/// Parse a `chr:pos` (or `chr:pos-end`) locus into `(contig, 0-based pos)`.
fn parse_locus(s: &str) -> Option<(String, u64)> {
    let (chrom, rest) = s.split_once(':')?;
    let pos_str = rest.split(['-', '–']).next()?.replace([',', '_'], "");
    let pos: u64 = pos_str.trim().parse().ok()?;
    Some((gx_core::normalize_contig(chrom), pos.saturating_sub(1)))
}

/// Human-readable byte size (B / KiB / MiB).
fn human_bytes(n: u64) -> String {
    const KIB: u64 = 1024;
    const MIB: u64 = KIB * 1024;
    if n >= MIB {
        format!("{:.1} MiB", n as f64 / MIB as f64)
    } else if n >= KIB {
        format!("{:.0} KiB", n as f64 / KIB as f64)
    } else {
        format!("{n} B")
    }
}

/// Coarse "n units ago" from a duration, for cache-freshness display.
fn humanize_ago(d: std::time::Duration) -> String {
    let s = d.as_secs();
    let (n, unit) = if s < 60 {
        return "just now".to_string();
    } else if s < 3600 {
        (s / 60, "minute")
    } else if s < 86_400 {
        (s / 3600, "hour")
    } else {
        (s / 86_400, "day")
    };
    format!("{n} {unit}{} ago", if n == 1 { "" } else { "s" })
}

fn cache_path() -> PathBuf {
    let dir = directories::ProjectDirs::from("io", "GenomeForge", "GenomeForge")
        .map(|d| d.cache_dir().to_path_buf())
        .unwrap_or_else(|| std::env::temp_dir().join("genomeforge"));
    let _ = std::fs::create_dir_all(&dir);
    dir.join("cache.redb")
}

/// Read a (possibly gzipped) PGS scoring file into text.
fn read_pgs_text(path: &std::path::Path) -> Result<String, String> {
    use std::io::Read;
    let bytes = std::fs::read(path).map_err(|e| e.to_string())?;
    let gz = path.extension().map(|e| e.eq_ignore_ascii_case("gz")).unwrap_or(false)
        || bytes.starts_with(&[0x1f, 0x8b]);
    if gz {
        let mut s = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut s)
            .map_err(|e| e.to_string())?;
        Ok(s)
    } else {
        String::from_utf8(bytes).map_err(|e| e.to_string())
    }
}
