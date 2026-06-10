//! CRISPR studio workspace: guide design, off-target search, edit simulation,
//! and AAV cargo planning. Tier-aware (Beginner/Intermediate/Expert).

use eframe::egui;
use egui_extras::{Column, TableBuilder};

use gx_crispr::{
    AavCargo, BaseEditor, Guide, Nuclease, OffTarget, PrimeEdit, base_edit, find_guides_for,
    find_matches, find_offtargets, hdr_replace, nhej_deletion, prime_edit, specificity,
};

use crate::theme::palette;
use crate::ux::{self, Tier};

/// Action the CRISPR sidebar asks the app to perform (file I/O lives in the app).
pub enum CrisprAction {
    LoadReference,
}

pub struct CrisprState {
    pub target: Vec<u8>,
    pub target_name: String,
    pub input_text: String,
    /// Nuclease / PAM the guides are enumerated for.
    pub nuclease: Nuclease,
    pub guides: Vec<Guide>,
    pub selected: Option<usize>,
    pub offtargets: Vec<OffTarget>,
    pub spec: f64,
    pub max_mm: usize,
    pub min_on: f32,
    pub del_len: usize,
    pub donor: String,
    /// Base the prime-edit demo installs at the cut site.
    pub pe_to: char,
    pub aav: AavCargo,
    last_key: Option<(usize, usize)>,
    /// A loaded reference (chromosome / genome FASTA) to search for off-targets.
    pub reference: Option<(String, Vec<u8>)>,
    pub ref_offtargets: Vec<OffTarget>,
    /// Guide index the `ref_offtargets` were computed for.
    ref_key: Option<usize>,
}

impl Default for CrisprState {
    fn default() -> Self {
        Self {
            target: Vec::new(),
            target_name: String::new(),
            input_text: String::new(),
            nuclease: Nuclease::SpCas9,
            guides: Vec::new(),
            selected: None,
            offtargets: Vec::new(),
            spec: 100.0,
            max_mm: 3,
            min_on: 0.0,
            del_len: 1,
            donor: String::new(),
            pe_to: 'A',
            aav: AavCargo::default(),
            last_key: None,
            reference: None,
            ref_offtargets: Vec::new(),
            ref_key: None,
        }
    }
}

impl CrisprState {
    pub fn set_target(&mut self, seq: Vec<u8>, name: String) {
        // Normalise to uppercase ACGTN so every downstream byte index is ASCII
        // (a loaded FASTA/GenBank can contain lowercase, IUPAC, or other bytes;
        // the paste path already filters, but loaded sequences bypass it).
        self.target = seq
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
            .collect();
        self.target_name = name;
        self.reanalyze();
        self.aav.transgene_bp = self.target.len() as u32;
    }

    /// Re-find guides for the current target + nuclease (used on a nuclease change).
    pub fn reanalyze(&mut self) {
        self.guides = find_guides_for(&self.target, self.nuclease);
        self.selected = None;
        self.offtargets.clear();
        self.spec = 100.0;
        self.last_key = None;
        self.ref_offtargets.clear();
        self.ref_key = None;
    }

    pub fn set_reference(&mut self, name: String, seq: Vec<u8>) {
        self.reference = Some((name, seq));
        self.ref_offtargets.clear();
        self.ref_key = None;
    }

    /// Search the loaded reference for off-targets of the selected guide.
    fn search_reference(&mut self) {
        if let (Some(i), Some((_, refseq))) = (self.selected, self.reference.as_ref())
            && let Some(g) = self.guides.get(i)
        {
            self.ref_offtargets = find_matches(g, refseq, self.max_mm);
            self.ref_key = Some(i);
        }
    }

    /// Specificity over target off-targets, plus reference hits when current.
    fn combined_specificity(&self) -> f64 {
        if self.ref_key == self.selected && !self.ref_offtargets.is_empty() {
            let mut all = self.offtargets.clone();
            all.extend(self.ref_offtargets.iter().cloned());
            specificity(&all)
        } else {
            self.spec
        }
    }

    fn select(&mut self, i: usize) {
        self.selected = Some(i);
        self.last_key = None;
        self.ref_offtargets.clear();
        self.ref_key = None;
    }

    fn ensure_offtargets(&mut self) {
        if let Some(i) = self.selected {
            let key = (i, self.max_mm);
            if self.last_key != Some(key) {
                if let Some(g) = self.guides.get(i) {
                    self.offtargets = find_offtargets(g, &self.target, self.max_mm);
                    self.spec = specificity(&self.offtargets);
                }
                self.last_key = Some(key);
            }
        }
    }
}

fn efficiency_word(score: f64) -> (&'static str, egui::Color32) {
    if score >= 0.6 {
        ("High", palette::BENIGN)
    } else if score >= 0.4 {
        ("Medium", palette::VUS)
    } else {
        ("Low", palette::PATHOGENIC)
    }
}

fn off_target_row(ui: &mut egui::Ui, o: &OffTarget, tier: Tier) {
    let (word, color) = if o.cfd >= 0.5 {
        ("high", palette::PATHOGENIC)
    } else if o.cfd >= 0.1 {
        ("med", palette::VUS)
    } else {
        ("low", palette::RULER_TEXT)
    };
    let detail = if tier.is_expert() {
        format!("CFD {:.3}", o.cfd)
    } else {
        format!("{word} risk")
    };
    ui.label(
        egui::RichText::new(format!(
            "{}{}  {} mm  ·  {}",
            o.strand.symbol(),
            o.start + 1,
            o.mismatches,
            detail
        ))
        .color(color)
        .size(11.0),
    );
}

// ============================ sidebar =======================================

pub fn sidebar(
    ui: &mut egui::Ui,
    state: &mut CrisprState,
    tier: Tier,
    loaded: Option<(&str, &[u8])>,
) -> Option<CrisprAction> {
    let mut action = None;
    ui.add_space(4.0);
    ui.heading("CRISPR studio");
    ux::explain_beginner(
        ui,
        tier,
        "CRISPR-Cas9 is molecular scissors. You pick a 20-letter 'guide' and Cas9 cuts the \
         matching DNA. Choose a target to begin.",
    );
    ui.separator();

    if let Some((name, seq)) = loaded
        && ui.button(format!("→ Use loaded sequence ({name})")).clicked() {
            state.set_target(seq.to_vec(), name.to_string());
        }
    ui.label("…or paste DNA and analyze:");
    ui.add(
        egui::TextEdit::multiline(&mut state.input_text)
            .desired_rows(3)
            .hint_text("Paste ACGT… sequence"),
    );
    if ui.button("Analyze pasted sequence").clicked() {
        let cleaned: Vec<u8> = state
            .input_text
            .bytes()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
            .collect();
        if cleaned.len() >= 23 {
            state.set_target(cleaned, "pasted".to_string());
        }
    }

    // Nuclease / PAM picker — re-finds guides for the current target on change.
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        ui.label("Nuclease:");
        let before = state.nuclease;
        egui::ComboBox::from_id_salt("nuclease")
            .selected_text(state.nuclease.label())
            .show_ui(ui, |ui| {
                for n in Nuclease::all() {
                    ui.selectable_value(&mut state.nuclease, n, n.label());
                }
            });
        if state.nuclease != before && !state.target.is_empty() {
            state.reanalyze();
        }
    });
    ux::explain_beginner(
        ui,
        tier,
        "Different nucleases recognise different short 'PAM' motifs next to the cut, which \
         changes where you can target. SpCas9 (NGG) is the classic one.",
    );

    ui.separator();
    if state.target.is_empty() {
        ui.label(
            egui::RichText::new("No target selected yet.")
                .color(palette::RULER_TEXT),
        );
    } else {
        ui.label(format!("Target: {} · {} bp", state.target_name, state.target.len()));
        ui.label(
            egui::RichText::new(format!("{} candidate guides", state.guides.len()))
                .size(11.0)
                .color(palette::RULER_TEXT),
        );
        ui.add(egui::Slider::new(&mut state.min_on, 0.0..=1.0).text("min on-target"));
        if tier.at_least(Tier::Intermediate) {
            ui.add(egui::Slider::new(&mut state.max_mm, 0..=4).text("off-target mismatches"));
        }
    }

    ui.separator();
    egui::CollapsingHeader::new("🌐 Off-target reference")
        .default_open(false)
        .show(ui, |ui| {
            ux::explain(
                ui,
                tier,
                "By default off-targets are searched only in the target above. Load a reference \
                 (a chromosome or whole-genome FASTA) to check off-targets across it.",
            );
            match &state.reference {
                Some((name, seq)) => {
                    ui.label(
                        egui::RichText::new(format!("Loaded: {name} · {} bp", seq.len()))
                            .size(11.0)
                            .color(palette::RULER_TEXT),
                    );
                }
                None => {
                    ui.label(egui::RichText::new("No reference loaded.").size(11.0).color(palette::RULER_TEXT));
                }
            }
            if ui.button("Load reference FASTA…").clicked() {
                action = Some(CrisprAction::LoadReference);
            }
        });

    ui.separator();
    egui::CollapsingHeader::new("AAV cargo planner")
        .default_open(tier.is_beginner())
        .show(ui, |ui| aav_ui(ui, &mut state.aav, tier));

    action
}

fn aav_ui(ui: &mut egui::Ui, aav: &mut AavCargo, tier: Tier) {
    ux::explain(
        ui,
        tier,
        "AAV is a common gene-delivery virus, but it can only carry ~4.7 kb of DNA. Check that \
         your promoter + gene + tail fit inside.",
    );
    ui.add(egui::Slider::new(&mut aav.promoter_bp, 0..=2000).text("promoter bp"));
    ui.add(egui::Slider::new(&mut aav.transgene_bp, 0..=6000).text("transgene bp"));
    ui.add(egui::Slider::new(&mut aav.polya_bp, 0..=600).text("polyA bp"));

    let frac = aav.fraction().min(1.5) as f32;
    let bar = egui::ProgressBar::new(frac.min(1.0))
        .text(format!("{} / {} bp", aav.total(), AavCargo::CAPACITY))
        .fill(if aav.fits() { palette::BENIGN } else { palette::PATHOGENIC });
    ui.add(bar);
    if aav.fits() {
        ui.label(
            egui::RichText::new(format!("✓ Fits — {} bp headroom", aav.headroom()))
                .color(palette::BENIGN),
        );
    } else {
        ui.label(
            egui::RichText::new(format!("✗ Over capacity by {} bp", -aav.headroom()))
                .color(palette::PATHOGENIC),
        );
    }
}

// ============================ central (guide table) =========================

pub fn central(ui: &mut egui::Ui, state: &mut CrisprState, tier: Tier) {
    if state.target.is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label(
                egui::RichText::new("Choose a target sequence in the sidebar to design guides.")
                    .color(palette::RULER_TEXT),
            );
        });
        return;
    }

    ux::explain(
        ui,
        tier,
        "Each row is a spot Cas9 could cut (next to an 'NGG' PAM). On-target estimates cutting \
         efficiency. Click a guide to inspect off-targets and simulate its edit.",
    );

    let filtered: Vec<usize> = state
        .guides
        .iter()
        .enumerate()
        .filter(|(_, g)| g.on_score >= state.min_on as f64)
        .map(|(i, _)| i)
        .collect();
    ui.label(
        egui::RichText::new(format!("{} guides shown (of {})", filtered.len(), state.guides.len()))
            .size(11.0)
            .color(palette::RULER_TEXT),
    );

    let expert = tier.is_expert();
    let mut clicked = None;
    let mut builder = TableBuilder::new(ui)
        .striped(true)
        .column(Column::auto().at_least(64.0)) // position
        .column(Column::auto().at_least(36.0)); // strand
    builder = builder.column(Column::remainder().at_least(200.0)); // guide + PAM
    builder = builder.column(Column::auto().at_least(54.0)); // GC
    builder = builder.column(Column::auto().at_least(if expert { 70.0 } else { 90.0 })); // on-target
    builder = builder.column(Column::auto().at_least(70.0)); // flags

    builder
        .header(18.0, |mut h| {
            h.col(|ui| { ui.strong("Pos"); });
            h.col(|ui| { ui.strong("Str"); });
            h.col(|ui| { ui.strong("Guide (5'→3') + PAM"); });
            h.col(|ui| { ui.strong("GC"); });
            h.col(|ui| { ui.strong(if expert { "On-tgt" } else { "Efficiency" }); });
            h.col(|ui| { ui.strong("Flags"); });
        })
        .body(|mut body| {
            for &i in &filtered {
                let g = &state.guides[i];
                body.row(18.0, |mut row| {
                    row.col(|ui| {
                        if ui.link((g.start + 1).to_string()).clicked() {
                            clicked = Some(i);
                        }
                    });
                    row.col(|ui| {
                        ui.label(g.strand.symbol().to_string());
                    });
                    row.col(|ui| {
                        ui.label(
                            egui::RichText::new(format!("{} {}", g.protospacer, g.pam))
                                .monospace()
                                .size(11.0),
                        );
                    });
                    row.col(|ui| {
                        ui.label(format!("{:.0}%", g.gc * 100.0));
                    });
                    row.col(|ui| {
                        let (word, color) = efficiency_word(g.on_score);
                        if expert {
                            ui.label(format!("{:.2}", g.on_score));
                        } else {
                            ui.label(egui::RichText::new(word).color(color));
                        }
                    });
                    row.col(|ui| {
                        if g.poly_t {
                            ui.label(egui::RichText::new("poly-T").color(palette::VUS))
                                .on_hover_text("Run of T's can terminate the guide RNA early.");
                        }
                    });
                });
            }
        });

    if let Some(i) = clicked {
        state.select(i);
    }
}

// ============================ detail (selected guide) =======================

pub fn detail(ui: &mut egui::Ui, state: &mut CrisprState, tier: Tier) {
    ui.add_space(4.0);
    ui.heading("Guide details");
    ui.separator();
    state.ensure_offtargets();

    let Some(i) = state.selected else {
        ui.label(
            egui::RichText::new("Click a guide in the table to inspect it.")
                .color(palette::RULER_TEXT),
        );
        return;
    };
    let g = state.guides[i].clone();

    egui::Grid::new("guide_detail").num_columns(2).striped(true).show(ui, |ui| {
        kv(ui, "Protospacer", &g.protospacer);
        kv(ui, "PAM", &g.pam);
        kv(ui, "Strand", g.strand.symbol().to_string().as_str());
        kv(ui, "Cut site", &format!("{}", g.cut_site + 1));
        kv(ui, "GC", &format!("{:.0}%", g.gc * 100.0));
        if tier.is_expert() {
            kv(ui, "On-target (Doench '14)", &format!("{:.3}", g.on_score));
        } else {
            let (w, _) = efficiency_word(g.on_score);
            kv(ui, "Efficiency", w);
        }
    });

    if ui
        .button("Copy guide")
        .on_hover_text("Copy the protospacer sequence to the clipboard")
        .clicked()
    {
        ui.ctx().copy_text(g.protospacer.clone());
    }

    ui.add_space(6.0);
    ux::explain(
        ui,
        tier,
        "Specificity shows how uniquely this guide matches. 100 means no near-matches were found \
         (in the target, plus any reference you've searched).",
    );
    let spec = state.combined_specificity();
    let spec_color = if spec >= 80.0 {
        palette::BENIGN
    } else if spec >= 50.0 {
        palette::VUS
    } else {
        palette::PATHOGENIC
    };
    ui.label(
        egui::RichText::new(format!("Specificity: {spec:.0}/100"))
            .strong()
            .color(spec_color),
    );

    let mut do_ref_search = false;
    egui::CollapsingHeader::new(format!("Off-targets in target ({})", state.offtargets.len()))
        .default_open(true)
        .show(ui, |ui| {
            if state.offtargets.is_empty() {
                ui.label(egui::RichText::new("None in the target sequence ✓").color(palette::BENIGN));
            } else {
                egui::ScrollArea::vertical().max_height(120.0).id_salt("ot_target").show(ui, |ui| {
                    for o in state.offtargets.iter().take(50) {
                        off_target_row(ui, o, tier);
                    }
                });
            }

            ui.add_space(4.0);
            ui.separator();
            match &state.reference {
                None => {
                    ui.label(
                        egui::RichText::new(
                            "Load a reference (sidebar ▸ Off-target reference) to also search a \
                             chromosome / whole genome.",
                        )
                        .size(11.0)
                        .color(palette::RULER_TEXT),
                    );
                }
                Some((name, seq)) => {
                    if ui.button(format!("🌐 Search {name} ({} bp)", seq.len())).clicked() {
                        do_ref_search = true;
                    }
                    if seq.len() > 50_000_000 {
                        ui.label(
                            egui::RichText::new("Large reference — search may take a few seconds.")
                                .size(11.0)
                                .color(palette::VUS),
                        );
                    }
                    if state.ref_key == state.selected {
                        if state.ref_offtargets.is_empty() {
                            ui.label(
                                egui::RichText::new(format!("No off-targets in {name} ✓"))
                                    .size(11.0)
                                    .color(palette::BENIGN),
                            );
                        } else {
                            ui.label(
                                egui::RichText::new(format!(
                                    "{} off-targets in {name}:",
                                    state.ref_offtargets.len()
                                ))
                                .size(11.0)
                                .color(palette::RULER_TEXT),
                            );
                            egui::ScrollArea::vertical().max_height(120.0).id_salt("ot_ref").show(ui, |ui| {
                                for o in state.ref_offtargets.iter().take(100) {
                                    off_target_row(ui, o, tier);
                                }
                            });
                        }
                    }
                }
            }
        });
    if do_ref_search {
        state.search_reference();
    }

    edit_sim(ui, state, &g, tier);

    ui.add_space(8.0);
    ux::disclaimer(ui);
}

fn edit_sim(ui: &mut egui::Ui, state: &mut CrisprState, g: &Guide, tier: Tier) {
    egui::CollapsingHeader::new("✂ Simulate edit")
        .default_open(true)
        .show(ui, |ui| {
            ux::explain(
                ui,
                tier,
                "After Cas9 cuts, the cell repairs the break. NHEJ often deletes a few letters \
                 (which can disable a gene); HDR can paste in a donor sequence precisely.",
            );

            ui.add(egui::Slider::new(&mut state.del_len, 1..=12).text("NHEJ deletion (bp)"));
            let nhej = nhej_deletion(&state.target, g.cut_site, state.del_len);
            show_edit(ui, &nhej.before, &nhej.after, nhej.frameshift);

            ui.add_space(6.0);
            ui.label("HDR donor (replaces the cut region):");
            ui.add(
                egui::TextEdit::singleline(&mut state.donor)
                    .hint_text("donor DNA, e.g. ACGT…")
                    .desired_width(f32::INFINITY),
            );
            if !state.donor.trim().is_empty() {
                let donor: Vec<u8> = state
                    .donor
                    .bytes()
                    .filter(|b| b.is_ascii_alphabetic())
                    .collect();
                let lo = g.cut_site.saturating_sub(3);
                let hi = (g.cut_site + 3).min(state.target.len());
                let hdr = hdr_replace(&state.target, lo, hi, &donor);
                show_edit(ui, &hdr.before, &hdr.after, hdr.frameshift);
            }
        });

    // ---- Base & prime editing (no double-strand break) --------------------
    ui.add_space(6.0);
    egui::CollapsingHeader::new("Base & prime editing")
        .default_open(false)
        .show(ui, |ui| {
            ux::explain(
                ui,
                tier,
                "Base editors change one letter within a small window of the guide, without \
                 cutting. Prime editing can install an arbitrary small edit at the nick.",
            );
            let mut shown_window = false;
            for editor in [BaseEditor::CytosineCBE, BaseEditor::AdenineABE] {
                let out = base_edit(g, editor);
                if !shown_window {
                    ui.label(
                        egui::RichText::new(format!(
                            "Edit window: protospacer positions {}–{}",
                            out.window.0 + 1,
                            out.window.1
                        ))
                        .size(10.0)
                        .color(palette::RULER_TEXT),
                    );
                    shown_window = true;
                }
                ui.add_space(2.0);
                ui.label(egui::RichText::new(out.editor.label()).strong());
                if out.edits.is_empty() {
                    ui.label(
                        egui::RichText::new("No editable base in the window.")
                            .size(11.0)
                            .color(palette::RULER_TEXT),
                    );
                } else {
                    ui.label(egui::RichText::new(format!("guide:  {}", g.protospacer)).monospace().size(11.0));
                    ui.label(
                        egui::RichText::new(format!("result: {}", out.edited_protospacer))
                            .monospace()
                            .size(11.0)
                            .color(palette::ACCENT),
                    );
                    let edits: Vec<String> = out
                        .edits
                        .iter()
                        .map(|e| format!("pos {} {}→{}", e.pos_in_protospacer + 1, e.from, e.to))
                        .collect();
                    ui.label(egui::RichText::new(edits.join(", ")).size(11.0).color(palette::RULER_TEXT));
                }
            }

            ui.add_space(6.0);
            ui.label(egui::RichText::new("Prime editing").strong());
            ui.horizontal(|ui| {
                ui.label(format!("Install at cut site ({}):", g.cut_site + 1));
                for b in ['A', 'C', 'G', 'T'] {
                    ui.selectable_value(&mut state.pe_to, b, b.to_string());
                }
            });
            let pe = prime_edit(
                &state.target,
                &PrimeEdit::Substitution { pos: g.cut_site, to: state.pe_to },
            );
            let lo = g.cut_site.saturating_sub(8);
            let hi = (g.cut_site + 9).min(pe.edited.len());
            show_edit(
                ui,
                &pe.original[lo..hi.min(pe.original.len())],
                &pe.edited[lo..hi],
                None,
            );
        });
}

fn show_edit(ui: &mut egui::Ui, before: &str, after: &str, frameshift: Option<bool>) {
    ui.label(egui::RichText::new(format!("before: {before}")).monospace().size(11.0));
    ui.label(
        egui::RichText::new(format!("after:  {after}"))
            .monospace()
            .size(11.0)
            .color(palette::ACCENT),
    );
    match frameshift {
        Some(true) => {
            ui.label(egui::RichText::new("⚠ frameshift — likely disrupts the protein").color(palette::PATHOGENIC));
        }
        Some(false) => {
            ui.label(egui::RichText::new("in-frame change").color(palette::BENIGN));
        }
        None => {}
    }
}

fn kv(ui: &mut egui::Ui, k: &str, v: &str) {
    ui.label(k);
    ui.label(egui::RichText::new(v).monospace());
    ui.end_row();
}
