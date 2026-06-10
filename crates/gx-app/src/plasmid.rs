//! Plasmid designer workspace: circular + linear sequence maps, restriction
//! mapping, ORF/translation and primer stats — rendered with egui's painter.

use std::collections::HashSet;
use std::f32::consts::{FRAC_PI_2, TAU};

use eframe::egui;
use egui::{Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, pos2, vec2};

use gx_core::FeatureKind;
use gx_io::SequenceRecord;
use gx_plasmid::{Orf, RestrictionSite, find_orfs, find_sites, site_counts, unique_cutters};

use crate::theme::{self, palette};
use crate::ux::{self, Tier};

const ORF_MIN_AA: usize = 40;

/// What the user currently has selected in the plasmid view.
#[derive(Clone, PartialEq)]
pub enum Selected {
    None,
    Feature(usize),
    Site(usize),
    Orf(usize),
}

pub struct PlasmidState {
    pub record_index: usize,
    pub sites: Vec<RestrictionSite>,
    pub orfs: Vec<Orf>,
    unique_names: HashSet<&'static str>,
    pub show_unique_only: bool,
    pub show_orfs: bool,
    pub selected: Selected,
    rotation: f32,
    lin_offset: f64,
    lin_bp_per_px: f64,
}

impl PlasmidState {
    pub fn new(record_index: usize, rec: &SequenceRecord) -> Self {
        let mut s = Self {
            record_index,
            sites: Vec::new(),
            orfs: Vec::new(),
            unique_names: HashSet::new(),
            show_unique_only: true,
            show_orfs: true,
            selected: Selected::None,
            rotation: 0.0,
            lin_offset: 0.0,
            lin_bp_per_px: (rec.len().max(1) as f64 / 900.0).max(0.05),
        };
        s.recompute(rec);
        s
    }

    pub fn recompute(&mut self, rec: &SequenceRecord) {
        self.sites = find_sites(&rec.seq, rec.circular);
        self.orfs = find_orfs(&rec.seq, ORF_MIN_AA);
        let counts = site_counts(&self.sites);
        self.unique_names = counts
            .into_iter()
            .filter(|(_, c)| *c == 1)
            .map(|(n, _)| n)
            .collect();
        self.selected = Selected::None;
    }

    fn visible_sites(&self) -> Vec<(usize, &RestrictionSite)> {
        self.sites
            .iter()
            .enumerate()
            .filter(|(_, s)| !self.show_unique_only || self.unique_names.contains(s.enzyme))
            .collect()
    }
}

fn feature_color(kind: &FeatureKind) -> Color32 {
    match kind {
        FeatureKind::Gene | FeatureKind::Cds => Color32::from_rgb(96, 142, 222),
        FeatureKind::Regulatory => Color32::from_rgb(96, 178, 110),
        FeatureKind::Other(s) => {
            let s = s.to_ascii_lowercase();
            if s.contains("promoter") {
                Color32::from_rgb(96, 178, 110)
            } else if s.contains("terminator") {
                Color32::from_rgb(224, 110, 96)
            } else if s.contains("origin") || s.contains("rep_origin") {
                Color32::from_rgb(168, 130, 220)
            } else if s.contains("primer") {
                Color32::from_rgb(226, 168, 90)
            } else {
                Color32::from_rgb(150, 154, 160)
            }
        }
        _ => Color32::from_rgb(150, 154, 160),
    }
}

// ============================ central maps ==================================

pub fn central(ui: &mut egui::Ui, rec: &SequenceRecord, state: &mut PlasmidState, tier: Tier) {
    ux::explain_beginner(
        ui,
        tier,
        "This is a map of the plasmid. Coloured arcs are features (genes, promoters); ticks \
         outside the ring are enzyme cut sites. Drag the circle to rotate; click anything to \
         inspect it.",
    );
    let avail = ui.available_size();
    if rec.circular {
        let circ_h = (avail.y * 0.6).min(avail.x).max(220.0);
        draw_circular(ui, rec, state, circ_h);
        ui.separator();
    }
    draw_linear(ui, rec, state);
}

fn draw_circular(ui: &mut egui::Ui, rec: &SequenceRecord, state: &mut PlasmidState, height: f32) {
    let (response, painter) =
        ui.allocate_painter(vec2(ui.available_width(), height), Sense::click_and_drag());
    let rect = response.rect;
    let center = rect.center();
    let radius = (rect.width().min(rect.height()) * 0.5 - 70.0).max(40.0);
    let len = rec.len().max(1) as f32;

    // rotate on drag (before computing angles, so the closure can copy rotation)
    if response.dragged() {
        let d = response.drag_delta();
        state.rotation += (d.x - d.y) * 0.004;
    }
    let rotation = state.rotation;
    let angle_of = move |pos: f32| -> f32 { rotation - FRAC_PI_2 + TAU * (pos / len) };
    let point = move |r: f32, ang: f32| -> Pos2 { center + vec2(r * ang.cos(), r * ang.sin()) };

    // backbone
    painter.circle_stroke(center, radius, Stroke::new(2.0, palette::RULER));

    // position ticks
    let step = theme::nice_step(rec.len() as f64 / 12.0) as f32;
    let mut p = 0.0;
    while p < len {
        let ang = angle_of(p);
        painter.line_segment(
            [point(radius - 5.0, ang), point(radius + 5.0, ang)],
            Stroke::new(1.0, palette::RULER.gamma_multiply(0.7)),
        );
        painter.text(
            point(radius + 16.0, ang),
            Align2::CENTER_CENTER,
            theme::format_pos_compact(p as u64),
            FontId::monospace(9.0),
            palette::RULER_TEXT,
        );
        p += step;
    }

    // feature arcs
    let click = response.clicked().then(|| response.interact_pointer_pos()).flatten();
    let mut new_sel: Option<Selected> = None;
    let band = radius - 16.0;
    for (i, f) in rec.features.iter().enumerate() {
        let a0 = angle_of(f.range.start as f32);
        let a1 = angle_of(f.range.end as f32);
        let color = feature_color(&f.kind);
        let selected = state.selected == Selected::Feature(i);
        let pts = arc_points(center, band, a0, a1);
        painter.add(egui::Shape::line(
            pts.clone(),
            Stroke::new(if selected { 12.0 } else { 8.0 }, color),
        ));
        // label at midpoint
        let mid = angle_of((f.range.start + f.range.end) as f32 * 0.5);
        painter.text(
            point(band - 16.0, mid),
            Align2::CENTER_CENTER,
            f.label(),
            FontId::proportional(10.0),
            color,
        );
        if let Some(cp) = click
            && pts.iter().any(|pp| pp.distance(cp) < 9.0) {
                new_sel = Some(Selected::Feature(i));
            }
    }

    // restriction sites (radial ticks outside)
    let mut last_label_ang = f32::NEG_INFINITY;
    for (i, s) in state.visible_sites() {
        let ang = angle_of(s.start as f32);
        let selected = state.selected == Selected::Site(i);
        let col = if state.unique_names.contains(s.enzyme) {
            palette::ACCENT
        } else {
            palette::RULER_TEXT
        };
        painter.line_segment(
            [point(radius + 6.0, ang), point(radius + 26.0, ang)],
            Stroke::new(if selected { 2.5 } else { 1.0 }, col),
        );
        // De-clutter: label only well-separated (or selected) sites so labels
        // don't pile up when many enzymes cut.
        if selected || (ang - last_label_ang).abs() > 0.16 {
            painter.text(
                point(radius + 42.0, ang),
                Align2::CENTER_CENTER,
                format!("{} ({})", s.enzyme, s.start + 1),
                FontId::proportional(10.0),
                col,
            );
            last_label_ang = ang;
        }
        if let Some(cp) = click
            && point(radius + 16.0, ang).distance(cp) < 12.0 {
                new_sel = Some(Selected::Site(i));
            }
    }

    // centre label
    painter.text(
        center,
        Align2::CENTER_CENTER,
        format!("{}\n{} bp\ncircular", rec.name, rec.len()),
        FontId::proportional(13.0),
        palette::RULER_TEXT,
    );

    if let Some(sel) = new_sel {
        state.selected = sel;
    }
}

fn draw_linear(ui: &mut egui::Ui, rec: &SequenceRecord, state: &mut PlasmidState) {
    let (response, painter) =
        ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let rect = response.rect;
    painter.rect_filled(rect, 0.0, ui.visuals().extreme_bg_color);
    let len = rec.len() as f64;

    // pan / zoom
    if response.dragged() {
        state.lin_offset -= response.drag_delta().x as f64 * state.lin_bp_per_px;
    }
    let scroll = ui.input(|i| i.smooth_scroll_delta.y);
    if scroll != 0.0
        && let Some(pp) = response.hover_pos() {
            let anchor = state.lin_offset + (pp.x - rect.left()) as f64 * state.lin_bp_per_px;
            let factor = if scroll > 0.0 { 0.85 } else { 1.0 / 0.85 };
            state.lin_bp_per_px = (state.lin_bp_per_px * factor).clamp(0.05, len / 50.0 + 1.0);
            state.lin_offset = anchor - (pp.x - rect.left()) as f64 * state.lin_bp_per_px;
        }
    let view_bp = rect.width() as f64 * state.lin_bp_per_px;
    state.lin_offset = state.lin_offset.clamp(0.0, (len - view_bp).max(0.0));

    let x_of = |base: f64| -> f32 {
        rect.left() + ((base - state.lin_offset) / state.lin_bp_per_px) as f32
    };

    // ruler
    let ruler_y = rect.top() + 16.0;
    painter.line_segment(
        [pos2(rect.left(), ruler_y), pos2(rect.right(), ruler_y)],
        Stroke::new(1.0, palette::RULER),
    );
    let step = theme::nice_step(view_bp / 8.0);
    let mut tick = (state.lin_offset / step).floor() * step;
    while tick < state.lin_offset + view_bp {
        if tick >= 0.0 {
            let x = x_of(tick);
            painter.line_segment([pos2(x, ruler_y - 4.0), pos2(x, ruler_y + 4.0)], Stroke::new(1.0, palette::RULER));
            painter.text(pos2(x, rect.top()), Align2::CENTER_TOP, theme::format_pos_compact(tick as u64), FontId::monospace(9.0), palette::RULER_TEXT);
        }
        tick += step;
    }

    let click = response.clicked().then(|| response.interact_pointer_pos()).flatten();
    let mut new_sel: Option<Selected> = None;

    // features as arrows
    let feat_y = ruler_y + 22.0;
    for (i, f) in rec.features.iter().enumerate() {
        let x0 = x_of(f.range.start as f64).max(rect.left());
        let x1 = x_of(f.range.end as f64).min(rect.right());
        if x1 < x0 {
            continue;
        }
        let color = feature_color(&f.kind);
        let selected = state.selected == Selected::Feature(i);
        let r = Rect::from_min_max(pos2(x0, feat_y - 7.0), pos2(x1.max(x0 + 2.0), feat_y + 7.0));
        painter.rect_filled(r, 2.0, color.gamma_multiply(if selected { 1.0 } else { 0.85 }));
        if r.width() > 36.0 {
            painter.text(r.center(), Align2::CENTER_CENTER, f.label(), FontId::proportional(10.0), Color32::WHITE);
        }
        if let Some(cp) = click
            && r.expand(2.0).contains(cp) {
                new_sel = Some(Selected::Feature(i));
            }
    }

    // ORFs (thin bars below features)
    if state.show_orfs {
        let orf_y = feat_y + 18.0;
        for (i, o) in state.orfs.iter().enumerate() {
            let x0 = x_of(o.start as f64).max(rect.left());
            let x1 = x_of(o.end as f64).min(rect.right());
            if x1 < x0 {
                continue;
            }
            let selected = state.selected == Selected::Orf(i);
            let yy = if o.strand == gx_core::Strand::Reverse { orf_y + 8.0 } else { orf_y };
            let r = Rect::from_min_max(pos2(x0, yy - 3.0), pos2(x1.max(x0 + 2.0), yy + 3.0));
            painter.rect_filled(r, 1.0, Color32::from_rgb(120, 150, 120).gamma_multiply(if selected { 1.0 } else { 0.7 }));
            if let Some(cp) = click
                && r.expand(3.0).contains(cp) {
                    new_sel = Some(Selected::Orf(i));
                }
        }
    }

    // restriction site ticks
    let site_y0 = rect.bottom() - 34.0;
    for (i, s) in state.visible_sites() {
        let x = x_of(s.start as f64);
        if x < rect.left() || x > rect.right() {
            continue;
        }
        let unique = state.unique_names.contains(s.enzyme);
        let col = if unique { palette::ACCENT } else { palette::RULER_TEXT };
        let selected = state.selected == Selected::Site(i);
        painter.line_segment([pos2(x, site_y0), pos2(x, site_y0 + 14.0)], Stroke::new(if selected { 2.5 } else { 1.0 }, col));
        painter.text(pos2(x, site_y0 - 1.0), Align2::CENTER_BOTTOM, s.enzyme, FontId::proportional(9.0), col);
        if let Some(cp) = click
            && (cp.x - x).abs() < 5.0 && cp.y > site_y0 - 14.0 {
                new_sel = Some(Selected::Site(i));
            }
    }

    // bases when zoomed in
    if state.lin_bp_per_px <= 0.125 {
        let seq = &rec.seq;
        let first = state.lin_offset.floor().max(0.0) as usize;
        let last = ((state.lin_offset + view_bp).ceil() as usize).min(seq.len());
        for (pos, &b) in seq.iter().enumerate().take(last).skip(first) {
            let x = x_of(pos as f64 + 0.5);
            painter.text(
                pos2(x, rect.bottom() - 16.0),
                Align2::CENTER_TOP,
                (b as char).to_string(),
                FontId::monospace(12.0),
                theme::base_color(b),
            );
        }
    }

    if let Some(sel) = new_sel {
        state.selected = sel;
    }
}

// ============================ side panels ===================================

/// Returns a newly chosen record index if the user switched records.
pub fn sidebar(
    ui: &mut egui::Ui,
    records: &[SequenceRecord],
    state: &mut PlasmidState,
    tier: Tier,
) -> Option<usize> {
    let rec = &records[state.record_index];
    let mut switch = None;

    ui.add_space(4.0);
    ui.heading("Plasmid designer");
    ui.separator();

    if records.len() > 1 {
        egui::ComboBox::from_label("Record")
            .selected_text(rec.name.clone())
            .show_ui(ui, |ui| {
                for (i, r) in records.iter().enumerate() {
                    if ui.selectable_label(i == state.record_index, &r.name).clicked() {
                        switch = Some(i);
                    }
                }
            });
        ui.separator();
    }

    let gc = gx_plasmid::gc_content(&rec.seq) * 100.0;
    ui.label(format!("{} · {} bp", rec.name, rec.len()));
    ui.label(
        egui::RichText::new(format!(
            "{} · GC {:.1}% · {} features",
            if rec.circular { "circular" } else { "linear" },
            gc,
            rec.features.len()
        ))
        .size(11.0)
        .color(palette::RULER_TEXT),
    );
    if let Some(d) = &rec.description {
        ui.label(egui::RichText::new(d).size(11.0).color(palette::RULER_TEXT));
    }

    ui.separator();
    ui.checkbox(&mut state.show_unique_only, "Unique cutters only");
    ui.checkbox(&mut state.show_orfs, "Show ORFs");

    ui.add_space(4.0);
    ux::explain(
        ui,
        tier,
        "Unique cutters cut the plasmid exactly once — the convenient places to insert new DNA.",
    );
    let uniq = unique_cutters(&state.sites);
    ui.label(egui::RichText::new(format!("Unique cutters ({})", uniq.len())).strong());
    egui::ScrollArea::vertical()
        .max_height(150.0)
        .id_salt("uniq")
        .show(ui, |ui| {
            for (i, s) in state.sites.iter().enumerate() {
                if state.unique_names.contains(s.enzyme)
                    && ui
                        .button(format!("{} @ {}", s.enzyme, s.start + 1))
                        .clicked()
                {
                    state.selected = Selected::Site(i);
                }
            }
        });

    ui.add_space(4.0);
    ui.label(egui::RichText::new(format!("ORFs ≥ {ORF_MIN_AA} aa ({})", state.orfs.len())).strong());
    egui::ScrollArea::vertical()
        .max_height(150.0)
        .id_salt("orfs")
        .show(ui, |ui| {
            for (i, o) in state.orfs.iter().enumerate() {
                if ui
                    .button(format!(
                        "{} {}..{} ({} aa)",
                        o.strand.symbol(),
                        o.start + 1,
                        o.end,
                        o.aa_len()
                    ))
                    .clicked()
                {
                    state.selected = Selected::Orf(i);
                }
            }
        });

    ui.add_space(4.0);
    egui::CollapsingHeader::new("✂ Digest / cloning")
        .default_open(false)
        .show(ui, |ui| cloning_ui(ui, rec, state, tier));

    switch
}

/// Restriction-digest preview: fragment sizes (as a mini gel) and sticky ends,
/// using whichever sites are currently shown (respecting "unique cutters only").
fn cloning_ui(ui: &mut egui::Ui, rec: &SequenceRecord, state: &PlasmidState, tier: Tier) {
    ux::explain_beginner(
        ui,
        tier,
        "A digest cuts the DNA with the enzymes shown, yielding fragments. This previews their \
         sizes (like a gel) and the sticky ends left behind — the basis of cloning.",
    );
    let sites: Vec<gx_plasmid::RestrictionSite> =
        state.visible_sites().iter().map(|(_, s)| (*s).clone()).collect();
    let mut frags = gx_plasmid::digest(rec.len(), rec.circular, &sites);
    frags.sort_by_key(|f| std::cmp::Reverse(f.len));

    ui.label(egui::RichText::new(format!("{} fragment(s)", frags.len())).strong());
    gel(ui, &frags);
    egui::ScrollArea::vertical()
        .max_height(150.0)
        .id_salt("digest")
        .show(ui, |ui| {
            for f in &frags {
                ui.label(
                    egui::RichText::new(format!(
                        "{} bp   {}  —  {}",
                        f.len,
                        end_label(&f.left),
                        end_label(&f.right)
                    ))
                    .size(11.0)
                    .color(palette::RULER_TEXT),
                );
            }
        });
}

fn end_label(e: &gx_plasmid::FragmentEnd) -> String {
    match e.enzyme {
        Some(name) => format!("{name} {}", e.overhang.label()),
        None => "free".to_string(),
    }
}

/// A tiny agarose-gel sketch: fragment bands placed by log(size) (larger = higher).
fn gel(ui: &mut egui::Ui, frags: &[gx_plasmid::Fragment]) {
    if frags.is_empty() {
        return;
    }
    let (rect, _) = ui.allocate_exact_size(vec2(ui.available_width().min(210.0), 70.0), Sense::hover());
    let painter = ui.painter();
    painter.rect_filled(rect, 3.0, Color32::from_gray(28));
    let max_len = frags.iter().map(|f| f.len).max().unwrap_or(1).max(1) as f64;
    let lo = (50.0_f64).ln();
    let hi = (max_len.max(60.0)).ln();
    for f in frags {
        let t = (((f.len.max(1)) as f64).ln() - lo) / (hi - lo + 1e-9);
        let y = rect.bottom() - 6.0 - t.clamp(0.0, 1.0) as f32 * (rect.height() - 12.0);
        painter.line_segment(
            [pos2(rect.left() + 12.0, y), pos2(rect.right() - 12.0, y)],
            Stroke::new(2.5, palette::EXON.gamma_multiply(0.9)),
        );
    }
}

pub fn detail(ui: &mut egui::Ui, rec: &SequenceRecord, state: &PlasmidState, tier: Tier) {
    ui.add_space(4.0);
    ui.heading("Details");
    ui.separator();
    ux::explain_beginner(
        ui,
        tier,
        "Selecting a feature shows its DNA, GC content, and the protein it codes for. Selecting \
         an enzyme site shows exactly where it cuts.",
    );

    match &state.selected {
        Selected::None => {
            ui.label(
                egui::RichText::new("Click a feature, restriction site, or ORF.")
                    .color(palette::RULER_TEXT),
            );
        }
        Selected::Feature(i) => {
            if let Some(f) = rec.features.get(*i) {
                let s = (f.range.start as usize).min(rec.seq.len());
                let e = (f.range.end as usize).min(rec.seq.len()).max(s);
                let sub = &rec.seq[s..e];
                detail_grid(ui, &[
                    ("Feature", f.label().to_string()),
                    ("Type", format!("{:?}", f.kind)),
                    ("Location", format!("{}..{}", f.range.start + 1, f.range.end)),
                    ("Length", format!("{} bp", sub.len())),
                    ("GC", format!("{:.1}%", gx_plasmid::gc_content(sub) * 100.0)),
                ]);
                sequence_tools(ui, sub);
            }
        }
        Selected::Site(i) => {
            if let Some(s) = state.sites.get(*i) {
                let unique = state.unique_names.contains(s.enzyme);
                detail_grid(ui, &[
                    ("Enzyme", s.enzyme.to_string()),
                    ("Recognition", s.site.to_string()),
                    ("Site start", (s.start + 1).to_string()),
                    ("Cut position", (s.cut + 1).to_string()),
                    ("Cutter", if unique { "unique".into() } else { "multiple sites".into() }),
                ]);
            }
        }
        Selected::Orf(i) => {
            if let Some(o) = state.orfs.get(*i) {
                detail_grid(ui, &[
                    ("ORF", format!("{}..{}", o.start + 1, o.end)),
                    ("Strand", o.strand.symbol().to_string()),
                    ("Frame", o.frame.to_string()),
                    ("Length", format!("{} aa", o.aa_len())),
                ]);
                ui.label("Protein:");
                egui::ScrollArea::vertical().max_height(160.0).show(ui, |ui| {
                    ui.add(egui::Label::new(egui::RichText::new(&o.protein).monospace().size(11.0)).wrap());
                });
            }
        }
    }

    ui.add_space(8.0);
    ux::disclaimer(ui);
}

fn sequence_tools(ui: &mut egui::Ui, sub: &[u8]) {
    if sub.is_empty() {
        return;
    }
    let stats = gx_plasmid::analyze(sub);
    // The basic formula is only meaningful for ~14–50 nt; show "—" otherwise.
    let tm_basic = if sub.len() >= 14 && stats.tm_basic > 0.0 {
        format!("{:.1} °C", stats.tm_basic)
    } else {
        "—".to_string()
    };
    let tm_nn = if stats.tm_nn.is_finite() {
        format!("{:.1} °C", stats.tm_nn)
    } else {
        "—".to_string()
    };
    ui.add_space(4.0);
    detail_grid(ui, &[("Tm (basic)", tm_basic), ("Tm (NN)", tm_nn)]);
    ui.collapsing("Translation (frame +1)", |ui| {
        let protein = gx_plasmid::translate(sub);
        egui::ScrollArea::vertical().max_height(140.0).show(ui, |ui| {
            ui.add(egui::Label::new(egui::RichText::new(protein).monospace().size(11.0)).wrap());
        });
    });
}

fn detail_grid(ui: &mut egui::Ui, rows: &[(&str, String)]) {
    egui::Grid::new("plasmid_detail").num_columns(2).striped(true).show(ui, |ui| {
        for (k, v) in rows {
            ui.label(*k);
            ui.label(v);
            ui.end_row();
        }
    });
}

// ============================ small helpers =================================

fn arc_points(center: Pos2, radius: f32, a0: f32, a1: f32) -> Vec<Pos2> {
    let span = (a1 - a0).abs().max(0.02);
    let steps = (span / 0.08).ceil().max(2.0) as usize;
    (0..=steps)
        .map(|k| {
            let t = a0 + (a1 - a0) * (k as f32 / steps as f32);
            center + vec2(radius * t.cos(), radius * t.sin())
        })
        .collect()
}

