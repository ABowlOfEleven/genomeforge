//! Ancestry & lineage workspace: predicted maternal (mtDNA) and paternal (Y-DNA)
//! haplogroups from the loaded genome, computed locally from the sample's SNP
//! calls (no network). Population-ancestry / PCA estimation also lives here.

use eframe::egui;
use egui::{Color32, Pos2, RichText, Sense, Stroke, pos2, vec2};

use gx_ancestry::{AncestryEstimate, SUPERPOPS, SUPERPOP_NAMES};
use gx_haplo::{Confidence, HaploCall};

use crate::theme::palette;
use crate::ux::{self, Tier};

/// Per-super-population colours, in SUPERPOPS order [AFR, AMR, EAS, EUR, SAS].
const POP_COLORS: [Color32; 5] = [
    Color32::from_rgb(0xE8, 0x8A, 0x3C), // African, orange
    Color32::from_rgb(0xD1, 0x49, 0x49), // Admixed American, red
    Color32::from_rgb(0x3F, 0xA9, 0x5B), // East Asian, green
    Color32::from_rgb(0x3C, 0x7E, 0xD1), // European, blue
    Color32::from_rgb(0x9B, 0x5D, 0xE5), // South Asian, purple
];

/// Per-document ancestry results, recomputed on each import.
#[derive(Default)]
pub struct AncestryState {
    pub maternal: Option<HaploCall>,
    pub paternal: Option<HaploCall>,
    pub composition: Option<AncestryEstimate>,
    /// Whether ancestry classification has already run for the current document.
    pub done: bool,
}

pub fn sidebar(ui: &mut egui::Ui, tier: Tier) {
    ui.add_space(4.0);
    ui.heading("Ancestry & lineage");
    ux::explain_beginner(
        ui,
        tier,
        "Your mitochondrial DNA traces a single maternal line (mother, her mother, \
         and so on); the Y chromosome traces a single paternal line in people who \
         have one. The deep branch each falls on is its haplogroup.",
    );
    ui.separator();
    ui.label(
        RichText::new(
            "Haplogroups are predicted entirely on your computer from the SNPs in \
             your file. Consumer chips read only a fraction of mtDNA / Y positions, \
             so the result is a major branch, not a fine subclade.",
        )
        .size(11.0)
        .color(palette::RULER_TEXT),
    );
    ui.add_space(6.0);
    ux::disclaimer(ui);
}

pub fn central(ui: &mut egui::Ui, state: &AncestryState, has_doc: bool, tier: Tier) {
    ui.add_space(8.0);
    ui.heading("Your deep ancestry");
    ui.separator();

    if !has_doc {
        ui.label(
            RichText::new("Load a genome (File ▸ Open) to estimate your ancestry.")
                .color(palette::RULER_TEXT),
        );
        return;
    }

    composition_section(ui, state, tier);

    ui.add_space(12.0);
    ui.label(RichText::new("Deep lineages").strong());
    ui.add_space(4.0);
    lineage_card(
        ui,
        "Maternal line (mtDNA)",
        &state.maternal,
        "No usable mitochondrial markers were genotyped in this file.",
        tier,
    );
    ui.add_space(8.0);
    lineage_card(
        ui,
        "Paternal line (Y-DNA)",
        &state.paternal,
        "No Y-chromosome data found. This is expected for samples from people \
         without a Y chromosome (typically XX).",
        tier,
    );

    ui.add_space(10.0);
    ui.group(|ui| {
        ui.label(
            RichText::new(
                "A haplogroup is a branch on the human family tree, shared by everyone \
                 descended from one ancestor along a single line. It says where that one \
                 line came from deep in the past; it is not your whole ancestry, and it \
                 says nothing about health.",
            )
            .size(11.0)
            .color(palette::RULER_TEXT),
        );
    });
}

pub fn detail(ui: &mut egui::Ui, state: &AncestryState, tier: Tier) {
    ui.add_space(4.0);
    ui.heading("How this was determined");
    ui.separator();

    let mut any = false;
    for call in [&state.maternal, &state.paternal].into_iter().flatten() {
        any = true;
        ui.label(RichText::new(call.lineage.label()).strong());
        ui.label(
            RichText::new(format!(
                "Assigned {} from {} matched marker(s) ({} relevant markers in your data).",
                call.haplogroup,
                call.supporting.len(),
                call.tested
            ))
            .size(11.0),
        );
        if !call.supporting.is_empty() && tier.at_least(Tier::Intermediate) {
            ui.label(
                RichText::new(format!("Markers: {}", call.supporting.join(", ")))
                    .size(10.0)
                    .color(palette::RULER_TEXT),
            );
        }
        ui.add_space(6.0);
    }

    if !any {
        ui.label(
            RichText::new("Load a genome to see which markers placed each lineage.")
                .color(palette::RULER_TEXT),
        );
        return;
    }

    ui.separator();
    ux::explain_beginner(
        ui,
        tier,
        "Each lineage is placed by reading 'defining mutations', positions where a \
         branch carries a specific letter. The more of a branch's mutations your chip \
         happens to cover, the deeper and more confident the call.",
    );
    ui.label(
        RichText::new(
            "Reference data: PhyloTree mtDNA Build 17 (maternal) and the ISOGG Y-DNA \
             tree (paternal).",
        )
        .size(10.0)
        .color(palette::RULER_TEXT),
    );
    ui.add_space(6.0);
    ux::disclaimer(ui);
}

fn composition_section(ui: &mut egui::Ui, state: &AncestryState, tier: Tier) {
    ui.label(RichText::new("Ancestry composition").strong());
    ux::explain_beginner(
        ui,
        tier,
        "This compares the ancestry-informative markers in your file against five \
         broad reference groups and estimates how much your genome resembles each. \
         It is genetic similarity at a continental scale, not an identity test.",
    );

    let Some(est) = &state.composition else {
        ui.label(
            RichText::new(
                "Not enough ancestry-informative markers in this file to estimate \
                 composition.",
            )
            .size(11.0)
            .color(palette::RULER_TEXT),
        );
        return;
    };

    ui.add_space(4.0);
    stacked_bar(ui, est);
    ui.add_space(6.0);

    // Ranked proportions with colour swatches (anything at least 1%).
    for (name, frac) in est.ranked() {
        if frac < 0.01 {
            continue;
        }
        let idx = SUPERPOP_NAMES.iter().position(|n| *n == name).unwrap_or(0);
        ui.horizontal(|ui| {
            let (sw, _) = ui.allocate_exact_size(vec2(12.0, 12.0), Sense::hover());
            ui.painter().rect_filled(sw, 2.0, POP_COLORS[idx]);
            ui.label(format!("{name}: {:.0}%", frac * 100.0));
        });
    }

    ui.add_space(8.0);
    simplex_plot(ui, est);

    ui.label(
        RichText::new(format!(
            "{} of {} markers used.",
            est.markers_used, est.markers_total
        ))
        .size(10.0)
        .color(palette::RULER_TEXT),
    );
    ui.label(RichText::new(&est.note).size(11.0).color(palette::RULER_TEXT));
}

/// A horizontal stacked bar of the five proportions, in SUPERPOPS order.
#[allow(clippy::needless_range_loop)] // parallel indexing of proportions + colours
fn stacked_bar(ui: &mut egui::Ui, est: &AncestryEstimate) {
    let w = ui.available_width().min(440.0);
    let (rect, _) = ui.allocate_exact_size(vec2(w, 20.0), Sense::hover());
    let mut x = rect.left();
    for k in 0..5 {
        let seg_w = rect.width() * est.proportions[k] as f32;
        let seg = egui::Rect::from_min_size(pos2(x, rect.top()), vec2(seg_w, rect.height()));
        ui.painter().rect_filled(seg, 0.0, POP_COLORS[k]);
        x += seg_w;
    }
}

/// Plot the five reference groups at the vertices of a pentagon and the sample at
/// the proportion-weighted centroid, showing where it sits among the references.
#[allow(clippy::needless_range_loop)] // vertices indexed with modular wraparound
fn simplex_plot(ui: &mut egui::Ui, est: &AncestryEstimate) {
    let size = 210.0;
    let (rect, _) = ui.allocate_exact_size(vec2(size, size), Sense::hover());
    let center = rect.center();
    let radius = size * 0.34;
    let mut verts = [center; 5];
    for (k, v) in verts.iter_mut().enumerate() {
        let ang = -std::f32::consts::FRAC_PI_2 + k as f32 * std::f32::consts::TAU / 5.0;
        *v = center + vec2(ang.cos() * radius, ang.sin() * radius);
    }
    let painter = ui.painter();
    // Faint pentagon outline.
    for k in 0..5 {
        painter.line_segment(
            [verts[k], verts[(k + 1) % 5]],
            Stroke::new(1.0, palette::RULER_TEXT.gamma_multiply(0.4)),
        );
    }
    // Reference anchors with short labels just outside each vertex.
    for k in 0..5 {
        painter.circle_filled(verts[k], 5.0, POP_COLORS[k]);
        let outward = (verts[k] - center).normalized() * 14.0;
        painter.text(
            verts[k] + outward,
            egui::Align2::CENTER_CENTER,
            SUPERPOPS[k],
            egui::FontId::proportional(10.0),
            POP_COLORS[k],
        );
    }
    // Sample point: proportion-weighted centroid of the vertices.
    let mut ix = 0.0;
    let mut iy = 0.0;
    for k in 0..5 {
        ix += verts[k].x * est.proportions[k] as f32;
        iy += verts[k].y * est.proportions[k] as f32;
    }
    let sample = Pos2::new(ix, iy);
    painter.circle_filled(sample, 6.0, palette::ACCENT);
    painter.circle_stroke(sample, 6.0, Stroke::new(1.5, Color32::WHITE));
}

fn lineage_card(
    ui: &mut egui::Ui,
    title: &str,
    call: &Option<HaploCall>,
    empty_msg: &str,
    _tier: Tier,
) {
    egui::Frame::group(ui.style()).show(ui, |ui| {
        ui.set_width(ui.available_width());
        ui.label(RichText::new(title).strong());
        match call {
            None => {
                ui.label(RichText::new(empty_msg).size(12.0).color(palette::RULER_TEXT));
            }
            Some(c) => {
                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new(&c.haplogroup)
                            .size(24.0)
                            .strong()
                            .color(palette::ACCENT),
                    );
                    confidence_badge(ui, c.confidence);
                });
                ui.label(
                    RichText::new(format!(
                        "{} defining marker(s) matched, {} tested",
                        c.supporting.len(),
                        c.tested
                    ))
                    .size(11.0)
                    .color(palette::RULER_TEXT),
                );
                ui.label(RichText::new(&c.note).size(11.0).color(palette::RULER_TEXT));
            }
        }
    });
}

fn confidence_badge(ui: &mut egui::Ui, c: Confidence) {
    let (txt, col) = match c {
        Confidence::High => ("High confidence", palette::BENIGN),
        Confidence::Moderate => ("Moderate", palette::VUS),
        Confidence::Low => ("Low / best guess", palette::RULER_TEXT),
    };
    egui::Frame::new()
        .fill(col.gamma_multiply(0.18))
        .corner_radius(egui::CornerRadius::same(4))
        .inner_margin(egui::Margin::symmetric(6, 2))
        .show(ui, |ui| {
            ui.label(RichText::new(txt).size(10.0).strong().color(col));
        });
}

#[cfg(test)]
mod tests {
    // End-to-end: the bundled demo genome must parse and classify to the
    // haplogroups its MT/Y markers encode. Guards both the import pipeline and
    // the demo file from regressions.
    #[test]
    fn demo_genome_yields_expected_haplogroups() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../assets/sample_genome_23andme.txt"
        );
        let store = match gx_io::import_path(std::path::Path::new(path)).expect("load demo genome") {
            gx_io::Imported::Variants(v) => v.store,
            gx_io::Imported::Sequences(_) => panic!("demo genome should be variants"),
        };
        let maternal = gx_haplo::classify_maternal(&store).expect("maternal call");
        assert!(
            maternal.haplogroup.starts_with("U5"),
            "expected U5, got {}",
            maternal.haplogroup
        );
        let paternal = gx_haplo::classify_paternal(&store).expect("paternal call");
        assert_eq!(paternal.haplogroup, "R1b", "got {}", paternal.haplogroup);

        // The demo carries a European-typical AIM profile, so the ancestry
        // estimate should be European-dominant.
        let est = gx_ancestry::estimate(&store).expect("ancestry estimate");
        let (top, frac) = est.ranked()[0];
        assert_eq!(top, "European", "top ancestry was {top} ({frac:.2})");
    }
}
