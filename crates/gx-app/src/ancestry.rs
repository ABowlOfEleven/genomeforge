//! Ancestry & lineage workspace: predicted maternal (mtDNA) and paternal (Y-DNA)
//! haplogroups from the loaded genome, computed locally from the sample's SNP
//! calls (no network). Population-ancestry / PCA estimation also lives here.

use eframe::egui;
use egui::RichText;

use gx_haplo::{Confidence, HaploCall};

use crate::theme::palette;
use crate::ux::{self, Tier};

/// Per-document ancestry results, recomputed on each import.
#[derive(Default)]
pub struct AncestryState {
    pub maternal: Option<HaploCall>,
    pub paternal: Option<HaploCall>,
    /// Whether haplogroup classification has already run for the current document.
    pub haplo_done: bool,
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
            RichText::new("Load a genome (File ▸ Open) to predict your haplogroups.")
                .color(palette::RULER_TEXT),
        );
        return;
    }

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
    }
}
