//! Pharmacogenomics workspace: how the loaded genome's variants relate to drug
//! response, via CPIC-derived star-allele / diplotype calling (see `gx_pgx`).
//! Computed locally; informational only and emphatically not medical advice.

use eframe::egui;
use egui::RichText;

use gx_pgx::{Confidence, GeneResult};

use crate::theme::palette;
use crate::ux::{self, Tier};

/// Per-document pharmacogenomics report, recomputed on each import.
#[derive(Default)]
pub struct PgxState {
    pub results: Vec<GeneResult>,
    pub selected: Option<usize>,
    /// Whether the report has been computed for the current document.
    pub done: bool,
}

pub fn sidebar(ui: &mut egui::Ui, state: &mut PgxState, has_doc: bool, tier: Tier) {
    ui.add_space(4.0);
    ui.heading("Pharmacogenomics");
    ux::explain_beginner(
        ui,
        tier,
        "Some genes change how your body processes certain medicines. This reads the \
         relevant variants in your file and shows what pharmacogenomics guidelines say \
         about them. It is background information, never a reason to change a medication.",
    );
    medical_warning(ui);
    ui.separator();

    if !has_doc {
        ui.label(
            RichText::new("Load a genome (File ▸ Open) to build a report.")
                .size(11.0)
                .color(palette::RULER_TEXT),
        );
        return;
    }
    if state.results.is_empty() {
        ui.label(
            RichText::new(
                "No pharmacogene markers were found in this file. Consumer chips cover \
                 only some of them.",
            )
            .size(11.0)
            .color(palette::RULER_TEXT),
        );
        return;
    }

    ui.label(RichText::new(format!("Genes found ({})", state.results.len())).strong());
    for (i, r) in state.results.iter().enumerate() {
        let selected = state.selected == Some(i);
        let label = format!("{}  {}", r.gene, r.diplotype);
        if ui.selectable_label(selected, label).clicked() {
            state.selected = Some(i);
        }
    }
}

pub fn central(ui: &mut egui::Ui, state: &PgxState, has_doc: bool, tier: Tier) {
    ui.add_space(8.0);
    ui.heading("Pharmacogenomics report");
    ui.separator();

    if !has_doc {
        ui.label(
            RichText::new("Load a genome to see how your variants relate to drug response.")
                .color(palette::RULER_TEXT),
        );
        return;
    }
    if state.results.is_empty() {
        ui.label(
            RichText::new(
                "No pharmacogene markers were genotyped in this file. This does not mean \
                 anything is normal or abnormal; the variants simply were not measured.",
            )
            .color(palette::RULER_TEXT),
        );
        return;
    }

    let Some(r) = state.selected.and_then(|i| state.results.get(i)) else {
        ui.label(
            RichText::new("Select a gene on the left to see its details.")
                .color(palette::RULER_TEXT),
        );
        ui.add_space(8.0);
        ui.label(RichText::new("Summary").strong());
        for r in &state.results {
            ui.label(format!("{}: {} ({})", r.gene, r.phenotype, r.diplotype));
        }
        return;
    };

    ui.horizontal(|ui| {
        ui.heading(&r.gene);
        ui.label(RichText::new(&r.diplotype).color(palette::RULER_TEXT));
    });
    ui.horizontal(|ui| {
        ui.label(RichText::new(&r.phenotype).size(18.0).strong());
        confidence_badge(ui, r.confidence);
    });

    if let Some(limit) = &r.limitation {
        ui.add_space(4.0);
        egui::Frame::new()
            .fill(palette::VUS.gamma_multiply(0.14))
            .stroke(egui::Stroke::new(1.0, palette::VUS))
            .corner_radius(egui::CornerRadius::same(6))
            .inner_margin(egui::Margin::same(8))
            .show(ui, |ui| {
                ui.label(
                    RichText::new(format!("⚠ Limitation: {limit}"))
                        .size(11.0)
                        .color(palette::VUS),
                );
            });
    }

    ui.add_space(10.0);
    ui.label(RichText::new("Medication notes (informational only)").strong());
    if r.drugs.is_empty() {
        ui.label(
            RichText::new("No specific drug notes for this result.")
                .size(11.0)
                .color(palette::RULER_TEXT),
        );
    }
    for d in &r.drugs {
        ui.add_space(4.0);
        ui.group(|ui| {
            ui.set_width(ui.available_width());
            ui.horizontal(|ui| {
                ui.label(RichText::new(&d.drug).strong());
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(RichText::new(&d.source).size(10.0).color(palette::RULER_TEXT));
                });
            });
            ui.label(RichText::new(&d.guidance).size(12.0));
        });
    }

    if tier.at_least(Tier::Intermediate) && !r.markers_used.is_empty() {
        ui.add_space(6.0);
        ui.label(
            RichText::new(format!("Based on: {}", r.markers_used.join(", ")))
                .size(10.0)
                .color(palette::RULER_TEXT),
        );
    }

    ui.add_space(8.0);
    ux::disclaimer(ui);
}

pub fn detail(ui: &mut egui::Ui, tier: Tier) {
    ui.add_space(4.0);
    ui.heading("About this report");
    ui.separator();
    ux::explain_beginner(
        ui,
        tier,
        "A 'metabolizer' status describes how fast a gene's enzyme works on certain \
         drugs. 'Poor' or 'rapid' is about processing speed, not about being unhealthy. \
         The same status can mean 'consider a different drug' for one medicine and \
         'no change' for another.",
    );
    ui.label(
        RichText::new(
            "Guidance is summarised from the Clinical Pharmacogenetics Implementation \
             Consortium (CPIC) and PharmGKB. Allele orientation is checked against dbSNP.",
        )
        .size(11.0)
        .color(palette::RULER_TEXT),
    );
    ui.add_space(6.0);
    medical_warning(ui);
    ui.add_space(4.0);
    ui.label(
        RichText::new(
            "Consumer-array data is not clinical grade and misses structural variants \
             (gene deletions/duplications) that matter for some genes, especially CYP2D6. \
             A 'normal' result never rules out untested variants.",
        )
        .size(11.0)
        .color(palette::RULER_TEXT),
    );
}

/// The prominent, repeated safety banner this workspace leans on.
fn medical_warning(ui: &mut egui::Ui) {
    egui::Frame::new()
        .fill(palette::VUS.gamma_multiply(0.12))
        .stroke(egui::Stroke::new(1.0, palette::VUS))
        .corner_radius(egui::CornerRadius::same(6))
        .inner_margin(egui::Margin::same(8))
        .show(ui, |ui| {
            ui.label(
                RichText::new(
                    "Not medical advice. Do not start, stop, or change any medication or \
                     dose based on this. Talk to a clinician or pharmacist, who can order \
                     validated testing.",
                )
                .size(11.0)
                .strong()
                .color(palette::VUS),
            );
        });
}

fn confidence_badge(ui: &mut egui::Ui, c: Confidence) {
    let (txt, col) = match c {
        Confidence::High => ("Well covered", palette::BENIGN),
        Confidence::Moderate => ("Partial coverage", palette::VUS),
        Confidence::Low => ("Limited", palette::RULER_TEXT),
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
    // End-to-end: the bundled demo genome's pharmacogene markers must produce the
    // expected report. Guards the import pipeline, the demo file, and gx-pgx.
    #[test]
    fn demo_genome_pgx_report() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../assets/sample_genome_23andme.txt"
        );
        let store = match gx_io::import_path(std::path::Path::new(path)).expect("load demo genome") {
            gx_io::Imported::Variants(v) => v.store,
            gx_io::Imported::Sequences(_) => panic!("demo genome should be variants"),
        };
        let report = gx_pgx::report(&store);
        let by_gene = |g: &str| report.iter().find(|r| r.gene == g);

        let cyp2c19 = by_gene("CYP2C19").expect("CYP2C19 present");
        assert_eq!(cyp2c19.phenotype, "Intermediate metabolizer");
        assert!(cyp2c19.drugs.iter().any(|d| d.drug == "Clopidogrel"));

        let slco = by_gene("SLCO1B1").expect("SLCO1B1 present");
        assert_eq!(slco.diplotype, "*1/*5");

        // CYP2D6 must always carry its structural-variant limitation.
        let cyp2d6 = by_gene("CYP2D6").expect("CYP2D6 present");
        assert!(cyp2d6.limitation.is_some());
    }
}
