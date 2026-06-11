//! Phenotype & risk workspace: polygenic risk scores (PGS Catalog) applied to
//! the loaded genome, plus Mendelian/ClinVar findings. Heavy on uncertainty
//! framing — a PRS is a probabilistic tendency, not a diagnosis.

use eframe::egui;
use egui::{Align2, Color32, FontId, Rect, Sense, Stroke, pos2, vec2};

use gx_annotate::VariantAnnotation;
use gx_pgs::PrsResult;

use crate::document::VariantDoc;
use crate::theme::{self, palette};
use crate::ux::{self, Tier};

pub struct CuratedPgs {
    pub id: &'static str,
    pub label: &'static str,
}

/// A few well-known scores; the custom box covers everything else.
pub const CURATED: &[CuratedPgs] = &[
    CuratedPgs { id: "PGS000001", label: "Breast cancer (PRS77)" },
    CuratedPgs { id: "PGS000004", label: "Breast cancer (313-variant)" },
    CuratedPgs { id: "PGS000018", label: "Coronary artery disease (metaGRS)" },
];

#[derive(Default)]
pub struct PhenotypeState {
    pub custom_id: String,
    pub results: Vec<PrsResult>,
    pub selected: Option<usize>,
    pub status: String,
}

/// What the sidebar wants the app to do (the app owns the worker + file I/O).
pub enum PhenoRequest {
    Fetch(String),
    LocalFile,
}

pub fn sidebar(ui: &mut egui::Ui, state: &mut PhenotypeState, tier: Tier) -> Option<PhenoRequest> {
    ui.add_space(4.0);
    ui.heading("Phenotype & risk");
    ux::explain_beginner(
        ui,
        tier,
        "A polygenic score adds up thousands of tiny genetic effects to estimate a tendency \
         toward a trait. It is a probability, not a diagnosis.",
    );
    ui.separator();

    let mut req = None;

    ui.label(egui::RichText::new("Example scores").strong());
    for c in CURATED {
        if ui.button(c.label).clicked() {
            req = Some(PhenoRequest::Fetch(c.id.to_string()));
        }
    }
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        ui.add(
            egui::TextEdit::singleline(&mut state.custom_id)
                .hint_text("PGS000123")
                .desired_width(96.0),
        );
        if ui.button("Fetch").clicked() && !state.custom_id.trim().is_empty() {
            req = Some(PhenoRequest::Fetch(state.custom_id.trim().to_string()));
        }
    });
    if ui.button("Load scoring file…").clicked() {
        req = Some(PhenoRequest::LocalFile);
    }
    if tier.at_least(Tier::Intermediate) {
        ui.label(
            egui::RichText::new("Scores come from the PGS Catalog (pgscatalog.org).")
                .size(11.0)
                .color(palette::RULER_TEXT),
        );
    }

    ui.separator();
    ui.label(egui::RichText::new("Computed scores").strong());
    if state.results.is_empty() {
        ui.label(
            egui::RichText::new("None yet. Fetch one above.")
                .size(11.0)
                .color(palette::RULER_TEXT),
        );
    }
    for (i, r) in state.results.iter().enumerate() {
        let selected = state.selected == Some(i);
        if ui
            .selectable_label(selected, format!("{} · {}", r.trait_name, r.id))
            .clicked()
        {
            state.selected = Some(i);
        }
    }
    if !state.status.is_empty() {
        ui.add_space(4.0);
        ui.label(egui::RichText::new(&state.status).size(11.0).color(palette::RULER_TEXT));
    }

    req
}

pub fn central(ui: &mut egui::Ui, state: &PhenotypeState, tier: Tier) {
    let Some(i) = state.selected else {
        ui.centered_and_justified(|ui| {
            ui.label(
                egui::RichText::new(
                    "Fetch a polygenic score in the sidebar to see how your genome compares.",
                )
                .color(palette::RULER_TEXT),
            );
        });
        return;
    };
    let r = &state.results[i];

    ui.add_space(8.0);
    ui.heading(&r.trait_name);
    ui.label(egui::RichText::new(&r.id).color(palette::RULER_TEXT));
    ui.separator();

    let cov = r.coverage() * 100.0;
    ui.label(format!(
        "Used {} of {} score variants ({:.0}% of the score is covered by your data).",
        r.matched, r.total, cov
    ));
    if cov < 50.0 {
        ui.label(
            egui::RichText::new(
                "⚠ Low coverage: your file is missing many of this score's variants, so the \
                 result is unreliable. Whole-genome data covers far more.",
            )
            .color(palette::VUS),
        );
    }
    if r.ambiguous > 0 && tier.at_least(Tier::Intermediate) {
        ui.label(
            egui::RichText::new(format!("{} variants skipped (allele/strand ambiguity).", r.ambiguous))
                .size(11.0)
                .color(palette::RULER_TEXT),
        );
    }

    ui.add_space(10.0);
    match r.percentile {
        Some(p) => {
            ui.label(
                egui::RichText::new(format!("Estimated population percentile: {p:.0}"))
                    .strong()
                    .size(16.0),
            );
            gauge(ui, p);
            ui.label(tier.pick(
                "In plain terms: out of 100 random people, your genetics for this trait rank \
                 roughly here. This is a rough estimate, not destiny.",
                "Your score sits near the {p}th percentile of the modelled distribution \
                 (normal approximation from effect-allele frequencies).",
                "Normal-approx percentile from Σ2·AF·w mean and Σ2·AF·(1−AF)·w² variance.",
            ).replace("{p}", &format!("{p:.0}")));
            if tier.at_least(Tier::Expert)
                && let (Some(z), Some(raw)) = (r.z, Some(r.raw_score)) {
                    ui.label(format!("raw = {raw:.4}, z = {z:.3}"));
                }
        }
        None => {
            ui.label(
                egui::RichText::new(format!("Raw weighted score: {:.4}", r.raw_score))
                    .strong()
                    .size(16.0),
            );
            ux::explain(
                ui,
                tier,
                "This scoring file didn't include population allele frequencies, so we can't place \
                 you on a percentile, only the raw number, which is only meaningful compared to a \
                 reference population.",
            );
        }
    }

    ui.add_space(12.0);
    ui.group(|ui| {
        ui.label(
            egui::RichText::new(
                "Polygenic scores explain only part of most traits; environment, lifestyle, \
                 rare variants, and ancestry all matter. Most published scores were trained in \
                 European-ancestry cohorts and transfer poorly to others. Not medical advice.",
            )
            .size(11.0)
            .color(palette::VUS),
        );
    });
}

pub fn detail(ui: &mut egui::Ui, doc: Option<&VariantDoc>, tier: Tier) {
    ui.add_space(4.0);
    ui.heading("Mendelian findings");
    ui.separator();
    ux::explain_beginner(
        ui,
        tier,
        "These are single variants ClinVar flags as clinically significant. Unlike a polygenic \
         score, one of these can matter on its own; discuss any with a clinician.",
    );

    let Some(doc) = doc else {
        ui.label(egui::RichText::new("Load your genome to see findings.").color(palette::RULER_TEXT));
        return;
    };

    let mut findings: Vec<&VariantAnnotation> = doc
        .annotations
        .values()
        .filter(|a| a.is_clinically_notable())
        .collect();
    findings.sort_by(|a, b| a.gene.cmp(&b.gene));

    if findings.is_empty() {
        ui.label(
            egui::RichText::new(
                "No notable findings among annotated variants yet. Browse the genome browser to \
                 annotate more of your variants.",
            )
            .size(11.0)
            .color(palette::RULER_TEXT),
        );
        return;
    }

    egui::ScrollArea::vertical().show(ui, |ui| {
        for a in findings {
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new(&a.rsid).monospace().strong());
                    if let Some(g) = &a.gene {
                        ui.label(egui::RichText::new(g).color(palette::GENE));
                    }
                });
                if let Some(sig) = &a.clinical_significance {
                    ui.label(egui::RichText::new(sig).color(theme::sig_color(sig)));
                }
                if !a.conditions.is_empty() {
                    ui.label(
                        egui::RichText::new(a.conditions.join("; "))
                            .size(11.0)
                            .color(palette::RULER_TEXT),
                    );
                }
            });
        }
    });
}

fn gauge(ui: &mut egui::Ui, pct: f64) {
    let width = ui.available_width().min(380.0);
    let (rect, _) = ui.allocate_exact_size(vec2(width, 30.0), Sense::hover());
    let painter = ui.painter();
    let track = Rect::from_min_size(pos2(rect.left(), rect.center().y - 5.0), vec2(rect.width(), 10.0));
    painter.rect_filled(track, 5.0, Color32::from_gray(54));
    let x = rect.left() + (pct.clamp(0.0, 100.0) / 100.0) as f32 * rect.width();
    painter.rect_filled(
        Rect::from_min_max(track.left_top(), pos2(x, track.bottom())),
        5.0,
        palette::ACCENT.gamma_multiply(0.55),
    );
    painter.line_segment([pos2(x, rect.top()), pos2(x, rect.bottom())], Stroke::new(2.5, palette::ACCENT));
    for (frac, label) in [(0.0, "0"), (0.5, "50"), (1.0, "100")] {
        let lx = rect.left() + frac * rect.width();
        let anchor = if frac == 0.0 {
            Align2::LEFT_TOP
        } else if frac == 1.0 {
            Align2::RIGHT_TOP
        } else {
            Align2::CENTER_TOP
        };
        painter.text(pos2(lx, rect.bottom()), anchor, label, FontId::proportional(9.0), palette::RULER_TEXT);
    }
}

