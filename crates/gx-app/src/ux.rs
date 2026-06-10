//! Experience tiers and small UI helpers that adapt the interface to the
//! user's expertise.
//!
//! - **Beginner** — written for a high-school / college student: plain language,
//!   lots of inline explanation, advanced numbers hidden behind summaries.
//! - **Intermediate** — builds on Beginner: shows the numbers *with* explanation.
//! - **Expert** — for biologists / genomics professionals: dense, all values,
//!   minimal hand-holding.

use crate::theme::{palette, tokens};
use eframe::egui;
use egui::RichText;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Tier {
    #[default]
    Beginner,
    Intermediate,
    Expert,
}

impl Tier {
    pub fn label(self) -> &'static str {
        match self {
            Tier::Beginner => "Beginner",
            Tier::Intermediate => "Intermediate",
            Tier::Expert => "Expert",
        }
    }

    pub fn blurb(self) -> &'static str {
        match self {
            Tier::Beginner => "Plain-language, guided — for students and the curious.",
            Tier::Intermediate => "More detail, with explanations of the numbers.",
            Tier::Expert => "Dense and complete — for biologists & professionals.",
        }
    }

    pub fn all() -> [Tier; 3] {
        [Tier::Beginner, Tier::Intermediate, Tier::Expert]
    }

    fn rank(self) -> u8 {
        match self {
            Tier::Beginner => 0,
            Tier::Intermediate => 1,
            Tier::Expert => 2,
        }
    }

    /// True if this tier is at least as advanced as `other`.
    pub fn at_least(self, other: Tier) -> bool {
        self.rank() >= other.rank()
    }

    pub fn is_beginner(self) -> bool {
        self == Tier::Beginner
    }

    pub fn is_expert(self) -> bool {
        self == Tier::Expert
    }

    /// Pick the wording appropriate to the tier.
    #[allow(dead_code)] // available for panels that want fully tier-specific copy
    pub fn pick<'a>(self, beginner: &'a str, intermediate: &'a str, expert: &'a str) -> &'a str {
        match self {
            Tier::Beginner => beginner,
            Tier::Intermediate => intermediate,
            Tier::Expert => expert,
        }
    }
}

/// An explanatory callout shown only to Beginner + Intermediate users.
///
/// Rendered as a subtle accent-tinted card with a coloured left rule — a
/// consistent "here's a hint" visual language distinct from chrome.
pub fn explain(ui: &mut egui::Ui, tier: Tier, text: &str) {
    if tier.is_expert() {
        return;
    }
    egui::Frame::new()
        .fill(tokens::ACCENT.gamma_multiply(0.09))
        .stroke(egui::Stroke::new(1.0, tokens::ACCENT.gamma_multiply(0.32)))
        .corner_radius(egui::CornerRadius::same(6))
        .inner_margin(egui::Margin::symmetric(9, 7))
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.spacing_mut().item_spacing.x = 6.0;
                ui.label(
                    RichText::new("TIP")
                        .size(10.0)
                        .strong()
                        .color(tokens::ACCENT_HOVER),
                );
                ui.label(RichText::new(text).size(11.5).color(tokens::TEXT));
            });
        });
}

/// A beginner-only callout (hidden for Intermediate and Expert).
pub fn explain_beginner(ui: &mut egui::Ui, tier: Tier, text: &str) {
    if tier.is_beginner() {
        explain(ui, tier, text);
    }
}

/// The standard "not medical advice" footer, used across workspaces.
pub fn disclaimer(ui: &mut egui::Ui) {
    ui.label(
        RichText::new("⚠ Research / educational use only — not medical advice.")
            .italics()
            .size(11.0)
            .color(palette::VUS),
    );
}

/// Render a tier picker (used in the top bar). Returns true if it changed.
pub fn tier_picker(ui: &mut egui::Ui, tier: &mut Tier) -> bool {
    let before = *tier;
    egui::ComboBox::from_id_salt("experience_tier")
        .selected_text(format!("Experience: {}", tier.label()))
        .show_ui(ui, |ui| {
            for t in Tier::all() {
                ui.selectable_value(tier, t, format!("{} — {}", t.label(), t.blurb()));
            }
        });
    *tier != before
}
