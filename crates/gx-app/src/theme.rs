//! GenomeForge's visual identity: a modern, opinionated dark theme built from a
//! small set of surface/text/accent tokens, a deliberate type scale, and
//! consistent spacing & rounding. Plus the shared track colour palette.
//!
//! Design intent — *modern, clean, but opinionated*:
//! - one confident brand accent (a luminous cyan-blue, echoing the helix icon),
//! - a near-black blue-slate surface ramp (not flat grey), layered by elevation,
//! - text contrast and a small-text floor tuned for accessibility (WCAG-ish),
//! - colour-blind-safe semantics: red/amber/green are *backed by shape channels*
//!   elsewhere (rings, glyphs) so meaning never rides on hue alone.

use egui::Color32;

/// Surface, text, and accent tokens — one source of truth for the whole app, so
/// a rebrand is a few-line change. Tracks/variants live in [`palette`].
pub mod tokens {
    use egui::Color32;

    // Elevation ramp — deepest (window) to highest (active widget). Blue-slate,
    // not neutral grey: gives the dark UI a cooler, more intentional feel.
    pub const BG: Color32 = Color32::from_rgb(15, 18, 23); // app / canvas base
    pub const PANEL: Color32 = Color32::from_rgb(22, 26, 32); // side & top panels
    pub const SURFACE: Color32 = Color32::from_rgb(30, 35, 43); // cards, buttons at rest
    pub const SURFACE_HOVER: Color32 = Color32::from_rgb(41, 47, 56);
    pub const SURFACE_ACTIVE: Color32 = Color32::from_rgb(51, 59, 70);
    pub const FAINT: Color32 = Color32::from_rgb(26, 30, 37); // striped rows
    pub const BORDER: Color32 = Color32::from_rgb(46, 53, 62);
    pub const BORDER_STRONG: Color32 = Color32::from_rgb(64, 73, 85);

    // Text — primary is high-contrast; muted clears ~4.5:1 on PANEL.
    pub const TEXT: Color32 = Color32::from_rgb(228, 232, 238);
    pub const TEXT_MUTED: Color32 = Color32::from_rgb(166, 174, 185);
    pub const TEXT_FAINT: Color32 = Color32::from_rgb(126, 134, 145);

    // Brand accent — luminous cyan-blue from the helix mark.
    pub const ACCENT: Color32 = Color32::from_rgb(86, 182, 224);
    pub const ACCENT_HOVER: Color32 = Color32::from_rgb(122, 206, 242);
}

/// Colours used across the genome tracks and variant views.
pub mod palette {
    use super::tokens;
    use egui::Color32;

    pub use tokens::ACCENT;
    /// Muted body/caption text on dark surfaces (re-exported for call sites).
    pub const RULER_TEXT: Color32 = tokens::TEXT_MUTED;
    pub const RULER: Color32 = Color32::from_rgb(74, 82, 92);
    pub const GENE: Color32 = Color32::from_rgb(104, 156, 226);
    pub const EXON: Color32 = Color32::from_rgb(146, 186, 255);
    pub const VARIANT: Color32 = Color32::from_rgb(132, 198, 236);

    // ClinVar significance colouring. Hue is reinforced with shape channels
    // (rings / glyphs) at the call sites, so these stay colour-blind-safe.
    pub const PATHOGENIC: Color32 = Color32::from_rgb(232, 92, 86);
    pub const LIKELY_PATH: Color32 = Color32::from_rgb(232, 142, 78);
    pub const VUS: Color32 = Color32::from_rgb(224, 188, 92);
    pub const BENIGN: Color32 = Color32::from_rgb(92, 192, 124);
    pub const UNKNOWN: Color32 = Color32::from_rgb(140, 148, 158);
}

/// Apply the theme to a freshly created context.
pub fn apply(ctx: &egui::Context) {
    use egui::{CornerRadius, FontFamily, FontId, Stroke, TextStyle};
    use tokens::*;

    let mut style = (*ctx.global_style()).clone();
    let mut v = egui::Visuals::dark();

    // --- surfaces ---------------------------------------------------------
    v.panel_fill = PANEL;
    v.window_fill = PANEL;
    v.extreme_bg_color = BG; // canvas, text-edit, combo backgrounds
    v.faint_bg_color = FAINT; // striped table rows
    v.code_bg_color = BG;
    v.window_stroke = Stroke::new(1.0, BORDER);
    v.window_corner_radius = CornerRadius::same(10);
    v.menu_corner_radius = CornerRadius::same(8);

    // --- accent & selection ----------------------------------------------
    v.hyperlink_color = ACCENT;
    v.selection.bg_fill = ACCENT.gamma_multiply(0.36);
    v.selection.stroke = Stroke::new(1.0, ACCENT);

    // --- widget states (rest / hover / pressed / open) -------------------
    let radius = CornerRadius::same(6);
    let w = &mut v.widgets;
    w.noninteractive.bg_fill = PANEL;
    w.noninteractive.weak_bg_fill = PANEL;
    w.noninteractive.bg_stroke = Stroke::new(1.0, BORDER);
    w.noninteractive.fg_stroke = Stroke::new(1.0, TEXT);
    w.noninteractive.corner_radius = radius;

    w.inactive.bg_fill = SURFACE;
    w.inactive.weak_bg_fill = SURFACE;
    w.inactive.bg_stroke = Stroke::new(1.0, BORDER);
    w.inactive.fg_stroke = Stroke::new(1.0, TEXT);
    w.inactive.corner_radius = radius;

    w.hovered.bg_fill = SURFACE_HOVER;
    w.hovered.weak_bg_fill = SURFACE_HOVER;
    w.hovered.bg_stroke = Stroke::new(1.0, ACCENT_HOVER.gamma_multiply(0.7));
    w.hovered.fg_stroke = Stroke::new(1.5, TEXT);
    w.hovered.corner_radius = radius;
    w.hovered.expansion = 1.0;

    w.active.bg_fill = SURFACE_ACTIVE;
    w.active.weak_bg_fill = SURFACE_ACTIVE;
    w.active.bg_stroke = Stroke::new(1.0, ACCENT);
    w.active.fg_stroke = Stroke::new(1.5, TEXT);
    w.active.corner_radius = radius;
    w.active.expansion = 1.0;

    w.open.bg_fill = SURFACE;
    w.open.weak_bg_fill = SURFACE;
    w.open.bg_stroke = Stroke::new(1.0, BORDER_STRONG);
    w.open.fg_stroke = Stroke::new(1.0, TEXT);
    w.open.corner_radius = radius;

    style.visuals = v;

    // --- typography: a deliberate scale; Small floored at 11px (a11y) ----
    style.text_styles = [
        (TextStyle::Heading, FontId::new(18.5, FontFamily::Proportional)),
        (TextStyle::Body, FontId::new(13.5, FontFamily::Proportional)),
        (TextStyle::Button, FontId::new(13.5, FontFamily::Proportional)),
        (TextStyle::Small, FontId::new(11.0, FontFamily::Proportional)),
        (TextStyle::Monospace, FontId::new(12.5, FontFamily::Monospace)),
    ]
    .into();

    // --- spacing & rhythm -------------------------------------------------
    style.spacing.item_spacing = egui::vec2(8.0, 7.0);
    style.spacing.button_padding = egui::vec2(11.0, 5.0);
    style.spacing.window_margin = egui::Margin::same(12);
    style.spacing.menu_margin = egui::Margin::same(8);
    style.spacing.indent = 18.0;
    style.spacing.interact_size.y = 26.0;
    style.spacing.scroll.bar_width = 9.0;

    ctx.set_global_style(style);
}

// ---- shared visualization helpers (one definition, used across workspaces) --

/// Map a ClinVar clinical-significance string to a colour.
pub fn sig_color(sig: &str) -> Color32 {
    let s = sig.to_ascii_lowercase();
    if s.contains("pathogenic") && s.contains("likely") {
        palette::LIKELY_PATH
    } else if s.contains("pathogenic") {
        palette::PATHOGENIC
    } else if s.contains("benign") {
        palette::BENIGN
    } else if s.contains("uncertain") || s.contains("conflicting") {
        palette::VUS
    } else {
        palette::UNKNOWN
    }
}

/// Colour for a DNA base; non-ACGT use the muted ruler colour.
pub fn base_color(b: u8) -> Color32 {
    match b.to_ascii_uppercase() {
        b'A' => Color32::from_rgb(96, 200, 120),
        b'C' => Color32::from_rgb(96, 150, 230),
        b'G' => Color32::from_rgb(230, 180, 80),
        b'T' => Color32::from_rgb(220, 100, 100),
        _ => palette::RULER_TEXT,
    }
}

/// Round a raw spacing up to a "nice" 1/2/5 × 10ⁿ value (ruler ticks).
pub fn nice_step(raw: f64) -> f64 {
    if raw <= 0.0 {
        return 1.0;
    }
    let pow = 10f64.powf(raw.log10().floor());
    let n = raw / pow;
    let mult = if n < 1.5 {
        1.0
    } else if n < 3.0 {
        2.0
    } else if n < 7.0 {
        5.0
    } else {
        10.0
    };
    (mult * pow).max(1.0)
}

/// Human-readable genomic position, e.g. `1.23 Mb`, `500 kb`, `123 bp`.
pub fn format_pos(p: u64) -> String {
    if p >= 1_000_000 {
        format!("{:.2} Mb", p as f64 / 1e6)
    } else if p >= 1_000 {
        format!("{:.0} kb", p as f64 / 1e3)
    } else {
        format!("{p} bp")
    }
}

/// Compact genomic position for tight labels (circular map), e.g. `1.2M`.
pub fn format_pos_compact(p: u64) -> String {
    if p >= 1_000_000 {
        format!("{:.2}M", p as f64 / 1e6)
    } else if p >= 1_000 {
        format!("{:.1}k", p as f64 / 1e3)
    } else {
        p.to_string()
    }
}
