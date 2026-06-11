//! The genome-browser canvas: a stack of tracks (ruler, ideogram, gene models,
//! the sample's variants, and — when zoomed in — the reference sequence) drawn
//! with egui's painter, plus pan / zoom / click-to-select interaction.

use eframe::egui;
use egui::{Align2, Color32, FontId, Rect, Sense, Stroke, pos2, vec2};

use gx_core::{Assembly, Feature, FeatureKind};
use gx_annotate::VariantAnnotation;

use crate::document::VariantDoc;
use crate::theme::{self, palette};
use crate::view::GenomeView;

/// Everything the browser needs to render one frame. The app pre-filters the
/// reference features and sequence to the visible region.
pub struct BrowserData<'a> {
    pub doc: Option<&'a VariantDoc>,
    pub features: Vec<&'a Feature>,
    /// `(region_start_0based, bases)` covering (at least) the visible window.
    pub sequence: Option<(u64, &'a str)>,
    pub selected: Option<(String, u64)>,
}

pub struct ClickedVariant {
    pub contig: String,
    pub pos: u64,
    pub rsid: Option<String>,
}

/// Colour for a variant lollipop: ClinVar significance when annotated, else the
/// neutral "variant" blue.
fn variant_color(ann: Option<&VariantAnnotation>) -> Color32 {
    match ann.and_then(|a| a.clinical_significance.as_deref()) {
        Some(sig) => theme::sig_color(sig),
        None => palette::VARIANT,
    }
}

pub fn show(
    ui: &mut egui::Ui,
    view: &mut GenomeView,
    assembly: Assembly,
    data: &BrowserData,
) -> Option<ClickedVariant> {
    let size = ui.available_size();
    let (response, painter) = ui.allocate_painter(size, Sense::click_and_drag());
    let rect = response.rect;
    painter.rect_filled(rect, 0.0, ui.visuals().extreme_bg_color);

    // ---- interaction ----
    if response.dragged() {
        view.pan_px(response.drag_delta().x);
    }
    let scroll_y = ui.input(|i| i.smooth_scroll_delta.y);
    if scroll_y != 0.0
        && let Some(p) = response.hover_pos() {
            let anchor = view.pos_of(p.x, rect.left(), rect.width());
            let factor = if scroll_y > 0.0 { 0.85 } else { 1.0 / 0.85 };
            view.zoom(factor, anchor);
        }
    let contig_len = gx_core::chrom_length(assembly, &view.contig).unwrap_or(250_000_000) as f64;
    view.clamp_to(contig_len);

    let (lo, hi) = view.visible_range(rect.width());

    // ---- track bands ----
    let mut y = rect.top() + 4.0;
    draw_ruler(&painter, rect, view, y);
    y += 30.0;
    draw_ideogram(&painter, rect, view, contig_len, y);
    y += 22.0;
    let gene_top = y;
    draw_genes(&painter, rect, view, &data.features, gene_top);
    y += 48.0;

    let seq_h = if view.shows_bases() && data.sequence.is_some() {
        20.0
    } else {
        0.0
    };
    let var_top = y;
    let var_bottom = rect.bottom() - seq_h - 6.0;
    let clicked = draw_variants(
        &painter,
        &response,
        rect,
        view,
        data,
        lo,
        hi,
        var_top,
        var_bottom,
    );

    if seq_h > 0.0 {
        draw_sequence(&painter, rect, view, data, rect.bottom() - seq_h);
    }

    // selection guide line spanning all tracks
    if let Some((c, pos)) = &data.selected
        && *c == view.contig && (*pos as f64) >= lo && (*pos as f64) <= hi {
            let x = view.x_of(*pos as f64, rect.left(), rect.width());
            painter.line_segment(
                [pos2(x, rect.top() + 52.0), pos2(x, var_bottom)],
                Stroke::new(1.0, palette::ACCENT.gamma_multiply(0.8)),
            );
        }

    clicked
}

fn draw_ruler(painter: &egui::Painter, rect: Rect, view: &GenomeView, top: f32) {
    let baseline = top + 22.0;
    painter.line_segment(
        [pos2(rect.left(), baseline), pos2(rect.right(), baseline)],
        Stroke::new(1.0, palette::RULER),
    );
    let (lo, hi) = view.visible_range(rect.width());
    let step = theme::nice_step((hi - lo) / 8.0);
    let mut tick = (lo / step).ceil() * step;
    while tick < hi {
        if tick >= 0.0 {
            let x = view.x_of(tick, rect.left(), rect.width());
            painter.line_segment(
                [pos2(x, baseline - 5.0), pos2(x, baseline + 4.0)],
                Stroke::new(1.0, palette::RULER),
            );
            painter.text(
                pos2(x, top),
                Align2::CENTER_TOP,
                theme::format_pos(tick as u64),
                FontId::monospace(10.0),
                palette::RULER_TEXT,
            );
        }
        tick += step;
    }
}

fn draw_ideogram(painter: &egui::Painter, rect: Rect, view: &GenomeView, contig_len: f64, top: f32) {
    let h = 10.0;
    let bar = Rect::from_min_size(pos2(rect.left() + 8.0, top), vec2(rect.width() - 16.0, h));
    painter.rect_filled(bar, 4.0, Color32::from_gray(54));
    painter.rect_stroke(
        bar,
        4.0,
        Stroke::new(1.0, palette::RULER),
        egui::StrokeKind::Inside,
    );
    // highlight the visible window
    let (lo, hi) = view.visible_range(rect.width());
    let x0 = bar.left() + (lo.max(0.0) / contig_len) as f32 * bar.width();
    let x1 = bar.left() + (hi.min(contig_len) / contig_len) as f32 * bar.width();
    let win = Rect::from_min_max(
        pos2(x0.max(bar.left()), top - 1.0),
        pos2(x1.min(bar.right()), top + h + 1.0),
    );
    painter.rect_filled(win, 3.0, palette::ACCENT.gamma_multiply(0.7));
    painter.text(
        pos2(rect.right() - 10.0, top + h + 2.0),
        Align2::RIGHT_TOP,
        format!("chr{}", view.contig),
        FontId::proportional(10.0),
        palette::RULER_TEXT,
    );
}

fn draw_genes(painter: &egui::Painter, rect: Rect, view: &GenomeView, features: &[&Feature], top: f32) {
    let gene_y = top + 6.0;
    let exon_y = top + 22.0;
    for f in features {
        let x0 = view
            .x_of(f.range.start as f64, rect.left(), rect.width())
            .max(rect.left());
        let x1 = view
            .x_of(f.range.end as f64, rect.left(), rect.width())
            .min(rect.right());
        if x1 < x0 {
            continue;
        }
        match f.kind {
            FeatureKind::Gene => {
                painter.line_segment(
                    [pos2(x0, gene_y), pos2(x1, gene_y)],
                    Stroke::new(2.0, palette::GENE),
                );
                painter.text(
                    pos2(x0, gene_y - 12.0),
                    Align2::LEFT_BOTTOM,
                    f.label(),
                    FontId::proportional(11.0),
                    palette::GENE,
                );
            }
            FeatureKind::Exon | FeatureKind::Cds => {
                let r = Rect::from_min_max(pos2(x0, exon_y - 5.0), pos2(x1.max(x0 + 1.5), exon_y + 5.0));
                painter.rect_filled(r, 1.0, palette::EXON);
            }
            _ => {}
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_variants(
    painter: &egui::Painter,
    response: &egui::Response,
    rect: Rect,
    view: &GenomeView,
    data: &BrowserData,
    lo: f64,
    hi: f64,
    top: f32,
    bottom: f32,
) -> Option<ClickedVariant> {
    let axis_y = bottom;
    painter.line_segment(
        [pos2(rect.left(), axis_y), pos2(rect.right(), axis_y)],
        Stroke::new(1.0, palette::RULER.gamma_multiply(0.6)),
    );

    let Some(doc) = data.doc else {
        painter.text(
            rect.center(),
            Align2::CENTER_CENTER,
            "Open a genome file (File ▸ Open). Drag to pan, scroll to zoom",
            FontId::proportional(14.0),
            palette::RULER_TEXT,
        );
        return None;
    };

    let start = lo.max(0.0) as u64;
    let end = hi.max(0.0) as u64;
    let visible = doc.store.variants_in(&view.contig, start, end);

    // Too many to draw individually: render a density histogram instead.
    if visible.len() > 3000 {
        draw_density(painter, rect, view, visible, top, bottom);
        painter.text(
            pos2(rect.left() + 8.0, top),
            Align2::LEFT_TOP,
            format!("{} variants, zoom in to inspect", visible.len()),
            FontId::proportional(11.0),
            palette::RULER_TEXT,
        );
        return None;
    }

    let head_y = top + 8.0;
    let mut clicked = None;
    let click_pos = response
        .clicked()
        .then(|| response.interact_pointer_pos())
        .flatten();
    let mut best: Option<(f32, &gx_core::Variant)> = None;

    for v in visible {
        let x = view.x_of(v.pos as f64, rect.left(), rect.width());
        let ann = doc.annotation_for(v);
        let color = variant_color(ann);
        painter.line_segment([pos2(x, axis_y), pos2(x, head_y)], Stroke::new(1.0, color));
        if ann.map(|a| a.is_clinically_notable()).unwrap_or(false) {
            // Emphasise clinically-significant variants beyond colour alone
            // (a white ring + larger head reads even for colour-blind users).
            painter.circle_filled(pos2(x, head_y), 5.0, color);
            painter.circle_stroke(pos2(x, head_y), 5.0, Stroke::new(1.5, Color32::WHITE));
        } else {
            painter.circle_filled(pos2(x, head_y), 3.5, color);
        }

        if let Some(cp) = click_pos {
            let d = (cp.x - x).abs();
            if d < 6.0 && best.map(|(bd, _)| d < bd).unwrap_or(true) {
                best = Some((d, v));
            }
        }
    }
    if let Some((_, v)) = best {
        clicked = Some(ClickedVariant {
            contig: v.contig.clone(),
            pos: v.pos,
            rsid: v.rsid.clone(),
        });
    }
    clicked
}

fn draw_density(
    painter: &egui::Painter,
    rect: Rect,
    view: &GenomeView,
    visible: &[gx_core::Variant],
    top: f32,
    bottom: f32,
) {
    let cols = (rect.width() as usize / 3).max(1);
    let mut bins = vec![0u32; cols];
    let (lo, hi) = view.visible_range(rect.width());
    let span = (hi - lo).max(1.0);
    for v in visible {
        let frac = (v.pos as f64 - lo) / span;
        if (0.0..1.0).contains(&frac) {
            bins[(frac * cols as f64) as usize] += 1;
        }
    }
    let max = bins.iter().copied().max().unwrap_or(1).max(1) as f32;
    let col_w = rect.width() / cols as f32;
    for (i, &c) in bins.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let h = (c as f32 / max) * (bottom - top - 4.0);
        let x = rect.left() + i as f32 * col_w;
        let r = Rect::from_min_max(pos2(x, bottom - h), pos2(x + col_w.max(1.0), bottom));
        painter.rect_filled(r, 0.0, palette::VARIANT.gamma_multiply(0.8));
    }
}

fn draw_sequence(painter: &egui::Painter, rect: Rect, view: &GenomeView, data: &BrowserData, top: f32) {
    let Some((region_start, seq)) = data.sequence else {
        return;
    };
    let (lo, hi) = view.visible_range(rect.width());
    let seq_bytes = seq.as_bytes();
    let first = lo.floor().max(region_start as f64) as u64;
    let last = (hi.ceil() as u64).min(region_start + seq_bytes.len() as u64);
    for pos in first..last {
        let idx = (pos - region_start) as usize;
        let Some(&b) = seq_bytes.get(idx) else { continue };
        let x = view.x_of(pos as f64 + 0.5, rect.left(), rect.width());
        let color = theme::base_color(b);
        painter.text(
            pos2(x, top + 2.0),
            Align2::CENTER_TOP,
            (b as char).to_string(),
            FontId::monospace(13.0),
            color,
        );
    }
}

