//! [`GenomeView`] — the genome browser viewport and its coordinate math.
//!
//! The viewport is described by a contig, a centre position (in bases), and a
//! zoom expressed as **bases per pixel**. All screen<->genome mapping goes
//! through here so the canvas, ruler, and tracks stay perfectly aligned.

/// Most zoomed-in: ~20 px per base (sequence letters are legible).
const MIN_BP_PER_PX: f64 = 0.05;
/// Most zoomed-out: a megabase every couple of pixels.
const MAX_BP_PER_PX: f64 = 2_000_000.0;

#[derive(Debug, Clone)]
pub struct GenomeView {
    pub contig: String,
    /// Genomic position (0-based, fractional) at the horizontal centre.
    pub center: f64,
    /// Zoom: bases per screen pixel. Larger = more zoomed out.
    pub bp_per_px: f64,
}

impl Default for GenomeView {
    fn default() -> Self {
        Self {
            contig: "1".to_string(),
            center: 1_000_000.0,
            bp_per_px: 2_000.0,
        }
    }
}

impl GenomeView {
    /// The genomic `[start, end)` span visible across a viewport `width_px` wide.
    pub fn visible_range(&self, width_px: f32) -> (f64, f64) {
        let half = width_px as f64 * 0.5 * self.bp_per_px;
        (self.center - half, self.center + half)
    }

    /// Screen x for a genomic position, given the viewport's left edge + width.
    pub fn x_of(&self, pos: f64, left_px: f32, width_px: f32) -> f32 {
        let (lo, _) = self.visible_range(width_px);
        left_px + ((pos - lo) / self.bp_per_px) as f32
    }

    /// Genomic position under a screen x.
    pub fn pos_of(&self, x: f32, left_px: f32, width_px: f32) -> f64 {
        let (lo, _) = self.visible_range(width_px);
        lo + (x - left_px) as f64 * self.bp_per_px
    }

    /// Pan by a horizontal pixel delta (dragging right moves the view left).
    pub fn pan_px(&mut self, dx_px: f32) {
        self.center -= dx_px as f64 * self.bp_per_px;
    }

    /// Zoom by `factor` (<1 zooms in) keeping `anchor_pos` fixed under the cursor.
    pub fn zoom(&mut self, factor: f64, anchor_pos: f64) {
        let new_bp = (self.bp_per_px * factor).clamp(MIN_BP_PER_PX, MAX_BP_PER_PX);
        // Preserve the anchor's pixel offset from centre across the zoom.
        let px_offset = (anchor_pos - self.center) / self.bp_per_px;
        self.bp_per_px = new_bp;
        self.center = anchor_pos - px_offset * new_bp;
    }

    /// Keep the centre within `[0, contig_len)` after pans/zooms.
    pub fn clamp_to(&mut self, contig_len: f64) {
        self.center = self.center.clamp(0.0, contig_len.max(1.0));
    }

    /// Jump the view to a span, framing it with a little padding.
    pub fn go_to(&mut self, contig: impl Into<String>, start: f64, end: f64, width_px: f32) {
        self.contig = contig.into();
        self.center = (start + end) * 0.5;
        let span = (end - start).max(1.0) * 1.2;
        self.bp_per_px = (span / width_px.max(1.0) as f64).clamp(MIN_BP_PER_PX, MAX_BP_PER_PX);
    }

    /// True when zoomed in far enough that each base is ≳8 px (letters legible).
    pub fn shows_bases(&self) -> bool {
        self.bp_per_px <= 0.125
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x_pos_round_trip() {
        let v = GenomeView {
            contig: "1".into(),
            center: 1000.0,
            bp_per_px: 10.0,
        };
        let width = 800.0;
        for &pos in &[500.0, 1000.0, 1500.0] {
            let x = v.x_of(pos, 0.0, width);
            let back = v.pos_of(x, 0.0, width);
            assert!((back - pos).abs() < 0.5, "pos {pos} -> {back}");
        }
        // Centre maps to the middle of the viewport.
        assert!((v.x_of(1000.0, 0.0, width) - 400.0).abs() < 0.5);
    }

    #[test]
    fn zoom_keeps_anchor_fixed() {
        let mut v = GenomeView {
            contig: "1".into(),
            center: 1000.0,
            bp_per_px: 10.0,
        };
        let width = 800.0;
        let anchor = v.pos_of(600.0, 0.0, width); // a point right of centre
        v.zoom(0.5, anchor);
        // The anchor should still sit at x=600 after zooming in.
        assert!((v.x_of(anchor, 0.0, width) - 600.0).abs() < 0.5);
        assert!((v.bp_per_px - 5.0).abs() < 1e-9);
    }

    #[test]
    fn zoom_is_clamped() {
        let mut v = GenomeView::default();
        for _ in 0..200 {
            v.zoom(0.5, v.center);
        }
        assert!(v.bp_per_px >= MIN_BP_PER_PX);
        for _ in 0..200 {
            v.zoom(2.0, v.center);
        }
        assert!(v.bp_per_px <= MAX_BP_PER_PX);
    }
}
