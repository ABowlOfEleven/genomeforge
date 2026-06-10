//! Base-editing and prime-editing outcome simulation.
//!
//! These are **transparent, applied-edit models**, not predictors of editing
//! *efficiency*. Base editing reports which bases a deaminase would convert
//! within its canonical activity window; prime editing reports the sequence that
//! results when a pegRNA installs exactly the specified edit. Neither estimates
//! the probability that the edit occurs — they are decision-support views of the
//! intended outcome, not a clinical pipeline.

use crate::guides::Guide;

/// A base editor, defined by the single-base substitution its deaminase performs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BaseEditor {
    /// Cytosine base editor (CBE): deaminates C → T (on the protospacer/target strand).
    CytosineCBE,
    /// Adenine base editor (ABE): deaminates A → G.
    AdenineABE,
}

impl BaseEditor {
    /// The base this editor acts on (upper-case).
    pub fn target(self) -> char {
        match self {
            BaseEditor::CytosineCBE => 'C',
            BaseEditor::AdenineABE => 'A',
        }
    }

    /// The base it is converted to (upper-case).
    pub fn result(self) -> char {
        match self {
            BaseEditor::CytosineCBE => 'T',
            BaseEditor::AdenineABE => 'G',
        }
    }

    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            BaseEditor::CytosineCBE => "Cytosine base editor (C→T)",
            BaseEditor::AdenineABE => "Adenine base editor (A→G)",
        }
    }
}

/// A single base substitution applied by a base editor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseEdit {
    /// 0-based position within the protospacer (counted from the 5' end).
    pub pos_in_protospacer: usize,
    pub from: char,
    pub to: char,
}

/// The predicted outcome of applying a base editor to a guide's protospacer.
#[derive(Debug, Clone)]
pub struct BaseEditOutcome {
    pub editor: BaseEditor,
    /// 0-based, half-open `[lo, hi)` editing window within the protospacer.
    pub window: (usize, usize),
    pub edits: Vec<BaseEdit>,
    /// The protospacer string after applying every edit.
    pub edited_protospacer: String,
}

/// The canonical base-editing activity window, 1-based and inclusive, measured
/// from the 5' (PAM-distal) end of a 20-nt protospacer. Positions 4–8 are the
/// widely cited "edit window" for typical BE3/ABE7.10-class editors.
const WINDOW_1BASED: (usize, usize) = (4, 8);

/// Apply `editor` to `guide`'s protospacer over the canonical editing window.
///
/// The window is positions 4..=8 (1-based, from the 5'/PAM-distal end), per the
/// commonly reported BE3/ABE activity window; it is clamped to the protospacer
/// length so shorter spacers (or longer Cas12a spacers) stay in bounds. Every
/// occurrence of the editor's target base within the window is converted and
/// recorded. This models *where* a deaminase can act, not the per-base
/// efficiency of editing.
pub fn base_edit(guide: &Guide, editor: BaseEditor) -> BaseEditOutcome {
    let proto: Vec<char> = guide.protospacer.chars().collect();
    let len = proto.len();

    // Convert the 1-based inclusive window to a 0-based half-open range, clamped.
    let lo = WINDOW_1BASED.0.saturating_sub(1).min(len);
    let hi = WINDOW_1BASED.1.min(len); // inclusive end (1-based) == exclusive end (0-based)

    let target = editor.target();
    let result = editor.result();

    let mut edited = proto.clone();
    let mut edits = Vec::new();
    for (i, edited_ch) in edited.iter_mut().enumerate().take(hi).skip(lo) {
        if edited_ch.eq_ignore_ascii_case(&target) {
            edits.push(BaseEdit {
                pos_in_protospacer: i,
                from: *edited_ch,
                to: result,
            });
            *edited_ch = result;
        }
    }

    BaseEditOutcome {
        editor,
        window: (lo, hi),
        edits,
        edited_protospacer: edited.into_iter().collect(),
    }
}

/// A prime-editing operation to install into a template via a pegRNA.
///
/// `pos` is 0-based on the template. For [`PrimeEdit::Insertion`] the new
/// sequence is inserted *before* `pos`; for [`PrimeEdit::Deletion`] `len` bases
/// starting at `pos` are removed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrimeEdit {
    /// Replace the single base at `pos` with `to`.
    Substitution { pos: usize, to: char },
    /// Insert `seq` immediately before `pos`.
    Insertion { pos: usize, seq: String },
    /// Delete `len` bases starting at `pos`.
    Deletion { pos: usize, len: usize },
}

/// The result of applying a prime edit to a template.
#[derive(Debug, Clone)]
pub struct PrimeEditOutcome {
    pub original: String,
    pub edited: String,
    pub description: String,
}

/// Apply `edit` to `template`, returning the original and edited sequences plus a
/// short human description.
///
/// This is a deterministic "the pegRNA installs exactly this edit" model: it
/// computes the intended product sequence. It does **not** model prime-editing
/// efficiency, nick repair, or indel by-products. Positions are 0-based and
/// clamped to the template length, so out-of-range inputs degrade gracefully.
pub fn prime_edit(template: &[u8], edit: &PrimeEdit) -> PrimeEditOutcome {
    let original = String::from_utf8_lossy(template).into_owned();
    let n = template.len();

    let (edited, description) = match edit {
        PrimeEdit::Substitution { pos, to } => {
            let p = (*pos).min(n.saturating_sub(1));
            let mut out: Vec<char> = original.chars().collect();
            let from = out.get(p).copied().unwrap_or('?');
            let desc = if p < out.len() {
                let upper = to.to_ascii_uppercase();
                out[p] = upper;
                format!("Substitution {from}→{upper} at position {}", p + 1)
            } else {
                format!("Substitution at position {} out of range", pos + 1)
            };
            (out.into_iter().collect::<String>(), desc)
        }
        PrimeEdit::Insertion { pos, seq } => {
            let p = (*pos).min(n);
            let ins: String = seq.to_ascii_uppercase();
            let mut out = String::with_capacity(original.len() + ins.len());
            out.push_str(&original[..byte_index(&original, p)]);
            out.push_str(&ins);
            out.push_str(&original[byte_index(&original, p)..]);
            let desc = format!("Insertion of {} ({} bp) before position {}", ins, ins.len(), p + 1);
            (out, desc)
        }
        PrimeEdit::Deletion { pos, len } => {
            let p = (*pos).min(n);
            let end = (p + *len).min(n);
            let removed = end - p;
            let mut out = String::with_capacity(n - removed);
            out.push_str(&original[..byte_index(&original, p)]);
            out.push_str(&original[byte_index(&original, end)..]);
            let desc = format!("Deletion of {} bp starting at position {}", removed, p + 1);
            (out, desc)
        }
    };

    PrimeEditOutcome { original, edited, description }
}

/// Byte index of the `ch`-th character. Sequences here are ASCII nucleotides, so
/// this is the identity for valid input; it stays correct if a stray non-ASCII
/// byte sneaks in.
fn byte_index(s: &str, ch: usize) -> usize {
    s.char_indices().nth(ch).map(|(i, _)| i).unwrap_or(s.len())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nuclease::Nuclease;
    use gx_core::Strand;

    fn guide_with(protospacer: &str) -> Guide {
        Guide {
            start: 0,
            strand: Strand::Forward,
            protospacer: protospacer.to_string(),
            pam: "TGG".to_string(),
            gc: 0.0,
            on_score: 0.0,
            poly_t: false,
            cut_site: 0,
            nuclease: Nuclease::SpCas9,
        }
    }

    #[test]
    fn cbe_edits_cs_in_window_only() {
        // Window is 0-based [3, 8). Put C's inside and outside the window.
        // index:        0123456789...
        let proto = "AAACCCCCAAACCCAAAAAA"; // C at 3,4,5,6,7 (in window) and 11,12,13 (out)
        let g = guide_with(proto);
        let out = base_edit(&g, BaseEditor::CytosineCBE);
        assert_eq!(out.window, (3, 8));
        // Five C's at positions 3..=7 should be edited to T.
        assert_eq!(out.edits.len(), 5);
        assert!(out.edits.iter().all(|e| e.from == 'C' && e.to == 'T'));
        assert!(out.edits.iter().all(|e| (3..8).contains(&e.pos_in_protospacer)));
        // Window C's became T; the out-of-window C's at 11..14 are untouched.
        assert_eq!(&out.edited_protospacer[3..8], "TTTTT");
        assert_eq!(&out.edited_protospacer[11..14], "CCC");
    }

    #[test]
    fn abe_leaves_cs_untouched() {
        let proto = "AAACCCCCAAACCCAAAAAA";
        let g = guide_with(proto);
        let out = base_edit(&g, BaseEditor::AdenineABE);
        // No A's in the [3,8) window of this protospacer -> no edits.
        assert!(out.edits.is_empty());
        assert_eq!(out.edited_protospacer, proto);
    }

    #[test]
    fn abe_edits_as_in_window() {
        // A's at every position; window [3,8) holds five A's -> all become G.
        let proto = "AAAAAAAAAAAAAAAAAAAA";
        let g = guide_with(proto);
        let out = base_edit(&g, BaseEditor::AdenineABE);
        assert_eq!(out.edits.len(), 5);
        assert!(out.edits.iter().all(|e| e.from == 'A' && e.to == 'G'));
        assert_eq!(&out.edited_protospacer[3..8], "GGGGG");
    }

    #[test]
    fn prime_substitution() {
        let out = prime_edit(b"ACGTACGT", &PrimeEdit::Substitution { pos: 2, to: 'a' });
        assert_eq!(out.original, "ACGTACGT");
        assert_eq!(out.edited, "ACATACGT");
        assert!(out.description.contains("Substitution"));
    }

    #[test]
    fn prime_insertion() {
        let out = prime_edit(b"ACGT", &PrimeEdit::Insertion { pos: 2, seq: "TTT".into() });
        assert_eq!(out.edited, "ACTTTGT");
        // Insertion at the end.
        let end = prime_edit(b"ACGT", &PrimeEdit::Insertion { pos: 4, seq: "GG".into() });
        assert_eq!(end.edited, "ACGTGG");
    }

    #[test]
    fn prime_deletion() {
        let out = prime_edit(b"ACGTACGT", &PrimeEdit::Deletion { pos: 2, len: 3 });
        assert_eq!(out.edited, "ACCGT");
        // Over-long deletion clamps to the end.
        let clamp = prime_edit(b"ACGT", &PrimeEdit::Deletion { pos: 2, len: 10 });
        assert_eq!(clamp.edited, "AC");
    }
}
