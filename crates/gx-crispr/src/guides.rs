//! SpCas9 guide (sgRNA) enumeration: find protospacers adjacent to an NGG PAM
//! on both strands and score them.
//!
//! [`find_guides`] is the SpCas9-specific entry point and is preserved for
//! backward compatibility; it now delegates to
//! [`crate::nuclease::find_guides_for`] with [`Nuclease::SpCas9`], which
//! generalizes the same logic to other nucleases/PAMs.

use gx_core::Strand;
use gx_plasmid::{gc_content, reverse_complement};

use crate::nuclease::{Nuclease, find_guides_for};
use crate::score::{has_poly_t, on_target_estimate};

pub const PROTOSPACER_LEN: usize = 20;

#[derive(Debug, Clone)]
pub struct Guide {
    /// 0-based top-strand start of the protospacer footprint.
    pub start: usize,
    pub strand: Strand,
    /// Spacer, 5'→3' (the guide/target-strand sequence); length depends on the
    /// nuclease (20 nt for the Cas9 family, 23 nt for Cas12a).
    pub protospacer: String,
    /// PAM, 5'→3'.
    pub pam: String,
    pub gc: f64,
    /// Heuristic on-target efficiency estimate in `[0, 1]`.
    pub on_score: f64,
    pub poly_t: bool,
    /// 0-based top-strand position of the predicted blunt DSB (3 bp 5' of PAM).
    pub cut_site: usize,
    /// The nuclease this guide was enumerated for.
    pub nuclease: Nuclease,
}

/// Construct a [`Guide`] for an arbitrary nuclease. Crate-internal helper shared
/// with [`crate::nuclease`]; the on-target score is supplied by the caller so
/// each enumerator can choose the appropriate model.
#[allow(clippy::too_many_arguments)]
pub(crate) fn make_guide_for(
    start: usize,
    strand: Strand,
    proto: &[u8],
    pam: &[u8],
    cut: usize,
    on_score: f64,
    nuclease: Nuclease,
) -> Guide {
    Guide {
        start,
        strand,
        protospacer: String::from_utf8_lossy(proto).into_owned(),
        pam: String::from_utf8_lossy(pam).into_owned(),
        gc: gc_content(proto),
        on_score,
        poly_t: has_poly_t(proto),
        cut_site: cut,
        nuclease,
    }
}

/// SpCas9 on-target score: the exact Doench 2014 model computed from the 30-mer
/// context (4 nt + 20 nt + 3 nt PAM + 3 nt) when that context is in-bounds,
/// otherwise the documented heuristic (only happens near sequence ends).
///
/// Only valid for SpCas9 — the Doench model was trained on SpCas9 and is not
/// applied to other nucleases.
pub(crate) fn on_target_spcas9(up: &[u8], start: usize, strand: Strand, proto: &[u8]) -> f64 {
    let n = up.len();
    let mer: Option<Vec<u8>> = match strand {
        Strand::Reverse if start >= 6 && start + 24 <= n => {
            Some(reverse_complement(&up[start - 6..start + 24]))
        }
        Strand::Reverse => None,
        _ if start >= 4 && start + 26 <= n => Some(up[start - 4..start + 26].to_vec()),
        _ => None,
    };
    mer.as_deref()
        .and_then(crate::doench::score_30mer)
        .unwrap_or_else(|| on_target_estimate(proto))
}

/// Enumerate all SpCas9 (NGG) guides on both strands, sorted by on-target score.
///
/// Preserved for backward compatibility; equivalent to
/// `find_guides_for(seq, Nuclease::SpCas9)`.
pub fn find_guides(seq: &[u8]) -> Vec<Guide> {
    find_guides_for(seq, Nuclease::SpCas9)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_forward_guide() {
        // 20 nt + "TGG" PAM at the end.
        let seq = b"ACGTACGTACGTACGTACGTTGG";
        let guides = find_guides(seq);
        let fwd: Vec<_> = guides.iter().filter(|g| g.strand == Strand::Forward).collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].protospacer, "ACGTACGTACGTACGTACGT");
        assert_eq!(fwd[0].pam, "TGG");
        assert_eq!(fwd[0].start, 0);
        assert_eq!(fwd[0].cut_site, 17); // 3 bp 5' of the PAM at index 20
    }

    #[test]
    fn finds_reverse_guide() {
        // "CCA" PAM on top strand -> NGG on the reverse strand guide.
        let seq = b"CCAACGTACGTACGTACGTACGT";
        let guides = find_guides(seq);
        assert!(guides.iter().any(|g| g.strand == Strand::Reverse && g.pam == "TGG"));
    }
}
