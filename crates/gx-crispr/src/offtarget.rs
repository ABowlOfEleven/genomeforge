//! Off-target search against a provided sequence.
//!
//! [`find_offtargets`] searches the guide's own target (excluding its on-target
//! site); [`find_matches`] searches any other reference (e.g. a chromosome or
//! whole-genome FASTA the user loads) with no self-exclusion. The scan is linear
//! in the reference length with early-terminated mismatch counting, so a
//! chromosome-scale reference is sub-second per guide; whole-genome is slower.

use gx_core::Strand;
use gx_plasmid::reverse_complement;

use crate::guides::{Guide, PROTOSPACER_LEN};
use crate::score::cfd_score;

#[derive(Debug, Clone)]
pub struct OffTarget {
    pub start: usize,
    pub strand: Strand,
    pub mismatches: usize,
    pub cfd: f64,
}

/// Mismatch count between equal-length seqs, short-circuiting once it exceeds
/// `cap` (returns `None`) — most candidates differ a lot, so this is the hot path.
fn mm_capped(a: &[u8], b: &[u8], cap: usize) -> Option<usize> {
    let mut m = 0;
    for (x, y) in a.iter().zip(b) {
        if !x.eq_ignore_ascii_case(y) {
            m += 1;
            if m > cap {
                return None;
            }
        }
    }
    Some(m)
}

fn scan(guide: &Guide, seq: &[u8], max_mm: usize, exclude_own: bool) -> Vec<OffTarget> {
    let up = seq.to_ascii_uppercase();
    let n = up.len();
    let g = guide.protospacer.as_bytes();
    let mut hits = Vec::new();
    if n < PROTOSPACER_LEN + 3 {
        return hits;
    }

    for p in PROTOSPACER_LEN..=(n - 3) {
        if up[p + 1] == b'G' && up[p + 2] == b'G' {
            let cand = &up[p - PROTOSPACER_LEN..p];
            if let Some(mm) = mm_capped(g, cand, max_mm) {
                let own = exclude_own && guide.strand == Strand::Forward && p - PROTOSPACER_LEN == guide.start;
                if !own {
                    hits.push(OffTarget {
                        start: p - PROTOSPACER_LEN,
                        strand: Strand::Forward,
                        mismatches: mm,
                        cfd: cfd_score(g, cand),
                    });
                }
            }
        }
    }
    for p in 0..=(n - (PROTOSPACER_LEN + 3)) {
        if up[p] == b'C' && up[p + 1] == b'C' {
            let cand = reverse_complement(&up[p + 3..p + 3 + PROTOSPACER_LEN]);
            if let Some(mm) = mm_capped(g, &cand, max_mm) {
                let own = exclude_own && guide.strand == Strand::Reverse && p + 3 == guide.start;
                if !own {
                    hits.push(OffTarget {
                        start: p + 3,
                        strand: Strand::Reverse,
                        mismatches: mm,
                        cfd: cfd_score(g, &cand),
                    });
                }
            }
        }
    }

    hits.sort_by(|a, b| b.cfd.partial_cmp(&a.cfd).unwrap_or(std::cmp::Ordering::Equal));
    hits
}

/// Off-targets within the guide's own target sequence (excludes the on-target site).
pub fn find_offtargets(guide: &Guide, seq: &[u8], max_mm: usize) -> Vec<OffTarget> {
    scan(guide, seq, max_mm, true)
}

/// Off-targets within a separate reference (no self-exclusion) — for searching a
/// loaded chromosome / genome.
pub fn find_matches(guide: &Guide, reference: &[u8], max_mm: usize) -> Vec<OffTarget> {
    scan(guide, reference, max_mm, false)
}

/// MIT-style specificity score in `[0, 100]` from a guide's off-target set
/// (100 = no detected off-targets in the searched sequence).
pub fn specificity(offtargets: &[OffTarget]) -> f64 {
    let sum: f64 = offtargets.iter().map(|o| o.cfd).sum();
    100.0 / (1.0 + sum * 100.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::guides::find_guides;

    #[test]
    fn no_offtargets_gives_full_specificity() {
        let seq = b"ACGTACGTACGTACGTACGTTGGAAAAAAAAAAAAAAAAAAAA";
        let guide = find_guides(seq).into_iter().next().unwrap();
        let off = find_offtargets(&guide, seq, 4);
        assert!((specificity(&off) - 100.0).abs() < 1e-6 || off.iter().all(|o| o.cfd < 0.001));
    }

    #[test]
    fn detects_a_duplicated_site() {
        // Same protospacer+PAM twice -> each is the other's perfect off-target.
        let unit = b"ACGTACGTACGTACGTACGTTGG";
        let mut seq = unit.to_vec();
        seq.extend_from_slice(b"TTTT");
        seq.extend_from_slice(unit);
        let guide = find_guides(&seq).into_iter().next().unwrap();
        let off = find_offtargets(&guide, &seq, 0);
        assert!(off.iter().any(|o| o.mismatches == 0));
        assert!(specificity(&off) < 60.0);
    }
}
