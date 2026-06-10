//! SpCas9 guide (sgRNA) enumeration: find protospacers adjacent to an NGG PAM
//! on both strands and score them.

use gx_core::Strand;
use gx_plasmid::{gc_content, reverse_complement};

use crate::score::{has_poly_t, on_target_estimate};

pub const PROTOSPACER_LEN: usize = 20;

#[derive(Debug, Clone)]
pub struct Guide {
    /// 0-based top-strand start of the 20 nt protospacer footprint.
    pub start: usize,
    pub strand: Strand,
    /// 20 nt spacer, 5'→3' (the guide/target-strand sequence).
    pub protospacer: String,
    /// 3 nt PAM, 5'→3'.
    pub pam: String,
    pub gc: f64,
    /// Heuristic on-target efficiency estimate in `[0, 1]`.
    pub on_score: f64,
    pub poly_t: bool,
    /// 0-based top-strand position of the predicted blunt DSB (3 bp 5' of PAM).
    pub cut_site: usize,
}

fn is_acgt(s: &[u8]) -> bool {
    s.iter().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

fn make_guide(start: usize, strand: Strand, proto: &[u8], pam: &[u8], cut: usize, on_score: f64) -> Guide {
    Guide {
        start,
        strand,
        protospacer: String::from_utf8_lossy(proto).into_owned(),
        pam: String::from_utf8_lossy(pam).into_owned(),
        gc: gc_content(proto),
        on_score,
        poly_t: has_poly_t(proto),
        cut_site: cut,
    }
}

/// On-target score: the exact Doench 2014 model computed from the 30-mer
/// context (4 nt + 20 nt + 3 nt PAM + 3 nt) when that context is in-bounds,
/// otherwise the documented heuristic (only happens near sequence ends).
fn on_target(up: &[u8], start: usize, strand: Strand, proto: &[u8]) -> f64 {
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
pub fn find_guides(seq: &[u8]) -> Vec<Guide> {
    let up = seq.to_ascii_uppercase();
    let n = up.len();
    let mut guides = Vec::new();
    if n < PROTOSPACER_LEN + 3 {
        return guides;
    }

    // Forward strand: 5'-[20 nt protospacer]-NGG-3'. PAM occupies [p, p+2].
    for p in PROTOSPACER_LEN..=(n - 3) {
        if up[p + 1] == b'G' && up[p + 2] == b'G' {
            let proto = &up[p - PROTOSPACER_LEN..p];
            if !is_acgt(proto) {
                continue;
            }
            make_into(&mut guides, &up, p - PROTOSPACER_LEN, Strand::Forward, proto, &up[p..p + 3], p - 3);
        }
    }

    // Reverse strand: PAM reads CCN on the top strand at [p, p+2]; the
    // protospacer is the reverse complement of the 20 nt 3' of it.
    for p in 0..=(n - (PROTOSPACER_LEN + 3)) {
        if up[p] == b'C' && up[p + 1] == b'C' {
            let region = &up[p + 3..p + 3 + PROTOSPACER_LEN];
            if !is_acgt(region) {
                continue;
            }
            let proto = reverse_complement(region);
            let pam = reverse_complement(&up[p..p + 3]); // CCN -> NGG (5'->3')
            make_into(&mut guides, &up, p + 3, Strand::Reverse, &proto, &pam, p + 6);
        }
    }

    guides.sort_by(|a, b| {
        b.on_score
            .partial_cmp(&a.on_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    guides
}

#[allow(clippy::too_many_arguments)]
fn make_into(out: &mut Vec<Guide>, up: &[u8], start: usize, strand: Strand, proto: &[u8], pam: &[u8], cut: usize) {
    let on = on_target(up, start, strand, proto);
    out.push(make_guide(start, strand, proto, pam, cut, on));
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
