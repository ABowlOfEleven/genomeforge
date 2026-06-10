//! Guide scoring.
//!
//! **On-target** here is a transparent, documented *heuristic* efficiency
//! estimate — not a trained model. For production design use a model like
//! Azimuth / Doench Rule Set 2. **Off-target** uses a CFD-style, position-
//! weighted mismatch penalty (PAM-proximal "seed" mismatches hurt most),
//! following the spirit of Doench et al. 2016.

use gx_plasmid::gc_content;

/// Heuristic on-target efficiency estimate in `[0, 1]` for a 20 nt protospacer.
///
/// Combines: optimal GC (~50%), penalties for extreme GC and for a poly-T
/// (U6 terminator) run, and a small bonus for a PAM-proximal G.
pub fn on_target_estimate(protospacer: &[u8]) -> f64 {
    if protospacer.is_empty() {
        return 0.0;
    }
    let gc = gc_content(protospacer);
    let mut s = 0.65;
    s -= (gc - 0.5).abs() * 0.7;
    if !(0.2..=0.85).contains(&gc) {
        s -= 0.2;
    }
    if has_poly_t(protospacer) {
        s -= 0.3;
    }
    if protospacer
        .last()
        .map(|b| b.eq_ignore_ascii_case(&b'G'))
        .unwrap_or(false)
    {
        s += 0.05;
    }
    s.clamp(0.0, 1.0)
}

/// True if the sequence contains a run of ≥4 T's (RNA-Pol-III terminator).
pub fn has_poly_t(seq: &[u8]) -> bool {
    seq.windows(4).any(|w| w.eq_ignore_ascii_case(b"TTTT"))
}

/// CFD-style position weight by distance from the PAM (1 = PAM-adjacent seed,
/// 20 = PAM-distal). Seed mismatches drive the score down hardest.
fn position_weight(pos_from_pam: usize) -> f64 {
    match pos_from_pam {
        1..=4 => 0.10,
        5..=8 => 0.30,
        9..=12 => 0.55,
        13..=16 => 0.75,
        _ => 0.90,
    }
}

/// CFD-style off-target score in `[0, 1]` for a candidate vs the guide, both
/// 20 nt 5'→3' with the PAM at the 3' end (higher = more likely to be cut).
pub fn cfd_score(guide: &[u8], candidate: &[u8]) -> f64 {
    let len = guide.len().min(candidate.len());
    let mut score = 1.0;
    for i in 0..len {
        if !guide[i].eq_ignore_ascii_case(&candidate[i]) {
            let pos_from_pam = len - i; // i = len-1 is PAM-adjacent
            score *= position_weight(pos_from_pam);
        }
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn on_target_prefers_balanced_gc_no_polyt() {
        let balanced = b"ACGTACGTACGTACGTACGT"; // 50% GC
        let polyt = b"ACGTTTTACGTACGTACGTA";
        let at_rich = b"ATATATATATATATATATAT";
        assert!(on_target_estimate(balanced) > on_target_estimate(polyt));
        assert!(on_target_estimate(balanced) > on_target_estimate(at_rich));
        assert!((0.0..=1.0).contains(&on_target_estimate(balanced)));
    }

    #[test]
    fn cfd_perfect_is_one_seed_mismatch_worst() {
        let g = b"ACGTACGTACGTACGTACGT";
        assert!((cfd_score(g, g) - 1.0).abs() < 1e-9);
        // mismatch at PAM-adjacent position (last) should hurt more than at the 5' end.
        let mut seed = g.to_vec();
        *seed.last_mut().unwrap() = b'A'; // was T
        let mut distal = g.to_vec();
        distal[0] = b'T'; // was A
        assert!(cfd_score(g, &seed) < cfd_score(g, &distal));
    }
}
