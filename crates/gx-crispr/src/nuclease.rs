//! Nuclease / PAM definitions and generalized guide enumeration.
//!
//! The original [`crate::guides::find_guides`] is hardcoded to SpCas9 (a 3'-side
//! `NGG` PAM with a 20 nt protospacer). This module generalizes that to a small
//! panel of commonly used nucleases — engineered SpCas9 variants with relaxed
//! PAMs and the 5'-PAM Cas12a — via [`find_guides_for`].
//!
//! The PAM strings are the **canonical/approximate** recognition motifs in IUPAC
//! code, not empirical activity tables; SpRY in particular is near-PAMless and
//! `NRN` here is a deliberately simplified stand-in. These are decision-support
//! estimates, not a clinical pipeline. Only SpCas9 carries the Doench 2014
//! on-target model (it was trained on SpCas9); other nucleases fall back to the
//! documented heuristic in [`crate::score::on_target_estimate`].

use gx_core::Strand;
use gx_plasmid::reverse_complement;

use crate::guides::{Guide, make_guide_for, on_target_spcas9};
use crate::score::on_target_estimate;

/// A supported CRISPR nuclease, identified by its PAM and protospacer geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nuclease {
    /// *Streptococcus pyogenes* Cas9 — `NGG` PAM, 20 nt protospacer (the original).
    SpCas9,
    /// SpCas9-NG — engineered relaxed `NG` PAM, 20 nt protospacer.
    SpCas9NG,
    /// SpRY — near-PAMless SpCas9 variant, approximated here as `NRN`, 20 nt.
    SpRY,
    /// Cas12a (Cpf1) — 5'-side `TTTV` PAM, 23 nt protospacer.
    Cas12a,
}

impl Nuclease {
    /// The recognition PAM as an IUPAC string (5'→3').
    pub fn pam(self) -> &'static str {
        match self {
            Nuclease::SpCas9 => "NGG",
            Nuclease::SpCas9NG => "NG",
            Nuclease::SpRY => "NRN",
            Nuclease::Cas12a => "TTTV",
        }
    }

    /// `true` if the PAM sits 5' of the protospacer (Cas12a); `false` for the
    /// 3'-PAM Cas9 family.
    pub fn pam_is_5prime(self) -> bool {
        matches!(self, Nuclease::Cas12a)
    }

    /// Protospacer length in nucleotides: 20 for the Cas9 variants, 23 for Cas12a.
    pub fn protospacer_len(self) -> usize {
        match self {
            Nuclease::Cas12a => 23,
            _ => 20,
        }
    }

    /// Human-readable label, e.g. `"SpCas9 (NGG)"`.
    pub fn label(self) -> &'static str {
        match self {
            Nuclease::SpCas9 => "SpCas9 (NGG)",
            Nuclease::SpCas9NG => "SpCas9-NG (NG)",
            Nuclease::SpRY => "SpRY (NRN)",
            Nuclease::Cas12a => "Cas12a (TTTV)",
        }
    }

    /// All supported nucleases.
    pub fn all() -> [Nuclease; 4] {
        [Nuclease::SpCas9, Nuclease::SpCas9NG, Nuclease::SpRY, Nuclease::Cas12a]
    }
}

/// IUPAC nucleotide match: does ambiguity `code` accept concrete `base`?
///
/// Supports the literals `A/C/G/T` plus the degenerate codes used by the PAMs in
/// this module (`N`, `R`, `Y`, `V`). Both arguments are matched case-insensitively.
pub fn iupac_matches(code: u8, base: u8) -> bool {
    let b = base.to_ascii_uppercase();
    match code.to_ascii_uppercase() {
        b'N' => matches!(b, b'A' | b'C' | b'G' | b'T'),
        b'R' => matches!(b, b'A' | b'G'),
        b'Y' => matches!(b, b'C' | b'T'),
        b'V' => matches!(b, b'A' | b'C' | b'G'),
        c @ (b'A' | b'C' | b'G' | b'T') => b == c,
        _ => false,
    }
}

/// True if every base of `window` satisfies the corresponding IUPAC `pam` code.
fn pam_matches(pam: &[u8], window: &[u8]) -> bool {
    pam.len() == window.len()
        && pam.iter().zip(window).all(|(&c, &b)| iupac_matches(c, b))
}

fn is_acgt(s: &[u8]) -> bool {
    s.iter().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

/// Enumerate all valid guide sites for `nuclease` on both strands of `seq`,
/// sorted by on-target score (descending).
///
/// For 3'-PAM nucleases (the Cas9 family) the protospacer of length `L` lies
/// immediately 5' of the PAM. For the 5'-PAM Cas12a the PAM lies immediately 5'
/// of the protospacer. The reverse strand is searched by matching the
/// reverse-complemented PAM motif on the top strand. Only SpCas9 receives the
/// Doench 2014 on-target score; all others use the documented heuristic.
pub fn find_guides_for(seq: &[u8], nuclease: Nuclease) -> Vec<Guide> {
    let up = seq.to_ascii_uppercase();
    let mut guides = Vec::new();
    if nuclease.pam_is_5prime() {
        find_5prime(&up, nuclease, &mut guides);
    } else {
        find_3prime(&up, nuclease, &mut guides);
    }
    guides.sort_by(|a, b| {
        b.on_score
            .partial_cmp(&a.on_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    guides
}

/// On-target score for `nuclease`: Doench 2014 (from the 30-mer context) for
/// SpCas9, the documented heuristic otherwise.
fn on_target(up: &[u8], start: usize, strand: Strand, proto: &[u8], nuclease: Nuclease) -> f64 {
    if nuclease == Nuclease::SpCas9 {
        on_target_spcas9(up, start, strand, proto)
    } else {
        on_target_estimate(proto)
    }
}

/// 3'-PAM layout: `5'-[protospacer L]-[PAM]-3'`.
fn find_3prime(up: &[u8], nuclease: Nuclease, out: &mut Vec<Guide>) {
    let l = nuclease.protospacer_len();
    let pam = nuclease.pam().as_bytes();
    let plen = pam.len();
    let n = up.len();
    if n < l + plen {
        return;
    }
    let rc_pam = reverse_complement(pam); // top-strand motif for a reverse-strand guide

    // Forward strand: protospacer at [p-l, p), PAM at [p, p+plen).
    for p in l..=(n - plen) {
        let pam_win = &up[p..p + plen];
        if !pam_matches(pam, pam_win) {
            continue;
        }
        let proto = &up[p - l..p];
        if !is_acgt(proto) {
            continue;
        }
        let cut = (p - 3).min(n);
        push_guide(out, up, p - l, Strand::Forward, proto, pam_win, cut, nuclease);
    }

    // Reverse strand: the PAM's reverse complement appears on the top strand at
    // [p, p+plen); the protospacer is the reverse complement of the L bases 3' of it.
    for p in 0..=(n - (l + plen)) {
        let pam_win = &up[p..p + plen];
        if !pam_matches(&rc_pam, pam_win) {
            continue;
        }
        let region = &up[p + plen..p + plen + l];
        if !is_acgt(region) {
            continue;
        }
        let proto = reverse_complement(region);
        let pam_53 = reverse_complement(pam_win); // back to 5'→3' on the guide strand
        let cut = (p + plen + 3).min(up.len());
        push_guide(out, up, p + plen, Strand::Reverse, &proto, &pam_53, cut, nuclease);
    }
}

/// 5'-PAM layout (Cas12a): `5'-[PAM]-[protospacer L]-3'`.
fn find_5prime(up: &[u8], nuclease: Nuclease, out: &mut Vec<Guide>) {
    let l = nuclease.protospacer_len();
    let pam = nuclease.pam().as_bytes();
    let plen = pam.len();
    let n = up.len();
    if n < l + plen {
        return;
    }
    let rc_pam = reverse_complement(pam);

    // Forward strand: PAM at [p, p+plen), protospacer at [p+plen, p+plen+l).
    for p in 0..=(n - (l + plen)) {
        let pam_win = &up[p..p + plen];
        if !pam_matches(pam, pam_win) {
            continue;
        }
        let proto = &up[p + plen..p + plen + l];
        if !is_acgt(proto) {
            continue;
        }
        // Cas12a makes a staggered cut distal to the PAM; we report a nominal
        // cut near the 3' end of the protospacer (no blunt-DSB claim).
        let cut = (p + plen + l).min(n);
        push_guide(out, up, p + plen, Strand::Forward, proto, pam_win, cut, nuclease);
    }

    // Reverse strand: the reverse-complemented PAM appears on the top strand 3'
    // of the protospacer region. Top-strand PAM window at [p, p+plen); the
    // protospacer (guide 5'→3') is the reverse complement of the L bases 5' of it.
    for p in l..=(n - plen) {
        let pam_win = &up[p..p + plen];
        if !pam_matches(&rc_pam, pam_win) {
            continue;
        }
        let region = &up[p - l..p];
        if !is_acgt(region) {
            continue;
        }
        let proto = reverse_complement(region);
        let pam_53 = reverse_complement(pam_win);
        let cut = (p - l).min(up.len());
        push_guide(out, up, p - l, Strand::Reverse, &proto, &pam_53, cut, nuclease);
    }
}

#[allow(clippy::too_many_arguments)]
fn push_guide(
    out: &mut Vec<Guide>,
    up: &[u8],
    start: usize,
    strand: Strand,
    proto: &[u8],
    pam: &[u8],
    cut: usize,
    nuclease: Nuclease,
) {
    let on = on_target(up, start, strand, proto, nuclease);
    out.push(make_guide_for(start, strand, proto, pam, cut, on, nuclease));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_matches_codes() {
        // N matches everything.
        for &b in b"ACGTacgt" {
            assert!(iupac_matches(b'N', b));
        }
        // R = A/G, Y = C/T, V = A/C/G.
        assert!(iupac_matches(b'R', b'A') && iupac_matches(b'R', b'G'));
        assert!(!iupac_matches(b'R', b'C') && !iupac_matches(b'R', b'T'));
        assert!(iupac_matches(b'Y', b'C') && iupac_matches(b'Y', b'T'));
        assert!(!iupac_matches(b'Y', b'A'));
        assert!(iupac_matches(b'V', b'A') && iupac_matches(b'V', b'C') && iupac_matches(b'V', b'G'));
        assert!(!iupac_matches(b'V', b'T'));
        // Literals, case-insensitive.
        assert!(iupac_matches(b'A', b'a') && iupac_matches(b'a', b'A'));
        assert!(!iupac_matches(b'A', b'C'));
        // Unknown code matches nothing.
        assert!(!iupac_matches(b'Z', b'A'));
    }

    #[test]
    fn spcas9_via_find_guides_for_matches_legacy() {
        let seq = b"ACGTACGTACGTACGTACGTTGG";
        let g = find_guides_for(seq, Nuclease::SpCas9);
        let fwd: Vec<_> = g.iter().filter(|g| g.strand == Strand::Forward).collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].protospacer, "ACGTACGTACGTACGTACGT");
        assert_eq!(fwd[0].pam, "TGG");
        assert_eq!(fwd[0].start, 0);
        assert_eq!(fwd[0].nuclease, Nuclease::SpCas9);
    }

    #[test]
    fn finds_cas12a_5prime_site() {
        // TTTA (TTTV) PAM followed by a 23 nt protospacer on the forward strand.
        let proto = b"ACGTACGTACGTACGTACGTACG"; // 23 nt
        let mut seq = b"TTTA".to_vec();
        seq.extend_from_slice(proto);
        let g = find_guides_for(&seq, Nuclease::Cas12a);
        let fwd: Vec<_> = g.iter().filter(|g| g.strand == Strand::Forward).collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].pam, "TTTA");
        assert_eq!(fwd[0].protospacer, "ACGTACGTACGTACGTACGTACG");
        assert_eq!(fwd[0].start, 4);
        assert_eq!(fwd[0].nuclease, Nuclease::Cas12a);
    }

    #[test]
    fn constructed_seq_has_both_spcas9_and_cas12a_sites() {
        // Region 1: SpCas9 forward site (20 nt + TGG).
        // Region 2: Cas12a forward site (TTTC + 23 nt).
        let mut seq = b"ACGTACGTACGTACGTACGTTGG".to_vec(); // SpCas9 NGG at end
        seq.extend_from_slice(b"AAAA");
        seq.extend_from_slice(b"TTTC"); // Cas12a TTTV PAM
        seq.extend_from_slice(b"GCGCGCGCGCGCGCGCGCGCGCG"); // 23 nt protospacer

        let cas9 = find_guides_for(&seq, Nuclease::SpCas9);
        assert!(cas9.iter().any(|g| g.strand == Strand::Forward && g.pam == "TGG"));

        let cas12a = find_guides_for(&seq, Nuclease::Cas12a);
        let fwd: Vec<_> = cas12a.iter().filter(|g| g.strand == Strand::Forward && g.pam == "TTTC").collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].protospacer, "GCGCGCGCGCGCGCGCGCGCGCG");
    }

    #[test]
    fn spry_relaxed_pam_finds_more_than_spcas9() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let spcas9 = find_guides_for(seq, Nuclease::SpCas9).len();
        let spry = find_guides_for(seq, Nuclease::SpRY).len();
        // The near-PAMless SpRY (NRN) should find at least as many sites.
        assert!(spry >= spcas9);
    }

    #[test]
    fn nuclease_metadata_is_consistent() {
        assert_eq!(Nuclease::all().len(), 4);
        assert_eq!(Nuclease::SpCas9.protospacer_len(), 20);
        assert_eq!(Nuclease::Cas12a.protospacer_len(), 23);
        assert!(Nuclease::Cas12a.pam_is_5prime());
        assert!(!Nuclease::SpCas9.pam_is_5prime());
        assert_eq!(Nuclease::SpCas9.pam(), "NGG");
    }
}
