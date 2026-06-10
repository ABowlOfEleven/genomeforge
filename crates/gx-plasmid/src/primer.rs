//! Primer / oligo analysis: GC content and melting temperature.
//!
//! Three Tm estimates are provided because each suits a different regime:
//! the Wallace rule (very short oligos), a salt-adjusted basic formula, and the
//! SantaLucia (1998) nearest-neighbour model (the accurate one for design).

use crate::seq::gc_content;

#[derive(Debug, Clone)]
pub struct PrimerStats {
    pub length: usize,
    pub gc_percent: f64,
    pub tm_wallace: f64,
    pub tm_basic: f64,
    pub tm_nn: f64,
}

pub fn analyze(seq: &[u8]) -> PrimerStats {
    PrimerStats {
        length: seq.len(),
        gc_percent: gc_content(seq) * 100.0,
        tm_wallace: tm_wallace(seq),
        tm_basic: tm_basic(seq),
        tm_nn: tm_nearest_neighbor(seq),
    }
}

fn counts(seq: &[u8]) -> (usize, usize) {
    let mut at = 0;
    let mut gc = 0;
    for b in seq {
        match b.to_ascii_uppercase() {
            b'A' | b'T' | b'U' => at += 1,
            b'G' | b'C' => gc += 1,
            _ => {}
        }
    }
    (at, gc)
}

/// Wallace rule: 2(A+T) + 4(G+C). Good only for short oligos (< ~14 nt).
pub fn tm_wallace(seq: &[u8]) -> f64 {
    let (at, gc) = counts(seq);
    2.0 * at as f64 + 4.0 * gc as f64
}

/// Salt-adjusted basic formula: 64.9 + 41·(GC − 16.4)/N.
pub fn tm_basic(seq: &[u8]) -> f64 {
    let n = seq.len();
    if n == 0 {
        return 0.0;
    }
    let (_, gc) = counts(seq);
    64.9 + 41.0 * (gc as f64 - 16.4) / n as f64
}

/// SantaLucia (1998) nearest-neighbour Tm in °C, at 50 mM Na⁺ and 0.25 µM oligo.
/// Returns `NaN` if the sequence contains non-ACGT bases.
pub fn tm_nearest_neighbor(seq: &[u8]) -> f64 {
    let s: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    if s.len() < 2 {
        return f64::NAN;
    }
    let mut dh = 0.0f64;
    let mut ds = 0.0f64;
    for w in s.windows(2) {
        match nn_params(w[0], w[1]) {
            Some((h, sv)) => {
                dh += h;
                ds += sv;
            }
            None => return f64::NAN,
        }
    }
    let (h0, s0) = init_params(s[0]);
    let (hn, sn) = init_params(*s.last().unwrap());
    dh += h0 + hn;
    ds += s0 + sn;

    // Entropic salt correction (SantaLucia 1998), Na⁺ = 50 mM.
    let na = 0.05_f64;
    ds += 0.368 * (s.len() as f64 - 1.0) * na.ln();

    let r = 1.987_f64; // cal/(mol·K)
    let ct = 0.25e-6_f64; // total strand conc.
    let tm_k = (dh * 1000.0) / (ds + r * (ct / 4.0).ln());
    tm_k - 273.15
}

/// ΔH (kcal/mol), ΔS (cal/mol·K) for a 5'→3' nearest-neighbour dinucleotide.
fn nn_params(a: u8, b: u8) -> Option<(f64, f64)> {
    let pair = [a, b];
    let v = match &pair {
        b"AA" => (-7.9, -22.2),
        b"AC" => (-8.4, -22.4),
        b"AG" => (-7.8, -21.0),
        b"AT" => (-7.2, -20.4),
        b"CA" => (-8.5, -22.7),
        b"CC" => (-8.0, -19.9),
        b"CG" => (-10.6, -27.2),
        b"CT" => (-7.8, -21.0),
        b"GA" => (-8.2, -22.2),
        b"GC" => (-9.8, -24.4),
        b"GG" => (-8.0, -19.9),
        b"GT" => (-8.4, -22.4),
        b"TA" => (-7.2, -21.3),
        b"TC" => (-8.2, -22.2),
        b"TG" => (-8.5, -22.7),
        b"TT" => (-7.9, -22.2),
        _ => return None,
    };
    Some(v)
}

/// Initiation parameters depending on a terminal base being G·C or A·T.
fn init_params(end: u8) -> (f64, f64) {
    if end == b'G' || end == b'C' {
        (0.1, -2.8)
    } else {
        (2.3, 4.1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wallace_and_basic() {
        // 4 G/C + 4 A/T -> Wallace 2*4 + 4*4 = 24
        assert!((tm_wallace(b"ATGCATGC") - 24.0).abs() < 1e-9);
        assert!(tm_basic(b"ATGCATGCATGCATGCATGC") > 0.0);
    }

    #[test]
    fn nn_in_plausible_range() {
        // M13 -21 primer-ish 17-mer; NN Tm should land in a sane primer range.
        let tm = tm_nearest_neighbor(b"GTAAAACGACGGCCAGT");
        assert!(tm.is_finite());
        assert!((40.0..70.0).contains(&tm), "tm was {tm}");
    }

    #[test]
    fn nn_rejects_non_acgt() {
        assert!(tm_nearest_neighbor(b"ATGCN").is_nan());
    }
}
