//! Low-level DNA sequence helpers: complement, reverse-complement, GC content,
//! and IUPAC pattern matching (for restriction sites with ambiguous bases).

/// IUPAC-aware complement of a single base (case-insensitive, returns uppercase).
pub fn complement_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        other => other,
    }
}

/// Reverse complement of a sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// GC fraction in `[0, 1]` (counts G, C, and the IUPAC S = G/C).
pub fn gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C' | b'S'))
        .count();
    gc as f64 / seq.len() as f64
}

/// Does an IUPAC `pattern` base match a concrete `base` (case-insensitive)?
pub fn iupac_matches(pattern: u8, base: u8) -> bool {
    let base = base.to_ascii_uppercase();
    let set: &[u8] = match pattern.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"GC",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => return false,
    };
    set.contains(&base)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rev_comp() {
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
        assert_eq!(reverse_complement(b"gaattc"), b"GAATTC"); // palindrome, upper-cased
    }

    #[test]
    fn gc() {
        assert!((gc_content(b"GGCC") - 1.0).abs() < 1e-9);
        assert!((gc_content(b"ATAT") - 0.0).abs() < 1e-9);
        assert!((gc_content(b"ATGC") - 0.5).abs() < 1e-9);
    }

    #[test]
    fn iupac() {
        assert!(iupac_matches(b'N', b'a'));
        assert!(iupac_matches(b'R', b'G'));
        assert!(!iupac_matches(b'R', b'C'));
        assert!(iupac_matches(b'Y', b'T'));
    }
}
