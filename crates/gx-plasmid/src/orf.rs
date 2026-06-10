//! Translation (standard genetic code) and open-reading-frame finding.

use gx_core::Strand;

use crate::seq::reverse_complement;

/// Translate one codon (3 bases) to a single-letter amino acid; `*` = stop,
/// `X` = unknown/ambiguous.
pub fn codon_to_aa(c: &[u8]) -> char {
    if c.len() < 3 {
        return 'X';
    }
    let norm = |b: u8| {
        let u = b.to_ascii_uppercase();
        if u == b'U' { b'T' } else { u }
    };
    let k = [norm(c[0]), norm(c[1]), norm(c[2])];
    match &k {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}

/// Translate a DNA sequence to protein (standard code); trailing partial codon
/// is ignored. Stop codons render as `*`.
pub fn translate(dna: &[u8]) -> String {
    dna.chunks_exact(3).map(codon_to_aa).collect()
}

#[derive(Debug, Clone)]
pub struct Orf {
    /// 0-based start on the top-strand coordinate (inclusive).
    pub start: usize,
    /// 0-based end on the top-strand coordinate (exclusive; includes the stop).
    pub end: usize,
    pub strand: Strand,
    pub frame: u8,
    /// Protein sequence (without the trailing stop).
    pub protein: String,
}

impl Orf {
    pub fn aa_len(&self) -> usize {
        self.protein.len()
    }
}

/// Find ORFs (ATG…stop) of at least `min_aa` amino acids in all six frames.
/// Coordinates are reported on the top strand.
pub fn find_orfs(seq: &[u8], min_aa: usize) -> Vec<Orf> {
    let mut orfs = Vec::new();
    scan_strand(seq, Strand::Forward, min_aa, seq.len(), &mut orfs);
    let rc = reverse_complement(seq);
    scan_strand(&rc, Strand::Reverse, min_aa, seq.len(), &mut orfs);
    orfs.sort_by(|a, b| a.start.cmp(&b.start).then(b.aa_len().cmp(&a.aa_len())));
    orfs
}

fn scan_strand(seq: &[u8], strand: Strand, min_aa: usize, orig_len: usize, out: &mut Vec<Orf>) {
    let n = seq.len();
    for frame in 0..3usize {
        let mut i = frame;
        while i + 3 <= n {
            let is_atg = seq[i].eq_ignore_ascii_case(&b'A')
                && seq[i + 1].eq_ignore_ascii_case(&b'T')
                && seq[i + 2].eq_ignore_ascii_case(&b'G');
            if is_atg {
                let mut j = i;
                let mut protein = String::new();
                let mut stop = false;
                while j + 3 <= n {
                    let aa = codon_to_aa(&seq[j..j + 3]);
                    j += 3;
                    if aa == '*' {
                        stop = true;
                        break;
                    }
                    protein.push(aa);
                }
                if stop && protein.len() >= min_aa {
                    // Map the [i, j) span on this strand back to top-strand coords.
                    let (start, end) = match strand {
                        Strand::Reverse => (orig_len - j, orig_len - i),
                        _ => (i, j),
                    };
                    out.push(Orf {
                        start,
                        end,
                        strand,
                        frame: frame as u8,
                        protein,
                    });
                    i = j; // non-overlapping within a frame
                    continue;
                }
            }
            i += 3;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translation() {
        // ATG GCC TAA -> M A *
        assert_eq!(translate(b"ATGGCCTAA"), "MA*");
        assert_eq!(codon_to_aa(b"atg"), 'M'); // lower-case handled
        assert_eq!(codon_to_aa(b"TGA"), '*');
    }

    #[test]
    fn finds_forward_orf() {
        // ATG (M) AAA(K) ... TAA, 6 codons before stop.
        let seq = b"CCCATGAAAGCCGTTTTGTAAGGG";
        let orfs = find_orfs(seq, 3);
        let fwd: Vec<_> = orfs.iter().filter(|o| o.strand == Strand::Forward).collect();
        assert!(!fwd.is_empty());
        assert_eq!(fwd[0].start, 3);
        assert!(fwd[0].protein.starts_with('M'));
        assert!(!fwd[0].protein.contains('*'));
    }

    #[test]
    fn finds_reverse_orf() {
        // Reverse complement of an ATG..stop should be discovered on the - strand.
        let coding = b"ATGAAAGCCGTTTTGTAA";
        let rc = reverse_complement(coding);
        let orfs = find_orfs(&rc, 3);
        assert!(orfs.iter().any(|o| o.strand == Strand::Reverse));
    }
}
