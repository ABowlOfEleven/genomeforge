//! DNA assembly simulation — Golden-Gate / restriction-ligation and Gibson —
//! built on the digest/overhang model in [`crate::cloning`].
//!
//! The [`Fragment`] produced by a digest carries coordinates and typed ends but
//! not the nucleotides, so assembly works on a sequence-carrying [`Part`]
//! instead. [`fragments_to_parts`] derives parts from a digest of a known
//! sequence; [`golden_gate`] and [`gibson`] then search for a circular product.
//!
//! ## Simplifying assumptions
//!
//! * **Overhang compatibility is by type+length** (via
//!   [`ends_compatible`]): two non-blunt ends ligate when they share the same
//!   overhang class and length. Sequence identity of the overhang bases is
//!   *assumed*, because [`Overhang`] stores only a length, not the bases.
//! * The search is **bounded to a small number of parts** (`MAX_PARTS`, below)
//!   and explores orderings/orientations by backtracking, which is adequate for
//!   the hand-sized constructs a designer assembles interactively.

use crate::cloning::{Fragment, FragmentEnd, Overhang, ends_compatible};
use crate::seq::reverse_complement;

/// Upper bound on the number of parts an assembly search will consider. Beyond
/// this the permutation search is refused (returns `None`) to avoid factorial
/// blow-up; real Golden-Gate / Gibson reactions are well under this.
pub const MAX_PARTS: usize = 10;

/// A linear DNA piece with typed ends — the sequence-carrying counterpart of a
/// [`Fragment`]. `seq` is the top-strand nucleotides spanning the fragment,
/// inclusive of any single-stranded overhang at either end.
#[derive(Debug, Clone)]
pub struct Part {
    pub name: String,
    pub seq: Vec<u8>,
    pub left: FragmentEnd,
    pub right: FragmentEnd,
}

impl Part {
    /// This part in the reverse-complement orientation: the sequence is
    /// reverse-complemented and the two ends swap (left↔right), since flipping
    /// the molecule exchanges which physical end is which.
    fn reverse_complement(&self) -> Part {
        Part {
            name: format!("{}(rc)", self.name),
            seq: reverse_complement(&self.seq),
            left: self.right.clone(),
            right: self.left.clone(),
        }
    }
}

/// The result of an assembly: a circular product plus the order in which the
/// input parts were joined and a human-readable description of each junction.
#[derive(Debug, Clone)]
pub struct Assembly {
    pub circular: bool,
    pub seq: Vec<u8>,
    /// Indices into the input `parts` slice, in assembled order.
    pub order: Vec<usize>,
    /// One label per junction (same length as `order` for a circular product).
    pub junctions: Vec<String>,
}

/// Slice a digest into sequence-carrying [`Part`]s. Each fragment's nucleotides
/// are copied from `seq` per its `start`/`end` coordinates; the circular wrap
/// fragment (whose `start > end`) is stitched across the origin. Overhang ends
/// are carried over verbatim from each [`Fragment`].
pub fn fragments_to_parts(seq: &[u8], circular: bool, frags: &[Fragment]) -> Vec<Part> {
    let n = seq.len();
    frags
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let bytes = if circular && f.start > f.end {
                // Wrap fragment: tail of the sequence followed by the head.
                let mut v = seq[f.start..n].to_vec();
                v.extend_from_slice(&seq[0..f.end]);
                v
            } else if circular && f.start == f.end && f.len == n {
                // Single cutter linearises the whole molecule: rotate to the cut.
                let mut v = seq[f.start..n].to_vec();
                v.extend_from_slice(&seq[0..f.start]);
                v
            } else {
                seq[f.start..f.end].to_vec()
            };
            Part {
                name: format!("frag{i}"),
                seq: bytes,
                left: f.left.clone(),
                right: f.right.clone(),
            }
        })
        .collect()
}

/// Length of the single-stranded overhang an end contributes to a junction.
fn overhang_len(o: Overhang) -> usize {
    match o {
        Overhang::Blunt => 0,
        Overhang::FivePrime(n) | Overhang::ThreePrime(n) => n as usize,
    }
}

/// Golden-Gate / restriction-ligation assembly.
///
/// Searches for an ordering of **all** `parts` — each free to take its
/// reverse-complement orientation — such that consecutive parts have
/// ligation-compatible ends ([`ends_compatible`]) and the chain closes into a
/// loop (last part's right end compatible with first part's left end). The
/// returned [`Assembly`] is circular.
///
/// Sequences are concatenated with **shared overhang bases counted once**: at a
/// junction where both ends carry a length-`n` overhang, those `n` bases live in
/// the upstream part's slice, so the downstream part contributes its sequence
/// minus its leading `n` overhang bases.
///
/// Compatibility is by overhang type+length only (see the module docs); blunt
/// ends are *not* used to close a junction here, since a blunt circular ligation
/// has no defined orientation and any blunt end would match any other.
///
/// Returns `None` if no valid circular ordering exists, or if `parts` is empty
/// or exceeds [`MAX_PARTS`].
pub fn golden_gate(parts: &[Part]) -> Option<Assembly> {
    let k = parts.len();
    if k == 0 || k > MAX_PARTS {
        return None;
    }

    // Each part may be used forward or reverse-complemented.
    let oriented: Vec<[Part; 2]> = parts
        .iter()
        .map(|p| [p.clone(), p.reverse_complement()])
        .collect();

    let mut used = vec![false; k];
    let mut order: Vec<usize> = Vec::with_capacity(k);
    let mut flips: Vec<bool> = Vec::with_capacity(k);

    // Anchor the first part forward (index 0) to break rotational/reflective
    // symmetry of the loop, then backtrack over the rest.
    used[0] = true;
    order.push(0);
    flips.push(false);

    // Bound the backtracking so a pathological input (many mutually compatible
    // overhangs that never close into a loop) can't freeze a caller that runs
    // this in a paint loop. Real assemblies close well within the budget.
    let mut budget: u32 = 200_000;
    if let Some(chain) = gg_extend(&oriented, &mut used, &mut order, &mut flips, &mut budget) {
        return Some(build_gg_assembly(&oriented, &chain.0, &chain.1));
    }
    None
}

/// Backtracking extension for [`golden_gate`]. Returns `(order, flips)` of a
/// complete, closing chain if one exists. `budget` caps the nodes explored.
fn gg_extend(
    oriented: &[[Part; 2]],
    used: &mut [bool],
    order: &mut Vec<usize>,
    flips: &mut Vec<bool>,
    budget: &mut u32,
) -> Option<(Vec<usize>, Vec<bool>)> {
    if *budget == 0 {
        return None;
    }
    *budget -= 1;
    let k = oriented.len();
    if order.len() == k {
        // Close the loop: last part's right end ↔ first part's left end.
        let last = &oriented[*order.last().unwrap()][*flips.last().unwrap() as usize];
        let first = &oriented[order[0]][flips[0] as usize];
        if junction_ligates(&last.right, &first.left) {
            return Some((order.clone(), flips.clone()));
        }
        return None;
    }

    let prev = &oriented[*order.last().unwrap()][*flips.last().unwrap() as usize];
    for next in 0..k {
        if used[next] {
            continue;
        }
        for &flip in &[false, true] {
            let cand = &oriented[next][flip as usize];
            if junction_ligates(&prev.right, &cand.left) {
                used[next] = true;
                order.push(next);
                flips.push(flip);
                if let Some(done) = gg_extend(oriented, used, order, flips, budget) {
                    return Some(done);
                }
                order.pop();
                flips.pop();
                used[next] = false;
            }
        }
    }
    None
}

/// A junction ligates if the two ends are compatible and non-blunt (see the
/// note in [`golden_gate`] on why blunt junctions are excluded).
fn junction_ligates(a: &FragmentEnd, b: &FragmentEnd) -> bool {
    a.overhang != Overhang::Blunt && ends_compatible(a, b)
}

/// Stitch the chosen chain into a circular sequence, dropping each downstream
/// part's leading overhang bases (shared with the upstream junction).
fn build_gg_assembly(oriented: &[[Part; 2]], order: &[usize], flips: &[bool]) -> Assembly {
    let chain: Vec<&Part> = order
        .iter()
        .zip(flips)
        .map(|(&i, &f)| &oriented[i][f as usize])
        .collect();

    let mut seq = Vec::new();
    let mut junctions = Vec::with_capacity(chain.len());
    for (pos, part) in chain.iter().enumerate() {
        // The base each part overlaps with its *predecessor* is the shared
        // overhang of the incoming junction, which the predecessor already
        // emitted — so trim it from this part's leading bases.
        let shared = overhang_len(part.left.overhang).min(part.seq.len());
        seq.extend_from_slice(&part.seq[shared..]);

        // Junction from this part's right end into the next part (wrapping).
        let next = &chain[(pos + 1) % chain.len()];
        junctions.push(format!(
            "{} {} | {} {}",
            part.name,
            part.right.overhang.label(),
            next.left.overhang.label(),
            next.name,
        ));
    }

    Assembly { circular: true, seq, order: order.to_vec(), junctions }
}

/// Gibson assembly.
///
/// Joins linear `parts` (name + top-strand sequence) where the 3' end of one
/// shares an **exact** overlap of at least `min_overlap` bp with the 5' start of
/// another (both orientations are tried via [`reverse_complement`]). The chain
/// is closed into a circle (the last part overlaps back into the first) and each
/// overlapping region is merged once.
///
/// Returns `None` if the parts cannot be chained into a closed circle, or if
/// `parts` is empty/exceeds [`MAX_PARTS`], or if `min_overlap` is zero.
pub fn gibson(parts: &[(String, Vec<u8>)], min_overlap: usize) -> Option<Assembly> {
    let k = parts.len();
    if k == 0 || k > MAX_PARTS || min_overlap == 0 {
        return None;
    }

    // Both orientations of each part.
    let oriented: Vec<[(String, Vec<u8>); 2]> = parts
        .iter()
        .map(|(name, s)| {
            [
                (name.clone(), s.clone()),
                (format!("{name}(rc)"), reverse_complement(s)),
            ]
        })
        .collect();

    let mut used = vec![false; k];
    let mut order: Vec<usize> = vec![0];
    let mut flips: Vec<bool> = vec![false];
    let mut overlaps: Vec<usize> = Vec::new(); // overlap into each successor
    used[0] = true;

    if let Some((order, flips, overlaps)) =
        gibson_extend(&oriented, min_overlap, &mut used, &mut order, &mut flips, &mut overlaps)
    {
        return Some(build_gibson_assembly(&oriented, &order, &flips, &overlaps));
    }
    None
}

/// Backtracking extension for [`gibson`]. `overlaps[i]` is the overlap length
/// between chain position `i` and `i+1`; the final entry (added on close) is the
/// wrap overlap from the last part back to the first.
#[allow(clippy::too_many_arguments)]
fn gibson_extend(
    oriented: &[[(String, Vec<u8>); 2]],
    min_overlap: usize,
    used: &mut [bool],
    order: &mut Vec<usize>,
    flips: &mut Vec<bool>,
    overlaps: &mut Vec<usize>,
) -> Option<(Vec<usize>, Vec<bool>, Vec<usize>)> {
    let k = oriented.len();
    if order.len() == k {
        let last = &oriented[*order.last().unwrap()][*flips.last().unwrap() as usize].1;
        let first = &oriented[order[0]][flips[0] as usize].1;
        if let Some(ov) = overlap_len(last, first, min_overlap) {
            overlaps.push(ov);
            return Some((order.clone(), flips.clone(), overlaps.clone()));
        }
        return None;
    }

    let prev = oriented[*order.last().unwrap()][*flips.last().unwrap() as usize].1.clone();
    for next in 0..k {
        if used[next] {
            continue;
        }
        for &flip in &[false, true] {
            let cand = &oriented[next][flip as usize].1;
            if let Some(ov) = overlap_len(&prev, cand, min_overlap) {
                used[next] = true;
                order.push(next);
                flips.push(flip);
                overlaps.push(ov);
                if let Some(done) =
                    gibson_extend(oriented, min_overlap, used, order, flips, overlaps)
                {
                    return Some(done);
                }
                overlaps.pop();
                flips.pop();
                order.pop();
                used[next] = false;
            }
        }
    }
    None
}

/// Largest exact overlap (≥ `min`) where a suffix of `a` equals a prefix of `b`.
/// Returns `None` if no overlap of at least `min` exists.
fn overlap_len(a: &[u8], b: &[u8], min: usize) -> Option<usize> {
    let max = a.len().min(b.len());
    // Prefer the longest overlap to avoid spuriously short merges.
    for len in (min..=max).rev() {
        if a[a.len() - len..] == b[..len] {
            return Some(len);
        }
    }
    None
}

/// Concatenate the Gibson chain, merging each overlap once.
fn build_gibson_assembly(
    oriented: &[[(String, Vec<u8>); 2]],
    order: &[usize],
    flips: &[bool],
    overlaps: &[usize],
) -> Assembly {
    let chain: Vec<&(String, Vec<u8>)> = order
        .iter()
        .zip(flips)
        .map(|(&i, &f)| &oriented[i][f as usize])
        .collect();

    let mut seq = Vec::new();
    let mut junctions = Vec::with_capacity(chain.len());
    for (pos, part) in chain.iter().enumerate() {
        // Drop the leading bases shared with the predecessor's overlap (the wrap
        // overlap is the last entry, shared between the last and first parts).
        let incoming = if pos == 0 {
            *overlaps.last().unwrap()
        } else {
            overlaps[pos - 1]
        };
        seq.extend_from_slice(&part.1[incoming..]);

        let ov = overlaps[pos];
        let next = &chain[(pos + 1) % chain.len()];
        junctions.push(format!("{} →{ov}bp→ {}", part.0, next.0));
    }

    Assembly { circular: true, seq, order: order.to_vec(), junctions }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cloning::digest;
    use crate::enzymes::find_sites;

    #[test]
    fn golden_gate_does_not_panic_when_part_shorter_than_overhang() {
        // A part whose sequence is shorter than its 4-bp overhang must not panic
        // the assembler (it can arise from two cut sites closer than the overhang).
        let end = FragmentEnd { enzyme: Some("X"), overhang: Overhang::FivePrime(4) };
        let part = |name: &str| Part {
            name: name.to_string(),
            seq: b"AC".to_vec(), // 2 bp, < the 4-bp overhang
            left: end.clone(),
            right: end.clone(),
        };
        // Should return a result (or None) without panicking.
        let _ = golden_gate(&[part("a"), part("b")]);
    }

    #[test]
    fn parts_from_circular_digest_reconstruct_sequence() {
        // Two cutters on a circular molecule: EcoRI (GAATTC) and BamHI (GGATCC).
        let seq = b"GAATTCAAAAGGATCCTTTT".to_vec();
        let sites = find_sites(&seq, true);
        let frags = digest(seq.len(), true, &sites);
        let parts = fragments_to_parts(&seq, true, &frags);
        // The fragment slices partition the molecule, so their concatenation
        // (in cut order) is a rotation of the original of identical length.
        let total: usize = parts.iter().map(|p| p.seq.len()).sum();
        assert_eq!(total, seq.len());
    }

    #[test]
    fn golden_gate_round_trips_circular() {
        // Circular plasmid cut by two enzymes → 2 sticky-ended parts that must
        // re-ligate back into a circle of the original length.
        let seq = b"GAATTCAAAAGGATCCTTTT".to_vec();
        let sites = find_sites(&seq, true);
        let frags = digest(seq.len(), true, &sites);
        assert_eq!(frags.len(), 2);
        let parts = fragments_to_parts(&seq, true, &frags);

        let asm = golden_gate(&parts).expect("two compatible sticky fragments must re-circularise");
        assert!(asm.circular);
        assert_eq!(asm.order.len(), parts.len());

        // Overhang accounting: each of the two junctions shares a 5'+4 overhang
        // (EcoRI and BamHI both leave 4-base 5' overhangs), counted once. The
        // raw fragment slices sum to seq.len(); each closed junction removes the
        // duplicated overhang from the downstream part, so the product is
        // seq.len() − (sum of the two junction overhangs).
        let shared: usize = parts
            .iter()
            .map(|p| overhang_len(p.left.overhang))
            .sum();
        let expected = seq.len() - shared;
        assert_eq!(asm.seq.len(), expected);
        assert_eq!(asm.seq.len(), seq.len() - 8); // 2 junctions × 4 bp
    }

    #[test]
    fn golden_gate_empty_is_none() {
        assert!(golden_gate(&[]).is_none());
    }

    #[test]
    fn gibson_closes_three_overlapping_fragments() {
        // Build a circle of length 30, then split it into three linear pieces
        // that share 6-bp overlaps at every join (including the wrap).
        let circle: Vec<u8> = b"ATGCGTAAACCCGGGTTTAGCTAGCATCGGA".to_vec()[..30].to_vec();
        let ov = 6;
        let n = circle.len();
        // Cut points at 0, 10, 20; each piece extends `ov` past its cut into the
        // next region so neighbours overlap by `ov`.
        let piece = |start: usize, end: usize| -> Vec<u8> {
            let mut v = Vec::new();
            let mut i = start;
            let stop = end + ov;
            while i < stop {
                v.push(circle[i % n]);
                i += 1;
            }
            v
        };
        let parts = vec![
            ("a".to_string(), piece(0, 10)),
            ("b".to_string(), piece(10, 20)),
            ("c".to_string(), piece(20, 30)),
        ];

        let asm = gibson(&parts, ov).expect("overlapping fragments must close into a circle");
        assert!(asm.circular);
        assert_eq!(asm.order.len(), 3);
        // Three pieces of length 16 (10 + 6 overlap) sharing three 6-bp overlaps:
        // 3*16 − 3*6 = 30, the original circle length.
        assert_eq!(asm.seq.len(), n);
    }

    #[test]
    fn gibson_non_overlapping_is_none() {
        let parts = vec![
            ("x".to_string(), b"AAAAAAAAAA".to_vec()),
            ("y".to_string(), b"CCCCCCCCCC".to_vec()),
            ("z".to_string(), b"GGGGGGGGGG".to_vec()),
        ];
        // No 6-bp exact overlaps between distinct homopolymers → cannot close.
        assert!(gibson(&parts, 6).is_none());
    }

    #[test]
    fn gibson_zero_overlap_rejected() {
        let parts = vec![("x".to_string(), b"ACGT".to_vec())];
        assert!(gibson(&parts, 0).is_none());
    }
}
