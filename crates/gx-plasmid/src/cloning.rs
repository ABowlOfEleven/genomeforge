//! Restriction digest → fragments, and ligation compatibility — the core of a
//! cloning simulation. Overhang type/length is derived from each enzyme's cut
//! offset (palindromic sites: overhang = site_len − 2·cut).

use std::collections::HashMap;

use crate::enzymes::{ENZYMES, RestrictionSite};

/// The end left behind by a cut (or a free molecule end on a linear molecule).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Overhang {
    Blunt,
    /// 5' overhang of N bases.
    FivePrime(u32),
    /// 3' overhang of N bases.
    ThreePrime(u32),
}

impl Overhang {
    pub fn from_cut(site_len: usize, cut: usize) -> Overhang {
        let d = site_len as i32 - 2 * cut as i32;
        match d.cmp(&0) {
            std::cmp::Ordering::Greater => Overhang::FivePrime(d as u32),
            std::cmp::Ordering::Less => Overhang::ThreePrime((-d) as u32),
            std::cmp::Ordering::Equal => Overhang::Blunt,
        }
    }

    pub fn label(self) -> String {
        match self {
            Overhang::Blunt => "blunt".to_string(),
            Overhang::FivePrime(n) => format!("5' +{n}"),
            Overhang::ThreePrime(n) => format!("3' +{n}"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct FragmentEnd {
    /// The enzyme that cut here, or `None` for a free end of a linear molecule.
    pub enzyme: Option<&'static str>,
    pub overhang: Overhang,
}

#[derive(Debug, Clone)]
pub struct Fragment {
    pub start: usize,
    pub end: usize,
    pub len: usize,
    pub left: FragmentEnd,
    pub right: FragmentEnd,
}

/// Digest a molecule of length `seq_len` at the given sites (already filtered to
/// the chosen enzymes). Linear molecules yield `cuts + 1` fragments; circular
/// molecules yield `cuts` fragments (or one, linearised, for a single cutter).
pub fn digest(seq_len: usize, circular: bool, sites: &[RestrictionSite]) -> Vec<Fragment> {
    let overhang_by_name: HashMap<&str, Overhang> = ENZYMES
        .iter()
        .map(|e| (e.name, Overhang::from_cut(e.site.len(), e.cut)))
        .collect();
    let oh = |s: &RestrictionSite| overhang_by_name.get(s.enzyme).copied().unwrap_or(Overhang::Blunt);

    let mut cuts: Vec<(usize, &'static str, Overhang)> =
        sites.iter().map(|s| (s.cut, s.enzyme, oh(s))).collect();
    cuts.sort_by_key(|(c, _, _)| *c);

    let n = seq_len;
    let mut frags = Vec::new();

    if cuts.is_empty() {
        frags.push(Fragment {
            start: 0,
            end: n,
            len: n,
            left: FragmentEnd { enzyme: None, overhang: Overhang::Blunt },
            right: FragmentEnd { enzyme: None, overhang: Overhang::Blunt },
        });
        return frags;
    }

    if circular {
        let m = cuts.len();
        for i in 0..m {
            let (c, enz, ohl) = cuts[i];
            let (nc, nenz, ohr) = cuts[(i + 1) % m];
            let len = if m == 1 {
                n
            } else if (i + 1) % m == 0 {
                (n - c) + nc
            } else {
                nc - c
            };
            frags.push(Fragment {
                start: c,
                end: nc,
                len,
                left: FragmentEnd { enzyme: Some(enz), overhang: ohl },
                right: FragmentEnd { enzyme: Some(nenz), overhang: ohr },
            });
        }
    } else {
        let mut bounds = vec![0usize];
        bounds.extend(cuts.iter().map(|(c, _, _)| *c));
        bounds.push(n);
        for j in 0..bounds.len() - 1 {
            let start = bounds[j];
            let end = bounds[j + 1];
            let left = if j == 0 {
                FragmentEnd { enzyme: None, overhang: Overhang::Blunt }
            } else {
                let (_, enz, o) = cuts[j - 1];
                FragmentEnd { enzyme: Some(enz), overhang: o }
            };
            let right = if j == bounds.len() - 2 {
                FragmentEnd { enzyme: None, overhang: Overhang::Blunt }
            } else {
                let (_, enz, o) = cuts[j];
                FragmentEnd { enzyme: Some(enz), overhang: o }
            };
            frags.push(Fragment { start, end, len: end.saturating_sub(start), left, right });
        }
    }
    frags
}

/// Whether two fragment ends can be ligated. Approximation: matching overhang
/// class + length (blunt↔blunt always; sticky ends need complementary
/// sequence, which holds for same-enzyme ends).
pub fn ends_compatible(a: &FragmentEnd, b: &FragmentEnd) -> bool {
    a.overhang == b.overhang
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzymes::find_sites;

    #[test]
    fn overhang_classes() {
        assert_eq!(Overhang::from_cut(6, 1), Overhang::FivePrime(4)); // EcoRI G^AATTC
        assert_eq!(Overhang::from_cut(6, 5), Overhang::ThreePrime(4)); // PstI CTGCA^G
        assert_eq!(Overhang::from_cut(6, 3), Overhang::Blunt); // SmaI CCC^GGG
    }

    #[test]
    fn linear_digest_two_cutters() {
        // EcoRI then BamHI -> 3 fragments, free ends on the outer two.
        let seq = b"AAAGAATTCAAAAAAAGGATCCAAA";
        let sites = find_sites(seq, false);
        let frags = digest(seq.len(), false, &sites);
        assert_eq!(frags.len(), 3);
        assert!(frags[0].left.enzyme.is_none()); // 5' free end
        assert!(frags[2].right.enzyme.is_none()); // 3' free end
        assert_eq!(frags.iter().map(|f| f.len).sum::<usize>(), seq.len());
    }

    #[test]
    fn circular_single_cutter_linearises() {
        let seq = b"AAAGAATTCAAAAAAAAAAA";
        let sites = find_sites(seq, true);
        let frags = digest(seq.len(), true, &sites);
        assert_eq!(frags.len(), 1);
        assert_eq!(frags[0].len, seq.len());
        assert_eq!(frags[0].left.enzyme, Some("EcoRI"));
    }

    #[test]
    fn compatible_ends() {
        let a = FragmentEnd { enzyme: Some("EcoRI"), overhang: Overhang::FivePrime(4) };
        let b = FragmentEnd { enzyme: Some("EcoRI"), overhang: Overhang::FivePrime(4) };
        let c = FragmentEnd { enzyme: Some("SmaI"), overhang: Overhang::Blunt };
        assert!(ends_compatible(&a, &b));
        assert!(!ends_compatible(&a, &c));
    }
}
