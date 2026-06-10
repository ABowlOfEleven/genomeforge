//! Restriction enzymes and site finding.
//!
//! The built-in table is common cloning enzymes with **palindromic** recognition
//! sites, so scanning the top strand finds every double-strand site. Circular
//! molecules are scanned with origin wrap-around.

use std::collections::HashMap;

use crate::seq::iupac_matches;

#[derive(Debug, Clone, Copy)]
pub struct Enzyme {
    pub name: &'static str,
    pub site: &'static str,
    /// Top-strand cut offset within the site (e.g. EcoRI `G^AATTC` = 1).
    pub cut: usize,
}

/// A curated set of common, palindromic restriction enzymes.
pub const ENZYMES: &[Enzyme] = &[
    Enzyme { name: "EcoRI", site: "GAATTC", cut: 1 },
    Enzyme { name: "BamHI", site: "GGATCC", cut: 1 },
    Enzyme { name: "HindIII", site: "AAGCTT", cut: 1 },
    Enzyme { name: "NotI", site: "GCGGCCGC", cut: 2 },
    Enzyme { name: "XhoI", site: "CTCGAG", cut: 1 },
    Enzyme { name: "SalI", site: "GTCGAC", cut: 1 },
    Enzyme { name: "PstI", site: "CTGCAG", cut: 5 },
    Enzyme { name: "SmaI", site: "CCCGGG", cut: 3 },
    Enzyme { name: "KpnI", site: "GGTACC", cut: 5 },
    Enzyme { name: "SacI", site: "GAGCTC", cut: 5 },
    Enzyme { name: "SpeI", site: "ACTAGT", cut: 1 },
    Enzyme { name: "XbaI", site: "TCTAGA", cut: 1 },
    Enzyme { name: "NcoI", site: "CCATGG", cut: 1 },
    Enzyme { name: "NdeI", site: "CATATG", cut: 2 },
    Enzyme { name: "EcoRV", site: "GATATC", cut: 3 },
    Enzyme { name: "NheI", site: "GCTAGC", cut: 1 },
    Enzyme { name: "BglII", site: "AGATCT", cut: 1 },
    Enzyme { name: "ApaI", site: "GGGCCC", cut: 5 },
    Enzyme { name: "SphI", site: "GCATGC", cut: 5 },
    Enzyme { name: "AflII", site: "CTTAAG", cut: 1 },
    Enzyme { name: "ClaI", site: "ATCGAT", cut: 2 },
    Enzyme { name: "MluI", site: "ACGCGT", cut: 1 },
    Enzyme { name: "DraI", site: "TTTAAA", cut: 3 },
    Enzyme { name: "PvuII", site: "CAGCTG", cut: 3 },
];

#[derive(Debug, Clone)]
pub struct RestrictionSite {
    pub enzyme: &'static str,
    pub site: &'static str,
    /// 0-based start of the recognition site on the top strand.
    pub start: usize,
    /// 0-based top-strand cut position.
    pub cut: usize,
}

/// Find all sites for the built-in enzyme set.
pub fn find_sites(seq: &[u8], circular: bool) -> Vec<RestrictionSite> {
    find_sites_with(seq, circular, ENZYMES)
}

/// Find all sites for a given enzyme set.
pub fn find_sites_with(seq: &[u8], circular: bool, enzymes: &[Enzyme]) -> Vec<RestrictionSite> {
    let n = seq.len();
    let mut sites = Vec::new();
    if n == 0 {
        return sites;
    }
    for e in enzymes {
        let pat = e.site.as_bytes();
        let m = pat.len();
        if m == 0 || m > n {
            continue;
        }
        let scan_end = if circular { n } else { n - m + 1 };
        for i in 0..scan_end {
            let hit = (0..m).all(|k| iupac_matches(pat[k], seq[(i + k) % n]));
            if hit {
                sites.push(RestrictionSite {
                    enzyme: e.name,
                    site: e.site,
                    start: i,
                    cut: (i + e.cut) % n,
                });
            }
        }
    }
    sites.sort_by_key(|s| s.start);
    sites
}

/// Counts of sites per enzyme name.
pub fn site_counts(sites: &[RestrictionSite]) -> HashMap<&'static str, usize> {
    let mut counts = HashMap::new();
    for s in sites {
        *counts.entry(s.enzyme).or_insert(0) += 1;
    }
    counts
}

/// Sites whose enzyme cuts the molecule exactly once (the useful cloning sites).
pub fn unique_cutters(sites: &[RestrictionSite]) -> Vec<&RestrictionSite> {
    let counts = site_counts(sites);
    sites
        .iter()
        .filter(|s| counts.get(s.enzyme) == Some(&1))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_ecori_and_cut() {
        // GAATTC at index 3.
        let seq = b"AAAGAATTCAAA";
        let sites = find_sites(seq, false);
        let ecori: Vec<_> = sites.iter().filter(|s| s.enzyme == "EcoRI").collect();
        assert_eq!(ecori.len(), 1);
        assert_eq!(ecori[0].start, 3);
        assert_eq!(ecori[0].cut, 4); // G^AATTC -> cut after the G at index 3
    }

    #[test]
    fn circular_wraps_origin() {
        // Site GAATTC split across the origin: "ATTC...GA"
        let seq = b"ATTCAAAAAAGA"; // positions 10,11 = G,A ; 0,1,2,3 = A,T,T,C
        let linear = find_sites(seq, false);
        assert!(!linear.iter().any(|s| s.enzyme == "EcoRI"));
        let circ = find_sites(seq, true);
        assert!(circ.iter().any(|s| s.enzyme == "EcoRI" && s.start == 10));
    }

    #[test]
    fn unique_cutter_detection() {
        let seq = b"GAATTCAAAGAATTCAAAGGATCC"; // two EcoRI, one BamHI
        let sites = find_sites(seq, false);
        let unique: Vec<_> = unique_cutters(&sites).iter().map(|s| s.enzyme).collect();
        assert!(unique.contains(&"BamHI"));
        assert!(!unique.contains(&"EcoRI"));
    }
}
