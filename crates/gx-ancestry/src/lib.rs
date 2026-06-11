//! `gx-ancestry` estimates continental ancestry from a sample's autosomal
//! ancestry-informative-marker (AIM) genotypes, by supervised admixture against
//! 1000 Genomes phase 3 super-population allele frequencies. Pure library
//! (no UI, no network), depends only on `gx-core`.
//!
//! # Honesty and scope
//!
//! This is a coarse, continental-level estimate of genetic similarity to five
//! reference super-populations (African, Admixed American, East Asian, European,
//! South Asian). It is NOT an identity or ethnicity test. The panel is small
//! (tens of markers, embedded in `panel`), so fine population structure and some
//! recent mixed ancestry cannot be resolved. The reference frequencies come from
//! the 1000 Genomes Project (phase 3) via Ensembl; see `scripts/gen_aim_panel.py`.

mod panel;

pub use panel::{SUPERPOP_NAMES, SUPERPOPS};

use gx_core::{Assembly, Genotype, VariantStore};

/// Minimum genotyped AIMs required to attempt an estimate.
const MIN_MARKERS: usize = 8;

/// A continental admixture estimate.
pub struct AncestryEstimate {
    /// Admixture proportions summing to ~1.0, in [`SUPERPOPS`] order
    /// (`[AFR, AMR, EAS, EUR, SAS]`).
    pub proportions: [f64; 5],
    /// How many panel AIMs were genotyped and used.
    pub markers_used: usize,
    /// Total AIMs in the reference panel.
    pub markers_total: usize,
    /// Caveats (plain text; no em-dashes or en-dashes).
    pub note: String,
}

impl AncestryEstimate {
    /// Super-population display names paired with their proportions, largest first.
    pub fn ranked(&self) -> Vec<(&'static str, f64)> {
        let mut v: Vec<(&'static str, f64)> =
            SUPERPOP_NAMES.iter().copied().zip(self.proportions).collect();
        v.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        v
    }
}

/// Estimate continental admixture proportions from a sample's AIM genotypes.
///
/// Returns `None` when fewer than [`MIN_MARKERS`] panel SNPs were genotyped, so
/// the caller can show an appropriate "not enough coverage" message.
pub fn estimate(store: &VariantStore) -> Option<AncestryEstimate> {
    // (effect-allele dosage, ploidy, per-super-population effect-allele freqs)
    let mut obs: Vec<(u8, u8, [f64; 5])> = Vec::new();
    for aim in panel::PANEL {
        let Some(geno) = lookup(store, aim) else {
            continue;
        };
        if geno.is_no_call() {
            continue;
        }
        let ploidy = geno.alleles.len().min(2) as u8;
        if ploidy == 0 {
            continue;
        }
        let dosage = geno
            .alleles
            .iter()
            .take(2)
            .filter(|a| a.as_bytes().first().map(u8::to_ascii_uppercase) == Some(aim.effect))
            .count() as u8;
        obs.push((dosage, ploidy, aim.freq));
    }
    if obs.len() < MIN_MARKERS {
        return None;
    }

    let proportions = em(&obs);
    let note = format!(
        "Estimated genetic similarity to five 1000 Genomes reference \
         super-populations from {} ancestry-informative markers. This is a coarse, \
         continental-level estimate, not an identity or ethnicity test; a small marker \
         panel cannot resolve fine structure or every recent mixed ancestry.",
        obs.len()
    );
    Some(AncestryEstimate {
        proportions,
        markers_used: obs.len(),
        markers_total: panel::PANEL.len(),
        note,
    })
}

/// Resolve a panel marker's genotype in the sample: by rsID first (robust across
/// builds), then by build-appropriate position.
fn lookup<'a>(store: &'a VariantStore, aim: &panel::Aim) -> Option<&'a Genotype> {
    if let Some(v) = store.find_by_rsid(aim.rsid) {
        return v.genotype.as_ref();
    }
    let pos1 = match store.assembly() {
        Assembly::Grch37 => aim.pos37 as u64,
        Assembly::Grch38 => aim.pos38 as u64,
    };
    store
        .variants_in(aim.contig, pos1.saturating_sub(1), pos1)
        .first()
        .and_then(|v| v.genotype.as_ref())
}

/// Supervised-admixture EM: find the mixture `q` of reference super-populations
/// that maximises the likelihood of the observed effect-allele counts. This is
/// the standard ADMIXTURE-style update with the reference frequencies fixed.
fn em(obs: &[(u8, u8, [f64; 5])]) -> [f64; 5] {
    let mut q = [0.2f64; 5];
    for _ in 0..500 {
        let mut acc = [0.0f64; 5];
        let mut total = 0.0f64;
        for &(dosage, ploidy, f) in obs {
            let g = dosage as f64;
            let n = ploidy as f64;
            let p = (0..5).map(|k| q[k] * f[k]).sum::<f64>().clamp(1e-9, 1.0 - 1e-9);
            for k in 0..5 {
                let post_eff = q[k] * f[k] / p;
                let post_ref = q[k] * (1.0 - f[k]) / (1.0 - p);
                acc[k] += g * post_eff + (n - g) * post_ref;
            }
            total += n;
        }
        let mut max_change = 0.0f64;
        for k in 0..5 {
            let nk = acc[k] / total;
            max_change = max_change.max((nk - q[k]).abs());
            q[k] = nk;
        }
        if max_change < 1e-7 {
            break;
        }
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;
    use gx_core::Variant;

    /// Build a store whose AIM genotypes are the rounded typical dosages of
    /// super-population `sp`, for recovering that population in the EM.
    fn store_for(sp: usize) -> VariantStore {
        let mut s = VariantStore::new(Assembly::Grch37);
        for aim in panel::PANEL {
            let dose = (2.0 * aim.freq[sp]).round() as usize;
            let eff = (aim.effect as char).to_string();
            let other = if aim.effect == b'A' { "G" } else { "A" }.to_string();
            let alleles = match dose {
                2 => vec![eff.clone(), eff],
                1 => vec![eff, other],
                _ => vec![other.clone(), other],
            };
            s.push(Variant {
                rsid: Some(aim.rsid.to_string()),
                contig: aim.contig.to_string(),
                pos: aim.pos37 as u64 - 1,
                ref_allele: None,
                alt_alleles: vec![],
                genotype: Some(Genotype::from_alleles(alleles)),
            });
        }
        s.finalize();
        s
    }

    fn argmax(p: &[f64; 5]) -> usize {
        (0..5)
            .max_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap())
            .unwrap()
    }

    #[test]
    fn recovers_well_separated_populations() {
        // AFR(0), EAS(2), EUR(3), SAS(4) are continentally distinct; a typical
        // member of each must come back with that population dominant.
        for sp in [0usize, 2, 3, 4] {
            let est = estimate(&store_for(sp)).expect("enough markers");
            assert_eq!(
                argmax(&est.proportions),
                sp,
                "{} not dominant: {:?}",
                SUPERPOPS[sp],
                est.proportions
            );
            let sum: f64 = est.proportions.iter().sum();
            assert!((sum - 1.0).abs() < 1e-6, "proportions sum {sum}");
        }
    }

    #[test]
    fn admixed_american_is_a_plausible_mixture() {
        // AMR is itself admixed, so we only require it to rank in the top three
        // rather than be strictly dominant.
        let est = estimate(&store_for(1)).expect("enough markers");
        let mut idx = [0, 1, 2, 3, 4];
        idx.sort_by(|&a, &b| est.proportions[b].partial_cmp(&est.proportions[a]).unwrap());
        assert!(idx[..3].contains(&1), "AMR not in top 3: {:?}", est.proportions);
    }

    #[test]
    fn too_few_markers_is_none() {
        let store = VariantStore::new(Assembly::Grch37);
        assert!(estimate(&store).is_none());
    }
}
