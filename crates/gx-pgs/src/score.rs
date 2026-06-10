//! Apply a [`ScoreFile`] to a sample's genotypes to compute a polygenic score.
//!
//! The raw score is `Σ dosage_i · weight_i` over variants present in the sample.
//! When the scoring file carries effect-allele frequencies, we also estimate a
//! **population percentile** from a normal approximation: the score's expected
//! mean `Σ 2·AF·w` and variance `Σ 2·AF·(1−AF)·w²` under Hardy–Weinberg, giving
//! `z = (score − mean)/sd` and `Φ(z)`. This ignores LD and environment and
//! assumes the reference frequencies apply to the individual — surfaced as a
//! caveat in the UI.

use gx_core::VariantStore;

use crate::model::ScoreFile;

#[derive(Debug, Clone)]
pub struct PrsResult {
    pub trait_name: String,
    pub id: String,
    /// Variants from the score found (and unambiguously matched) in the sample.
    pub matched: usize,
    /// Total variants in the scoring file.
    pub total: usize,
    /// Variants skipped due to allele/strand ambiguity.
    pub ambiguous: usize,
    pub raw_score: f64,
    pub mean: Option<f64>,
    pub sd: Option<f64>,
    pub z: Option<f64>,
    /// Population percentile in `[0, 100]`, when AF data allowed an estimate.
    pub percentile: Option<f64>,
}

impl PrsResult {
    pub fn coverage(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.matched as f64 / self.total as f64
        }
    }
}

fn complement1(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        other => other,
    }
}

/// Complement of a single-base allele, or `None` for multi-base (indels).
fn complement(allele: &str) -> Option<String> {
    let bytes = allele.as_bytes();
    if bytes.len() == 1 {
        Some((complement1(bytes[0]) as char).to_string())
    } else {
        None
    }
}

/// Effect-allele dosage (0/1/2) for a genotype, handling strand flips. Returns
/// `None` if the sample alleles can't be reconciled with the score's alleles.
fn dosage(sample: &[String], effect: &str, other: Option<&str>) -> Option<u8> {
    let up: Vec<String> = sample.iter().map(|a| a.to_ascii_uppercase()).collect();
    if up.is_empty() {
        return None;
    }

    // Palindromic SNV (effect/other are complementary, e.g. A/T or C/G): the
    // strand is unresolvable from genotype alone, so the dosage can't be trusted.
    if let Some(o) = other {
        if complement(effect).as_deref() == Some(o) {
            return None;
        }
    }

    let count_ea = |ea: &str| up.iter().filter(|a| a.as_str() == ea).count() as u8;
    let direct = |ea: &str, oa: Option<&str>| -> Option<u8> {
        match oa {
            // With both alleles known, require every sample allele to be one of them.
            Some(o) if up.iter().all(|a| a == ea || a == o) => Some(count_ea(ea)),
            Some(_) => None,
            // Effect allele only: match if it appears at all.
            None if up.iter().any(|a| a == ea) => Some(count_ea(ea)),
            None => None,
        }
    };

    if let Some(d) = direct(effect, other) {
        return Some(d);
    }
    // Opposite-strand fallback — only when both alleles are known. With the
    // other allele unknown a flip is unresolvable and could fabricate a dosage.
    let o = other?;
    let ce = complement(effect)?;
    let co = complement(o)?;
    direct(&ce, Some(&co))
}

/// Apply a scoring file to a sample's [`VariantStore`].
pub fn apply(score: &ScoreFile, store: &VariantStore) -> PrsResult {
    let mut raw = 0.0;
    let mut mean = 0.0;
    let mut var = 0.0;
    let mut matched = 0usize;
    let mut ambiguous = 0usize;
    let mut af_count = 0usize;

    for sv in &score.variants {
        let Some(rsid) = &sv.rsid else { continue };
        let Some(variant) = store.find_by_rsid(rsid) else { continue };
        let Some(gt) = &variant.genotype else { continue };
        if gt.is_no_call() {
            continue;
        }
        match dosage(&gt.alleles, &sv.effect_allele, sv.other_allele.as_deref()) {
            Some(d) => {
                raw += d as f64 * sv.weight;
                matched += 1;
                if let Some(af) = sv.allele_freq {
                    mean += 2.0 * af * sv.weight;
                    var += 2.0 * af * (1.0 - af) * sv.weight * sv.weight;
                    af_count += 1;
                }
            }
            None => ambiguous += 1,
        }
    }

    // Only report a percentile if AF covered most matched variants.
    let (mean_o, sd_o, z_o, pct_o) = if matched > 0 && af_count as f64 >= 0.8 * matched as f64 && var > 0.0 {
        let sd = var.sqrt();
        let z = (raw - mean) / sd;
        (Some(mean), Some(sd), Some(z), Some(normal_cdf(z) * 100.0))
    } else {
        (None, None, None, None)
    };

    PrsResult {
        trait_name: score.trait_name.clone(),
        id: score.id.clone(),
        matched,
        total: score.variants.len(),
        ambiguous,
        raw_score: raw,
        mean: mean_o,
        sd: sd_o,
        z: z_o,
        percentile: pct_o,
    }
}

/// Standard normal CDF via an erf approximation (Abramowitz & Stegun 7.1.26).
pub fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2))
}

fn erf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-x * x).exp();
    if x >= 0.0 {
        y
    } else {
        -y
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::ScoreVariant;
    use gx_core::{Assembly, Genotype, Variant};

    fn variant(rsid: &str, gt: &str) -> Variant {
        Variant {
            rsid: Some(rsid.to_string()),
            contig: "1".into(),
            pos: 0,
            ref_allele: None,
            alt_alleles: vec![],
            genotype: Some(Genotype::from_call(gt)),
        }
    }

    fn sv(rsid: &str, effect: &str, other: &str, weight: f64, af: Option<f64>) -> ScoreVariant {
        ScoreVariant {
            rsid: Some(rsid.into()),
            chrom: None,
            pos: None,
            effect_allele: effect.into(),
            other_allele: Some(other.into()),
            weight,
            allele_freq: af,
        }
    }

    #[test]
    fn dosage_counts_effect_alleles() {
        assert_eq!(dosage(&["A".into(), "A".into()], "A", Some("G")), Some(2));
        assert_eq!(dosage(&["A".into(), "G".into()], "A", Some("G")), Some(1));
        assert_eq!(dosage(&["G".into(), "G".into()], "A", Some("G")), Some(0));
    }

    #[test]
    fn dosage_handles_strand_flip() {
        // Effect A/G but sample reported on the other strand as T/C.
        assert_eq!(dosage(&["T".into(), "T".into()], "A", Some("G")), Some(2));
        assert_eq!(dosage(&["T".into(), "C".into()], "A", Some("G")), Some(1));
    }

    #[test]
    fn palindromic_snp_is_ambiguous() {
        // A/T and C/G are reverse-complement pairs -> strand unresolvable -> skip.
        assert_eq!(dosage(&["A".into(), "T".into()], "A", Some("T")), None);
        assert_eq!(dosage(&["G".into(), "G".into()], "C", Some("G")), None);
    }

    #[test]
    fn effect_only_does_not_fabricate_via_flip() {
        // Other allele unknown and effect absent: must not flip-and-count.
        assert_eq!(dosage(&["G".into(), "G".into()], "C", None), None);
    }

    #[test]
    fn applies_score_and_coverage() {
        let mut store = VariantStore::new(Assembly::Grch38);
        store.push(variant("rs100", "AA")); // effect A -> dosage 2, weight 0.5 -> 1.0
        store.push(variant("rs200", "TC")); // effect T -> dosage 1, weight -0.2 -> -0.2
        store.finalize();

        let score = ScoreFile {
            id: "PGS1".into(),
            trait_name: "Test".into(),
            build: None,
            declared_count: Some(3),
            variants: vec![
                sv("rs100", "A", "G", 0.5, Some(0.3)),
                sv("rs200", "T", "C", -0.2, Some(0.1)),
                sv("rs999", "G", "A", 1.0, Some(0.2)), // not in sample
            ],
        };
        let r = apply(&score, &store);
        assert_eq!(r.matched, 2);
        assert_eq!(r.total, 3);
        assert!((r.raw_score - 0.8).abs() < 1e-9);
        assert!(r.percentile.is_some());
        let p = r.percentile.unwrap();
        assert!((0.0..=100.0).contains(&p));
    }

    #[test]
    fn no_af_means_no_percentile() {
        let mut store = VariantStore::new(Assembly::Grch38);
        store.push(variant("rs100", "AA"));
        store.finalize();
        let score = ScoreFile {
            id: "PGS1".into(),
            trait_name: "Test".into(),
            build: None,
            declared_count: None,
            variants: vec![sv("rs100", "A", "G", 0.5, None)],
        };
        let r = apply(&score, &store);
        assert_eq!(r.matched, 1);
        assert!(r.percentile.is_none());
    }
}
