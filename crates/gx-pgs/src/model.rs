use serde::{Deserialize, Serialize};

/// One weighted variant in a polygenic score.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoreVariant {
    pub rsid: Option<String>,
    pub chrom: Option<String>,
    /// 0-based position, if the file provided one.
    pub pos: Option<u64>,
    pub effect_allele: String,
    pub other_allele: Option<String>,
    pub weight: f64,
    /// Effect-allele frequency, if the file provided one (enables a percentile).
    pub allele_freq: Option<f64>,
}

/// A parsed PGS Catalog scoring file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoreFile {
    pub id: String,
    pub trait_name: String,
    pub build: Option<String>,
    /// `variants_number` declared in the header (may differ from parsed count).
    pub declared_count: Option<usize>,
    pub variants: Vec<ScoreVariant>,
}

impl ScoreFile {
    pub fn variant_count(&self) -> usize {
        self.variants.len()
    }

    /// Fraction of variants that carry an effect-allele frequency.
    pub fn af_coverage(&self) -> f64 {
        if self.variants.is_empty() {
            return 0.0;
        }
        let with = self.variants.iter().filter(|v| v.allele_freq.is_some()).count();
        with as f64 / self.variants.len() as f64
    }
}
