use serde::{Deserialize, Serialize};

use crate::range::Locus;

/// Zygosity of a genotype call.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Zygosity {
    /// No call (e.g. `--`).
    NoCall,
    /// A single allele observed (haploid contexts: X/Y/MT in XY individuals).
    Hemizygous,
    /// Two identical alleles.
    Homozygous,
    /// Two differing alleles.
    Heterozygous,
}

/// The sample's observed alleles at a locus.
///
/// Each allele is stored as a short string so the type covers SNVs (`"A"`),
/// consumer-array indel encodings (`"I"` / `"D"`), and short indels alike. An
/// empty allele list (or `--`) is a no-call.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
pub struct Genotype {
    pub alleles: Vec<String>,
}

impl Genotype {
    /// Parse a consumer-array genotype field such as `"AG"`, `"TT"`, `"--"`,
    /// `"I"`, or `"DI"`. Each character is one allele; `-` and `0` mean missing.
    pub fn from_call(raw: &str) -> Self {
        let alleles = raw
            .trim()
            .chars()
            .filter(|c| !matches!(c, '-' | '0' | '.'))
            .map(|c| c.to_ascii_uppercase().to_string())
            .collect();
        Self { alleles }
    }

    /// Build directly from explicit allele strings (e.g. from a VCF GT).
    pub fn from_alleles(alleles: Vec<String>) -> Self {
        Self { alleles }
    }

    pub fn is_no_call(&self) -> bool {
        self.alleles.is_empty()
    }

    pub fn zygosity(&self) -> Zygosity {
        match self.alleles.as_slice() {
            [] => Zygosity::NoCall,
            [_] => Zygosity::Hemizygous,
            [a, b] if a == b => Zygosity::Homozygous,
            _ => Zygosity::Heterozygous,
        }
    }

    /// Compact display, e.g. `"A/G"` or `"--"` for a no-call.
    pub fn display(&self) -> String {
        if self.alleles.is_empty() {
            "--".to_string()
        } else {
            self.alleles.join("/")
        }
    }
}

/// A variant locus together with the sample's call at it.
///
/// `rsid` is the dbSNP identifier when known (the build-independent key used for
/// annotation). For consumer arrays `ref_allele`/`alt_alleles` are often unknown
/// and the meaningful data is the `genotype`; for VCF imports they are populated.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Variant {
    pub rsid: Option<String>,
    pub contig: String,
    /// 0-based position.
    pub pos: u64,
    pub ref_allele: Option<String>,
    pub alt_alleles: Vec<String>,
    pub genotype: Option<Genotype>,
}

impl Variant {
    pub fn locus(&self) -> Locus {
        Locus::new(self.contig.clone(), self.pos)
    }

    /// A stable display label: the rsID if present, else `contig:pos` (1-based).
    pub fn display_id(&self) -> String {
        match &self.rsid {
            Some(id) => id.clone(),
            None => format!("{}:{}", self.contig, self.pos + 1),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_calls() {
        assert_eq!(Genotype::from_call("AG").alleles, ["A", "G"]);
        assert_eq!(Genotype::from_call("tt").alleles, ["T", "T"]);
        assert!(Genotype::from_call("--").is_no_call());
        assert!(Genotype::from_call("").is_no_call());
        assert_eq!(Genotype::from_call("DI").alleles, ["D", "I"]);
        assert_eq!(Genotype::from_call("A").alleles, ["A"]);
    }

    #[test]
    fn zygosity_classification() {
        assert_eq!(Genotype::from_call("--").zygosity(), Zygosity::NoCall);
        assert_eq!(Genotype::from_call("A").zygosity(), Zygosity::Hemizygous);
        assert_eq!(Genotype::from_call("AA").zygosity(), Zygosity::Homozygous);
        assert_eq!(Genotype::from_call("AG").zygosity(), Zygosity::Heterozygous);
    }

    #[test]
    fn display_id_falls_back_to_locus() {
        let v = Variant {
            rsid: None,
            contig: "1".into(),
            pos: 999,
            ref_allele: None,
            alt_alleles: vec![],
            genotype: None,
        };
        assert_eq!(v.display_id(), "1:1000");
    }
}
