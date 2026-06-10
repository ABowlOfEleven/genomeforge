use serde::{Deserialize, Serialize};

use crate::range::GenomicRange;

/// The kind of an annotation feature drawn on the gene-model track.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum FeatureKind {
    Gene,
    Transcript,
    Exon,
    Cds,
    FivePrimeUtr,
    ThreePrimeUtr,
    Regulatory,
    Other(String),
}

impl FeatureKind {
    /// Parse the `feature_type` / SO term Ensembl returns.
    pub fn from_ensembl(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "gene" => FeatureKind::Gene,
            "transcript" | "mrna" => FeatureKind::Transcript,
            "exon" => FeatureKind::Exon,
            "cds" => FeatureKind::Cds,
            "five_prime_utr" => FeatureKind::FivePrimeUtr,
            "three_prime_utr" => FeatureKind::ThreePrimeUtr,
            other if other.contains("regulatory") || other.contains("enhancer") => {
                FeatureKind::Regulatory
            }
            other => FeatureKind::Other(other.to_string()),
        }
    }
}

/// An annotation feature (gene, transcript, exon, …) with its genomic span.
///
/// `parent` links exons/CDS to their transcript and transcripts to their gene,
/// mirroring the GFF3 hierarchy Ensembl exposes.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Feature {
    pub id: String,
    pub name: Option<String>,
    pub kind: FeatureKind,
    pub range: GenomicRange,
    pub parent: Option<String>,
    pub biotype: Option<String>,
}

impl Feature {
    /// Best human-facing label: the symbol/name if present, else the id.
    pub fn label(&self) -> &str {
        self.name.as_deref().unwrap_or(&self.id)
    }
}
