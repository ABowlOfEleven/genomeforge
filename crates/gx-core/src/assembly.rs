use serde::{Deserialize, Serialize};

/// Human reference-genome assembly build.
///
/// rsID-based annotation (dbSNP, ClinVar, gnomAD via MyVariant.info) is
/// **build-independent**, so variant annotation is robust regardless of which
/// build a file uses. Genome-browser *coordinates* and reference-sequence
/// lookups are **not** build-independent, so we track the build of every
/// imported dataset and route Ensembl REST calls to the matching host.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum Assembly {
    /// GRCh38 / hg38 — the current human reference.
    #[default]
    Grch38,
    /// GRCh37 / hg19 — common in older consumer arrays (e.g. 23andMe v1–v4).
    Grch37,
}

impl Assembly {
    /// Ensembl REST base URL for this build (no trailing slash).
    pub fn ensembl_rest_host(self) -> &'static str {
        match self {
            Assembly::Grch38 => "https://rest.ensembl.org",
            Assembly::Grch37 => "https://grch37.rest.ensembl.org",
        }
    }

    /// Human-facing label, e.g. `"GRCh38 (hg38)"`.
    pub fn label(self) -> &'static str {
        match self {
            Assembly::Grch38 => "GRCh38 (hg38)",
            Assembly::Grch37 => "GRCh37 (hg19)",
        }
    }

    /// UCSC-style short name, e.g. `"hg38"`.
    pub fn ucsc_name(self) -> &'static str {
        match self {
            Assembly::Grch38 => "hg38",
            Assembly::Grch37 => "hg19",
        }
    }

    /// Ensembl assembly name, e.g. `"GRCh38"` (used by the REST assembly-map API).
    pub fn ensembl_name(self) -> &'static str {
        match self {
            Assembly::Grch38 => "GRCh38",
            Assembly::Grch37 => "GRCh37",
        }
    }

    /// All builds, for UI pickers.
    pub fn all() -> [Assembly; 2] {
        [Assembly::Grch38, Assembly::Grch37]
    }
}

impl std::fmt::Display for Assembly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.label())
    }
}
