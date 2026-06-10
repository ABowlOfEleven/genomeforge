//! The in-memory document: the imported dataset plus everything fetched for it.

use std::collections::{HashMap, HashSet};

use gx_core::{Assembly, Feature, Variant, VariantStore};
use gx_annotate::VariantAnnotation;
use gx_io::{Imported, SequenceRecord};

/// An imported variant dataset (consumer raw or VCF) plus its annotations.
pub struct VariantDoc {
    pub store: VariantStore,
    pub assembly: Assembly,
    pub format_label: String,
    pub warnings: Vec<String>,
    /// rsID → annotation, filled in as the worker resolves them.
    pub annotations: HashMap<String, VariantAnnotation>,
}

impl VariantDoc {
    pub fn annotation_for(&self, v: &Variant) -> Option<&VariantAnnotation> {
        v.rsid.as_ref().and_then(|id| self.annotations.get(id))
    }

    /// Count of variants whose ClinVar call is clinically notable.
    pub fn notable_count(&self) -> usize {
        self.annotations
            .values()
            .filter(|a| a.is_clinically_notable())
            .count()
    }
}

pub enum Document {
    Variants(VariantDoc),
    Sequences(Vec<SequenceRecord>),
}

impl Document {
    pub fn from_imported(imported: Imported, assembly_override: Option<Assembly>) -> Self {
        match imported {
            Imported::Variants(v) => {
                let assembly = assembly_override.unwrap_or(v.assembly);
                Document::Variants(VariantDoc {
                    store: v.store,
                    assembly,
                    format_label: v.format_label,
                    warnings: v.warnings,
                    annotations: HashMap::new(),
                })
            }
            Imported::Sequences(s) => Document::Sequences(s),
        }
    }

    pub fn variant_doc(&self) -> Option<&VariantDoc> {
        match self {
            Document::Variants(v) => Some(v),
            _ => None,
        }
    }

    pub fn variant_doc_mut(&mut self) -> Option<&mut VariantDoc> {
        match self {
            Document::Variants(v) => Some(v),
            _ => None,
        }
    }
}

/// Reference data fetched from Ensembl for the current build, accumulated as the
/// user browses. Lives separately from the document because it is assembly-level,
/// not sample-level.
#[derive(Default)]
pub struct RefData {
    features: Vec<Feature>,
    feature_ids: HashSet<String>,
    /// Tiles (by region string) already requested, to avoid duplicate fetches.
    pub fetched_tiles: HashSet<String>,
    /// Sequence regions already requested.
    pub requested_seq: HashSet<String>,
    /// The most recently fetched reference sequence: `(contig, start_0based, bases)`.
    pub current_seq: Option<(String, u64, String)>,
}

impl RefData {
    /// Merge newly fetched features, de-duplicating by feature id.
    pub fn add_features(&mut self, features: Vec<Feature>) {
        for f in features {
            let key = if f.id.is_empty() {
                format!("{}:{}-{}", f.range.contig, f.range.start, f.range.end)
            } else {
                f.id.clone()
            };
            if self.feature_ids.insert(key) {
                self.features.push(f);
            }
        }
    }

    /// Features overlapping `[start, end)` on `contig`.
    pub fn features_in(&self, contig: &str, start: u64, end: u64) -> Vec<&Feature> {
        self.features
            .iter()
            .filter(|f| f.range.contig == contig && f.range.start < end && start < f.range.end)
            .collect()
    }

    pub fn clear(&mut self) {
        *self = RefData::default();
    }
}
