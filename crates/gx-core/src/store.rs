use std::collections::{BTreeMap, HashMap};

use rust_lapper::{Interval, Lapper};

use crate::assembly::Assembly;
use crate::feature::Feature;
use crate::variant::Variant;

/// In-memory store of an imported sample's variants.
///
/// Variants are bucketed by contig and, once [`finalize`](Self::finalize)d,
/// sorted by position so the browser can pull "everything visible in
/// `[start, end)`" with a binary search rather than a linear scan — important
/// when a consumer file carries ~1M SNPs.
#[derive(Debug, Default)]
pub struct VariantStore {
    assembly: Assembly,
    by_contig: BTreeMap<String, Vec<Variant>>,
    rsid_index: HashMap<String, (String, usize)>,
    total: usize,
    finalized: bool,
}

impl VariantStore {
    pub fn new(assembly: Assembly) -> Self {
        Self {
            assembly,
            ..Default::default()
        }
    }

    pub fn assembly(&self) -> Assembly {
        self.assembly
    }

    /// Add a variant. Call [`finalize`](Self::finalize) once all are added.
    pub fn push(&mut self, variant: Variant) {
        self.by_contig
            .entry(variant.contig.clone())
            .or_default()
            .push(variant);
        self.total += 1;
        self.finalized = false;
    }

    /// Sort each contig by position and (re)build the rsID lookup index.
    pub fn finalize(&mut self) {
        self.rsid_index.clear();
        for (contig, variants) in &mut self.by_contig {
            variants.sort_by_key(|v| v.pos);
            for (idx, v) in variants.iter().enumerate() {
                if let Some(rsid) = &v.rsid {
                    self.rsid_index.insert(rsid.clone(), (contig.clone(), idx));
                }
            }
        }
        self.finalized = true;
    }

    pub fn len(&self) -> usize {
        self.total
    }

    pub fn is_empty(&self) -> bool {
        self.total == 0
    }

    /// Contig names present, in sorted order.
    pub fn contigs(&self) -> impl Iterator<Item = &str> {
        self.by_contig.keys().map(String::as_str)
    }

    /// Variants whose position falls in `[start, end)` on `contig`.
    ///
    /// Returns an empty slice if the store has not been finalized or the contig
    /// is unknown.
    pub fn variants_in(&self, contig: &str, start: u64, end: u64) -> &[Variant] {
        if !self.finalized {
            return &[];
        }
        let Some(variants) = self.by_contig.get(contig) else {
            return &[];
        };
        let lo = variants.partition_point(|v| v.pos < start);
        let hi = variants.partition_point(|v| v.pos < end);
        &variants[lo..hi]
    }

    /// Min/max position present on a contig (for browser extent).
    pub fn contig_bounds(&self, contig: &str) -> Option<(u64, u64)> {
        let variants = self.by_contig.get(contig)?;
        let first = variants.first()?.pos;
        let last = variants.last()?.pos;
        Some((first, last))
    }

    pub fn find_by_rsid(&self, rsid: &str) -> Option<&Variant> {
        let (contig, idx) = self.rsid_index.get(rsid)?;
        self.by_contig.get(contig)?.get(*idx)
    }

    /// All rsIDs present, for batch annotation.
    pub fn rsids(&self) -> Vec<String> {
        self.by_contig
            .values()
            .flatten()
            .filter_map(|v| v.rsid.clone())
            .collect()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Variant> {
        self.by_contig.values().flatten()
    }
}

/// Interval-indexed store of gene-model features for fast "what overlaps the
/// visible window?" queries, backed by one [`Lapper`] per contig.
#[derive(Debug, Default)]
pub struct FeatureStore {
    by_contig: HashMap<String, Lapper<u32, Feature>>,
}

impl FeatureStore {
    /// Build from a flat list of features, grouping by contig.
    pub fn from_features(features: Vec<Feature>) -> Self {
        let mut grouped: HashMap<String, Vec<Interval<u32, Feature>>> = HashMap::new();
        for f in features {
            let start = f.range.start.min(u32::MAX as u64) as u32;
            let stop = f.range.end.min(u32::MAX as u64) as u32;
            grouped
                .entry(f.range.contig.clone())
                .or_default()
                .push(Interval {
                    start,
                    stop,
                    val: f,
                });
        }
        let by_contig = grouped
            .into_iter()
            .map(|(contig, intervals)| (contig, Lapper::new(intervals)))
            .collect();
        Self { by_contig }
    }

    /// Features overlapping `[start, end)` on `contig`.
    pub fn overlapping(&self, contig: &str, start: u64, end: u64) -> Vec<&Feature> {
        let Some(lapper) = self.by_contig.get(contig) else {
            return Vec::new();
        };
        let start = start.min(u32::MAX as u64) as u32;
        let end = end.min(u32::MAX as u64) as u32;
        lapper.find(start, end).map(|iv| &iv.val).collect()
    }

    pub fn is_empty(&self) -> bool {
        self.by_contig.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::FeatureKind;
    use crate::range::{GenomicRange, Strand};

    fn snp(contig: &str, pos: u64, rsid: &str) -> Variant {
        Variant {
            rsid: Some(rsid.to_string()),
            contig: contig.to_string(),
            pos,
            ref_allele: None,
            alt_alleles: vec![],
            genotype: None,
        }
    }

    #[test]
    fn range_query_uses_sorted_order() {
        let mut store = VariantStore::new(Assembly::Grch38);
        // Push out of order to prove finalize sorts.
        store.push(snp("1", 500, "rs5"));
        store.push(snp("1", 100, "rs1"));
        store.push(snp("1", 300, "rs3"));
        store.push(snp("2", 100, "rs_other"));

        // Before finalize, queries return nothing (avoids unsorted slicing bugs).
        assert!(store.variants_in("1", 0, 1000).is_empty());

        store.finalize();
        let hits = store.variants_in("1", 150, 600);
        let ids: Vec<_> = hits.iter().map(|v| v.rsid.clone().unwrap()).collect();
        assert_eq!(ids, ["rs3", "rs5"]);
        assert_eq!(store.len(), 4);
        assert_eq!(store.contig_bounds("1"), Some((100, 500)));
    }

    #[test]
    fn rsid_lookup() {
        let mut store = VariantStore::new(Assembly::Grch38);
        store.push(snp("3", 42, "rs42"));
        store.finalize();
        assert_eq!(store.find_by_rsid("rs42").unwrap().pos, 42);
        assert!(store.find_by_rsid("rs999").is_none());
    }

    #[test]
    fn feature_overlap_query() {
        let gene = Feature {
            id: "ENSG1".into(),
            name: Some("BRCA1".into()),
            kind: FeatureKind::Gene,
            range: GenomicRange::new("17", 1000, 5000, Strand::Reverse),
            parent: None,
            biotype: Some("protein_coding".into()),
        };
        let store = FeatureStore::from_features(vec![gene]);
        assert_eq!(store.overlapping("17", 2000, 3000).len(), 1);
        assert_eq!(store.overlapping("17", 6000, 7000).len(), 0);
        assert_eq!(store.overlapping("1", 2000, 3000).len(), 0);
    }
}
