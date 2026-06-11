//! `gx-haplo` predicts a sample's maternal (mtDNA) and paternal (Y-DNA) major
//! haplogroups from consumer-genome SNP calls. It is a pure library (no UI) and
//! depends only on `gx-core`.
//!
//! # Honesty and scope
//!
//! This runs on real human DNA, so the design is deliberately conservative:
//!
//! * Marker positions and derived alleles are sourced from authoritative
//!   references, PhyloTree mtDNA Build 17 (van Oven & Kayser; see `mt`) and the
//!   ISOGG Y-DNA tree via the YBrowse SNP tables (see `ydna`), cross-checked
//!   against NCBI dbSNP placements where an authoritative rsID exists. Markers we
//!   could not source confidently are omitted rather than guessed.
//! * Output is framed as a "predicted major haplogroup (coarse resolution)". We
//!   prefer a correct coarse call (e.g. "U" or "R1b") over a precise wrong one.
//! * The rCRS reference sequence is itself mtDNA haplogroup H2a2a1, so H cannot be
//!   resolved by positive mutations; West-Eurasian R0 samples lacking U/J/T/etc.
//!   diagnostics are reported as "H / HV (R0 cluster)" with an explicit note that
//!   fine H subclades are not resolved.
//!
//! # Algorithm
//!
//! Each lineage is a rooted tree of haplogroup nodes carrying defining markers
//! (derived state at a position, with rsID where known). [`classify_maternal`]
//! and [`classify_paternal`] walk from the root, descending into a child whenever
//! that child has positive marker support, and report the deepest node reached.
//! See `tree::classify`.

mod mt;
mod tree;
mod ydna;

use gx_core::VariantStore;

use tree::Walk;

/// How much confidence to attach to a prediction.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Confidence {
    /// Several markers matched along the path and the terminal node had direct
    /// derived support.
    High,
    /// A coarse call with some support, but thin or stopping at an interior node.
    Moderate,
    /// Very little support; treat as a best guess at a broad cluster.
    Low,
}

/// Which parental lineage a call describes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Lineage {
    Maternal,
    Paternal,
}

impl Lineage {
    /// Human-facing label.
    pub fn label(self) -> &'static str {
        match self {
            Lineage::Maternal => "maternal (mtDNA)",
            Lineage::Paternal => "paternal (Y-DNA)",
        }
    }
}

/// A predicted haplogroup for one lineage.
pub struct HaploCall {
    /// Predicted major haplogroup, e.g. `"U5"`, `"H / HV (R0 cluster)"`, `"R1b"`.
    pub haplogroup: String,
    /// Which parental lineage this describes.
    pub lineage: Lineage,
    /// Confidence in the call.
    pub confidence: Confidence,
    /// Names of defining markers matched (derived) along the assigned path.
    pub supporting: Vec<String>,
    /// How many of this lineage's markers were actually genotyped in the data.
    pub tested: usize,
    /// Caveats and resolution note (plain text; no em-dashes or en-dashes).
    pub note: String,
}

/// Predict the maternal haplogroup from mtDNA SNPs.
///
/// Returns `None` when the data has no usable MT calls at any modelled marker.
pub fn classify_maternal(store: &VariantStore) -> Option<HaploCall> {
    let walk = tree::classify(store, mt::CONTIG, &mt::ROOT)?;
    Some(finish(walk, Lineage::Maternal, MATERNAL_BASE_NOTE))
}

/// Predict the paternal haplogroup from Y-chromosome SNPs.
///
/// Returns `None` when no usable Y data is present (typical for XX samples); the
/// integrator should show an appropriate message in that case.
pub fn classify_paternal(store: &VariantStore) -> Option<HaploCall> {
    let walk = tree::classify(store, ydna::CONTIG, &ydna::ROOT)?;
    Some(finish(walk, Lineage::Paternal, PATERNAL_BASE_NOTE))
}

const MATERNAL_BASE_NOTE: &str =
    "Predicted major maternal (mtDNA) haplogroup at coarse resolution from \
PhyloTree Build 17 markers; not a full mtDNA sequence classification.";

const PATERNAL_BASE_NOTE: &str =
    "Predicted major paternal (Y-DNA) haplogroup at coarse resolution from ISOGG \
tree markers; not a full Y sequence classification.";

/// Turn a tree walk into a public [`HaploCall`], deriving confidence and note.
fn finish(walk: Walk, lineage: Lineage, base_note: &str) -> HaploCall {
    let matched = walk.supporting.len();
    // A "cluster" assignment is one we deliberately keep coarse (e.g. the
    // West-Eurasian H / HV / R0 cluster), so it never earns High confidence.
    let is_cluster = walk.haplogroup.contains("cluster");

    // Confidence heuristic, driven by how many defining markers actually matched
    // along the assigned path (not by depth, since some path nodes are structural
    // pass-throughs whose own SNP we do not model):
    // * High: three or more markers matched and the call is a specific clade.
    // * Moderate: one or two markers matched (or a coarse cluster call).
    // * Low: nothing matched (broadest cluster consistent with the data).
    let confidence = if matched >= 3 && !is_cluster {
        Confidence::High
    } else if matched >= 1 {
        Confidence::Moderate
    } else {
        Confidence::Low
    };

    let mut note = String::from(base_note);
    if let Some(extra) = walk.note {
        note.push(' ');
        note.push_str(extra);
    }
    if matched == 0 {
        note.push(' ');
        note.push_str(
            "No defining mutation was directly matched, so this is the broadest \
cluster consistent with the available calls; treat it as a best guess.",
        );
    }

    HaploCall {
        haplogroup: walk.haplogroup,
        lineage,
        confidence,
        supporting: walk.supporting,
        tested: walk.tested_total,
        note,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gx_core::{Assembly, Genotype, Variant};

    /// Build an MT variant at a 1-based rCRS position with a single derived call.
    fn mt_var(pos1: u64, allele: &str) -> Variant {
        Variant {
            rsid: None,
            contig: "MT".into(),
            pos: pos1 - 1, // store is 0-based
            ref_allele: None,
            alt_alleles: vec![],
            genotype: Some(Genotype::from_alleles(vec![allele.into()])),
        }
    }

    /// Build a Y variant at a 1-based position with a single derived call.
    fn y_var(pos1: u64, allele: &str) -> Variant {
        Variant {
            rsid: None,
            contig: "Y".into(),
            pos: pos1 - 1,
            ref_allele: None,
            alt_alleles: vec![],
            genotype: Some(Genotype::from_alleles(vec![allele.into()])),
        }
    }

    fn finalized(asm: Assembly, vars: Vec<Variant>) -> VariantStore {
        let mut s = VariantStore::new(asm);
        for v in vars {
            s.push(v);
        }
        s.finalize();
        s
    }

    #[test]
    fn maternal_u_trio_classifies_as_u() {
        // The U-defining trio: rCRS-derived G/G/A at 11467/12308/12372.
        let store = finalized(
            Assembly::Grch37,
            vec![
                mt_var(11467, "G"),
                mt_var(12308, "G"),
                mt_var(12372, "A"),
            ],
        );
        let call = classify_maternal(&store).expect("MT data present");
        assert_eq!(call.lineage, Lineage::Maternal);
        assert!(
            call.haplogroup == "U" || call.haplogroup.starts_with('U') || call.haplogroup == "K",
            "expected a U-something, got {}",
            call.haplogroup
        );
        // All three U markers should be in the supporting set.
        for m in ["11467", "12308", "12372"] {
            assert!(call.supporting.iter().any(|s| s == m), "missing {m}");
        }
        assert!(call.tested >= 3);
    }

    #[test]
    fn maternal_u5_resolves_subclade() {
        let store = finalized(
            Assembly::Grch38,
            vec![
                mt_var(11467, "G"),
                mt_var(12308, "G"),
                mt_var(12372, "A"),
                mt_var(16192, "T"),
                mt_var(16270, "T"),
                mt_var(14793, "G"),
                mt_var(16256, "T"),
            ],
        );
        let call = classify_maternal(&store).expect("MT data present");
        assert!(
            call.haplogroup.starts_with("U5"),
            "expected U5x, got {}",
            call.haplogroup
        );
        assert_eq!(call.confidence, Confidence::High);
    }

    #[test]
    fn maternal_jt_shared_markers() {
        // 4216 (R2'JT) + 11251 (JT) are the shared J/T anchor; add T markers.
        let store = finalized(
            Assembly::Grch37,
            vec![
                mt_var(4216, "C"),
                mt_var(11251, "G"),
                mt_var(16126, "C"),
                mt_var(4917, "G"),
                mt_var(13368, "A"),
                mt_var(14905, "A"),
            ],
        );
        let call = classify_maternal(&store).expect("MT data present");
        assert!(
            call.haplogroup == "T" || call.haplogroup.starts_with('T') || call.haplogroup == "JT",
            "expected T/JT, got {}",
            call.haplogroup
        );
    }

    #[test]
    fn maternal_h_hv_conservative_cluster() {
        // West-Eurasian R sample: R markers present, HV marker present, but NO
        // U/J/T diagnostics. Must report the conservative H / HV cluster.
        let store = finalized(
            Assembly::Grch37,
            vec![
                mt_var(12705, "C"), // R
                mt_var(16223, "C"), // R
                mt_var(14766, "T"), // HV
            ],
        );
        let call = classify_maternal(&store).expect("MT data present");
        assert!(
            call.haplogroup.contains("HV") || call.haplogroup.contains('H'),
            "expected H / HV cluster, got {}",
            call.haplogroup
        );
        assert!(
            call.note.contains("H2a2a1") || call.note.contains("not resolved"),
            "note should explain coarse H resolution: {}",
            call.note
        );
        // No em-dashes / en-dashes anywhere in the note.
        assert!(!call.note.contains('\u{2014}') && !call.note.contains('\u{2013}'));
    }

    #[test]
    fn paternal_r1b_full_backbone() {
        // Walk M168 -> M89 -> M9 -> M45 -> M207 -> M269 (GRCh37 positions).
        let store = finalized(
            Assembly::Grch37,
            vec![
                y_var(14813991, "T"), // M168
                y_var(21917313, "T"), // M89
                y_var(21730257, "G"), // M9
                y_var(21867787, "A"), // M45
                y_var(15581983, "G"), // M207 (R)
                y_var(22739367, "C"), // M269 (R1b)
            ],
        );
        let call = classify_paternal(&store).expect("Y data present");
        assert_eq!(call.haplogroup, "R1b", "got {}", call.haplogroup);
        assert_eq!(call.lineage, Lineage::Paternal);
        assert_eq!(call.confidence, Confidence::High);
        for m in ["M168", "M89", "M9", "M45", "M207", "M269"] {
            assert!(call.supporting.iter().any(|s| s == m), "missing {m}");
        }
    }

    #[test]
    fn paternal_r1b_grch38_positions() {
        // Same lineage but using GRCh38 coordinates, proving build selection.
        let store = finalized(
            Assembly::Grch38,
            vec![
                y_var(12702062, "T"), // M168
                y_var(19755427, "T"), // M89
                y_var(19568371, "G"), // M9
                y_var(19705901, "A"), // M45
                y_var(13470103, "G"), // M207
                y_var(20577481, "C"), // M269
            ],
        );
        let call = classify_paternal(&store).expect("Y data present");
        assert_eq!(call.haplogroup, "R1b");
    }

    #[test]
    fn paternal_j2_classifies() {
        let store = finalized(
            Assembly::Grch37,
            vec![
                y_var(14813991, "T"), // M168
                y_var(21917313, "T"), // M89
                y_var(22749853, "C"), // M304 (J)
                y_var(14969634, "G"), // M172 (J2)
            ],
        );
        let call = classify_paternal(&store).expect("Y data present");
        assert_eq!(call.haplogroup, "J2", "got {}", call.haplogroup);
        assert!(call.supporting.iter().any(|s| s == "M172"));
    }

    #[test]
    fn paternal_shortcut_without_backbone() {
        // Array with only the terminal R1b SNP (no M168/M89/M9 backbone) should
        // still resolve to R1b via the root shortcut, at lower confidence.
        let store = finalized(Assembly::Grch37, vec![y_var(22739367, "C")]);
        let call = classify_paternal(&store).expect("Y data present");
        assert_eq!(call.haplogroup, "R1b");
    }

    #[test]
    fn empty_store_returns_none() {
        let store = finalized(Assembly::Grch37, vec![]);
        assert!(classify_maternal(&store).is_none());
        assert!(classify_paternal(&store).is_none());
    }

    #[test]
    fn unrelated_autosomal_data_returns_none() {
        // Variants on chr1 only: no MT or Y markers -> both None.
        let store = finalized(
            Assembly::Grch37,
            vec![Variant {
                rsid: Some("rs1".into()),
                contig: "1".into(),
                pos: 100,
                ref_allele: None,
                alt_alleles: vec![],
                genotype: Some(Genotype::from_alleles(vec!["A".into()])),
            }],
        );
        assert!(classify_maternal(&store).is_none());
        assert!(classify_paternal(&store).is_none());
    }

    #[test]
    fn no_call_at_marker_is_not_usable() {
        // A no-call genotype at a U position must not count as tested/derived.
        let store = finalized(
            Assembly::Grch37,
            vec![Variant {
                rsid: None,
                contig: "MT".into(),
                pos: 11467 - 1,
                ref_allele: None,
                alt_alleles: vec![],
                genotype: Some(Genotype::from_call("--")),
            }],
        );
        assert!(classify_maternal(&store).is_none());
    }

    #[test]
    fn lineage_labels() {
        assert_eq!(Lineage::Maternal.label(), "maternal (mtDNA)");
        assert_eq!(Lineage::Paternal.label(), "paternal (Y-DNA)");
    }
}
