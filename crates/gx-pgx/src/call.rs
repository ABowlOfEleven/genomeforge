//! Driver: resolve each gene's defining markers against the sample's variants,
//! count variant-allele copies (strand-aware), and assemble [`GeneResult`]s.
//!
//! Strand handling: consumer arrays report on the reference plus strand, but PGx
//! literature sometimes quotes the opposite strand (see crate docs and `genes.rs`).
//! Each [`Marker`] lists the accepted variant and reference allele letters for
//! BOTH strands, and we match a sample allele if it appears in either list. If a
//! sample allele matches neither, the marker is treated as not informative for
//! that genotype (it is dropped rather than guessed), which keeps an unexpected
//! strand or multiallelic call from producing a wrong copy count.

use gx_core::VariantStore;

use crate::genes::{GENES, GeneDef, Marker, MarkerCall};
use crate::model::GeneResult;

/// Resolve one marker against the store. Returns `None` if the marker is not
/// present or is a no-call; otherwise returns the variant-allele copy count and
/// ploidy.
fn resolve_marker(store: &VariantStore, marker: &Marker) -> Option<ResolvedMarker> {
    let variant = store.find_by_rsid(marker.rsid)?;
    let genotype = variant.genotype.as_ref()?;
    if genotype.is_no_call() {
        return None;
    }

    let mut variant_copies: u8 = 0;
    let mut counted: u8 = 0;
    for allele in &genotype.alleles {
        let a = allele.to_ascii_uppercase();
        if marker.variant_alleles.iter().any(|v| v.eq_ignore_ascii_case(&a)) {
            variant_copies += 1;
            counted += 1;
        } else if marker.reference_alleles.iter().any(|r| r.eq_ignore_ascii_case(&a)) {
            counted += 1;
        }
        // Alleles matching neither list (unexpected strand/base) are ignored.
    }

    if counted == 0 {
        // Genotyped but every allele was off-target; not informative.
        return None;
    }

    Some(ResolvedMarker {
        rsid: marker.rsid,
        variant_copies,
        ploidy: counted,
    })
}

struct ResolvedMarker {
    rsid: &'static str,
    variant_copies: u8,
    ploidy: u8,
}

/// Call one gene. Returns `None` if none of its markers were genotyped.
fn call_gene(store: &VariantStore, def: &GeneDef) -> Option<GeneResult> {
    let resolved: Vec<ResolvedMarker> = def
        .markers
        .iter()
        .filter_map(|m| resolve_marker(store, m))
        .collect();

    if resolved.is_empty() {
        return None;
    }

    let markers_used: Vec<String> = resolved.iter().map(|r| r.rsid.to_string()).collect();

    let calls: Vec<MarkerCall> = resolved
        .iter()
        .map(|r| MarkerCall {
            rsid: r.rsid,
            variant_copies: r.variant_copies,
            ploidy: r.ploidy,
        })
        .collect();

    let outcome = (def.call)(&calls);

    Some(GeneResult {
        gene: def.gene.to_string(),
        diplotype: outcome.diplotype,
        phenotype: outcome.phenotype,
        drugs: outcome.drugs,
        confidence: outcome.confidence,
        limitation: outcome.limitation,
        markers_used,
    })
}

/// Build the full report, one [`GeneResult`] per gene with at least one
/// genotyped marker, in `GENES` order.
pub fn report(store: &VariantStore) -> Vec<GeneResult> {
    GENES.iter().filter_map(|def| call_gene(store, def)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use gx_core::{Assembly, Genotype, Variant, VariantStore};

    /// Helper: build a store from (rsid, contig, genotype-letters) tuples.
    fn store_with(entries: &[(&str, &str, &str)]) -> VariantStore {
        let mut store = VariantStore::new(Assembly::Grch37);
        for (i, (rsid, contig, call)) in entries.iter().enumerate() {
            store.push(Variant {
                rsid: Some((*rsid).to_string()),
                contig: (*contig).to_string(),
                pos: i as u64 * 100,
                ref_allele: None,
                alt_alleles: vec![],
                genotype: Some(Genotype::from_call(call)),
            });
        }
        store.finalize();
        store
    }

    fn result_for<'a>(results: &'a [GeneResult], gene: &str) -> Option<&'a GeneResult> {
        results.iter().find(|r| r.gene == gene)
    }

    #[test]
    fn cyp2c19_homozygous_star2_is_poor_metabolizer_with_clopidogrel() {
        // rs4244285 A/A = *2/*2 (plus strand A is the *2 allele).
        let store = store_with(&[("rs4244285", "10", "AA")]);
        let results = report(&store);
        let r = result_for(&results, "CYP2C19").expect("CYP2C19 called");
        assert_eq!(r.diplotype, "*2/*2");
        assert_eq!(r.phenotype, "Poor metabolizer");
        assert!(
            r.drugs.iter().any(|d| d.drug == "Clopidogrel"),
            "expected a clopidogrel guidance entry"
        );
        assert_eq!(r.markers_used, vec!["rs4244285".to_string()]);
        assert_eq!(r.limitation, None);
    }

    #[test]
    fn cyp2c19_star1_star17_is_rapid() {
        let store = store_with(&[("rs12248560", "10", "CT")]);
        let r = report(&store);
        let r = result_for(&r, "CYP2C19").unwrap();
        assert_eq!(r.diplotype, "*1/*17");
        assert_eq!(r.phenotype, "Rapid metabolizer");
    }

    #[test]
    fn vkorc1_aa_is_warfarin_sensitive() {
        // Literature -1639A == plus-strand T. Test the literature "A" letter is
        // accepted (strand robustness) and flags high sensitivity.
        let store = store_with(&[("rs9923231", "16", "AA")]);
        let results = report(&store);
        let r = result_for(&results, "VKORC1").expect("VKORC1 called");
        assert!(
            r.phenotype.contains("High warfarin sensitivity"),
            "phenotype was: {}",
            r.phenotype
        );
        assert!(r.drugs.iter().any(|d| d.drug == "Warfarin"));
    }

    #[test]
    fn vkorc1_plus_strand_tt_also_warfarin_sensitive() {
        // Plus-strand T/T is the SAME genotype as literature A/A; both must flag.
        let store = store_with(&[("rs9923231", "16", "TT")]);
        let results = report(&store);
        let r = result_for(&results, "VKORC1").unwrap();
        assert!(r.phenotype.contains("High warfarin sensitivity"));
    }

    #[test]
    fn slco1b1_cc_is_low_function_simvastatin() {
        let store = store_with(&[("rs4149056", "12", "CC")]);
        let results = report(&store);
        let r = result_for(&results, "SLCO1B1").unwrap();
        assert_eq!(r.diplotype, "*5/*5");
        assert!(r.drugs.iter().any(|d| d.drug == "Simvastatin"));
    }

    #[test]
    fn tpmt_3a_heterozygote_is_intermediate() {
        // *3A = *3B (rs1800460 T) + *3C (rs1142345 C) in cis; one copy each het.
        let store = store_with(&[
            ("rs1800460", "6", "CT"),
            ("rs1142345", "6", "TC"),
        ]);
        let results = report(&store);
        let r = result_for(&results, "TPMT").unwrap();
        assert_eq!(r.diplotype, "*1/*3A");
        assert_eq!(r.phenotype, "Intermediate metabolizer");
    }

    #[test]
    fn cyp2d6_carries_limitation() {
        // Any CYP2D6 marker present must yield a non-None limitation.
        let store = store_with(&[("rs3892097", "22", "CT")]);
        let results = report(&store);
        let r = result_for(&results, "CYP2D6").expect("CYP2D6 called");
        assert!(r.limitation.is_some(), "CYP2D6 must carry a limitation");
        assert!(
            r.limitation.as_ref().unwrap().contains("deletion")
                || r.limitation.as_ref().unwrap().contains("duplication"),
            "limitation should mention structural variants"
        );
    }

    #[test]
    fn dpyd_poor_metabolizer_no_function_homozygous() {
        // rs3918290 (*2A) A/A = activity score 0 = poor metabolizer.
        let store = store_with(&[("rs3918290", "1", "AA")]);
        let results = report(&store);
        let r = result_for(&results, "DPYD").unwrap();
        assert_eq!(r.phenotype, "Poor metabolizer");
        assert!(r.limitation.is_some());
    }

    #[test]
    fn no_pgx_markers_returns_empty() {
        let store = store_with(&[
            ("rs9999999", "1", "AA"),
            ("rs1111111", "2", "GG"),
        ]);
        assert!(report(&store).is_empty());
    }

    #[test]
    fn empty_store_returns_empty() {
        let mut store = VariantStore::new(Assembly::Grch38);
        store.finalize();
        assert!(report(&store).is_empty());
    }

    #[test]
    fn no_call_marker_is_skipped() {
        // A genotyped-but-no-call marker should not produce a result.
        let store = store_with(&[("rs4244285", "10", "--")]);
        assert!(result_for(&report(&store), "CYP2C19").is_none());
    }

    #[test]
    fn g6pd_hemizygous_variant_is_deficient() {
        // Male hemizygous: single A allele (literature) at rs1050828.
        let store = store_with(&[("rs1050828", "X", "A")]);
        let results = report(&store);
        let r = result_for(&results, "G6PD").unwrap();
        assert_eq!(r.phenotype, "Deficient");
        assert!(r.drugs.iter().any(|d| d.drug == "Rasburicase"));
        assert!(r.limitation.is_some());
    }

    #[test]
    fn off_target_allele_is_not_informative() {
        // rs4244285 with an unexpected base (not A/G on either strand) -> skipped.
        let store = store_with(&[("rs4244285", "10", "CC")]);
        // C is not in *2 variant {A} nor reference {G}; marker is dropped, so
        // CYP2C19 has no usable markers and is omitted.
        assert!(result_for(&report(&store), "CYP2C19").is_none());
    }
}
