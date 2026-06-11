//! `gx-pgx` - pharmacogenomics (PGx) reporting from a consumer genome.
//!
//! Given a [`gx_core::VariantStore`] imported from a consumer array (23andMe,
//! AncestryDNA) or VCF, this crate performs SNP-based star-allele / diplotype
//! calling for a curated set of pharmacogenes, maps each diplotype to a CPIC
//! metabolizer phenotype, and attaches concise, CPIC-derived drug-response
//! guidance. It is a pure library: no UI, no I/O, no network.
//!
//! # This is NOT medical advice
//!
//! Every result this crate produces is informational only. It is NOT a clinical
//! test, it is NOT validated for diagnosis, and NOTHING here should be used to
//! start, stop, or change any medication or dose. Consumer-array genotypes are
//! not clinical-grade, cover only a small subset of the variants that define each
//! gene, and CANNOT detect the structural variants (gene deletions, duplications,
//! hybrids) that are clinically important for several pharmacogenes. Decisions
//! about medication must be made by a qualified clinician using validated testing.
//! The integrating UI must surface this framing prominently.
//!
//! # Data sources
//!
//! Allele-defining rsIDs, allele functions, diplotype-to-phenotype mappings, and
//! drug guidance are derived from the Clinical Pharmacogenetics Implementation
//! Consortium (CPIC) guidelines and allele-definition / functionality tables
//! (cpicpgx.org), PharmGKB, and dbSNP for strand/allele orientation. Specific
//! citations appear in `genes.rs` next to each gene's definition. Where a gene is
//! not reliably callable from SNP data it is either omitted or returned with an
//! explicit [`GeneResult::limitation`].
//!
//! # Strand handling
//!
//! Consumer arrays report alleles on the plus (forward) strand of the reference.
//! The allele letters used in PGx literature are sometimes quoted on the opposite
//! strand (most notably VKORC1 rs9923231, quoted as -1639G>A in the literature but
//! C>T on the plus strand). Each defining marker in `genes.rs` lists the plus-
//! strand allele AND, where a different convention is common, its complement, and
//! the caller matches either. See `call.rs` for the matching logic.

mod call;
mod genes;
mod model;

pub use model::{Confidence, DrugGuidance, GeneResult};

/// Produce a pharmacogenomics report from a sample's variants.
///
/// Returns one [`GeneResult`] per pharmacogene that had at least one defining
/// marker genotyped (and called, not a no-call) in `store`. Genes with no usable
/// markers are omitted entirely, so an empty input (or a sample with no PGx
/// markers) yields an empty `Vec`.
///
/// Results are returned in a stable, deterministic order (the order the genes are
/// defined in `genes.rs`).
pub fn report(store: &gx_core::VariantStore) -> Vec<GeneResult> {
    call::report(store)
}
