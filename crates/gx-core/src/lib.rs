//! `gx-core` — the shared domain vocabulary for GenomeForge.
//!
//! Everything else (`gx-io`, `gx-annotate`, `gx-app`) speaks in these types:
//! genome [`Assembly`] builds, [`GenomicRange`]s, [`Variant`]s and their
//! [`Genotype`]s, gene-model [`Feature`]s, and the in-memory [`VariantStore`] /
//! [`FeatureStore`] that back the browser's range queries.
//!
//! Coordinate convention: internally all positions are **0-based, half-open**
//! `[start, end)` (like BED/noodles). Conversion to the 1-based inclusive form
//! used by humans, VCF, and Ensembl happens only at the edges (display / REST).

mod assembly;
mod feature;
mod range;
mod sizes;
mod store;
mod variant;

pub use assembly::Assembly;
pub use feature::{Feature, FeatureKind};
pub use range::{GenomicRange, Locus, Strand, normalize_contig};
pub use sizes::{chrom_length, chrom_order};
pub use store::{FeatureStore, VariantStore};
pub use variant::{Genotype, Variant, Zygosity};
