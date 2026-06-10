//! `gx-plasmid` — sequence analysis for the plasmid designer.
//!
//! Pure functions over `&[u8]` DNA (depends only on `gx-core`):
//! reverse-complement / GC ([`seq`]), restriction mapping ([`enzymes`]),
//! ORF finding + translation ([`orf`]), primer Tm ([`primer`]), restriction
//! digest + ligation ([`cloning`]), and DNA assembly ([`assembly`]).

pub mod assembly;
pub mod cloning;
pub mod enzymes;
pub mod orf;
pub mod primer;
pub mod seq;

pub use assembly::{Assembly, Part, fragments_to_parts, gibson, golden_gate};
pub use cloning::{Fragment, FragmentEnd, Overhang, digest, ends_compatible};
pub use enzymes::{ENZYMES, Enzyme, RestrictionSite, find_sites, site_counts, unique_cutters};
pub use orf::{Orf, codon_to_aa, find_orfs, translate};
pub use primer::{PrimerStats, analyze};
pub use seq::{gc_content, reverse_complement};
