//! `gx-crispr` — CRISPR/Cas9 guide design, off-target search, edit simulation,
//! and AAV cargo planning. Pure functions over `&[u8]` (depends on `gx-core`
//! and `gx-plasmid`).
//!
//! Scoring methods are documented and intentionally transparent; see
//! [`score`]. They are decision-support estimates, not a clinical pipeline.

pub mod aav;
pub mod doench;
pub mod edit;
pub mod guides;
pub mod offtarget;
pub mod score;

pub use aav::AavCargo;
pub use edit::{EditPreview, hdr_replace, nhej_deletion};
pub use guides::{Guide, PROTOSPACER_LEN, find_guides};
pub use offtarget::{OffTarget, find_matches, find_offtargets, specificity};
pub use score::{cfd_score, on_target_estimate};
