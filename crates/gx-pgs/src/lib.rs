//! `gx-pgs` — parse PGS Catalog scoring files and apply them to a sample's
//! genotypes to compute polygenic scores (with an optional, clearly-caveated
//! population percentile). Depends only on `gx-core`.

mod model;
mod parse;
mod score;

pub use model::{ScoreFile, ScoreVariant};
pub use parse::parse;
pub use score::{PrsResult, apply, normal_cdf};
