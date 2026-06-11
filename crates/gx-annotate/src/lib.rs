//! `gx-annotate` — hybrid (online + cached) variant/gene annotation.
//!
//! [`AnnotationService`] is the entry point: cache-first lookups backed by
//! [`MyVariantClient`] (variant annotation) and [`EnsemblClient`] (sequence,
//! gene models, symbol lookup), persisted in a SQLite [`Cache`]. It is owned by
//! the app's worker thread; the UI talks to it over a channel.

mod annotation;
mod cache;
mod clients;
mod error;
mod literature;
mod pgs;
mod service;
mod update;

pub use annotation::{GeneLocation, VariantAnnotation, feature_from_ensembl};
pub use cache::{Cache, CacheStats};
pub use clients::{EnsemblClient, MyVariantClient};
pub use error::{AnnotateError, Result};
pub use literature::{Article, PubMedClient};
pub use pgs::PgsClient;
pub use service::AnnotationService;
pub use update::{UpdateInfo, is_newer, latest_release};
