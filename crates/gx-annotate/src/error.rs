use thiserror::Error;

/// Errors from the annotation layer (network, parsing, cache).
#[derive(Debug, Error)]
pub enum AnnotateError {
    #[error("HTTP request failed: {0}")]
    Http(String),

    #[error("unexpected response shape: {0}")]
    Response(String),

    #[error("JSON error: {0}")]
    Json(String),

    #[error("cache/database error: {0}")]
    Db(String),

    #[error("offline: {0} is not in the local cache")]
    OfflineMiss(String),
}

pub type Result<T> = std::result::Result<T, AnnotateError>;

/// Map any displayable backend error (redb's several error types) into [`AnnotateError::Db`].
pub(crate) fn db_err<E: std::fmt::Display>(e: E) -> AnnotateError {
    AnnotateError::Db(e.to_string())
}

impl From<serde_json::Error> for AnnotateError {
    fn from(e: serde_json::Error) -> Self {
        AnnotateError::Json(e.to_string())
    }
}

impl From<ureq::Error> for AnnotateError {
    fn from(e: ureq::Error) -> Self {
        AnnotateError::Http(e.to_string())
    }
}
