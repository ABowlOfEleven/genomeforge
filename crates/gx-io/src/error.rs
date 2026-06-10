use std::path::PathBuf;

use thiserror::Error;

/// Errors raised while importing genome / sequence files.
#[derive(Debug, Error)]
pub enum IoError {
    #[error("I/O error reading {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("could not determine a supported file type for {0}")]
    UnknownFormat(PathBuf),

    #[error("the file appears empty or has no parseable records: {0}")]
    Empty(PathBuf),

    #[error("parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("VCF error: {0}")]
    Vcf(String),

    #[error("FASTA error: {0}")]
    Fasta(String),

    #[error("GenBank error: {0}")]
    GenBank(String),
}

pub type Result<T> = std::result::Result<T, IoError>;
