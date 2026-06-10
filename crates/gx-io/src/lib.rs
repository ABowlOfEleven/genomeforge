//! `gx-io` — importers that turn genome / sequence files into `gx-core` types.
//!
//! * [`load_raw`] — 23andMe / AncestryDNA consumer raw data → [`VariantImport`]
//! * [`load_vcf`] — VCF (incl. bgzipped) → [`VariantImport`]
//! * [`load_fasta`] / [`load_genbank`] — sequences → [`SequenceRecord`]
//!
//! [`import_path`] dispatches by extension, falling back to content sniffing.

mod error;
mod raw;
mod sequence;
mod vcf;

use std::path::Path;

pub use error::{IoError, Result};
pub use raw::{RawFormat, VariantImport, load_raw, looks_like_raw};
pub use sequence::{SeqFormat, SequenceRecord, load_fasta, load_genbank};
pub use vcf::load_vcf;

/// What an import yielded: sample variants, or one/more sequences.
pub enum Imported {
    Variants(VariantImport),
    Sequences(Vec<SequenceRecord>),
}

/// Import any supported file, choosing the parser by extension and, when the
/// extension is ambiguous (`.txt`, none), by sniffing the leading bytes.
pub fn import_path(path: &Path) -> Result<Imported> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase)
        .unwrap_or_default();
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    match ext.as_str() {
        "vcf" | "bcf" => return Ok(Imported::Variants(load_vcf(path)?)),
        "fa" | "fasta" | "fna" | "ffn" | "faa" => {
            return Ok(Imported::Sequences(load_fasta(path)?));
        }
        "gb" | "gbk" | "gbff" | "genbank" => {
            return Ok(Imported::Sequences(load_genbank(path)?));
        }
        "gz" if name.ends_with(".vcf.gz") => return Ok(Imported::Variants(load_vcf(path)?)),
        _ => {}
    }

    // Ambiguous extension: sniff content.
    sniff_and_load(path)
}

fn sniff_and_load(path: &Path) -> Result<Imported> {
    let head = raw::peek_head(path);
    let trimmed = head.trim_start();
    if trimmed.starts_with('>') {
        Ok(Imported::Sequences(load_fasta(path)?))
    } else if head.to_ascii_lowercase().contains("##fileformat=vcf") {
        Ok(Imported::Variants(load_vcf(path)?))
    } else if trimmed.to_ascii_uppercase().starts_with("LOCUS") {
        Ok(Imported::Sequences(load_genbank(path)?))
    } else {
        // Default: consumer raw (handles the common `.txt` case).
        Ok(Imported::Variants(load_raw(path)?))
    }
}
