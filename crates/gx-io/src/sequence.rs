//! Sequence importers: FASTA (`noodles-fasta`) and GenBank (`gb-io`).
//!
//! Both produce [`SequenceRecord`]s — annotated linear or circular molecules
//! that feed the sequence/plasmid view. GenBank brings its feature table along;
//! FASTA is bare sequence.

use std::path::Path;

use gb_io::reader::SeqReader;
use gb_io::seq::Topology;
use noodles::fasta;

use gx_core::{Feature, FeatureKind, GenomicRange, Strand};

use crate::error::{IoError, Result};
use crate::raw::open_maybe_gzip;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFormat {
    Fasta,
    GenBank,
}

/// An imported nucleotide (or protein) sequence with optional annotations.
#[derive(Debug, Clone)]
pub struct SequenceRecord {
    pub name: String,
    pub description: Option<String>,
    pub seq: Vec<u8>,
    pub circular: bool,
    pub features: Vec<Feature>,
    pub source_format: SeqFormat,
}

impl SequenceRecord {
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// Load one or more sequences from a FASTA file (gzip transparently handled).
pub fn load_fasta(path: &Path) -> Result<Vec<SequenceRecord>> {
    let mut reader = fasta::io::Reader::new(open_maybe_gzip(path)?);

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result.map_err(|e| IoError::Fasta(e.to_string()))?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let seq = record.sequence().as_ref().to_vec();
        records.push(SequenceRecord {
            name,
            description: None,
            seq,
            circular: false,
            features: Vec::new(),
            source_format: SeqFormat::Fasta,
        });
    }

    if records.is_empty() {
        return Err(IoError::Empty(path.to_path_buf()));
    }
    Ok(records)
}

/// Load one or more sequences (with feature tables) from a GenBank file
/// (gzip transparently handled).
pub fn load_genbank(path: &Path) -> Result<Vec<SequenceRecord>> {
    let mut records = Vec::new();
    for result in SeqReader::new(open_maybe_gzip(path)?) {
        let seq = result.map_err(|e| IoError::GenBank(e.to_string()))?;
        let name = seq.name.clone().unwrap_or_else(|| "unnamed".to_string());
        let circular = matches!(seq.topology, Topology::Circular);

        let features = seq
            .features
            .iter()
            .filter_map(|f| genbank_feature(&name, f))
            .collect();

        records.push(SequenceRecord {
            name,
            description: seq.definition.clone(),
            seq: seq.seq,
            circular,
            features,
            source_format: SeqFormat::GenBank,
        });
    }

    if records.is_empty() {
        return Err(IoError::Empty(path.to_path_buf()));
    }
    Ok(records)
}

/// Convert a gb-io feature into a `gx-core` [`Feature`] on the record's own
/// contig. gb-io's `find_bounds` already returns 0-based half-open `(start, end)`
/// (GenBank `1..6` parses to `(0, 6)`), matching our internal convention.
fn genbank_feature(contig: &str, f: &gb_io::seq::Feature) -> Option<Feature> {
    let (start, end) = f.location.find_bounds().ok()?;
    let start = start.max(0) as u64;
    let end = end.max(0) as u64;
    let label = f
        .qualifier_values("label")
        .next()
        .or_else(|| f.qualifier_values("gene").next())
        .or_else(|| f.qualifier_values("product").next())
        .map(|s| s.to_string());

    Some(Feature {
        id: label.clone().unwrap_or_else(|| f.kind.to_string()),
        name: label,
        kind: FeatureKind::from_ensembl(f.kind.as_ref()),
        range: GenomicRange::new(contig.to_string(), start, end, Strand::Unknown),
        parent: None,
        biotype: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    fn write_temp(name: &str, contents: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!("gxio_test_{name}"));
        let mut f = File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    #[test]
    fn reads_fasta() {
        let path = write_temp("seq.fasta", ">seq1 a description\nACGTACGT\nACGT\n");
        let records = load_fasta(&path).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "seq1");
        assert_eq!(records[0].seq, b"ACGTACGTACGT");
        assert!(!records[0].circular);
    }

    const GENBANK: &str = "\
LOCUS       TESTSEQ                   12 bp    DNA     circular SYN 01-JAN-2020
DEFINITION  A tiny test plasmid.
FEATURES             Location/Qualifiers
     gene            1..6
                     /gene=\"abc\"
ORIGIN
        1 acgtacgtac gt
//
";

    #[test]
    fn reads_genbank() {
        let path = write_temp("plasmid.gb", GENBANK);
        let records = load_genbank(&path).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq.len(), 12);
        assert!(records[0].circular);
        assert!(!records[0].features.is_empty());
        assert_eq!(records[0].features[0].name.as_deref(), Some("abc"));
        // GenBank "1..6" -> 0-based half-open [0, 6) (6 bases), no off-by-one.
        assert_eq!(records[0].features[0].range.start, 0);
        assert_eq!(records[0].features[0].range.end, 6);
    }
}
