//! Importer for consumer DNA "raw data" files (23andMe, AncestryDNA).
//!
//! These are simple tab-delimited tables of genotyped SNPs. The two major
//! vendors differ in shape:
//!
//! * **23andMe** — 4 columns `rsid / chromosome / position / genotype`, with a
//!   commented (`#`) preamble; chromosomes are `1..22, X, Y, MT`.
//! * **AncestryDNA** — 5 columns `rsid / chromosome / position / allele1 /
//!   allele2`, with a commented preamble and a (non-comment) column-header row;
//!   chromosomes are numeric, where `23=X, 24=Y, 25=PAR(→X), 26=MT`.
//!
//! Positions in both formats are **1-based**; we convert to the crate-internal
//! 0-based convention on the way in. Consumer arrays are almost always GRCh37,
//! so that is the default unless the preamble says otherwise — and rsID-based
//! annotation is build-independent regardless.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use gx_core::{Assembly, Genotype, Variant, VariantStore, normalize_contig};

use crate::error::{IoError, Result};

/// Detected vendor format of a consumer raw file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RawFormat {
    TwentyThreeAndMe,
    AncestryDna,
    /// Recognised as tabular SNP data but vendor not positively identified.
    Generic,
}

impl RawFormat {
    pub fn label(self) -> &'static str {
        match self {
            RawFormat::TwentyThreeAndMe => "23andMe raw data",
            RawFormat::AncestryDna => "AncestryDNA raw data",
            RawFormat::Generic => "Consumer raw data",
        }
    }
}

/// Result of importing a variant-bearing file.
pub struct VariantImport {
    pub store: VariantStore,
    /// Build we believe the coordinates use (default GRCh37 for consumer raw).
    pub assembly: Assembly,
    pub format_label: String,
    /// Non-fatal notes for the user (e.g. assumed build, skipped lines).
    pub warnings: Vec<String>,
}

/// Open a path, transparently decompressing `*.gz`.
pub(crate) fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|source| IoError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let is_gz = path
        .extension()
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    if is_gz {
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)) as Box<dyn BufRead>)
    } else {
        Ok(Box::new(BufReader::new(file)) as Box<dyn BufRead>)
    }
}

/// Load a 23andMe / AncestryDNA raw file into a [`VariantStore`].
pub fn load_raw(path: &Path) -> Result<VariantImport> {
    let reader = open_maybe_gzip(path)?;
    let mut import = parse_raw(reader)?;
    if import.store.is_empty() {
        return Err(IoError::Empty(path.to_path_buf()));
    }
    import.store.finalize();
    Ok(import)
}

/// Map a chromosome token to a normalised contig name, handling AncestryDNA's
/// numeric encoding (`23=X, 24=Y, 25=PAR→X, 26=MT`).
fn map_chrom(token: &str) -> String {
    match token {
        "23" => "X".to_string(),
        "24" => "Y".to_string(),
        "25" => "X".to_string(), // pseudoautosomal region, reported on X
        "26" => "MT".to_string(),
        other => normalize_contig(other),
    }
}

/// Sniff vendor and build from the commented preamble.
fn classify_preamble(comments: &str) -> (Option<RawFormat>, Option<Assembly>) {
    let lower = comments.to_ascii_lowercase();
    let vendor = if lower.contains("ancestrydna") || lower.contains("ancestry.com") {
        Some(RawFormat::AncestryDna)
    } else if lower.contains("23andme") {
        Some(RawFormat::TwentyThreeAndMe)
    } else {
        None
    };
    let assembly = if lower.contains("grch38") || lower.contains("build 38") {
        Some(Assembly::Grch38)
    } else if lower.contains("grch37") || lower.contains("build 37") || lower.contains("build 36") {
        Some(Assembly::Grch37)
    } else {
        None
    };
    (vendor, assembly)
}

fn parse_raw<R: BufRead>(reader: R) -> Result<VariantImport> {
    let mut comments = String::new();
    let mut vendor: Option<RawFormat> = None;
    let mut warnings = Vec::new();
    let mut skipped = 0usize;

    // We don't know the build until the preamble is parsed, so collect rows
    // first and construct the store once.
    let mut rows: Vec<Variant> = Vec::new();

    for (idx, line) in reader.lines().enumerate() {
        let line = line.map_err(|source| IoError::Io {
            path: Path::new("<stream>").to_path_buf(),
            source,
        })?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        if let Some(comment) = line.strip_prefix('#') {
            comments.push_str(comment);
            comments.push('\n');
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        // AncestryDNA's column-header row is a non-comment line. Match case-
        // insensitively, and only treat it as a header (not data) when other
        // header tokens are present too.
        let lower = line.to_ascii_lowercase();
        if fields.first().map(|s| s.eq_ignore_ascii_case("rsid")).unwrap_or(false)
            && lower.contains("position")
        {
            if lower.contains("allele1") {
                vendor.get_or_insert(RawFormat::AncestryDna);
            }
            continue;
        }

        if fields.len() < 4 {
            skipped += 1;
            continue;
        }

        let rsid = fields[0].trim();
        let chrom = map_chrom(fields[1].trim());
        let Ok(pos1) = fields[2].trim().parse::<u64>() else {
            skipped += 1;
            continue;
        };

        let genotype = if fields.len() >= 5 {
            // AncestryDNA: separate allele columns.
            vendor.get_or_insert(RawFormat::AncestryDna);
            let combined = format!("{}{}", fields[3].trim(), fields[4].trim());
            Genotype::from_call(&combined)
        } else {
            vendor.get_or_insert(RawFormat::TwentyThreeAndMe);
            Genotype::from_call(fields[3].trim())
        };

        let rsid = if rsid.starts_with("rs") || rsid.starts_with("i") {
            Some(rsid.to_string())
        } else {
            None
        };

        rows.push(Variant {
            rsid,
            contig: chrom,
            pos: pos1.saturating_sub(1), // 1-based file -> 0-based internal
            ref_allele: None,
            alt_alleles: vec![],
            genotype: Some(genotype),
        });
        let _ = idx;
    }

    let (pre_vendor, pre_assembly) = classify_preamble(&comments);
    let format = pre_vendor.or(vendor).unwrap_or(RawFormat::Generic);
    let assembly_final = pre_assembly.unwrap_or(Assembly::Grch37);

    if pre_assembly.is_none() {
        warnings.push(format!(
            "Assembly build not stated in file; assuming {}. rsID annotation is \
             build-independent, but verify the build for correct browser coordinates.",
            assembly_final.label()
        ));
    }
    if skipped > 0 {
        warnings.push(format!("Skipped {skipped} unparseable line(s)."));
    }

    let mut store = VariantStore::new(assembly_final);
    for v in rows {
        store.push(v);
    }

    Ok(VariantImport {
        store,
        assembly: assembly_final,
        format_label: format.label().to_string(),
        warnings,
    })
}

/// Read up to the first 4 KiB of a (possibly gzipped) file as lossy UTF-8.
/// Cheap, for format sniffing during dispatch.
pub(crate) fn peek_head(path: &Path) -> String {
    let Ok(mut reader) = open_maybe_gzip(path) else {
        return String::new();
    };
    let mut head = [0u8; 4096];
    let n = reader.read(&mut head).unwrap_or(0);
    String::from_utf8_lossy(&head[..n]).into_owned()
}

/// Heuristic: does this look like a consumer raw file (vs VCF/FASTA/GenBank)?
pub fn looks_like_raw(path: &Path) -> bool {
    let text = peek_head(path).to_ascii_lowercase();
    let trimmed = text.trim_start();
    if trimmed.starts_with('>') || text.contains("##fileformat=vcf") || trimmed.starts_with("locus") {
        return false;
    }
    text.contains("23andme")
        || text.contains("ancestrydna")
        || text.contains("rsid\tchromosome\tposition")
        || text.contains("\trs")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const TWENTYTHREE: &str = "\
# This data file generated by 23andMe at: build 37
# rsid\tchromosome\tposition\tgenotype
rs4477212\t1\t82154\tAA
rs3094315\t1\t752566\tAG
rsIndel\t1\t800000\tDI
rsX\tX\t2700000\tT
i5000940\tMT\t16500\t--
";

    const ANCESTRY: &str = "\
#AncestryDNA raw data download
#This file was generated by AncestryDNA
rsid\tchromosome\tposition\tallele1\tallele2
rs4477212\t1\t82154\tA\tA
rs3094315\t1\t752566\tA\tG
rsY\t24\t2700000\tT\tT
rsMT\t26\t16500\tG\tG
";

    #[test]
    fn parses_23andme() {
        let imp = parse_raw(Cursor::new(TWENTYTHREE)).unwrap();
        let mut store = imp.store;
        store.finalize();
        assert_eq!(store.len(), 5);
        assert_eq!(imp.format_label, "23andMe raw data");
        assert_eq!(imp.assembly, Assembly::Grch37);
        // 1-based 82154 -> 0-based 82153
        let v = store.find_by_rsid("rs4477212").unwrap();
        assert_eq!(v.pos, 82153);
        assert_eq!(v.genotype.as_ref().unwrap().display(), "A/A");
        // Indel encoding survives.
        let indel = store.find_by_rsid("rsIndel").unwrap();
        assert_eq!(indel.genotype.as_ref().unwrap().alleles, ["D", "I"]);
        // No-call retained.
        assert!(
            store
                .find_by_rsid("i5000940")
                .unwrap()
                .genotype
                .as_ref()
                .unwrap()
                .is_no_call()
        );
    }

    #[test]
    fn parses_ancestry_with_numeric_chromosomes() {
        let imp = parse_raw(Cursor::new(ANCESTRY)).unwrap();
        let mut store = imp.store;
        store.finalize();
        assert_eq!(imp.format_label, "AncestryDNA raw data");
        // 24 -> Y, 26 -> MT
        assert_eq!(store.find_by_rsid("rsY").unwrap().contig, "Y");
        assert_eq!(store.find_by_rsid("rsMT").unwrap().contig, "MT");
        assert_eq!(
            store.find_by_rsid("rs3094315").unwrap().genotype.as_ref().unwrap().display(),
            "A/G"
        );
    }
}
