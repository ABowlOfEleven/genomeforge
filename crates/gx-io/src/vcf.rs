//! VCF importer (Variant Call Format), backed by `noodles-vcf`.
//!
//! Reads sites + rsIDs + REF/ALT into a [`VariantStore`]. Per-sample genotype
//! (the `GT` field) extraction is deliberately left for a later pass — Phase 1
//! annotates and browses by rsID / locus, which is build-independent.
//!
//! VCFs from modern sequencing are usually GRCh38, so that is the default build
//! (surfaced as a warning so the user can correct it).

use std::path::Path;

use noodles::vcf;
// Bring the `iter()` methods on the ID / alternate-allele fields into scope.
// These traits share names with the concrete field structs, so import as `_`.
use noodles::vcf::variant::record::{AlternateBases as _, Ids as _};

use gx_core::{Assembly, Variant, VariantStore, normalize_contig};

use crate::error::{IoError, Result};
use crate::raw::VariantImport;

/// Convert one parsed VCF record into a [`Variant`], or `None` if it lacks a
/// usable position. Kept reader-agnostic so it can be unit-tested directly.
fn record_to_variant(record: &vcf::Record) -> Option<Variant> {
    let contig = normalize_contig(record.reference_sequence_name());
    let pos = match record.variant_start()? {
        Ok(p) => (p.get() as u64).saturating_sub(1), // 1-based -> 0-based
        Err(_) => return None,
    };
    // The ID field can hold several `;`-separated identifiers; keep the first
    // genuine dbSNP rsID (others, e.g. COSMIC ids, aren't dbSNP-queryable).
    let rsid = record
        .ids()
        .iter()
        .find(|s| {
            s.strip_prefix("rs")
                .map(|r| !r.is_empty() && r.bytes().all(|b| b.is_ascii_digit()))
                .unwrap_or(false)
        })
        .map(str::to_string);
    let ref_bases = record.reference_bases();
    let ref_allele = (!ref_bases.is_empty()).then(|| ref_bases.to_string());
    let alt_alleles = record
        .alternate_bases()
        .iter()
        .filter_map(|r| r.ok())
        .map(|s| s.to_string())
        .collect();

    Some(Variant {
        rsid,
        contig,
        pos,
        ref_allele,
        alt_alleles,
        genotype: None,
    })
}

/// Load a `.vcf` / `.vcf.gz` file into a [`VariantImport`].
pub fn load_vcf(path: &Path) -> Result<VariantImport> {
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(path)
        .map_err(|source| IoError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    let _header = reader
        .read_header()
        .map_err(|e| IoError::Vcf(format!("header: {e}")))?;

    let assembly = Assembly::Grch38;
    let mut store = VariantStore::new(assembly);
    let mut skipped = 0usize;

    for result in reader.records() {
        let record = result.map_err(|e| IoError::Vcf(e.to_string()))?;
        match record_to_variant(&record) {
            Some(v) => store.push(v),
            None => skipped += 1,
        }
    }

    if store.is_empty() {
        return Err(IoError::Empty(path.to_path_buf()));
    }
    store.finalize();

    let mut warnings = vec![format!(
        "Assembly build not detected from VCF; assuming {}. Override it if your \
         calls were made against a different reference.",
        assembly.label()
    )];
    if skipped > 0 {
        warnings.push(format!("Skipped {skipped} record(s) without a position."));
    }

    Ok(VariantImport {
        store,
        assembly,
        format_label: "VCF".to_string(),
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const VCF: &str = "\
##fileformat=VCFv4.3
##contig=<ID=1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\trs123\tA\tG\t.\t.\t.
chr1\t200\t.\tC\tT,A\t.\t.\t.
";

    #[test]
    fn converts_records() {
        let mut reader = vcf::io::Reader::new(Cursor::new(VCF.as_bytes()));
        reader.read_header().unwrap();
        let variants: Vec<Variant> = reader
            .records()
            .map(|r| record_to_variant(&r.unwrap()).unwrap())
            .collect();

        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].rsid.as_deref(), Some("rs123"));
        assert_eq!(variants[0].contig, "1");
        assert_eq!(variants[0].pos, 99); // 1-based 100 -> 0-based 99
        assert_eq!(variants[0].ref_allele.as_deref(), Some("A"));
        assert_eq!(variants[0].alt_alleles, ["G"]);

        assert_eq!(variants[1].rsid, None);
        assert_eq!(variants[1].contig, "1"); // chr1 normalised
        assert_eq!(variants[1].alt_alleles, ["T", "A"]);
    }
}
