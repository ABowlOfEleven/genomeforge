//! Annotation data models and the defensive JSON extraction that turns
//! MyVariant.info / Ensembl REST payloads into them.
//!
//! Real-world MyVariant records are deeply nested and a field can be an object
//! *or* an array of objects depending on how many sources reported it. Rather
//! than model every shape, we keep the full payload in [`VariantAnnotation::raw`]
//! and pull the headline fields with the forgiving [`dig`] walker, which steps
//! into the first element whenever it meets an array.

use serde::{Deserialize, Serialize};
use serde_json::Value;

use gx_core::{Feature, FeatureKind, GenomicRange, Strand};

/// Headline annotation for one variant, suitable for the detail panel + table.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct VariantAnnotation {
    pub rsid: String,
    pub found: bool,
    pub gene: Option<String>,
    pub consequence: Option<String>,
    pub clinical_significance: Option<String>,
    pub conditions: Vec<String>,
    pub gnomad_af: Option<f64>,
    pub cadd_phred: Option<f64>,
    pub dbsnp_ref: Option<String>,
    pub dbsnp_alt: Vec<String>,
    pub gwas_traits: Vec<String>,
    /// The full MyVariant payload, retained for the "raw" view / future fields.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub raw: Option<Value>,
}

impl VariantAnnotation {
    /// A not-found stub for an rsID MyVariant had no record for.
    pub fn not_found(rsid: impl Into<String>) -> Self {
        Self {
            rsid: rsid.into(),
            found: false,
            ..Default::default()
        }
    }

    /// True if this variant has a clinically meaningful ClinVar call.
    pub fn is_clinically_notable(&self) -> bool {
        self.clinical_significance
            .as_deref()
            .map(|s| {
                let s = s.to_ascii_lowercase();
                s.contains("pathogenic") || s.contains("risk") || s.contains("drug")
            })
            .unwrap_or(false)
    }

    /// Build from one element of a MyVariant batch response.
    pub fn from_myvariant(query: &str, hit: &Value) -> Self {
        if hit.get("notfound").and_then(Value::as_bool) == Some(true) {
            return Self::not_found(query);
        }

        let gene = dig_str(hit, &["dbnsfp", "genename"])
            .or_else(|| dig_str(hit, &["dbsnp", "gene", "symbol"]))
            .or_else(|| dig_str(hit, &["snpeff", "ann", "genename"]))
            .or_else(|| dig_str(hit, &["clinvar", "gene", "symbol"]))
            .or_else(|| dig_str(hit, &["cadd", "gene", "genename"]));

        let consequence = dig_str(hit, &["snpeff", "ann", "effect"])
            .or_else(|| dig_str(hit, &["cadd", "consequence"]))
            .or_else(|| dig_str(hit, &["dbnsfp", "ensembl", "consequence"]));

        let clinical_significance = dig_str(hit, &["clinvar", "rcv", "clinical_significance"])
            .or_else(|| dig_str(hit, &["clinvar", "clinical_significance"]));

        let conditions = collect_strings(dig_raw(hit, &["clinvar", "rcv", "conditions", "name"]));

        let gnomad_af = dig_f64(hit, &["gnomad_genome", "af", "af"])
            .or_else(|| dig_f64(hit, &["gnomad_exome", "af", "af"]))
            .or_else(|| dig_f64(hit, &["dbnsfp", "gnomad_genomes", "af"]));

        let cadd_phred = dig_f64(hit, &["cadd", "phred"]);

        let dbsnp_ref = dig_str(hit, &["dbsnp", "ref"]).or_else(|| dig_str(hit, &["vcf", "ref"]));
        let dbsnp_alt = collect_strings(dig_raw(hit, &["dbsnp", "alt"]))
            .into_iter()
            .chain(collect_strings(dig_raw(hit, &["vcf", "alt"])))
            .collect();

        let gwas_traits = collect_strings(dig_raw(hit, &["gwassnps", "trait"]));

        Self {
            rsid: query.to_string(),
            found: true,
            gene,
            consequence,
            clinical_significance,
            conditions,
            gnomad_af,
            cadd_phred,
            dbsnp_ref,
            dbsnp_alt,
            gwas_traits,
            raw: Some(hit.clone()),
        }
    }
}

/// A gene's location, from an Ensembl `lookup/symbol` response.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneLocation {
    pub id: String,
    pub symbol: Option<String>,
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub biotype: Option<String>,
}

impl GeneLocation {
    pub fn range(&self) -> GenomicRange {
        GenomicRange::new(self.contig.clone(), self.start, self.end, self.strand)
    }

    /// Parse an Ensembl `lookup/symbol/homo_sapiens/{sym}` JSON object.
    pub fn from_ensembl(v: &Value) -> Option<Self> {
        let id = v.get("id")?.as_str()?.to_string();
        let contig = gx_core::normalize_contig(v.get("seq_region_name")?.as_str()?);
        // Ensembl coordinates are 1-based inclusive -> 0-based half-open.
        let start = v.get("start")?.as_u64()?.saturating_sub(1);
        let end = v.get("end")?.as_u64()?;
        let strand = Strand::from_sign(v.get("strand").and_then(Value::as_i64).unwrap_or(0));
        Some(Self {
            id,
            symbol: v.get("display_name").and_then(Value::as_str).map(str::to_string),
            contig,
            start,
            end,
            strand,
            biotype: v.get("biotype").and_then(Value::as_str).map(str::to_string),
        })
    }
}

/// Parse one feature object from an Ensembl `overlap/region` array.
pub fn feature_from_ensembl(v: &Value) -> Option<Feature> {
    let contig = gx_core::normalize_contig(v.get("seq_region_name")?.as_str()?);
    let start = v.get("start")?.as_u64()?.saturating_sub(1);
    let end = v.get("end")?.as_u64()?;
    let strand = Strand::from_sign(v.get("strand").and_then(Value::as_i64).unwrap_or(0));
    let kind = FeatureKind::from_ensembl(v.get("feature_type").and_then(Value::as_str).unwrap_or(""));
    let id = v
        .get("id")
        .or_else(|| v.get("gene_id"))
        .or_else(|| v.get("transcript_id"))
        .or_else(|| v.get("exon_id"))
        .and_then(Value::as_str)
        .unwrap_or("")
        .to_string();
    let name = v
        .get("external_name")
        .and_then(Value::as_str)
        .map(str::to_string);
    let parent = v.get("Parent").and_then(Value::as_str).map(str::to_string);
    let biotype = v.get("biotype").and_then(Value::as_str).map(str::to_string);

    Some(Feature {
        id,
        name,
        kind,
        range: GenomicRange::new(contig, start, end, strand),
        parent,
        biotype,
    })
}

// ---- defensive JSON helpers ------------------------------------------------

/// If `v` is an array, the first element; otherwise `v` itself.
fn first(v: &Value) -> &Value {
    match v.as_array() {
        Some(arr) => arr.first().unwrap_or(&Value::Null),
        None => v,
    }
}

/// Walk `path` from `v`, stepping into the first element of any array met.
fn dig<'a>(v: &'a Value, path: &[&str]) -> Option<&'a Value> {
    let mut cur = v;
    for key in path {
        cur = first(cur).get(*key)?;
    }
    Some(first(cur))
}

/// Like [`dig`] but returns the terminal node as-is (no final `first()`), so a
/// list-valued leaf is preserved for [`collect_strings`] (e.g. multiple ClinVar
/// condition names).
fn dig_raw<'a>(v: &'a Value, path: &[&str]) -> Option<&'a Value> {
    let mut cur = v;
    for key in path {
        cur = first(cur).get(*key)?;
    }
    Some(cur)
}

fn dig_str(v: &Value, path: &[&str]) -> Option<String> {
    dig(v, path).and_then(|x| x.as_str().map(str::to_string))
}

fn dig_f64(v: &Value, path: &[&str]) -> Option<f64> {
    dig(v, path).and_then(Value::as_f64)
}

/// Collect string value(s) at a node that may be a string or array of strings.
fn collect_strings(v: Option<&Value>) -> Vec<String> {
    match v {
        Some(Value::String(s)) => vec![s.clone()],
        Some(Value::Array(arr)) => arr
            .iter()
            .filter_map(|x| x.as_str().map(str::to_string))
            .collect(),
        _ => Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn extracts_clinvar_and_gnomad_from_nested_shapes() {
        // `rcv` as an array; conditions name as object — the messy real shape.
        let hit = json!({
            "_id": "rs397507444",
            "dbsnp": { "gene": { "symbol": "CFTR" }, "rsid": "rs397507444" },
            "clinvar": {
                "rcv": [
                    { "clinical_significance": "Pathogenic",
                      "conditions": { "name": "Cystic fibrosis" } }
                ]
            },
            "gnomad_genome": { "af": { "af": 0.0001 } },
            "cadd": { "phred": 25.3, "consequence": "NON_SYNONYMOUS" }
        });
        let ann = VariantAnnotation::from_myvariant("rs397507444", &hit);
        assert!(ann.found);
        assert_eq!(ann.gene.as_deref(), Some("CFTR"));
        assert_eq!(ann.clinical_significance.as_deref(), Some("Pathogenic"));
        assert_eq!(ann.conditions, ["Cystic fibrosis"]);
        assert_eq!(ann.gnomad_af, Some(0.0001));
        assert_eq!(ann.cadd_phred, Some(25.3));
        assert!(ann.is_clinically_notable());
    }

    #[test]
    fn handles_not_found() {
        let hit = json!({ "query": "rsZZZ", "notfound": true });
        let ann = VariantAnnotation::from_myvariant("rsZZZ", &hit);
        assert!(!ann.found);
        assert!(ann.gene.is_none());
    }

    #[test]
    fn parses_ensembl_gene_lookup() {
        let v = json!({
            "id": "ENSG00000012048",
            "display_name": "BRCA1",
            "seq_region_name": "17",
            "start": 43044295,
            "end": 43125483,
            "strand": -1,
            "biotype": "protein_coding"
        });
        let g = GeneLocation::from_ensembl(&v).unwrap();
        assert_eq!(g.symbol.as_deref(), Some("BRCA1"));
        assert_eq!(g.contig, "17");
        assert_eq!(g.start, 43044294); // 1-based -> 0-based
        assert_eq!(g.strand, Strand::Reverse);
    }
}
