//! Parser for PGS Catalog scoring files.
//!
//! Format: `#key=value` metadata lines, then a tab-separated column header, then
//! data rows. Column sets vary; we locate columns by name and prefer the
//! harmonized (`hm_*`) rsID / position when present.

use std::collections::HashMap;

use crate::model::{ScoreFile, ScoreVariant};

/// Treat empty / `.` fields as absent.
fn nonempty(s: Option<&str>) -> Option<&str> {
    s.filter(|x| !x.is_empty() && *x != ".")
}

pub fn parse(text: &str) -> Result<ScoreFile, String> {
    let mut id = String::new();
    let mut trait_name = String::new();
    let mut build = None;
    let mut declared_count = None;
    let mut cols: Option<HashMap<String, usize>> = None;
    let mut variants = Vec::new();

    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('#') {
            // metadata `key=value` (also handles ## / ### section lines, ignored)
            let rest = rest.trim_start_matches('#');
            if let Some((k, v)) = rest.split_once('=') {
                match k.trim() {
                    "pgs_id" => id = v.trim().to_string(),
                    "trait_reported" => trait_name = v.trim().to_string(),
                    "trait_mapped" if trait_name.is_empty() => trait_name = v.trim().to_string(),
                    "variants_number" => declared_count = v.trim().parse().ok(),
                    "HmPOS_build" | "genome_build" => {
                        let b = v.trim();
                        if b != "NR" && !b.is_empty() {
                            build = Some(b.to_string());
                        }
                    }
                    _ => {}
                }
            }
            continue;
        }

        // First non-comment line is the column header.
        if cols.is_none() {
            let map: HashMap<String, usize> = line
                .split('\t')
                .enumerate()
                .map(|(i, name)| (name.trim().to_string(), i))
                .collect();
            cols = Some(map);
            continue;
        }

        let cols = cols.as_ref().unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |name: &str| -> Option<&str> {
            cols.get(name).and_then(|&i| fields.get(i)).copied().map(str::trim)
        };

        let Some(weight) = nonempty(get("effect_weight")).and_then(|s| s.parse::<f64>().ok()) else {
            continue;
        };
        let Some(effect_allele) = nonempty(get("effect_allele")).map(|s| s.to_ascii_uppercase()) else {
            continue;
        };
        let other_allele = nonempty(get("other_allele")).map(|s| s.to_ascii_uppercase());
        let rsid = nonempty(get("hm_rsID"))
            .or_else(|| nonempty(get("rsID")))
            .map(str::to_string);
        let chrom = nonempty(get("hm_chr")).or_else(|| nonempty(get("chr_name"))).map(str::to_string);
        let pos = nonempty(get("hm_pos"))
            .or_else(|| nonempty(get("chr_position")))
            .and_then(|s| s.parse::<u64>().ok())
            .map(|p| p.saturating_sub(1));
        let allele_freq = nonempty(get("allelefrequency_effect")).and_then(|s| s.parse::<f64>().ok());

        variants.push(ScoreVariant {
            rsid,
            chrom,
            pos,
            effect_allele,
            other_allele,
            weight,
            allele_freq,
        });
    }

    if cols.is_none() {
        return Err("no column header found — not a PGS scoring file".to_string());
    }
    if variants.is_empty() {
        return Err("scoring file contained no usable variants".to_string());
    }
    if id.is_empty() {
        id = "PGS".to_string();
    }
    if trait_name.is_empty() {
        trait_name = "polygenic score".to_string();
    }

    Ok(ScoreFile {
        id,
        trait_name,
        build,
        declared_count,
        variants,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
###PGS CATALOG SCORING FILE
#pgs_id=PGS999999
#trait_reported=Test Trait
#variants_number=3
#HmPOS_build=GRCh37
rsID\tchr_name\teffect_allele\tother_allele\teffect_weight\tallelefrequency_effect\thm_rsID\thm_chr\thm_pos
rs100\t1\tA\tG\t0.5\t0.30\trs100\t1\t1001
rs200\t2\tT\tC\t-0.2\t0.10\trs200\t2\t2001
rs300\t3\tG\tA\t0.1\t\trs300\t3\t3001
";

    #[test]
    fn parses_metadata_and_variants() {
        let sf = parse(SAMPLE).unwrap();
        assert_eq!(sf.id, "PGS999999");
        assert_eq!(sf.trait_name, "Test Trait");
        assert_eq!(sf.build.as_deref(), Some("GRCh37"));
        assert_eq!(sf.declared_count, Some(3));
        assert_eq!(sf.variant_count(), 3);

        let v = &sf.variants[0];
        assert_eq!(v.rsid.as_deref(), Some("rs100"));
        assert_eq!(v.effect_allele, "A");
        assert_eq!(v.other_allele.as_deref(), Some("G"));
        assert!((v.weight - 0.5).abs() < 1e-9);
        assert_eq!(v.allele_freq, Some(0.30));
        assert_eq!(v.pos, Some(1000)); // 1-based 1001 -> 0-based 1000
        // third variant has no AF
        assert_eq!(sf.variants[2].allele_freq, None);
    }
}
