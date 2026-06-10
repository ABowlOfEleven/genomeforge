//! Blocking HTTP clients for the two annotation sources. These run on the
//! worker thread, never the UI thread.
//!
//! * [`MyVariantClient`] — batch rsID → variant annotation (dbSNP, ClinVar,
//!   gnomAD, CADD, SnpEff, …) via <https://myvariant.info>.
//! * [`EnsemblClient`] — reference sequence, gene models, and gene-symbol
//!   lookup via the build-matched Ensembl REST host.

use std::time::Duration;

use serde_json::Value;

use gx_core::{Assembly, Feature, GenomicRange};

use crate::annotation::{GeneLocation, VariantAnnotation, feature_from_ensembl};
use crate::error::{AnnotateError, Result};

/// Fields requested from MyVariant — the union backing the detail panel.
const MYVARIANT_FIELDS: &str =
    "dbsnp,clinvar,gnomad_genome,gnomad_exome,cadd,snpeff,dbnsfp,gwassnps,vcf";

fn build_agent(timeout_secs: u64) -> ureq::Agent {
    ureq::Agent::config_builder()
        .timeout_global(Some(Duration::from_secs(timeout_secs)))
        .build()
        .into()
}

/// Client for the MyVariant.info batch query endpoint.
pub struct MyVariantClient {
    agent: ureq::Agent,
    base: String,
    api_key: Option<String>,
}

impl MyVariantClient {
    pub fn new(api_key: Option<String>) -> Self {
        Self {
            agent: build_agent(30),
            base: "https://myvariant.info/v1".to_string(),
            api_key,
        }
    }

    /// Annotate up to ~1000 rsIDs in a single POST. Returns one
    /// [`VariantAnnotation`] per response element (found or not-found).
    pub fn annotate_batch(&self, rsids: &[String]) -> Result<Vec<VariantAnnotation>> {
        if rsids.is_empty() {
            return Ok(Vec::new());
        }
        let joined = rsids.join(",");
        let url = format!("{}/query", self.base);

        let mut form: Vec<(&str, &str)> = vec![
            ("q", joined.as_str()),
            ("scopes", "dbsnp.rsid"),
            ("fields", MYVARIANT_FIELDS),
        ];
        if let Some(key) = &self.api_key {
            form.push(("api_key", key.as_str()));
        }

        let body: Value = self
            .agent
            .post(&url)
            .send_form(form)? // ureq wants an IntoIterator of owned (K, V) tuples
            .body_mut()
            .with_config()
            .limit(128 * 1024 * 1024) // large batches blow past ureq's 10 MB default
            .read_json()?;

        let arr = body
            .as_array()
            .ok_or_else(|| AnnotateError::Response("MyVariant: expected a JSON array".into()))?;

        Ok(arr
            .iter()
            .map(|hit| {
                let query = hit
                    .get("query")
                    .and_then(Value::as_str)
                    .unwrap_or_default()
                    .to_string();
                VariantAnnotation::from_myvariant(&query, hit)
            })
            .collect())
    }
}

/// Client for Ensembl REST, build-aware (routes to GRCh38 or GRCh37 host).
pub struct EnsemblClient {
    agent: ureq::Agent,
}

impl EnsemblClient {
    pub fn new() -> Self {
        Self {
            agent: build_agent(30),
        }
    }

    fn get_json(&self, url: &str) -> Result<Value> {
        let body = self
            .agent
            .get(url)
            .header("Accept", "application/json")
            .call()?
            .body_mut()
            .with_config()
            .limit(128 * 1024 * 1024)
            .read_json::<Value>()?;
        Ok(body)
    }

    /// Resolve a gene symbol (e.g. `BRCA1`) to its genomic location. A genuine
    /// 400/404 means "no such symbol" (`None`); other errors (rate limit, 5xx,
    /// transport) propagate so the UI can distinguish "absent" from "offline".
    pub fn gene_location(&self, assembly: Assembly, symbol: &str) -> Result<Option<GeneLocation>> {
        let url = format!(
            "{}/lookup/symbol/homo_sapiens/{}?expand=0;content-type=application/json",
            assembly.ensembl_rest_host(),
            symbol,
        );
        match self
            .agent
            .get(&url)
            .header("Accept", "application/json")
            .call()
        {
            Ok(mut resp) => {
                let v = resp
                    .body_mut()
                    .with_config()
                    .limit(8 * 1024 * 1024)
                    .read_json::<Value>()?;
                Ok(GeneLocation::from_ensembl(&v))
            }
            Err(ureq::Error::StatusCode(400 | 404)) => Ok(None),
            Err(e) => Err(AnnotateError::Http(e.to_string())),
        }
    }

    /// Gene / transcript / exon features overlapping a region.
    pub fn features(&self, assembly: Assembly, range: &GenomicRange) -> Result<Vec<Feature>> {
        let url = format!(
            "{}/overlap/region/human/{}?feature=gene;feature=transcript;feature=exon;content-type=application/json",
            assembly.ensembl_rest_host(),
            range.ensembl_region(),
        );
        let body = self.get_json(&url)?;
        let arr = body
            .as_array()
            .ok_or_else(|| AnnotateError::Response("Ensembl overlap: expected array".into()))?;
        Ok(arr.iter().filter_map(feature_from_ensembl).collect())
    }

    /// Reference sequence for a region (uppercase ACGT…).
    pub fn sequence(&self, assembly: Assembly, range: &GenomicRange) -> Result<String> {
        let url = format!(
            "{}/sequence/region/human/{}?content-type=application/json",
            assembly.ensembl_rest_host(),
            range.ensembl_region(),
        );
        let body = self.get_json(&url)?;
        body.get("seq")
            .and_then(Value::as_str)
            .map(str::to_string)
            .ok_or_else(|| AnnotateError::Response("Ensembl sequence: no `seq` field".into()))
    }
}

impl Default for EnsemblClient {
    fn default() -> Self {
        Self::new()
    }
}
