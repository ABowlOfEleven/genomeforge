//! NCBI PubMed literature lookup via E-utilities (esearch + esummary).
//!
//! Powers the variant detail panel's "related publications" drill-down. The
//! search term is the variant's rsID — high precision, since rsIDs are indexed
//! verbatim in many abstracts — and we return a handful of the most relevant
//! articles. Runs on the worker thread; never on the paint loop.

use std::time::Duration;

use serde_json::Value;

use crate::error::{AnnotateError, Result};

/// A single PubMed article summary (enough to show a clickable reference).
#[derive(Debug, Clone)]
pub struct Article {
    pub pmid: String,
    pub title: String,
    pub journal: String,
    pub year: String,
}

impl Article {
    /// Canonical PubMed URL for this article.
    pub fn url(&self) -> String {
        format!("https://pubmed.ncbi.nlm.nih.gov/{}/", self.pmid)
    }
}

/// Client for the NCBI E-utilities PubMed endpoints.
pub struct PubMedClient {
    agent: ureq::Agent,
    base: String,
}

impl PubMedClient {
    pub fn new() -> Self {
        Self {
            agent: ureq::Agent::config_builder()
                .timeout_global(Some(Duration::from_secs(20)))
                .build()
                .into(),
            base: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils".to_string(),
        }
    }

    /// Search PubMed for `term`, returning up to `retmax` article summaries,
    /// most relevant first. An empty result is normal (nothing indexed).
    pub fn search(&self, term: &str, retmax: usize) -> Result<Vec<Article>> {
        let ids = self.esearch(term, retmax)?;
        if ids.is_empty() {
            return Ok(Vec::new());
        }
        self.esummary(&ids)
    }

    /// esearch → ranked list of PMIDs.
    fn esearch(&self, term: &str, retmax: usize) -> Result<Vec<String>> {
        let url = format!(
            "{}/esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax={}&tool=GenomeForge&term={}",
            self.base,
            retmax,
            encode(term),
        );
        let v: Value = self.agent.get(&url).call()?.body_mut().read_json()?;
        Ok(v.get("esearchresult")
            .and_then(|e| e.get("idlist"))
            .and_then(Value::as_array)
            .map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect())
            .unwrap_or_default())
    }

    /// esummary → article titles / journals / years, preserving PMID order.
    fn esummary(&self, ids: &[String]) -> Result<Vec<Article>> {
        let url = format!(
            "{}/esummary.fcgi?db=pubmed&retmode=json&tool=GenomeForge&id={}",
            self.base,
            ids.join(","),
        );
        let v: Value = self.agent.get(&url).call()?.body_mut().read_json()?;
        let result = v
            .get("result")
            .ok_or_else(|| AnnotateError::Response("PubMed esummary: no result object".into()))?;

        let mut out = Vec::new();
        for id in ids {
            let Some(item) = result.get(id) else { continue };
            let title = item
                .get("title")
                .and_then(Value::as_str)
                .unwrap_or("")
                .trim()
                .trim_end_matches('.')
                .to_string();
            if title.is_empty() {
                continue;
            }
            let journal = item.get("source").and_then(Value::as_str).unwrap_or("").to_string();
            let year = item
                .get("pubdate")
                .and_then(Value::as_str)
                .unwrap_or("")
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            out.push(Article { pmid: id.clone(), title, journal, year });
        }
        Ok(out)
    }
}

impl Default for PubMedClient {
    fn default() -> Self {
        Self::new()
    }
}

/// Minimal percent-encoding for a query term (alnum and `-_.~` pass through,
/// spaces become `+`, everything else is percent-escaped).
fn encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => out.push(b as char),
            b' ' => out.push('+'),
            _ => out.push_str(&format!("%{b:02X}")),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_query_terms() {
        assert_eq!(encode("rs6025"), "rs6025");
        assert_eq!(encode("Factor V Leiden"), "Factor+V+Leiden");
        assert_eq!(encode("a/b"), "a%2Fb");
    }

    #[test]
    fn article_url_is_canonical() {
        let a = Article {
            pmid: "12345".into(),
            title: "x".into(),
            journal: "Nature".into(),
            year: "2020".into(),
        };
        assert_eq!(a.url(), "https://pubmed.ncbi.nlm.nih.gov/12345/");
    }
}
