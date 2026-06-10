//! Downloader for PGS Catalog scoring files (EBI FTP over HTTPS).

use std::io::Read;
use std::time::Duration;

use flate2::read::MultiGzDecoder;

use crate::error::{AnnotateError, Result};

pub struct PgsClient {
    agent: ureq::Agent,
}

impl PgsClient {
    pub fn new() -> Self {
        let agent = ureq::Agent::config_builder()
            .timeout_global(Some(Duration::from_secs(90)))
            .build()
            .into();
        Self { agent }
    }

    /// Download and decompress a harmonized scoring file by PGS id, returning
    /// the raw text. Tries the GRCh37 then GRCh38 harmonized files, then the
    /// un-harmonized file.
    pub fn fetch_scoring_text(&self, id: &str) -> Result<String> {
        let id = id.trim().to_ascii_uppercase();
        if id.is_empty() {
            return Err(AnnotateError::Response("empty PGS id".into()));
        }
        let base = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores";
        let candidates = [
            format!("{base}/{id}/ScoringFiles/Harmonized/{id}_hmPOS_GRCh37.txt.gz"),
            format!("{base}/{id}/ScoringFiles/Harmonized/{id}_hmPOS_GRCh38.txt.gz"),
            format!("{base}/{id}/ScoringFiles/{id}.txt.gz"),
        ];
        let mut last = AnnotateError::Response(format!("no scoring file found for {id}"));
        for url in candidates {
            match self.download_gz(&url) {
                Ok(text) => return Ok(text),
                Err(e) => last = e,
            }
        }
        Err(last)
    }

    fn download_gz(&self, url: &str) -> Result<String> {
        let mut resp = self.agent.get(url).call()?;
        // PGS scoring files routinely exceed ureq's default 10 MB body cap.
        let bytes = resp
            .body_mut()
            .with_config()
            .limit(512 * 1024 * 1024)
            .read_to_vec()
            .map_err(|e| AnnotateError::Http(format!("download read: {e}")))?;
        let mut text = String::new();
        MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .map_err(|e| AnnotateError::Response(format!("gunzip: {e}")))?;
        Ok(text)
    }
}

impl Default for PgsClient {
    fn default() -> Self {
        Self::new()
    }
}
