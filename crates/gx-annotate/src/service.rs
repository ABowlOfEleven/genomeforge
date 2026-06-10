//! [`AnnotationService`] — the hybrid orchestrator the worker thread owns.
//!
//! Every lookup is cache-first. When online, a cache miss is fetched and then
//! written back; when offline, a miss yields empty/`None` rather than an error
//! so the UI degrades gracefully to whatever has been cached.

use std::collections::HashMap;

use gx_core::{Assembly, Feature, GenomicRange};

use crate::annotation::{GeneLocation, VariantAnnotation};
use crate::cache::Cache;
use crate::clients::{EnsemblClient, MyVariantClient};
use crate::error::{AnnotateError, Result};
use crate::literature::{Article, PubMedClient};
use crate::pgs::PgsClient;

/// MyVariant accepts up to 1000 ids/batch; stay just under.
const BATCH_SIZE: usize = 900;

pub struct AnnotationService {
    cache: Cache,
    myvariant: MyVariantClient,
    ensembl: EnsemblClient,
    pgs: PgsClient,
    pubmed: PubMedClient,
    online: bool,
    /// Build of the currently loaded dataset (routes Ensembl host).
    assembly: Assembly,
}

impl AnnotationService {
    pub fn new(cache: Cache, online: bool, assembly: Assembly, api_key: Option<String>) -> Self {
        Self {
            cache,
            myvariant: MyVariantClient::new(api_key),
            ensembl: EnsemblClient::new(),
            pgs: PgsClient::new(),
            pubmed: PubMedClient::new(),
            online,
            assembly,
        }
    }

    /// Up to `limit` PubMed articles relevant to a variant (searched by rsID —
    /// the most precise term). Errors when offline; an empty list means nothing
    /// is indexed. Not cached (user-initiated, and results evolve).
    pub fn literature(&self, rsid: &str, limit: usize) -> Result<Vec<Article>> {
        if !self.online {
            return Err(AnnotateError::OfflineMiss(format!("PubMed {rsid}")));
        }
        self.pubmed.search(rsid, limit)
    }

    /// Download a PGS Catalog scoring file's raw text (no caching — files vary
    /// in size). Errors when offline.
    pub fn fetch_pgs_text(&self, id: &str) -> Result<String> {
        if !self.online {
            return Err(AnnotateError::OfflineMiss(format!("PGS {id}")));
        }
        self.pgs.fetch_scoring_text(id)
    }

    pub fn set_online(&mut self, online: bool) {
        self.online = online;
    }

    pub fn is_online(&self) -> bool {
        self.online
    }

    pub fn set_assembly(&mut self, assembly: Assembly) {
        self.assembly = assembly;
    }

    /// Annotate a set of rsIDs, cache-first. Returns one entry per rsID that
    /// resolved (cached or freshly fetched).
    pub fn annotate(&self, rsids: &[String]) -> Result<HashMap<String, VariantAnnotation>> {
        let (mut resolved, missing) = self.cache.split_cached(rsids)?;
        if self.online {
            for chunk in missing.chunks(BATCH_SIZE) {
                // A single failed chunk must not discard the cached hits and the
                // chunks that already succeeded — degrade gracefully to partial.
                match self.myvariant.annotate_batch(chunk) {
                    Ok(fetched) => {
                        for ann in fetched {
                            if !ann.rsid.is_empty() {
                                self.cache.put_annotation(&ann)?;
                                resolved.insert(ann.rsid.clone(), ann);
                            }
                        }
                    }
                    Err(e) => log::warn!("annotation batch failed ({} ids): {e}", chunk.len()),
                }
            }
        }
        Ok(resolved)
    }

    /// Gene / transcript / exon features overlapping `range`, cache-first.
    pub fn features(&self, range: &GenomicRange) -> Result<Vec<Feature>> {
        let region = range.ensembl_region();
        if let Some(cached) = self.cache.get_features(self.assembly, &region)? {
            return Ok(cached);
        }
        if !self.online {
            return Ok(Vec::new());
        }
        let features = self.ensembl.features(self.assembly, range)?;
        self.cache.put_features(self.assembly, &region, &features)?;
        Ok(features)
    }

    /// Reference sequence for `range`, cache-first.
    pub fn sequence(&self, range: &GenomicRange) -> Result<Option<String>> {
        let region = range.ensembl_region();
        if let Some(seq) = self.cache.get_sequence(self.assembly, &region)? {
            return Ok(Some(seq));
        }
        if !self.online {
            return Ok(None);
        }
        let seq = self.ensembl.sequence(self.assembly, range)?;
        self.cache.put_sequence(self.assembly, &region, &seq)?;
        Ok(Some(seq))
    }

    /// Resolve a gene symbol to a location, cache-first.
    pub fn gene_location(&self, symbol: &str) -> Result<Option<GeneLocation>> {
        if let Some(gene) = self.cache.get_gene(self.assembly, symbol)? {
            return Ok(Some(gene));
        }
        if !self.online {
            return Ok(None);
        }
        let gene = self.ensembl.gene_location(self.assembly, symbol)?;
        if let Some(g) = &gene {
            self.cache.put_gene(self.assembly, symbol, g)?;
        }
        Ok(gene)
    }
}
