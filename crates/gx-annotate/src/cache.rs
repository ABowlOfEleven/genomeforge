//! Local `redb` cache for the hybrid annotation model.
//!
//! Pure-Rust embedded storage (no C, no build script). Every network result is
//! stored so repeat lookups are instant and work offline. Keys are only rsIDs
//! and region strings — raw genotypes never leave the user's machine. Values
//! are JSON blobs of the corresponding `gx-core` / annotation types.

use std::collections::HashMap;
use std::path::Path;

use redb::{Database, ReadableDatabase, ReadableTableMetadata, TableDefinition};

use gx_core::{Assembly, Feature};

use crate::annotation::{GeneLocation, VariantAnnotation};
use crate::error::{Result, db_err};

const VARIANT_ANN: TableDefinition<&str, &str> = TableDefinition::new("variant_annotation");
const REGION_FEATURES: TableDefinition<&str, &str> = TableDefinition::new("region_features");
const REGION_SEQUENCE: TableDefinition<&str, &str> = TableDefinition::new("region_sequence");
const GENE_LOOKUP: TableDefinition<&str, &str> = TableDefinition::new("gene_lookup");

const ALL_TABLES: [TableDefinition<&str, &str>; 4] =
    [VARIANT_ANN, REGION_FEATURES, REGION_SEQUENCE, GENE_LOOKUP];

/// Per-table entry counts, for the Settings cache panel.
#[derive(Debug, Clone, Copy, Default)]
pub struct CacheStats {
    pub annotations: u64,
    pub feature_regions: u64,
    pub sequence_regions: u64,
    pub genes: u64,
}

impl CacheStats {
    /// Total cached entries across all tables.
    pub fn total(&self) -> u64 {
        self.annotations + self.feature_regions + self.sequence_regions + self.genes
    }
}

pub struct Cache {
    db: Database,
}

impl Cache {
    /// Open (creating if needed) a cache database at `path`.
    pub fn open(path: &Path) -> Result<Self> {
        let db = Database::create(path).map_err(db_err)?;
        let cache = Self { db };
        cache.ensure_tables()?;
        Ok(cache)
    }

    /// An ephemeral in-memory cache (tests, "no persistence" mode).
    pub fn open_in_memory() -> Result<Self> {
        let db = Database::builder()
            .create_with_backend(redb::backends::InMemoryBackend::new())
            .map_err(db_err)?;
        let cache = Self { db };
        cache.ensure_tables()?;
        Ok(cache)
    }

    /// Materialise all tables up front so reads on a fresh DB don't error.
    fn ensure_tables(&self) -> Result<()> {
        let w = self.db.begin_write().map_err(db_err)?;
        for table in ALL_TABLES {
            w.open_table(table).map_err(db_err)?;
        }
        w.commit().map_err(db_err)?;
        Ok(())
    }

    fn put(&self, table: TableDefinition<&str, &str>, key: &str, value: &str) -> Result<()> {
        let w = self.db.begin_write().map_err(db_err)?;
        {
            let mut t = w.open_table(table).map_err(db_err)?;
            t.insert(key, value).map_err(db_err)?;
        }
        w.commit().map_err(db_err)?;
        Ok(())
    }

    fn get(&self, table: TableDefinition<&str, &str>, key: &str) -> Result<Option<String>> {
        let r = self.db.begin_read().map_err(db_err)?;
        let t = r.open_table(table).map_err(db_err)?;
        let value = t.get(key).map_err(db_err)?;
        Ok(value.map(|g| g.value().to_string()))
    }

    fn region_key(assembly: Assembly, region: &str) -> String {
        format!("{}:{}", assembly.ucsc_name(), region)
    }

    // ---- variant annotations ----------------------------------------------

    pub fn get_annotation(&self, rsid: &str) -> Result<Option<VariantAnnotation>> {
        Ok(match self.get(VARIANT_ANN, rsid)? {
            Some(json) => Some(serde_json::from_str(&json)?),
            None => None,
        })
    }

    pub fn put_annotation(&self, ann: &VariantAnnotation) -> Result<()> {
        let json = serde_json::to_string(ann)?;
        self.put(VARIANT_ANN, &ann.rsid, &json)
    }

    /// Partition `rsids` into the ones already cached and the ones still missing.
    pub fn split_cached(
        &self,
        rsids: &[String],
    ) -> Result<(HashMap<String, VariantAnnotation>, Vec<String>)> {
        let mut cached = HashMap::new();
        let mut missing = Vec::new();
        for id in rsids {
            match self.get_annotation(id)? {
                Some(ann) => {
                    cached.insert(id.clone(), ann);
                }
                None => missing.push(id.clone()),
            }
        }
        Ok((cached, missing))
    }

    // ---- region features ---------------------------------------------------

    pub fn get_features(&self, assembly: Assembly, region: &str) -> Result<Option<Vec<Feature>>> {
        let key = Self::region_key(assembly, region);
        Ok(match self.get(REGION_FEATURES, &key)? {
            Some(json) => Some(serde_json::from_str(&json)?),
            None => None,
        })
    }

    pub fn put_features(&self, assembly: Assembly, region: &str, features: &[Feature]) -> Result<()> {
        let key = Self::region_key(assembly, region);
        let json = serde_json::to_string(features)?;
        self.put(REGION_FEATURES, &key, &json)
    }

    // ---- region sequence ---------------------------------------------------

    pub fn get_sequence(&self, assembly: Assembly, region: &str) -> Result<Option<String>> {
        let key = Self::region_key(assembly, region);
        self.get(REGION_SEQUENCE, &key)
    }

    pub fn put_sequence(&self, assembly: Assembly, region: &str, seq: &str) -> Result<()> {
        let key = Self::region_key(assembly, region);
        self.put(REGION_SEQUENCE, &key, seq)
    }

    // ---- gene lookup -------------------------------------------------------

    pub fn get_gene(&self, assembly: Assembly, symbol: &str) -> Result<Option<GeneLocation>> {
        let key = Self::region_key(assembly, &symbol.to_ascii_uppercase());
        Ok(match self.get(GENE_LOOKUP, &key)? {
            Some(json) => Some(serde_json::from_str(&json)?),
            None => None,
        })
    }

    pub fn put_gene(&self, assembly: Assembly, symbol: &str, gene: &GeneLocation) -> Result<()> {
        let key = Self::region_key(assembly, &symbol.to_ascii_uppercase());
        let json = serde_json::to_string(gene)?;
        self.put(GENE_LOOKUP, &key, &json)
    }

    // ---- maintenance -------------------------------------------------------

    /// Entry counts per table.
    pub fn stats(&self) -> Result<CacheStats> {
        let r = self.db.begin_read().map_err(db_err)?;
        let count = |table: TableDefinition<&str, &str>| -> Result<u64> {
            r.open_table(table).map_err(db_err)?.len().map_err(db_err)
        };
        Ok(CacheStats {
            annotations: count(VARIANT_ANN)?,
            feature_regions: count(REGION_FEATURES)?,
            sequence_regions: count(REGION_SEQUENCE)?,
            genes: count(GENE_LOOKUP)?,
        })
    }

    /// Drop every cached entry, leaving the (empty) tables in place.
    pub fn clear(&self) -> Result<()> {
        let w = self.db.begin_write().map_err(db_err)?;
        for table in ALL_TABLES {
            w.delete_table(table).map_err(db_err)?;
        }
        w.commit().map_err(db_err)?;
        self.ensure_tables()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotation_round_trip_and_split() {
        let cache = Cache::open_in_memory().unwrap();
        let mut ann = VariantAnnotation::not_found("rs1");
        ann.found = true;
        ann.gene = Some("BRCA1".into());
        cache.put_annotation(&ann).unwrap();

        let got = cache.get_annotation("rs1").unwrap().unwrap();
        assert_eq!(got.gene.as_deref(), Some("BRCA1"));

        let (cached, missing) = cache.split_cached(&["rs1".into(), "rs2".into()]).unwrap();
        assert!(cached.contains_key("rs1"));
        assert_eq!(missing, ["rs2"]);
    }

    #[test]
    fn sequence_round_trip_is_build_namespaced() {
        let cache = Cache::open_in_memory().unwrap();
        cache.put_sequence(Assembly::Grch38, "1:100-200", "ACGT").unwrap();
        assert_eq!(
            cache.get_sequence(Assembly::Grch38, "1:100-200").unwrap().as_deref(),
            Some("ACGT")
        );
        assert!(cache.get_sequence(Assembly::Grch37, "1:100-200").unwrap().is_none());
    }

    #[test]
    fn missing_keys_return_none_not_error() {
        let cache = Cache::open_in_memory().unwrap();
        assert!(cache.get_annotation("rsNONE").unwrap().is_none());
        assert!(cache.get_gene(Assembly::Grch38, "NOPE").unwrap().is_none());
    }
}
