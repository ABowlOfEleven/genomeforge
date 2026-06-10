//! User settings, persisted via eframe's storage.

use serde::{Deserialize, Serialize};

use gx_core::Assembly;

use crate::ux::Tier;

#[derive(Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct Settings {
    /// When false, the app never touches the network and uses only the cache.
    pub online: bool,
    /// Optional MyVariant.info API key (lifts the anonymous rate limit).
    pub api_key: String,
    /// Force a particular build instead of trusting the imported file's.
    pub assembly_override: Option<Assembly>,
    /// UI experience tier (beginner / intermediate / expert).
    pub tier: Tier,
    /// Whether each workspace's first-run tutorial has been shown.
    pub seen_genome: bool,
    pub seen_plasmid: bool,
    pub seen_crispr: bool,
    pub seen_phenotype: bool,
    /// Recently opened files (most recent first), for File ▸ Open Recent.
    #[serde(default)]
    pub recent_files: Vec<String>,
}

impl Settings {
    /// Record `path` as the most-recently-opened file (deduped, capped at 8).
    pub fn push_recent(&mut self, path: &std::path::Path) {
        let p = path.to_string_lossy().to_string();
        self.recent_files.retain(|x| x != &p);
        self.recent_files.insert(0, p);
        self.recent_files.truncate(8);
    }
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            online: true,
            api_key: String::new(),
            assembly_override: None,
            tier: Tier::Beginner,
            seen_genome: false,
            seen_plasmid: false,
            seen_crispr: false,
            seen_phenotype: false,
            recent_files: Vec::new(),
        }
    }
}

impl Settings {
    /// The MyVariant key as an `Option`, empty string treated as unset.
    pub fn api_key_opt(&self) -> Option<String> {
        let key = self.api_key.trim();
        (!key.is_empty()).then(|| key.to_string())
    }
}
