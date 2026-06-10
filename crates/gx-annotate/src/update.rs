//! Check GitHub for the latest released version of the app.
//!
//! Used by the in-app "check for updates" feature. Network-only (no cache); the
//! caller decides when to run it and whether the result is newer than the
//! running build.

use std::time::Duration;

use serde_json::Value;

use crate::error::{AnnotateError, Result};

/// The latest GitHub Release for a repository.
#[derive(Debug, Clone)]
pub struct UpdateInfo {
    /// The release tag, e.g. `v0.1.1`.
    pub tag: String,
    /// The release title (falls back to the tag).
    pub name: String,
    /// The release web page.
    pub url: String,
}

/// Fetch the latest published GitHub Release for `owner/repo`
/// (e.g. `"ABowlOfEleven/genomeforge"`). Returns `Ok(None)` when the repository
/// has no releases yet (404). Other transport/HTTP errors propagate.
pub fn latest_release(repo: &str) -> Result<Option<UpdateInfo>> {
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .timeout_global(Some(Duration::from_secs(15)))
        .build()
        .into();
    let url = format!("https://api.github.com/repos/{repo}/releases/latest");
    match agent
        .get(&url)
        .header("User-Agent", "GenomeForge")
        .header("Accept", "application/vnd.github+json")
        .call()
    {
        Ok(mut resp) => {
            let v: Value = resp.body_mut().with_config().limit(4 * 1024 * 1024).read_json()?;
            Ok(parse_release(&v))
        }
        // No releases published yet.
        Err(ureq::Error::StatusCode(404)) => Ok(None),
        Err(e) => Err(AnnotateError::Http(e.to_string())),
    }
}

/// Extract an [`UpdateInfo`] from a GitHub release JSON object (network-free, so
/// it is unit-testable).
fn parse_release(v: &Value) -> Option<UpdateInfo> {
    let tag = v.get("tag_name").and_then(Value::as_str)?.to_string();
    if tag.is_empty() {
        return None;
    }
    let name = v.get("name").and_then(Value::as_str).filter(|s| !s.is_empty()).unwrap_or(&tag).to_string();
    let url = v.get("html_url").and_then(Value::as_str).unwrap_or_default().to_string();
    Some(UpdateInfo { tag, name, url })
}

/// Compare two dotted version strings (leading `v` ignored). Returns true when
/// `latest` is strictly newer than `current`. Non-numeric components compare as 0,
/// so this is a pragmatic semver-lite check, not a full pre-release-aware compare.
pub fn is_newer(latest: &str, current: &str) -> bool {
    fn parts(s: &str) -> Vec<u64> {
        s.trim_start_matches(['v', 'V'])
            .split(['.', '-', '+'])
            .map(|p| p.parse::<u64>().unwrap_or(0))
            .collect()
    }
    let (a, b) = (parts(latest), parts(current));
    for i in 0..a.len().max(b.len()) {
        let (x, y) = (a.get(i).copied().unwrap_or(0), b.get(i).copied().unwrap_or(0));
        if x != y {
            return x > y;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_release_json() {
        let v: Value = serde_json::from_str(
            r#"{"tag_name":"v0.1.1","name":"GenomeForge v0.1.1","html_url":"https://example/releases/v0.1.1"}"#,
        )
        .unwrap();
        let u = parse_release(&v).unwrap();
        assert_eq!(u.tag, "v0.1.1");
        assert_eq!(u.name, "GenomeForge v0.1.1");
        assert!(u.url.ends_with("v0.1.1"));
    }

    #[test]
    fn name_falls_back_to_tag() {
        let v: Value = serde_json::from_str(r#"{"tag_name":"v0.2.0","name":""}"#).unwrap();
        assert_eq!(parse_release(&v).unwrap().name, "v0.2.0");
    }

    #[test]
    fn version_ordering() {
        assert!(is_newer("v0.1.1", "0.1.0"));
        assert!(is_newer("0.2.0", "v0.1.9"));
        assert!(is_newer("1.0.0", "0.9.9"));
        assert!(!is_newer("0.1.0", "0.1.0"));
        assert!(!is_newer("0.1.0", "0.1.1"));
    }
}
