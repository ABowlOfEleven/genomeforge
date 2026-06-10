use serde::{Deserialize, Serialize};

/// DNA strand. `Unknown` is used when a source does not specify one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum Strand {
    #[default]
    Unknown,
    Forward,
    Reverse,
}

impl Strand {
    /// Interpret a signed strand indicator (`+1` / `-1`, as Ensembl returns).
    pub fn from_sign(sign: i64) -> Self {
        match sign {
            s if s > 0 => Strand::Forward,
            s if s < 0 => Strand::Reverse,
            _ => Strand::Unknown,
        }
    }

    pub fn symbol(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
            Strand::Unknown => '.',
        }
    }
}

/// Normalise a chromosome / contig name to the form Ensembl uses.
///
/// Strips a leading `chr` (any case) and canonicalises the mitochondrial
/// contig to `MT`. Autosome numbers pass through unchanged; `X`/`Y` are
/// upper-cased.
///
/// ```
/// # use gx_core::normalize_contig;
/// assert_eq!(normalize_contig("chr1"), "1");
/// assert_eq!(normalize_contig("chrX"), "X");
/// assert_eq!(normalize_contig("chrM"), "MT");
/// assert_eq!(normalize_contig("MT"), "MT");
/// ```
pub fn normalize_contig(name: &str) -> String {
    let trimmed = name.trim();
    let stripped = if trimmed.len() > 3 && trimmed[..3].eq_ignore_ascii_case("chr") {
        &trimmed[3..]
    } else {
        trimmed
    };
    match stripped.to_ascii_uppercase().as_str() {
        "M" | "MT" => "MT".to_string(),
        upper => upper.to_string(),
    }
}

/// A half-open genomic interval `[start, end)` on a single contig, 0-based.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GenomicRange {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

impl GenomicRange {
    pub fn new(contig: impl Into<String>, start: u64, end: u64, strand: Strand) -> Self {
        Self {
            contig: contig.into(),
            start,
            end,
            strand,
        }
    }

    /// Length in bases.
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }

    /// Centre position (0-based), useful for label placement.
    pub fn midpoint(&self) -> u64 {
        self.start + self.len() / 2
    }

    /// Does this range contain a 0-based position on the same contig?
    pub fn contains(&self, contig: &str, pos: u64) -> bool {
        self.contig == contig && pos >= self.start && pos < self.end
    }

    /// Do two ranges overlap (same contig, non-empty intersection)?
    pub fn overlaps(&self, other: &GenomicRange) -> bool {
        self.contig == other.contig && self.start < other.end && other.start < self.end
    }

    /// 1-based inclusive, thousands-separated display, e.g. `1:1,001-2,000`.
    pub fn display_region(&self) -> String {
        format!(
            "{}:{}-{}",
            self.contig,
            group_thousands(self.start + 1),
            group_thousands(self.end)
        )
    }

    /// Ensembl REST `region` form `contig:start-end` (1-based inclusive).
    pub fn ensembl_region(&self) -> String {
        format!("{}:{}-{}", self.contig, self.start + 1, self.end)
    }
}

/// A single 0-based genomic position.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Locus {
    pub contig: String,
    pub pos: u64,
}

impl Locus {
    pub fn new(contig: impl Into<String>, pos: u64) -> Self {
        Self {
            contig: contig.into(),
            pos,
        }
    }

    /// 1-based display, e.g. `1:1,001`.
    pub fn display(&self) -> String {
        format!("{}:{}", self.contig, group_thousands(self.pos + 1))
    }
}

/// Format an integer with `,` thousands separators.
pub(crate) fn group_thousands(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    let first = bytes.len() % 3;
    for (i, &b) in bytes.iter().enumerate() {
        // `i >= first` must be checked first: it both guards the `i - first`
        // subtraction against underflow and is the correct grouping condition.
        if i >= first && i != 0 && (i - first).is_multiple_of(3) {
            out.push(',');
        }
        out.push(b as char);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contig_normalisation() {
        assert_eq!(normalize_contig("chr1"), "1");
        assert_eq!(normalize_contig("CHR12"), "12");
        assert_eq!(normalize_contig("chrX"), "X");
        assert_eq!(normalize_contig("chrM"), "MT");
        assert_eq!(normalize_contig("MT"), "MT");
        assert_eq!(normalize_contig(" 7 "), "7");
    }

    #[test]
    fn overlap_and_contains() {
        let a = GenomicRange::new("1", 100, 200, Strand::Forward);
        let b = GenomicRange::new("1", 150, 250, Strand::Reverse);
        let c = GenomicRange::new("1", 200, 300, Strand::Forward);
        let d = GenomicRange::new("2", 150, 250, Strand::Forward);
        assert!(a.overlaps(&b));
        assert!(!a.overlaps(&c)); // half-open: 200 is not in [100,200)
        assert!(!a.overlaps(&d)); // different contig
        assert!(a.contains("1", 199));
        assert!(!a.contains("1", 200));
        assert_eq!(a.len(), 100);
    }

    #[test]
    fn thousands() {
        assert_eq!(group_thousands(0), "0");
        assert_eq!(group_thousands(999), "999");
        assert_eq!(group_thousands(1000), "1,000");
        assert_eq!(group_thousands(1234567), "1,234,567");
    }

    #[test]
    fn display_forms() {
        let r = GenomicRange::new("1", 1000, 2000, Strand::Unknown);
        assert_eq!(r.display_region(), "1:1,001-2,000");
        assert_eq!(r.ensembl_region(), "1:1001-2000");
    }
}
