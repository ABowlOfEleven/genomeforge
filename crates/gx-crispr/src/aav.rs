//! AAV cargo capacity planning.
//!
//! Single-stranded AAV packages roughly **4.7 kb** between the two ITRs
//! (ITRs ~145 bp each). This helps a user check whether a promoter + transgene
//! + polyA cassette will fit.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AavCargo {
    pub itr_bp: u32,
    pub promoter_bp: u32,
    pub transgene_bp: u32,
    pub polya_bp: u32,
}

impl Default for AavCargo {
    fn default() -> Self {
        Self {
            itr_bp: 145,
            promoter_bp: 600, // ~CMV
            transgene_bp: 720,
            polya_bp: 250, // ~SV40 late polyA
        }
    }
}

impl AavCargo {
    /// Practical ssAAV packaging limit in bp (ITR-to-ITR).
    pub const CAPACITY: u32 = 4700;

    /// Total cargo size (both ITRs + cassette).
    pub fn total(&self) -> u32 {
        self.itr_bp.saturating_mul(2) + self.promoter_bp + self.transgene_bp + self.polya_bp
    }

    /// Fraction of capacity used.
    pub fn fraction(&self) -> f64 {
        self.total() as f64 / Self::CAPACITY as f64
    }

    pub fn fits(&self) -> bool {
        self.total() <= Self::CAPACITY
    }

    /// Headroom (negative if over capacity).
    pub fn headroom(&self) -> i64 {
        Self::CAPACITY as i64 - self.total() as i64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_fits() {
        let c = AavCargo::default();
        assert!(c.fits());
        assert!(c.fraction() < 1.0);
        assert_eq!(c.total(), 2 * 145 + 600 + 720 + 250);
    }

    #[test]
    fn oversize_transgene_does_not_fit() {
        let c = AavCargo {
            transgene_bp: 5000,
            ..AavCargo::default()
        };
        assert!(!c.fits());
        assert!(c.headroom() < 0);
    }
}
