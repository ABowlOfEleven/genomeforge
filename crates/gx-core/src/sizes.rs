//! Human chromosome lengths, used for browser extent / ideogram scaling.
//!
//! Values are GRCh38 primary-assembly lengths. GRCh37 differs by well under a
//! percent — negligible for framing the view — so both builds share this table
//! for now; precise per-build cytobands can come later.

use crate::assembly::Assembly;

/// `(contig, length_bp)` for the 24 nuclear chromosomes + MT, GRCh38.
const GRCH38: &[(&str, u64)] = &[
    ("1", 248_956_422),
    ("2", 242_193_529),
    ("3", 198_295_559),
    ("4", 190_214_555),
    ("5", 181_538_259),
    ("6", 170_805_979),
    ("7", 159_345_973),
    ("8", 145_138_636),
    ("9", 138_394_717),
    ("10", 133_797_422),
    ("11", 135_086_622),
    ("12", 133_275_309),
    ("13", 114_364_328),
    ("14", 107_043_718),
    ("15", 101_991_189),
    ("16", 90_338_345),
    ("17", 83_257_441),
    ("18", 80_373_285),
    ("19", 58_617_616),
    ("20", 64_444_167),
    ("21", 46_709_983),
    ("22", 50_818_468),
    ("X", 156_040_895),
    ("Y", 57_227_415),
    ("MT", 16_569),
];

/// Length of a chromosome in bases, or `None` for unknown contigs.
pub fn chrom_length(_assembly: Assembly, contig: &str) -> Option<u64> {
    GRCH38
        .iter()
        .find(|(name, _)| *name == contig)
        .map(|(_, len)| *len)
}

/// The canonical chromosome ordering (`1..22, X, Y, MT`).
pub fn chrom_order() -> impl Iterator<Item = &'static str> {
    GRCH38.iter().map(|(name, _)| *name)
}
