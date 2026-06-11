//! Public result types returned by [`crate::report`].

/// How much trust to place in a single gene's call, given consumer-array data.
///
/// This reflects how completely the genotyped markers capture the gene's
/// clinically relevant variation, NOT the strength of the underlying CPIC
/// evidence. A gene whose phenotype-defining alleles are all simple SNPs that the
/// array reads directly can be `High`; a gene where the array sees only a tag SNP
/// or a fraction of the defining alleles is `Moderate` or `Low`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Confidence {
    High,
    Moderate,
    Low,
}

/// One drug's CPIC-derived implication for a given diplotype/phenotype.
///
/// `guidance` is plain text, intentionally concise, and carries an implication
/// (e.g. increased risk, consider an alternative) rather than a prescription. It
/// is informational only and is not a recommendation to change therapy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DrugGuidance {
    /// Drug or drug-class name, e.g. `"Clopidogrel"`.
    pub drug: String,
    /// Concise CPIC-derived implication, plain text.
    pub guidance: String,
    /// Attribution, e.g. `"CPIC"`.
    pub source: String,
}

/// The call for a single pharmacogene.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneResult {
    /// Gene symbol, e.g. `"CYP2C19"`.
    pub gene: String,
    /// Diplotype (`"*1/*2"`) for star-allele genes, or a genotype string
    /// (`"rs9923231 A/A"`) for genes reported by a single defining SNP.
    pub diplotype: String,
    /// CPIC metabolizer phenotype or equivalent functional status,
    /// e.g. `"Intermediate metabolizer"`.
    pub phenotype: String,
    /// Per-drug guidance entries derived from CPIC for this phenotype.
    pub drugs: Vec<DrugGuidance>,
    /// Confidence in the call given consumer-array coverage.
    pub confidence: Confidence,
    /// Caveat the UI MUST surface (e.g. structural variants not callable);
    /// `None` only when the genotyped markers fully capture the gene.
    pub limitation: Option<String>,
    /// rsIDs that were actually found and used to make this call.
    pub markers_used: Vec<String>,
}
