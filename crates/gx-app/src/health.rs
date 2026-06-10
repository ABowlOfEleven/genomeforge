//! A small curated set of well-known health- and trait-associated SNPs, offered
//! as one-click jumps in the genome browser (the analogue of the Phenotype
//! workspace's "Example scores").
//!
//! Selection criteria: each is widely genotyped on consumer chips (23andMe /
//! AncestryDNA) so it is likely to be present in a loaded genome, and each has a
//! clear, well-established trait association. This is an educational starting
//! point for exploration — NOT a clinical panel, and not exhaustive. Full
//! interpretation (e.g. APOE ε-status or thrombophilia risk) needs both sites of
//! a haplotype, zygosity, and clinical context; clicking simply jumps to the SNP
//! and pulls its public annotation.

pub struct HealthSnp {
    pub rsid: &'static str,
    pub gene: &'static str,
    /// Short, layman-friendly trait label.
    pub trait_: &'static str,
}

/// Common health/trait SNPs, roughly grouped (clotting, metabolism, traits).
pub const HEALTH_SNPS: &[HealthSnp] = &[
    HealthSnp { rsid: "rs429358", gene: "APOE", trait_: "Alzheimer's & heart risk (ε4 site)" },
    HealthSnp { rsid: "rs7412", gene: "APOE", trait_: "Alzheimer's & heart risk (ε2 site)" },
    HealthSnp { rsid: "rs6025", gene: "F5", trait_: "Factor V Leiden — clotting" },
    HealthSnp { rsid: "rs1799963", gene: "F2", trait_: "Prothrombin G20210A — clotting" },
    HealthSnp { rsid: "rs1801133", gene: "MTHFR", trait_: "MTHFR C677T — folate metabolism" },
    HealthSnp { rsid: "rs1800562", gene: "HFE", trait_: "Hemochromatosis C282Y — iron overload" },
    HealthSnp { rsid: "rs4988235", gene: "LCT", trait_: "Lactase persistence (milk digestion)" },
    HealthSnp { rsid: "rs671", gene: "ALDH2", trait_: "Alcohol flush reaction" },
    HealthSnp { rsid: "rs1815739", gene: "ACTN3", trait_: "Muscle fibre type (R577X)" },
    HealthSnp { rsid: "rs762551", gene: "CYP1A2", trait_: "Caffeine metabolism speed" },
    HealthSnp { rsid: "rs9939609", gene: "FTO", trait_: "Body-weight / BMI tendency" },
    HealthSnp { rsid: "rs334", gene: "HBB", trait_: "Sickle-cell / malaria resistance" },
];
