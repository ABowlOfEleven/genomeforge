//! Curated pharmacogene definitions: allele-defining rsIDs, the copy-count logic
//! that turns genotypes into a diplotype/phenotype, and CPIC-derived guidance.
//!
//! Every gene below cites its authoritative source(s). Allele letters are given on
//! the reference PLUS (forward) strand (the strand consumer arrays report), with
//! the complement listed where the PGx literature commonly quotes the opposite
//! strand. Where a gene's clinically relevant variation cannot be captured from
//! SNP data, the gene is omitted, or returned with an explicit limitation.
//!
//! Sources (retrieved 2026-06):
//! - CPIC guidelines and gene allele-definition/functionality tables, cpicpgx.org
//! - PharmGKB clinical annotations, pharmgkb.org
//! - dbSNP (NCBI) for plus-strand reference/alternate alleles
//!
//! NOT MEDICAL ADVICE. Informational only; do not change any medication.

use crate::model::{Confidence, DrugGuidance};

/// A defining marker: an rsID plus the variant (non-reference, function-changing)
/// allele letters that count toward the variant, on the plus strand AND any
/// commonly-quoted complement.
pub struct Marker {
    pub rsid: &'static str,
    /// Allele letters that count as "the variant allele" for this marker.
    /// Includes the plus-strand letter and, where the literature quotes the
    /// opposite strand, its complement, so a chip reporting either is matched.
    pub variant_alleles: &'static [&'static str],
    /// Allele letters that count as the reference / non-variant allele
    /// (plus strand and complement). Used to distinguish a true reference call
    /// from an off-target/unexpected allele.
    pub reference_alleles: &'static [&'static str],
}

/// The genotype observed for one marker after resolution against the store.
pub struct MarkerCall<'a> {
    pub rsid: &'a str,
    /// Number of copies of the variant allele observed: 0, 1, or 2.
    /// For hemizygous (single-allele) calls on X/MT this is 0 or 1.
    pub variant_copies: u8,
    /// Total called alleles (1 for hemizygous, 2 for diploid).
    pub ploidy: u8,
}

/// Outcome of calling one gene from its resolved markers.
pub struct GeneCall {
    pub diplotype: String,
    pub phenotype: String,
    pub drugs: Vec<DrugGuidance>,
    pub confidence: Confidence,
    pub limitation: Option<String>,
}

/// A pharmacogene definition.
pub struct GeneDef {
    pub gene: &'static str,
    pub markers: &'static [Marker],
    /// Given the called markers (only those actually genotyped are passed),
    /// produce the gene call. `calls` is non-empty (genes with no usable markers
    /// are skipped before this is invoked).
    pub call: fn(calls: &[MarkerCall]) -> GeneCall,
}

fn g(drug: &str, guidance: &str) -> DrugGuidance {
    DrugGuidance {
        drug: drug.to_string(),
        guidance: guidance.to_string(),
        source: "CPIC".to_string(),
    }
}

/// Look up a called marker by rsID.
fn find<'a>(calls: &'a [MarkerCall], rsid: &str) -> Option<&'a MarkerCall<'a>> {
    calls.iter().find(|c| c.rsid == rsid)
}

/// Copies of the variant allele for `rsid`, or 0 if that marker was not genotyped.
fn copies(calls: &[MarkerCall], rsid: &str) -> u8 {
    find(calls, rsid).map(|c| c.variant_copies).unwrap_or(0)
}

// ---------------------------------------------------------------------------
// CYP2C19
// Source: CPIC clopidogrel guideline 2022 update (PMID 35034351); CPIC
// CYP2C19 allele definition & frequency tables; PharmGKB. dbSNP: rs4244285
// fwd G>A (*2 = A), rs4986893 fwd G>A (*3 = A), rs12248560 fwd C>T (*17 = T).
// Phenotype table: NBK379740 (CPIC 2017 assignment).
// *2 and *3 are no function; *17 is increased function.
// ---------------------------------------------------------------------------
const CYP2C19_MARKERS: &[Marker] = &[
    // *2 c.681G>A; plus strand G>A, so variant = A (literature also A).
    Marker {
        rsid: "rs4244285",
        variant_alleles: &["A"],
        reference_alleles: &["G"],
    },
    // *3 c.636G>A; plus strand G>A, variant = A.
    Marker {
        rsid: "rs4986893",
        variant_alleles: &["A"],
        reference_alleles: &["G"],
    },
    // *17 c.-806C>T; plus strand C>T, variant = T.
    Marker {
        rsid: "rs12248560",
        variant_alleles: &["T"],
        reference_alleles: &["C"],
    },
];

fn call_cyp2c19(calls: &[MarkerCall]) -> GeneCall {
    let no_fn = copies(calls, "rs4244285") + copies(calls, "rs4986893"); // *2 + *3 copies
    let no_fn = no_fn.min(2);
    let inc_fn = copies(calls, "rs12248560").min(2); // *17 copies

    // Build a star-allele diplotype label from copy counts. We cannot phase, so
    // we report the simplest consistent diplotype (loss-of-function alleles take
    // priority over *17 in labeling).
    let mut alleles: Vec<&str> = Vec::new();
    alleles.extend(std::iter::repeat_n(
        "*2",
        copies(calls, "rs4244285").min(2) as usize,
    ));
    alleles.extend(std::iter::repeat_n(
        "*3",
        copies(calls, "rs4986893").min(2) as usize,
    ));
    for _ in 0..inc_fn {
        if alleles.len() < 2 {
            alleles.push("*17");
        }
    }
    while alleles.len() < 2 {
        alleles.insert(0, "*1");
    }
    alleles.truncate(2);
    let diplotype = format!("{}/{}", alleles[0], alleles[1]);

    // CPIC phenotype assignment (NBK379740).
    let phenotype = match (no_fn, inc_fn) {
        (2, _) => "Poor metabolizer",
        (1, 0) => "Intermediate metabolizer",
        (1, _) => "Intermediate metabolizer", // *2/*17 provisional IM
        (0, 2) => "Ultrarapid metabolizer",
        (0, 1) => "Rapid metabolizer",
        _ => "Normal metabolizer",
    }
    .to_string();

    let drugs = match (no_fn, inc_fn) {
        (2, _) => vec![
            g(
                "Clopidogrel",
                "Poor metabolizer: significantly reduced formation of the active metabolite and reduced antiplatelet response, with increased risk of poor outcomes after a cardiovascular event. CPIC notes an alternative antiplatelet agent (for example prasugrel or ticagrelor) may be considered if not contraindicated. Informational only; do not change therapy.",
            ),
            g(
                "Citalopram/Escitalopram (SSRIs)",
                "Poor metabolizer: higher drug exposure is expected for SSRIs metabolized by CYP2C19, which CPIC associates with a possible need for a lower starting dose. Informational only; do not change therapy.",
            ),
            g(
                "Proton pump inhibitors",
                "Poor metabolizer: increased PPI exposure, which CPIC notes can improve efficacy but may warrant dose consideration. Informational only; do not change therapy.",
            ),
        ],
        (1, _) => vec![
            g(
                "Clopidogrel",
                "Intermediate metabolizer: reduced active-metabolite formation and reduced antiplatelet response; CPIC notes increased risk of poor cardiovascular outcomes and that an alternative agent may be considered for higher-risk indications. Informational only; do not change therapy.",
            ),
            g(
                "Proton pump inhibitors",
                "Intermediate metabolizer: modestly increased PPI exposure. Informational only; do not change therapy.",
            ),
        ],
        (0, n) if n >= 1 => vec![
            g(
                "Clopidogrel",
                "Rapid/Ultrarapid metabolizer: normal-to-increased active-metabolite formation; CPIC labels clopidogrel response as expected to be normal for these phenotypes. Informational only; do not change therapy.",
            ),
            g(
                "Citalopram/Escitalopram (SSRIs)",
                "Ultrarapid metabolizer (if *17/*17): increased metabolism may lower exposure of CYP2C19-dependent SSRIs; CPIC notes a possible reduced response. Informational only; do not change therapy.",
            ),
        ],
        _ => vec![
            g(
                "Clopidogrel",
                "Normal metabolizer: normal active-metabolite formation and expected antiplatelet response. Informational only; do not change therapy.",
            ),
        ],
    };

    GeneCall {
        diplotype,
        phenotype,
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// CYP2C9
// Source: CPIC warfarin 2017 update (PMID 28198005) and CYP2C9 allele
// functionality / activity-score tables; PharmGKB. dbSNP: rs1799853 (*2) fwd
// C>T, rs1057910 (*3) fwd A>C. *2 = decreased function (AS 0.5), *3 = no
// function (AS 0). *1 = 1.0. Phenotype by activity score.
// ---------------------------------------------------------------------------
const CYP2C9_MARKERS: &[Marker] = &[
    // *2 c.430C>T; plus strand C>T, variant = T.
    Marker {
        rsid: "rs1799853",
        variant_alleles: &["T"],
        reference_alleles: &["C"],
    },
    // *3 c.1075A>C; plus strand A>C, variant = C.
    Marker {
        rsid: "rs1057910",
        variant_alleles: &["C"],
        reference_alleles: &["A"],
    },
];

fn call_cyp2c9(calls: &[MarkerCall]) -> GeneCall {
    let star2 = copies(calls, "rs1799853").min(2); // AS 0.5 each
    let star3 = copies(calls, "rs1057910").min(2); // AS 0.0 each

    let mut alleles: Vec<&str> = Vec::new();
    alleles.extend(std::iter::repeat_n("*3", star3 as usize));
    for _ in 0..star2 {
        if alleles.len() < 2 {
            alleles.push("*2");
        }
    }
    while alleles.len() < 2 {
        alleles.insert(0, "*1");
    }
    alleles.truncate(2);
    let diplotype = format!("{}/{}", alleles[0], alleles[1]);

    // Activity score over the (up to) two reported alleles.
    let reduced = (star2 + star3).min(2);
    let normal = 2 - reduced;
    // *1 = 1.0, *2 = 0.5, *3 = 0.0.
    let n_star3 = star3.min(reduced);
    let n_star2 = reduced - n_star3;
    let score = normal as f32 * 1.0 + n_star2 as f32 * 0.5 + n_star3 as f32 * 0.0;

    let phenotype = if score >= 2.0 {
        "Normal metabolizer"
    } else if score >= 1.0 {
        "Intermediate metabolizer"
    } else {
        "Poor metabolizer"
    }
    .to_string();

    let drugs = vec![
        g(
            "Warfarin",
            "Reduced CYP2C9 function lowers S-warfarin clearance and can substantially lower the dose required; CPIC's pharmacogenetics-guided warfarin algorithm (with VKORC1) accounts for this and CPIC notes increased over-anticoagulation/bleeding risk for decreased-function carriers. Informational only; do not change therapy.",
        ),
        g(
            "Phenytoin",
            "Reduced CYP2C9 function lowers phenytoin clearance; CPIC notes a possible need for a lower maintenance dose to avoid toxicity. Informational only; do not change therapy.",
        ),
        g(
            "NSAIDs",
            "Reduced CYP2C9 function increases exposure of several NSAIDs (for example celecoxib, ibuprofen); CPIC notes considering the lowest effective dose. Informational only; do not change therapy.",
        ),
    ];

    GeneCall {
        diplotype,
        phenotype,
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// VKORC1
// Source: CPIC warfarin 2017 update (PMID 28198005); PharmGKB. dbSNP rs9923231
// is fwd-strand C>T (ref C, variant T). The PGx literature quotes it as
// -1639G>A on the OPPOSITE strand, so the sensitivity allele is "T" on the plus
// strand == "A" in literature notation. We match BOTH letters for robustness to
// 23andMe strand reporting.
// ---------------------------------------------------------------------------
const VKORC1_MARKERS: &[Marker] = &[Marker {
    rsid: "rs9923231",
    // Plus-strand variant T; literature -1639A is the same base on the minus
    // strand, so accept A too.
    variant_alleles: &["T", "A"],
    reference_alleles: &["C", "G"],
}];

fn call_vkorc1(calls: &[MarkerCall]) -> GeneCall {
    let sens = copies(calls, "rs9923231").min(2);
    // Report using literature -1639G>A notation, which is what users see.
    let geno = match sens {
        2 => "rs9923231 A/A (-1639 A/A)",
        1 => "rs9923231 G/A (-1639 G/A)",
        _ => "rs9923231 G/G (-1639 G/G)",
    }
    .to_string();
    let phenotype = match sens {
        2 => "High warfarin sensitivity (low expression)",
        1 => "Intermediate warfarin sensitivity",
        _ => "Normal warfarin sensitivity",
    }
    .to_string();

    let drugs = vec![g(
        "Warfarin",
        match sens {
            2 => "The -1639A/A genotype lowers VKORC1 expression and is associated with markedly lower warfarin dose requirements and higher early over-anticoagulation risk; CPIC's genotype-guided dosing algorithm incorporates this. Informational only; do not change therapy.",
            1 => "The -1639G/A genotype is associated with somewhat lower warfarin dose requirements; CPIC's genotype-guided dosing algorithm incorporates this. Informational only; do not change therapy.",
            _ => "The -1639G/G genotype is associated with typical warfarin dose requirements. Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype: geno,
        phenotype,
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// SLCO1B1
// Source: CPIC statins 2022 update (PMID 35152405) and simvastatin 2014 update;
// PharmGKB. dbSNP rs4149056 fwd T>C; *5 = C (decreased function). Reporting the
// single defining SNP; full star-allele set needs additional markers.
// ---------------------------------------------------------------------------
const SLCO1B1_MARKERS: &[Marker] = &[Marker {
    rsid: "rs4149056",
    variant_alleles: &["C"],
    reference_alleles: &["T"],
}];

fn call_slco1b1(calls: &[MarkerCall]) -> GeneCall {
    let dec = copies(calls, "rs4149056").min(2);
    let (diplotype, phenotype, func) = match dec {
        2 => ("*5/*5", "Low/poor function", "low"),
        1 => ("*1/*5", "Decreased function", "decreased"),
        _ => ("*1/*1", "Normal function", "normal"),
    };

    let drugs = vec![g(
        "Simvastatin",
        match func {
            "low" => "Poor SLCO1B1 transporter function (two decreased-function copies) raises simvastatin exposure and is associated with substantially higher risk of statin-associated muscle symptoms/myopathy; CPIC notes considering a lower dose or an alternative statin. Informational only; do not change therapy.",
            "decreased" => "Decreased SLCO1B1 transporter function raises simvastatin exposure and is associated with increased risk of statin-associated muscle symptoms/myopathy; CPIC notes considering a lower dose or an alternative statin. Informational only; do not change therapy.",
            _ => "Normal SLCO1B1 transporter function; typical simvastatin myopathy risk from this marker. Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype: diplotype.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// TPMT
// Source: CPIC thiopurines (TPMT/NUDT15) 2018 update (PMID 30447069); PharmGKB.
// *2 = rs1800462 (fwd C>G), *3B = rs1800460 (fwd C>T), *3C = rs1142345
// (fwd T>C). *3A = *3B + *3C in cis (both variants present). All are no
// function. Phenotype by number of no-function alleles.
// ---------------------------------------------------------------------------
const TPMT_MARKERS: &[Marker] = &[
    // *2 c.238G>C on the gene's coding (minus) strand; dbSNP fwd is C>G, variant G.
    Marker {
        rsid: "rs1800462",
        variant_alleles: &["G", "C"],
        reference_alleles: &["C", "G"],
    },
    // *3B c.460G>A; dbSNP fwd C>T, variant T.
    Marker {
        rsid: "rs1800460",
        variant_alleles: &["T", "A"],
        reference_alleles: &["C", "G"],
    },
    // *3C c.719A>G; dbSNP fwd T>C, variant C.
    Marker {
        rsid: "rs1142345",
        variant_alleles: &["C", "G"],
        reference_alleles: &["T", "A"],
    },
];

fn call_tpmt(calls: &[MarkerCall]) -> GeneCall {
    let star2 = copies(calls, "rs1800462").min(2);
    let star3b = copies(calls, "rs1800460").min(2); // *3B SNP
    let star3c = copies(calls, "rs1142345").min(2); // *3C SNP

    // *3A carries both *3B and *3C in cis. When equal copies of both are present
    // we count them as *3A (one no-function allele per copy). Any *3C copies not
    // matched by a *3B copy are *3C; unmatched *3B copies are *3B.
    let star3a = star3b.min(star3c);
    let extra_3b = star3b - star3a;
    let extra_3c = star3c - star3a;

    let no_fn = (star2 + star3a + extra_3b + extra_3c).min(2);

    let mut alleles: Vec<&str> = Vec::new();
    for _ in 0..star2 {
        if alleles.len() < 2 {
            alleles.push("*2");
        }
    }
    for _ in 0..star3a {
        if alleles.len() < 2 {
            alleles.push("*3A");
        }
    }
    for _ in 0..extra_3b {
        if alleles.len() < 2 {
            alleles.push("*3B");
        }
    }
    for _ in 0..extra_3c {
        if alleles.len() < 2 {
            alleles.push("*3C");
        }
    }
    while alleles.len() < 2 {
        alleles.insert(0, "*1");
    }
    alleles.truncate(2);
    let diplotype = format!("{}/{}", alleles[0], alleles[1]);

    let phenotype = match no_fn {
        2 => "Poor metabolizer",
        1 => "Intermediate metabolizer",
        _ => "Normal metabolizer",
    }
    .to_string();

    let drugs = vec![g(
        "Thiopurines (azathioprine, mercaptopurine)",
        match no_fn {
            2 => "Poor metabolizer: very high risk of severe, life-threatening myelosuppression at standard thiopurine doses; CPIC notes drastically reduced dosing or an alternative agent. Informational only; do not change therapy.",
            1 => "Intermediate metabolizer: increased risk of thiopurine-induced myelosuppression; CPIC notes a reduced starting dose may be considered. Informational only; do not change therapy.",
            _ => "Normal metabolizer: typical thiopurine tolerance from TPMT (NUDT15 should also be considered). Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype,
        phenotype,
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// NUDT15
// Source: CPIC thiopurines 2018 update (PMID 30447069); PharmGKB. *3 defining
// SNP rs116855232 (fwd C>T), no function. Phenotype by no-function copies.
// ---------------------------------------------------------------------------
const NUDT15_MARKERS: &[Marker] = &[Marker {
    rsid: "rs116855232",
    variant_alleles: &["T", "A"],
    reference_alleles: &["C", "G"],
}];

fn call_nudt15(calls: &[MarkerCall]) -> GeneCall {
    let no_fn = copies(calls, "rs116855232").min(2);
    let (diplotype, phenotype) = match no_fn {
        2 => ("*3/*3", "Poor metabolizer"),
        1 => ("*1/*3", "Intermediate metabolizer"),
        _ => ("*1/*1", "Normal metabolizer"),
    };

    let drugs = vec![g(
        "Thiopurines (azathioprine, mercaptopurine)",
        match no_fn {
            2 => "Poor metabolizer: very high risk of severe thiopurine-induced myelosuppression; CPIC notes substantial dose reduction or an alternative agent. Informational only; do not change therapy.",
            1 => "Intermediate metabolizer: increased risk of thiopurine-induced myelosuppression; CPIC notes a reduced starting dose may be considered. Informational only; do not change therapy.",
            _ => "Normal metabolizer: typical thiopurine tolerance from NUDT15 (TPMT should also be considered). Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype: diplotype.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// DPYD
// Source: CPIC fluoropyrimidines 2017 update (PMID 29152729); PharmGKB. Uses an
// activity-score system. Defining decreased/no-function variants captured here:
//   *2A   rs3918290   (fwd G>A) no function   (AS contribution 0)
//   *13   rs55886062  (fwd A>C) no function   (AS contribution 0)
//   D949V rs67376798  (fwd T>A) decreased     (AS contribution 0.5)
//   HapB3 rs56038477  (fwd C>T) decreased     (AS contribution 0.5)
// Normal-function allele AS = 1.0; diplotype AS in [0, 2].
// ---------------------------------------------------------------------------
const DPYD_MARKERS: &[Marker] = &[
    Marker {
        rsid: "rs3918290",
        variant_alleles: &["A", "T"],
        reference_alleles: &["G", "C"],
    },
    Marker {
        rsid: "rs55886062",
        variant_alleles: &["C", "G"],
        reference_alleles: &["A", "T"],
    },
    Marker {
        rsid: "rs67376798",
        variant_alleles: &["A", "T"],
        reference_alleles: &["T", "A"],
    },
    Marker {
        rsid: "rs56038477",
        variant_alleles: &["T", "A"],
        reference_alleles: &["C", "G"],
    },
];

fn call_dpyd(calls: &[MarkerCall]) -> GeneCall {
    // Per-allele functional value, smallest (most deleterious) first.
    // Collect one value per variant copy observed (cap total at 2 reported).
    let mut variant_vals: Vec<f32> = Vec::new();
    // *2A and *13 are no function (0.0); D949V and HapB3 are decreased (0.5).
    variant_vals.extend(std::iter::repeat_n(
        0.0,
        copies(calls, "rs3918290") as usize,
    ));
    variant_vals.extend(std::iter::repeat_n(
        0.0,
        copies(calls, "rs55886062") as usize,
    ));
    variant_vals.extend(std::iter::repeat_n(
        0.5,
        copies(calls, "rs67376798") as usize,
    ));
    variant_vals.extend(std::iter::repeat_n(
        0.5,
        copies(calls, "rs56038477") as usize,
    ));
    // Keep the most deleterious up to 2 alleles; remaining alleles are normal (1.0).
    variant_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    variant_vals.truncate(2);
    let mut allele_vals = variant_vals.clone();
    while allele_vals.len() < 2 {
        allele_vals.push(1.0);
    }
    let score: f32 = allele_vals.iter().sum();

    let phenotype = if score >= 2.0 {
        "Normal metabolizer"
    } else if score >= 1.0 {
        "Intermediate metabolizer"
    } else {
        "Poor metabolizer"
    }
    .to_string();

    let diplotype = format!("DPYD activity score {score:.1}");

    let drugs = vec![g(
        "Fluoropyrimidines (5-fluorouracil, capecitabine)",
        match phenotype.as_str() {
            "Poor metabolizer" => "Greatly reduced DPD activity (activity score 0 to <1): high risk of severe, potentially fatal fluoropyrimidine toxicity; CPIC notes avoiding these drugs or using a strongly reduced dose with monitoring. Informational only; do not change therapy.",
            "Intermediate metabolizer" => "Reduced DPD activity (activity score 1 to <2): increased risk of severe fluoropyrimidine toxicity; CPIC notes a reduced starting dose with titration. Informational only; do not change therapy.",
            _ => "Normal DPD activity for the markers tested; note this SNP panel does not capture all DPYD decreased-function variants. Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype,
        phenotype,
        drugs,
        confidence: Confidence::Moderate,
        limitation: Some(
            "DPYD has many rare decreased-function variants; only a few common ones are genotyped here, so a normal result does not rule out DPD deficiency. Informational only; do not change therapy.".to_string(),
        ),
    }
}

// ---------------------------------------------------------------------------
// G6PD (X-linked)
// Source: CPIC rasburicase 2014 guideline (PMID 24787449); PharmGKB. The G6PD
// A- deficiency (African) is tagged by rs1050828 (G6PD c.202G>A, V68M). dbSNP
// fwd C>T on the genomic plus strand (gene on minus strand), variant T == A in
// gene/literature notation. X-linked: males are hemizygous (one allele).
// ---------------------------------------------------------------------------
const G6PD_MARKERS: &[Marker] = &[Marker {
    rsid: "rs1050828",
    variant_alleles: &["T", "A"],
    reference_alleles: &["C", "G"],
}];

fn call_g6pd(calls: &[MarkerCall]) -> GeneCall {
    let call = find(calls, "rs1050828");
    let variant = call.map(|c| c.variant_copies).unwrap_or(0);
    let ploidy = call.map(|c| c.ploidy).unwrap_or(2);

    // X-linked: hemizygous (male) variant = deficient; female heterozygous =
    // variable/intermediate; female homozygous = deficient.
    let (geno, phenotype, status) = if ploidy <= 1 {
        if variant >= 1 {
            (
                "rs1050828 (A-) hemizygous variant",
                "Deficient",
                "deficient",
            )
        } else {
            ("rs1050828 hemizygous reference", "Normal", "normal")
        }
    } else {
        match variant {
            2 => ("rs1050828 A-/A- homozygous", "Deficient", "deficient"),
            1 => (
                "rs1050828 A-/B heterozygous",
                "Variable (heterozygous; X-inactivation dependent)",
                "variable",
            ),
            _ => ("rs1050828 B/B reference", "Normal", "normal"),
        }
    };

    let drugs = vec![g(
        "Rasburicase",
        match status {
            "deficient" => "Rasburicase is contraindicated in G6PD deficiency because of the risk of acute hemolytic anemia and methemoglobinemia; CPIC notes avoiding rasburicase and using an alternative. Informational only; do not change therapy.",
            "variable" => "Heterozygous females can have a deficient phenotype depending on X-inactivation; CPIC notes a quantitative G6PD enzyme activity test is recommended before rasburicase. Informational only; do not change therapy.",
            _ => "No G6PD deficiency detected at this marker; note this tags only the common A- variant and does not rule out other deficiency variants. Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype: geno.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::Moderate,
        limitation: Some(
            "G6PD is X-linked and has 200+ deficiency variants; only the common A- variant (rs1050828) is tested here, so a normal result does not rule out deficiency, and heterozygous females are not reliably classified from genotype alone. CPIC recommends a quantitative enzyme assay. Informational only; do not change therapy.".to_string(),
        ),
    }
}

// ---------------------------------------------------------------------------
// CYP3A5
// Source: CPIC tacrolimus 2015 guideline (PMID 25801146); PharmGKB. *3 =
// rs776746 (c.219-237G>A; dbSNP fwd C>T), the non-expresser allele. *1 is the
// expresser (functional) allele. Carrying any *1 makes an expresser.
// ---------------------------------------------------------------------------
const CYP3A5_MARKERS: &[Marker] = &[Marker {
    rsid: "rs776746",
    // *3 (non-expresser) is the variant; dbSNP fwd C>T, also reported as A on the
    // gene strand. We treat T (or A) as the *3 non-expresser allele.
    variant_alleles: &["T", "A"],
    reference_alleles: &["C", "G"],
}];

fn call_cyp3a5(calls: &[MarkerCall]) -> GeneCall {
    let star3 = copies(calls, "rs776746").min(2); // non-expresser copies
    let expresser_copies = 2 - star3; // *1 copies (from this single marker)
    let (diplotype, phenotype, expresser) = match star3 {
        2 => ("*3/*3", "Non-expresser (poor metabolizer)", false),
        1 => ("*1/*3", "Expresser (intermediate metabolizer)", true),
        _ => ("*1/*1", "Expresser (normal metabolizer)", true),
    };
    let _ = expresser_copies;

    let drugs = vec![g(
        "Tacrolimus",
        if expresser {
            "CYP3A5 expresser: metabolizes tacrolimus faster and typically reaches lower blood concentrations at a standard dose; CPIC notes a higher starting dose may be needed to reach target troughs. Informational only; do not change therapy."
        } else {
            "CYP3A5 non-expresser: typical tacrolimus metabolism; CPIC notes the standard recommended starting dose. Informational only; do not change therapy."
        },
    )];

    GeneCall {
        diplotype: diplotype.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// CYP4F2
// Source: CPIC warfarin 2017 update (PMID 28198005); PharmGKB. *3 = rs2108622
// (c.1297G>A, V433M; dbSNP fwd C>T). The *3/T allele modestly raises warfarin
// dose requirement. Single-SNP report.
// ---------------------------------------------------------------------------
const CYP4F2_MARKERS: &[Marker] = &[Marker {
    rsid: "rs2108622",
    variant_alleles: &["T", "A"],
    reference_alleles: &["C", "G"],
}];

fn call_cyp4f2(calls: &[MarkerCall]) -> GeneCall {
    let star3 = copies(calls, "rs2108622").min(2);
    let (diplotype, phenotype) = match star3 {
        2 => ("*3/*3", "Increased warfarin requirement"),
        1 => ("*1/*3", "Slightly increased warfarin requirement"),
        _ => ("*1/*1", "Typical warfarin requirement"),
    };

    let drugs = vec![g(
        "Warfarin",
        match star3 {
            0 => "No CYP4F2*3 allele; this marker does not raise the warfarin dose estimate. Informational only; do not change therapy.",
            _ => "CYP4F2*3 is associated with a modestly higher warfarin dose requirement; CPIC's genotype-guided dosing algorithm adds a small dose increment for carriers. Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype: diplotype.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// IFNL3/IFNL4 (near IL28B)
// Source: CPIC PEG-interferon-alpha / IFNL3 2014 guideline (PMID 24096968);
// PharmGKB. rs12979860: dbSNP fwd C>T. The C allele is the favorable-response
// allele; CC = favorable. Report by C-allele count.
// ---------------------------------------------------------------------------
const IFNL3_MARKERS: &[Marker] = &[Marker {
    // Here the "variant" we count is the favorable C allele (so variant_copies =
    // number of C alleles).
    rsid: "rs12979860",
    variant_alleles: &["C", "G"],
    reference_alleles: &["T", "A"],
}];

fn call_ifnl3(calls: &[MarkerCall]) -> GeneCall {
    let c_copies = copies(calls, "rs12979860").min(2);
    let (geno, phenotype) = match c_copies {
        2 => ("rs12979860 C/C", "Favorable response genotype"),
        1 => ("rs12979860 C/T", "Less favorable response genotype"),
        _ => ("rs12979860 T/T", "Unfavorable response genotype"),
    };

    let drugs = vec![g(
        "Peginterferon-alpha + ribavirin (HCV)",
        match c_copies {
            2 => "rs12979860 C/C is associated with a roughly two-fold higher chance of sustained virologic response to peginterferon-alpha-based HCV therapy; CPIC notes this favors a PEG-interferon-containing regimen where still used. Informational only; do not change therapy.",
            _ => "A non-C/C genotype at rs12979860 is associated with a lower chance of sustained virologic response to peginterferon-alpha-based HCV therapy; CPIC notes this is one factor in regimen selection (largely superseded by direct-acting antivirals). Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype: geno.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::High,
        limitation: None,
    }
}

// ---------------------------------------------------------------------------
// MT-RNR1 (mitochondrial)
// Source: CPIC aminoglycosides / MT-RNR1 2021 guideline (PMID 34032273);
// PharmGKB. m.1555A>G = rs267606617. Mitochondrial (haploid); presence of G is
// the increased-risk allele.
// ---------------------------------------------------------------------------
const MTRNR1_MARKERS: &[Marker] = &[Marker {
    rsid: "rs267606617",
    variant_alleles: &["G", "C"],
    reference_alleles: &["A", "T"],
}];

fn call_mtrnr1(calls: &[MarkerCall]) -> GeneCall {
    let variant = find(calls, "rs267606617")
        .map(|c| c.variant_copies)
        .unwrap_or(0);
    let risk = variant >= 1;
    let (geno, phenotype) = if risk {
        ("m.1555A>G present", "Increased ototoxicity risk")
    } else {
        ("m.1555A reference", "Normal risk")
    };

    let drugs = vec![g(
        "Aminoglycosides (gentamicin, etc.)",
        if risk {
            "The m.1555A>G variant strongly predisposes to aminoglycoside-induced hearing loss, which is usually bilateral and irreversible; CPIC notes avoiding aminoglycosides unless the infection is life-threatening and no alternative exists. Informational only; do not change therapy."
        } else {
            "The m.1555A>G risk variant was not detected; note other rarer MT-RNR1 risk variants are not tested here. Informational only; do not change therapy."
        },
    )];

    GeneCall {
        diplotype: geno.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::Moderate,
        limitation: Some(
            "MT-RNR1 is mitochondrial; heteroplasmy (mix of variant and normal mitochondria) is not resolved from array data, and only the m.1555A>G variant is tested. Informational only; do not change therapy.".to_string(),
        ),
    }
}

// ---------------------------------------------------------------------------
// CYP2D6 (SNP-only, INCOMPLETE)
// Source: CPIC CYP2D6 guidelines and gene resource; PharmGKB. *4 = rs3892097
// (c.506-1G>A; dbSNP fwd C>T) is a common no-function allele, but CYP2D6 is
// dominated by gene deletions (*5), duplications, and hybrid alleles that are
// NOT detectable from SNP/array data. We report ONLY a partial *4-based call
// with a mandatory limitation.
// ---------------------------------------------------------------------------
const CYP2D6_MARKERS: &[Marker] = &[Marker {
    rsid: "rs3892097",
    variant_alleles: &["T", "A"],
    reference_alleles: &["C", "G"],
}];

fn call_cyp2d6(calls: &[MarkerCall]) -> GeneCall {
    let star4 = copies(calls, "rs3892097").min(2);
    let (diplotype, phenotype) = match star4 {
        2 => ("*4/*4", "Poor metabolizer (partial, *4-based)"),
        1 => ("*1/*4", "Intermediate metabolizer (partial, *4-based)"),
        _ => ("*1/*1", "Presumed normal metabolizer (partial, *4-based)"),
    };

    let drugs = vec![
        g(
            "CYP2D6-metabolized drugs (e.g. codeine, tamoxifen, many antidepressants)",
            match star4 {
                2 => "Two CYP2D6*4 no-function alleles suggest a poor-metabolizer phenotype, which CPIC associates with reduced activation of prodrugs (for example codeine) and higher exposure of active drugs; HOWEVER this SNP-only call is incomplete (see limitation). Informational only; do not change therapy.",
                1 => "One CYP2D6*4 no-function allele; phenotype depends on the second allele, which array data cannot fully resolve (see limitation). Informational only; do not change therapy.",
                _ => "No CYP2D6*4 allele detected, but array data cannot rule out other reduced-function alleles, deletions, or duplications, so a normal phenotype is NOT confirmed (see limitation). Informational only; do not change therapy.",
            },
        ),
    ];

    GeneCall {
        diplotype: diplotype.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::Low,
        limitation: Some(
            "CYP2D6 cannot be accurately called from array/SNP data: clinically important alleles include whole-gene deletions, duplications (which change activity score), and hybrid genes that SNP genotyping does not detect. This result reflects only the *4 SNP (rs3892097) and may be wrong; a validated CYP2D6 test (with copy-number analysis) is required. Informational only; do not change therapy.".to_string(),
        ),
    }
}

// ---------------------------------------------------------------------------
// UGT1A1 (tag SNP, INCOMPLETE)
// Source: CPIC atazanavil / UGT1A1 guideline; PharmGKB. *28 is a TA-repeat
// (A(TA)7TAA) promoter insertion that SNP arrays cannot read directly. rs887829
// (fwd C>T) is a well-established tag in linkage disequilibrium with *28; T is
// the reduced-function tag. Report with a mandatory limitation.
// ---------------------------------------------------------------------------
const UGT1A1_MARKERS: &[Marker] = &[Marker {
    rsid: "rs887829",
    variant_alleles: &["T", "A"],
    reference_alleles: &["C", "G"],
}];

fn call_ugt1a1(calls: &[MarkerCall]) -> GeneCall {
    let tag = copies(calls, "rs887829").min(2);
    let (geno, phenotype) = match tag {
        2 => (
            "rs887829 T/T (tag for *28/*28)",
            "Likely reduced function (poor metabolizer)",
        ),
        1 => (
            "rs887829 C/T (tag for *1/*28)",
            "Likely intermediate function",
        ),
        _ => ("rs887829 C/C (tag for *1/*1)", "Likely normal function"),
    };

    let drugs = vec![g(
        "Atazanavir; irinotecan",
        match tag {
            0 => "The reduced-function tag was not detected; reduced UGT1A1 function (and elevated bilirubin / Gilbert syndrome) is less likely, but the *28 repeat is not directly genotyped here. Informational only; do not change therapy.",
            _ => "The rs887829 T tag is associated with reduced UGT1A1 activity (Gilbert-type), which CPIC links to a higher chance of atazanavir-associated hyperbilirubinemia and, for irinotecan, increased toxicity risk; the actual *28 repeat is not directly genotyped here. Informational only; do not change therapy.",
        },
    )];

    GeneCall {
        diplotype: geno.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::Low,
        limitation: Some(
            "UGT1A1*28 is a TA-repeat (promoter) variant that SNP arrays cannot read directly; this call uses the rs887829 tag SNP, which is correlated with but not identical to *28, and other UGT1A1 alleles (for example *6) are not tested. Informational only; do not change therapy.".to_string(),
        ),
    }
}

// ---------------------------------------------------------------------------
// HLA-B*57:01 (tag SNP, INCOMPLETE)
// Source: CPIC abacavir / HLA-B 2014 guideline (PMID 24561393); PharmGKB.
// HLA-B*57:01 status determines abacavir hypersensitivity risk. rs2395029
// (HCP5 G>T, A on the gene strand) is a tag in strong LD with *57:01 but is NOT
// the allele itself. Report with a mandatory limitation.
// ---------------------------------------------------------------------------
const HLAB_MARKERS: &[Marker] = &[Marker {
    rsid: "rs2395029",
    // The minor allele (G on dbSNP fwd; reported as G/T) tags HLA-B*57:01.
    variant_alleles: &["G", "C"],
    reference_alleles: &["T", "A"],
}];

fn call_hlab(calls: &[MarkerCall]) -> GeneCall {
    let tag = copies(calls, "rs2395029").min(2);
    let positive = tag >= 1;
    let (geno, phenotype) = if positive {
        (
            "rs2395029 tag positive (suggests HLA-B*57:01)",
            "Possible HLA-B*57:01 positive",
        )
    } else {
        (
            "rs2395029 tag negative",
            "HLA-B*57:01 not suggested",
        )
    };

    let drugs = vec![g(
        "Abacavir",
        if positive {
            "The rs2395029 tag suggests possible HLA-B*57:01 carriage, which CPIC associates with high risk of a serious abacavir hypersensitivity reaction; CPIC notes abacavir is contraindicated in confirmed HLA-B*57:01-positive individuals. The tag is NOT a confirmatory HLA test (see limitation). Informational only; do not change therapy."
        } else {
            "The rs2395029 tag was not detected, which makes HLA-B*57:01 carriage less likely, but the tag does not rule it out; a validated HLA-B*57:01 test is required before abacavir (see limitation). Informational only; do not change therapy."
        },
    )];

    GeneCall {
        diplotype: geno.to_string(),
        phenotype: phenotype.to_string(),
        drugs,
        confidence: Confidence::Low,
        limitation: Some(
            "HLA-B*57:01 cannot be determined from this array marker: rs2395029 is only a tag SNP in linkage disequilibrium with the allele (LD is incomplete and population-dependent), not the HLA allele itself. A validated HLA-B*57:01 typing test is required before abacavir. Informational only; do not change therapy.".to_string(),
        ),
    }
}

/// All pharmacogene definitions, in a stable reporting order.
pub const GENES: &[GeneDef] = &[
    GeneDef {
        gene: "CYP2C19",
        markers: CYP2C19_MARKERS,
        call: call_cyp2c19,
    },
    GeneDef {
        gene: "CYP2C9",
        markers: CYP2C9_MARKERS,
        call: call_cyp2c9,
    },
    GeneDef {
        gene: "VKORC1",
        markers: VKORC1_MARKERS,
        call: call_vkorc1,
    },
    GeneDef {
        gene: "CYP4F2",
        markers: CYP4F2_MARKERS,
        call: call_cyp4f2,
    },
    GeneDef {
        gene: "SLCO1B1",
        markers: SLCO1B1_MARKERS,
        call: call_slco1b1,
    },
    GeneDef {
        gene: "TPMT",
        markers: TPMT_MARKERS,
        call: call_tpmt,
    },
    GeneDef {
        gene: "NUDT15",
        markers: NUDT15_MARKERS,
        call: call_nudt15,
    },
    GeneDef {
        gene: "DPYD",
        markers: DPYD_MARKERS,
        call: call_dpyd,
    },
    GeneDef {
        gene: "G6PD",
        markers: G6PD_MARKERS,
        call: call_g6pd,
    },
    GeneDef {
        gene: "CYP3A5",
        markers: CYP3A5_MARKERS,
        call: call_cyp3a5,
    },
    GeneDef {
        gene: "IFNL3",
        markers: IFNL3_MARKERS,
        call: call_ifnl3,
    },
    GeneDef {
        gene: "MT-RNR1",
        markers: MTRNR1_MARKERS,
        call: call_mtrnr1,
    },
    GeneDef {
        gene: "CYP2D6",
        markers: CYP2D6_MARKERS,
        call: call_cyp2d6,
    },
    GeneDef {
        gene: "UGT1A1",
        markers: UGT1A1_MARKERS,
        call: call_ugt1a1,
    },
    GeneDef {
        gene: "HLA-B",
        markers: HLAB_MARKERS,
        call: call_hlab,
    },
];
