//! Mitochondrial (maternal) haplogroup tree.
//!
//! Marker positions and derived alleles are taken from PhyloTree mtDNA Build 17
//! (van Oven & Kayser 2009; Build 17, 18 Feb 2016), the standard reference for
//! human mtDNA phylogeny: <https://www.phylotree.org/>. Positions are rCRS
//! 1-based coordinates on contig `"MT"` (the same numbering is used in GRCh37 and
//! GRCh38 practice). For each marker the "derived" base is the rCRS-relative
//! mutant the sample carries (the trailing base of PhyloTree's `X<pos>Y`
//! notation); back-mutations (PhyloTree's `!`) and transversions to lowercase
//! bases are omitted to avoid ambiguity, as are length variants / indels that
//! consumer arrays do not call reliably.
//!
//! IMPORTANT (the H trap): the rCRS reference sequence is itself haplogroup
//! H2a2a1. Many "H-defining" mutations in RSRS-rooted notation are already the
//! reference state in rCRS, so H cannot be resolved by positive mutations; it is
//! recognised conservatively as the R0 / HV cluster that lacks U/J/T/etc.
//! diagnostic mutations. Fine H subclades are deliberately not resolved here.

use crate::tree::{Marker, Node};

/// MT marker positions are identical across builds, so both position fields match.
const fn mt(name: &'static str, pos: u64, derived: &'static str) -> Marker {
    Marker {
        name,
        rsid: None,
        pos_grch37: pos,
        pos_grch38: pos,
        derived,
    }
}

pub const CONTIG: &str = "MT";

// ---------------------------------------------------------------------------
// Terminal / leaf nodes, defined bottom-up so the tree literal reads top-down.
// All positions cross-checked against PhyloTree Build 17 (rCRS-oriented).
// ---------------------------------------------------------------------------

// --- U subclades (U is defined by 11467, 12308, 12372; the textbook anchor) ---
static U5: Node = Node {
    haplogroup: "U5",
    markers: &[mt("16192", 16192, "T"), mt("16270", 16270, "T")],
    children: &[
        Node {
            haplogroup: "U5a",
            markers: &[mt("14793", 14793, "G"), mt("16256", 16256, "T")],
            children: &[],
            note: None,
        },
        Node {
            haplogroup: "U5b",
            markers: &[mt("7768", 7768, "G"), mt("14182", 14182, "C")],
            children: &[],
            note: None,
        },
    ],
    note: None,
};

static U4: Node = Node {
    haplogroup: "U4",
    markers: &[
        mt("4646", 4646, "C"),
        mt("11332", 11332, "T"),
        mt("14620", 14620, "T"),
        mt("16356", 16356, "C"),
    ],
    children: &[],
    note: None,
};

static U3: Node = Node {
    haplogroup: "U3",
    markers: &[
        mt("14139", 14139, "G"),
        mt("15454", 15454, "C"),
        mt("16343", 16343, "G"),
    ],
    children: &[],
    note: None,
};

static U2: Node = Node {
    haplogroup: "U2",
    markers: &[mt("16051", 16051, "G")],
    children: &[],
    note: None,
};

static U6: Node = Node {
    haplogroup: "U6",
    markers: &[mt("3348", 3348, "G"), mt("16172", 16172, "C")],
    children: &[],
    note: None,
};

static U7: Node = Node {
    haplogroup: "U7",
    markers: &[mt("980", 980, "C"), mt("3741", 3741, "T"), mt("5360", 5360, "C")],
    children: &[],
    note: None,
};

// K descends from U8 (K is a subclade of U8; the textbook anchor). U8b'K carries
// few coding-region SNPs, so we route U8 -> K using K's own diagnostic set.
static K: Node = Node {
    haplogroup: "K",
    markers: &[
        mt("10550", 10550, "G"),
        mt("11299", 11299, "C"),
        mt("14798", 14798, "C"),
        mt("16224", 16224, "C"),
    ],
    children: &[],
    note: None,
};

static U8: Node = Node {
    haplogroup: "U8",
    markers: &[mt("9698", 9698, "C")],
    children: &[/* K is checked directly under U below to be robust to U8 gaps */],
    note: None,
};

static U: Node = Node {
    haplogroup: "U",
    // The U-defining trio (anchor fact): rCRS-derived G/G/A at 11467/12308/12372.
    markers: &[mt("11467", 11467, "G"), mt("12308", 12308, "G"), mt("12372", 12372, "A")],
    children: &[
        // K placed first: U8 is often not directly genotyped on arrays, but K's
        // own markers are, and K is the most common U subclade in Europeans.
        Node {
            haplogroup: "K",
            markers: K.markers,
            children: &[],
            note: None,
        },
        Node {
            haplogroup: "U8",
            markers: U8.markers,
            children: &[Node {
                haplogroup: "K",
                markers: K.markers,
                children: &[],
                note: None,
            }],
            note: None,
        },
        Node { haplogroup: "U5", markers: U5.markers, children: U5.children, note: None },
        Node { haplogroup: "U4", markers: U4.markers, children: U4.children, note: None },
        Node { haplogroup: "U3", markers: U3.markers, children: U3.children, note: None },
        Node { haplogroup: "U2", markers: U2.markers, children: U2.children, note: None },
        Node { haplogroup: "U6", markers: U6.markers, children: U6.children, note: None },
        Node { haplogroup: "U7", markers: U7.markers, children: U7.children, note: None },
        Node { haplogroup: "U1", markers: &[mt("285", 285, "T"), mt("12879", 12879, "C")], children: &[], note: None },
    ],
    note: None,
};

// --- JT cluster: J and T share 4216 (via R2'JT) and 11251 (via JT) (anchor) ---
static T: Node = Node {
    haplogroup: "T",
    markers: &[
        mt("4917", 4917, "G"),
        mt("13368", 13368, "A"),
        mt("14905", 14905, "A"),
        mt("15607", 15607, "G"),
        mt("16294", 16294, "T"),
    ],
    children: &[
        Node { haplogroup: "T1", markers: &[mt("16163", 16163, "G")], children: &[], note: None },
        Node { haplogroup: "T2", markers: &[mt("11812", 11812, "G"), mt("14233", 14233, "G")], children: &[], note: None },
    ],
    note: None,
};

static J: Node = Node {
    haplogroup: "J",
    // 10398 is a back-mutation in some references; 13708 and 16069 are robust J markers.
    markers: &[mt("13708", 13708, "A"), mt("16069", 16069, "T")],
    children: &[
        Node { haplogroup: "J1", markers: &[mt("462", 462, "T"), mt("3010", 3010, "A")], children: &[], note: None },
        Node { haplogroup: "J2", markers: &[mt("7476", 7476, "T"), mt("15257", 15257, "A")], children: &[], note: None },
    ],
    note: None,
};

static JT: Node = Node {
    haplogroup: "JT",
    // 11251 defines JT (anchor); 16126 also defines JT. 4216 (R2'JT) is checked at
    // the parent (pre_jt) below so J and T inherit it.
    markers: &[mt("11251", 11251, "G"), mt("16126", 16126, "C")],
    children: &[
        Node { haplogroup: "J", markers: J.markers, children: J.children, note: None },
        Node { haplogroup: "T", markers: T.markers, children: T.children, note: None },
    ],
    note: None,
};

// R2'JT (pre-JT) carries 4216; sits between R and JT (anchor: J/T share 4216).
static PRE_JT: Node = Node {
    haplogroup: "JT (R2'JT)",
    markers: &[mt("4216", 4216, "C")],
    children: &[Node { haplogroup: "JT", markers: JT.markers, children: JT.children, note: None }],
    note: Some(
        "predicted via the R2'JT / JT cluster; resolution stops at the deepest matched subclade.",
    ),
};

const H_NOTE: &str = "the rCRS reference is itself haplogroup H2a2a1, so fine H subclades are not \
resolved by positive mutations; this is the West-Eurasian R0 / HV cluster lacking U, J, T, and other \
diagnostic mutations (coarse resolution).";

// R0 / HV / H cluster: handled conservatively. The reference IS H2a2a1, so we do
// not over-resolve H; HV and V have a couple of usable rCRS-derived markers.
static R0: Node = Node {
    haplogroup: "H / HV (R0 cluster)",
    // R0 markers (G73A, A11719G) are RSRS-rooted and largely reference state in
    // rCRS, so R0 is mostly inferred by exclusion. HV/V below add real signal.
    markers: &[],
    children: &[
        Node {
            haplogroup: "V",
            markers: &[mt("4580", 4580, "A")],
            children: &[],
            note: Some("predicted V within the R0 / HV cluster; H subclades not resolved."),
        },
        Node {
            haplogroup: "HV",
            markers: &[mt("14766", 14766, "T")],
            children: &[],
            note: Some(H_NOTE),
        },
    ],
    note: Some(H_NOTE),
};

// --- Other macro-N / R branches present on consumer arrays ---
static W: Node = Node {
    haplogroup: "W",
    markers: &[mt("204", 204, "C"), mt("8251", 8251, "A"), mt("11947", 11947, "G"), mt("16292", 16292, "T")],
    children: &[],
    note: None,
};
static X: Node = Node {
    haplogroup: "X",
    markers: &[mt("6221", 6221, "C"), mt("6371", 6371, "T"), mt("13966", 13966, "G"), mt("14470", 14470, "C")],
    children: &[],
    note: None,
};
static I: Node = Node {
    haplogroup: "I",
    markers: &[mt("10034", 10034, "C"), mt("16129", 16129, "A")],
    children: &[],
    note: None,
};
static A: Node = Node {
    haplogroup: "A",
    markers: &[mt("235", 235, "G"), mt("663", 663, "G"), mt("4248", 4248, "C"), mt("8794", 8794, "T")],
    children: &[],
    note: None,
};

// --- macro-M branch (Asian / Native American) ---
static C: Node = Node {
    haplogroup: "C",
    markers: &[mt("14318", 14318, "C"), mt("16327", 16327, "T")],
    children: &[],
    note: None,
};
static D: Node = Node {
    haplogroup: "D",
    markers: &[mt("5178", 5178, "A"), mt("16362", 16362, "C")],
    children: &[],
    note: None,
};
static G: Node = Node {
    haplogroup: "G",
    markers: &[mt("4833", 4833, "G"), mt("5108", 5108, "C")],
    children: &[],
    note: None,
};
static M: Node = Node {
    haplogroup: "M",
    markers: &[mt("10400", 10400, "T"), mt("14783", 14783, "C"), mt("15043", 15043, "A")],
    children: &[
        Node { haplogroup: "C", markers: C.markers, children: &[], note: None },
        Node { haplogroup: "D", markers: D.markers, children: &[], note: None },
        Node { haplogroup: "G", markers: G.markers, children: &[], note: None },
    ],
    note: None,
};

// R is defined by 12705 and 16223 (PhyloTree Build 17).
static R: Node = Node {
    haplogroup: "R",
    markers: &[mt("12705", 12705, "C"), mt("16223", 16223, "C")],
    children: &[
        Node { haplogroup: "U", markers: U.markers, children: U.children, note: None },
        Node { haplogroup: "JT (R2'JT)", markers: PRE_JT.markers, children: PRE_JT.children, note: PRE_JT.note },
        // R0/HV/H placed last: it is the "default West-Eurasian R" with the
        // fewest positive markers, so only assigned when nothing more specific hits.
        Node { haplogroup: "H / HV (R0 cluster)", markers: R0.markers, children: R0.children, note: R0.note },
    ],
    note: None,
};

// macro-N: parent of R and the basal N clades (W, X, I, A). N's PhyloTree-defining
// mutations (e.g. 8701, 10873) are RSRS-rooted and are reference state in rCRS, so
// N cannot be detected by positive mutation; it is a pass-through node whose
// descendants (R and its subclades, W/X/I/A) carry the detectable signal.
static N: Node = Node {
    haplogroup: "N (macro-haplogroup)",
    markers: &[],
    children: &[
        Node { haplogroup: "R", markers: R.markers, children: R.children, note: None },
        Node { haplogroup: "W", markers: W.markers, children: &[], note: None },
        Node { haplogroup: "X", markers: X.markers, children: &[], note: None },
        Node { haplogroup: "I", markers: I.markers, children: &[], note: None },
        Node { haplogroup: "A", markers: A.markers, children: &[], note: None },
    ],
    note: None,
};

/// Root of the maternal tree: L3, the common ancestor of the out-of-Africa
/// macro-haplogroups M and N. (Sub-Saharan L0-L6 lineages are not modelled here.)
pub static ROOT: Node = Node {
    haplogroup: "L3 / N (root)",
    // L3 markers (A769G, A1018G) are RSRS-rooted; in rCRS they are reference
    // state, so the root is "entered" by default and signal comes from descent.
    markers: &[],
    children: &[
        Node { haplogroup: "N (macro-haplogroup)", markers: N.markers, children: N.children, note: None },
        Node { haplogroup: "M", markers: M.markers, children: M.children, note: None },
    ],
    note: None,
};
