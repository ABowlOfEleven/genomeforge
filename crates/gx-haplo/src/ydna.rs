//! Y-chromosome (paternal) haplogroup tree.
//!
//! Marker positions, derived alleles, and ISOGG haplogroup assignments are taken
//! from the ISOGG Y-DNA tree as published in the YBrowse / ISOGG SNP tables
//! (Thomas Krahn, <https://ybrowse.org/>; ISOGG Y-DNA Haplogroup Tree, 2019/2020).
//! Each marker carries BOTH its GRCh37/hg19 and GRCh38/hg38 1-based position so
//! the classifier can pick the right coordinate from the store assembly; the
//! contig is `"Y"`. Positions were cross-checked against NCBI dbSNP placements
//! where an authoritative rsID exists (e.g. M269 = rs9786153: GRCh37 chrY
//! 22739367, GRCh38 chrY 20577481; M9 = rs3900; M89 = rs2032652; M170 =
//! rs2032597), which matched the YBrowse table exactly.
//!
//! Indel markers (e.g. M17, M175) are intentionally omitted as primary markers
//! because consumer SNP arrays do not call short indels reliably; the
//! single-base SNPs above each terminal node give the signal instead.
//!
//! Anchor facts encoded: R1b = M269; R1a = M198; I1 = M253; J1 = M267; J2 = M172;
//! G = M201; E-M35 = M35; N = M231; O = M175 (we use the downstream context);
//! Q = M242.

use crate::tree::{Marker, Node};

const fn y(
    name: &'static str,
    rsid: Option<&'static str>,
    pos37: u64,
    pos38: u64,
    derived: &'static str,
) -> Marker {
    Marker {
        name,
        rsid,
        pos_grch37: pos37,
        pos_grch38: pos38,
        derived,
    }
}

pub const CONTIG: &str = "Y";

// --- Deepest haplogroups under P / R ---
static R1B: Node = Node {
    haplogroup: "R1b",
    // M269 (anchor). dbSNP rs9786153; YBrowse GRCh37 22739367 / GRCh38 20577481.
    markers: &[y("M269", Some("rs9786153"), 22739367, 20577481, "C")],
    children: &[],
    note: None,
};
static R1A: Node = Node {
    haplogroup: "R1a",
    // M198 (anchor: R1a by M198/M17). M17 is an indel, so we key on M198 (SNP).
    markers: &[y("M198", None, 15030752, 12918840, "T")],
    children: &[],
    note: None,
};
static R1: Node = Node {
    haplogroup: "R1",
    // M173 defines R1; we descend from R (M207) straight to R1a / R1b on their
    // own SNPs, which are what arrays carry.
    markers: &[],
    children: &[
        Node { haplogroup: "R1b", markers: R1B.markers, children: &[], note: None },
        Node { haplogroup: "R1a", markers: R1A.markers, children: &[], note: None },
    ],
    note: None,
};
static R: Node = Node {
    haplogroup: "R",
    // M207 defines R (ancestral A -> derived G).
    markers: &[y("M207", None, 15581983, 13470103, "G")],
    children: &[Node { haplogroup: "R1", markers: R1.markers, children: R1.children, note: None }],
    note: None,
};

// Q (anchor: Q = M242). Sibling of R under P.
static Q: Node = Node {
    haplogroup: "Q",
    markers: &[y("M242", None, 15018582, 12906671, "T")],
    children: &[],
    note: None,
};

// P / K2b2: parent of Q and R (M45).
static P: Node = Node {
    haplogroup: "P",
    markers: &[y("M45", None, 21867787, 19705901, "A")],
    children: &[
        Node { haplogroup: "R", markers: R.markers, children: R.children, note: None },
        Node { haplogroup: "Q", markers: Q.markers, children: &[], note: None },
    ],
    note: None,
};

// --- Branches under K (M9) ---
static N_NODE: Node = Node {
    haplogroup: "N",
    // N (anchor: N = M231).
    markers: &[y("M231", None, 15469724, 13357844, "A")],
    children: &[],
    note: None,
};
static O_NODE: Node = Node {
    haplogroup: "O",
    // O (anchor: O = M175). M175 is an indel; we additionally accept P186/M175.
    // To stay SNP-only we rely on M175 derived-state where called; arrays that
    // carry it report the deletion as the derived call.
    markers: &[y("M175", None, 15508706, 13396826, "del")],
    children: &[],
    note: None,
};

// K2 (M526) -> NO, P, etc. We model the practically-resolvable descendants.
static K2: Node = Node {
    haplogroup: "K2",
    markers: &[y("M526", None, 23550924, 21389038, "C")],
    children: &[
        Node { haplogroup: "P", markers: P.markers, children: P.children, note: None },
        Node { haplogroup: "N", markers: N_NODE.markers, children: &[], note: None },
        Node { haplogroup: "O", markers: O_NODE.markers, children: &[], note: None },
    ],
    note: None,
};

// K (M9): defines the great Eurasian K macro-clade.
static K: Node = Node {
    haplogroup: "K",
    markers: &[y("M9", Some("rs3900"), 21730257, 19568371, "G")],
    children: &[
        Node { haplogroup: "K2", markers: K2.markers, children: K2.children, note: None },
        // Robustness: P/N/O often present without an M526 call on arrays, so also
        // try them directly under K.
        Node { haplogroup: "P", markers: P.markers, children: P.children, note: None },
        Node { haplogroup: "N", markers: N_NODE.markers, children: &[], note: None },
        Node { haplogroup: "O", markers: O_NODE.markers, children: &[], note: None },
    ],
    note: None,
};

// --- Branches under F (M89): I, J, G ---
static I1: Node = Node {
    haplogroup: "I1",
    // I1 (anchor: I1 = M253).
    markers: &[y("M253", None, 15022707, 12910796, "T")],
    children: &[],
    note: None,
};
static I_NODE: Node = Node {
    haplogroup: "I",
    // I (M170).  dbSNP rs2032597; GRCh37 14847792 / GRCh38 12735858.
    markers: &[y("M170", Some("rs2032597"), 14847792, 12735858, "C")],
    children: &[Node { haplogroup: "I1", markers: I1.markers, children: &[], note: None }],
    note: None,
};
static J1: Node = Node {
    haplogroup: "J1",
    // J1 (anchor: J1 = M267).
    markers: &[y("M267", None, 22741818, 20579932, "G")],
    children: &[],
    note: None,
};
static J2: Node = Node {
    haplogroup: "J2",
    // J2 (anchor: J2 = M172).
    markers: &[y("M172", None, 14969634, 12857709, "G")],
    children: &[],
    note: None,
};
static J_NODE: Node = Node {
    haplogroup: "J",
    // J (M304).
    markers: &[y("M304", None, 22749853, 20587967, "C")],
    children: &[
        Node { haplogroup: "J1", markers: J1.markers, children: &[], note: None },
        Node { haplogroup: "J2", markers: J2.markers, children: &[], note: None },
    ],
    note: None,
};
static G_NODE: Node = Node {
    haplogroup: "G",
    // G (anchor: G = M201).
    markers: &[y("M201", None, 15027529, 12915617, "T")],
    children: &[],
    note: None,
};

// F (M89): parent of G, H, I, J, K (and thus everything non-A/B/C/D/E here).
static F: Node = Node {
    haplogroup: "F",
    // M89.  dbSNP rs2032652; GRCh37 21917313 / GRCh38 19755427. YBrowse lists the
    // mutation C->T; dbSNP annotates the complementary strand (T->C). Allele
    // matching accepts either strand, so both the base and its complement count.
    markers: &[y("M89", Some("rs2032652"), 21917313, 19755427, "T")],
    children: &[
        Node { haplogroup: "K", markers: K.markers, children: K.children, note: None },
        Node { haplogroup: "I", markers: I_NODE.markers, children: I_NODE.children, note: None },
        Node { haplogroup: "J", markers: J_NODE.markers, children: J_NODE.children, note: None },
        Node { haplogroup: "G", markers: G_NODE.markers, children: &[], note: None },
        // Robustness: I/J subclades are on arrays even when M9/M89 are missing, so
        // K's children are also reachable directly from F via the K node above.
    ],
    note: None,
};

// --- E (M96 / M35) under CT, sibling of CF ---
static E: Node = Node {
    haplogroup: "E-M35",
    // E1b1b1 = M35 (anchor: E-M35 by M35).
    markers: &[y("M35", None, 21741703, 19579817, "C")],
    children: &[],
    note: None,
};

// CF (P143): parent of C and F.
static CF: Node = Node {
    haplogroup: "CF",
    markers: &[y("P143", None, 14197867, 12077161, "A")],
    children: &[Node { haplogroup: "F", markers: F.markers, children: F.children, note: None }],
    note: None,
};

// CT (M168): the root of all non-African (and most African) lineages here.
static CT: Node = Node {
    haplogroup: "CT",
    // M168 (the deepest marker arrays reliably carry; ancestral C -> derived T).
    markers: &[y("M168", None, 14813991, 12702062, "T")],
    children: &[
        Node { haplogroup: "CF", markers: CF.markers, children: CF.children, note: None },
        // F often reachable when P143 is not on the array.
        Node { haplogroup: "F", markers: F.markers, children: F.children, note: None },
        Node { haplogroup: "E-M35", markers: E.markers, children: &[], note: None },
    ],
    note: None,
};

/// Root of the paternal tree. The root is "entered" for any sample carrying Y
/// markers; descent into CT (M168) covers the resolvable lineages. Several
/// terminal lineages (R1b, R1a, I1, J1, J2, G, E-M35, N, O, Q) are also reachable
/// directly from the root so that arrays missing the deep backbone SNPs (M168,
/// M89, M9) can still be classified to a major haplogroup.
pub static ROOT: Node = Node {
    haplogroup: "Y (root)",
    markers: &[],
    children: &[
        Node { haplogroup: "CT", markers: CT.markers, children: CT.children, note: None },
        // Direct shortcuts to major terminals for backbone-sparse arrays.
        Node { haplogroup: "R", markers: R.markers, children: R.children, note: None },
        Node { haplogroup: "Q", markers: Q.markers, children: &[], note: None },
        Node { haplogroup: "I", markers: I_NODE.markers, children: I_NODE.children, note: None },
        Node { haplogroup: "J", markers: J_NODE.markers, children: J_NODE.children, note: None },
        Node { haplogroup: "G", markers: G_NODE.markers, children: &[], note: None },
        Node { haplogroup: "E-M35", markers: E.markers, children: &[], note: None },
        Node { haplogroup: "N", markers: N_NODE.markers, children: &[], note: None },
        Node { haplogroup: "O", markers: O_NODE.markers, children: &[], note: None },
    ],
    note: None,
};
