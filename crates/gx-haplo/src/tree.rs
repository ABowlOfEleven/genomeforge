//! Shared haplogroup-tree model and the tree-walking classifier used by both the
//! mitochondrial (`mt`) and Y-chromosome (`ydna`) modules.
//!
//! A lineage is a rooted tree of [`Node`]s. Each node names a haplogroup and
//! carries the [`Marker`]s that, when found in their derived (mutant) state,
//! support descending into that node. Classification walks from the root,
//! descending into a child whenever that child has positive marker support, and
//! reports the deepest supported node. See [`classify`].

use gx_core::VariantStore;

/// A single haplogroup-defining variant: the sample is "derived" at this marker
/// when its genotype carries `derived` at the reference position.
///
/// Positions are 1-based reference coordinates (the file's natural numbering);
/// the matching code converts to the store's 0-based, half-open convention.
#[derive(Clone, Copy, Debug)]
pub struct Marker {
    /// Short marker name, e.g. `"M269"` (Y) or a position label like `"11467"` (MT).
    pub name: &'static str,
    /// dbSNP rsID when an authoritative one is known, else `None`. Used only as a
    /// fallback lookup; primary matching is by position + derived allele.
    pub rsid: Option<&'static str>,
    /// 1-based position on the relevant contig for GRCh37/hg19. For MT this is the
    /// rCRS position (also used for GRCh38).
    pub pos_grch37: u64,
    /// 1-based position for GRCh38/hg38. For MT, equal to `pos_grch37`.
    pub pos_grch38: u64,
    /// The derived (mutant) base the sample must carry to support this marker,
    /// uppercase single base. Indel markers are omitted from the tables because
    /// consumer arrays do not call them reliably.
    pub derived: &'static str,
}

impl Marker {
    /// 1-based position for the given assembly.
    fn pos_for(&self, asm: gx_core::Assembly) -> u64 {
        match asm {
            gx_core::Assembly::Grch37 => self.pos_grch37,
            gx_core::Assembly::Grch38 => self.pos_grch38,
        }
    }
}

/// A haplogroup node in a lineage tree.
pub struct Node {
    /// Reported haplogroup label, e.g. `"U5"`, `"R1b"`, `"H / HV (R0 cluster)"`.
    pub haplogroup: &'static str,
    /// Markers whose derived state supports descending into this node. A node may
    /// have several; any one matching is treated as support (consumer arrays have
    /// gaps), but more matches raise confidence.
    pub markers: &'static [Marker],
    /// Child nodes, tried in order; the first child with positive support is taken.
    pub children: &'static [Node],
    /// Optional note appended when this node is the terminal assignment (used for
    /// the conservative H / HV cluster caveat).
    pub note: Option<&'static str>,
}

/// Whether a marker was found in the data and, if so, whether it was derived.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MarkerState {
    /// Marker position/rsID was genotyped and carries the derived allele.
    Derived,
    /// Marker was genotyped but carries the ancestral (or some other) allele.
    Ancestral,
    /// Marker was not present in the data (not on the array / no-call).
    Missing,
}

/// Test a single marker against the store: is the derived allele present?
pub fn marker_state(store: &VariantStore, contig: &str, m: &Marker) -> MarkerState {
    let asm = store.assembly();
    let pos1 = m.pos_for(asm);

    // Primary match: by position. variants_in uses 0-based half-open [start,end);
    // a 1-based position P maps to the half-open window [P-1, P).
    let start = pos1.saturating_sub(1);
    let mut genotyped = false;
    for v in store.variants_in(contig, start, pos1) {
        if let Some(gt) = &v.genotype {
            if gt.is_no_call() {
                continue;
            }
            genotyped = true;
            if allele_matches(gt, m.derived) {
                return MarkerState::Derived;
            }
        }
    }

    // Fallback: rsID lookup (covers files whose positions differ, e.g. unusual
    // contig naming) only when we have an authoritative rsID for this marker.
    if let Some(rsid) = m.rsid {
        if let Some(v) = store.find_by_rsid(rsid) {
            if let Some(gt) = &v.genotype {
                if !gt.is_no_call() {
                    genotyped = true;
                    if allele_matches(gt, m.derived) {
                        return MarkerState::Derived;
                    }
                }
            }
        }
    }

    if genotyped {
        MarkerState::Ancestral
    } else {
        MarkerState::Missing
    }
}

/// Does a genotype carry the derived base, on either reported strand?
///
/// Consumer files may report a locus on the opposite strand, so we accept the
/// derived base or its complement. (Both lineages here are effectively haploid,
/// so any allele equal to the target counts.)
fn allele_matches(gt: &gx_core::Genotype, derived: &str) -> bool {
    let want = derived.to_ascii_uppercase();
    let want_comp = complement(&want);
    gt.alleles.iter().any(|a| {
        let a = a.to_ascii_uppercase();
        a == want || want_comp.as_deref() == Some(a.as_str())
    })
}

/// Single-base complement, or `None` for multi-base / non-ACGT alleles.
fn complement(allele: &str) -> Option<String> {
    if allele.len() != 1 {
        return None;
    }
    let c = match allele.as_bytes()[0].to_ascii_uppercase() {
        b'A' => 'T',
        b'T' => 'A',
        b'G' => 'C',
        b'C' => 'G',
        _ => return None,
    };
    Some(c.to_string())
}

/// Outcome of walking a lineage tree.
pub struct Walk {
    /// Label of the deepest supported node.
    pub haplogroup: String,
    /// Marker names matched (derived) along the assigned path, root to leaf.
    pub supporting: Vec<String>,
    /// Total markers across the whole lineage tree that were genotyped.
    pub tested_total: usize,
    /// Note from the terminal node, if any.
    pub note: Option<&'static str>,
}

/// Walk `root` against the store, descending into supported children.
///
/// Returns `None` only when nothing in the entire tree was genotyped (so we
/// cannot say anything at all). Otherwise returns the deepest supported node;
/// if even the root lacks support but some markers were tested, the root is
/// returned with empty `supporting` (the caller decides confidence).
pub fn classify(store: &VariantStore, contig: &str, root: &'static Node) -> Option<Walk> {
    let tested_total = count_tested(store, contig, root);
    if tested_total == 0 {
        return None;
    }

    let mut supporting: Vec<String> = Vec::new();
    // Report the deepest *directly supported* node. We may pass through
    // intermediate nodes that have no direct support (their own SNP is absent
    // from the data) to reach a supported terminal, but those pass-through nodes
    // must not be reported as the assignment, only the supported ones.
    let mut label = root.haplogroup;
    let mut note = root.note;

    // Tally the root's own markers (the root is always "entered").
    tally(store, contig, root, &mut supporting);

    let mut current = root;
    loop {
        let mut advanced = false;
        for child in current.children {
            if node_supported(store, contig, child) {
                tally(store, contig, child, &mut supporting);
                if directly_supported(store, contig, child) {
                    // A real clade assignment: adopt its label and note.
                    label = child.haplogroup;
                    note = child.note;
                }
                current = child;
                advanced = true;
                break;
            }
        }
        if !advanced {
            break;
        }
    }

    Some(Walk {
        haplogroup: label.to_string(),
        supporting,
        tested_total,
        note,
    })
}

/// Does this node have its own direct support (a derived marker)?
fn directly_supported(store: &VariantStore, contig: &str, node: &Node) -> bool {
    node.markers
        .iter()
        .any(|m| marker_state(store, contig, m) == MarkerState::Derived)
}

/// Is there a reason to descend into this node, anywhere in its subtree?
///
/// True when the node itself is directly supported, or when any descendant is.
/// This lets the walk pass through nodes whose own defining SNP is absent from
/// the data (structural pass-throughs like R1, or backbone SNPs not on the array,
/// or upstream clades whose mutations are reference state in rCRS) so that a deep
/// terminal marker (e.g. the U trio, or M269) still pulls the assignment down to
/// its clade. The deepest *directly supported* node is what gets reported.
fn node_supported(store: &VariantStore, contig: &str, node: &Node) -> bool {
    directly_supported(store, contig, node)
        || node
            .children
            .iter()
            .any(|c| node_supported(store, contig, c))
}

/// Record derived-marker names for a single node along the assigned path.
fn tally(store: &VariantStore, contig: &str, node: &Node, supporting: &mut Vec<String>) {
    for m in node.markers {
        if marker_state(store, contig, m) == MarkerState::Derived {
            supporting.push(m.name.to_string());
        }
    }
}

/// Count, across the whole subtree, how many distinct markers were genotyped
/// (derived or ancestral) so we can tell "no data" from "data, no support".
fn count_tested(store: &VariantStore, contig: &str, node: &Node) -> usize {
    let mut n = 0;
    for m in node.markers {
        if marker_state(store, contig, m) != MarkerState::Missing {
            n += 1;
        }
    }
    for c in node.children {
        n += count_tested(store, contig, c);
    }
    n
}
