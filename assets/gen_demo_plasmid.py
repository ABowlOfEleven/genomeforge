"""Generate a small synthetic circular demo plasmid in GenBank format.
Not a real construct — curated so GenomeForge has features, unique cutters,
and real ORFs to display."""
import random

random.seed(7)
SAFE_CODONS = ["GCT", "GAA", "AAA", "CTG", "GGT", "ATT", "ACC", "CGT",
               "GTT", "TCT", "CCG", "CAG", "GAT", "TAT", "AAC", "TGG"]

seq = ""
features = []  # (kind, start0, end0_excl, {qualifiers})


def add(s: str) -> int:
    global seq
    start = len(seq)
    seq += s
    return start


def feat(kind, start, end, **quals):
    features.append((kind, start, end, quals))


def orf(n_codons):
    return "ATG" + "".join(random.choice(SAFE_CODONS) for _ in range(n_codons)) + "TAA"


add("ATCATCATCATC")  # spacer (no restriction sites)

# pUC-style MCS: each enzyme exactly once -> unique cutters
mcs = "GAATTC GAGCTC GGTACC GGATCC TCTAGA GTCGAC CTGCAG AAGCTT".replace(" ", "")
s = add(mcs)
feat("misc_feature", s, s + len(mcs), label="MCS")

add("ATCATCATC")

prom = "TTGACA" + "ATAT" * 5 + "TATAAT" + "ATAT" * 3
s = add(prom)
feat("promoter", s, s + len(prom), label="Plac")

add("ATCATC")

cds1 = orf(96)
s = add(cds1)
feat("CDS", s, s + len(cds1), label="demoFP", gene="dfp")

add("ATCATCATC")

cds2 = orf(132)
s = add(cds2)
feat("CDS", s, s + len(cds2), label="AmpR", gene="bla")

add("ATCATCATC")

ori = "GCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCC"
s = add(ori)
feat("rep_origin", s, s + len(ori), label="ori")

add("ATCATCATCATC")

n = len(seq)
out = []
out.append(f"LOCUS       pGFDEMO        {n:>6} bp    DNA     circular SYN 09-JUN-2026")
out.append("DEFINITION  Synthetic demo plasmid for GenomeForge (not a real construct).")
out.append("ACCESSION   pGFDEMO")
out.append("VERSION     pGFDEMO.1")
out.append("KEYWORDS    .")
out.append("SOURCE      synthetic DNA construct")
out.append("  ORGANISM  synthetic DNA construct")
out.append("FEATURES             Location/Qualifiers")
out.append(f"     source          1..{n}")
out.append('                     /organism="synthetic DNA construct"')
for kind, s0, e0, quals in features:
    out.append(f"     {kind:<15} {s0 + 1}..{e0}")
    for k, v in quals.items():
        out.append(f'                     /{k}="{v}"')
out.append("ORIGIN")
low = seq.lower()
for i in range(0, len(low), 60):
    chunk = low[i:i + 60]
    groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
    out.append(f"{i + 1:>9} {groups}")
out.append("//")

with open("E:/Repos/genomeforge/assets/demo_plasmid.gb", "w") as fh:
    fh.write("\n".join(out) + "\n")
print(f"wrote {n} bp, {len(features)} features")
