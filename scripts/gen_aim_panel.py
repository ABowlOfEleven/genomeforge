#!/usr/bin/env python3
"""Generate the ancestry-informative-marker (AIM) panel for gx-ancestry.

Fetches 1000 Genomes phase 3 super-population allele frequencies (Ensembl REST,
pops=1) and GRCh37/GRCh38 positions for a curated set of continental AIM rsIDs,
then writes a self-contained Rust data module `crates/gx-ancestry/src/panel.rs`.
Also prints EUR-typical demo genotype rows for the bundled sample genome.

This is a developer/build-time script. The app does not run it; the generated
panel.rs is committed. Frequencies are 1000G phase 3 (Ensembl), build-independent.
"""
import json
import sys
import time
import urllib.request

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

# Curated continental ancestry-informative SNPs: pigmentation/selection loci and
# published AISNP-panel markers (Kidd et al. 55-AISNP and similar). Any rsID that
# does not return clean biallelic 5-super-population data with both positions is
# dropped automatically below.
CANDIDATES = [
    "rs2814778", "rs1426654", "rs16891982", "rs3827760", "rs12913832",
    "rs1042602", "rs1800407", "rs1800414", "rs2402130", "rs12498138",
    "rs2238151", "rs7997709", "rs1834640", "rs6058017", "rs2378249",
    "rs260690", "rs3811801", "rs1229984", "rs671", "rs3814134",
    "rs4833103", "rs7421394", "rs2070586", "rs459920", "rs6437783",
    "rs9522149", "rs1572018", "rs1369290", "rs735480", "rs10108270",
    "rs4918664", "rs174570", "rs1079597", "rs2196051", "rs1871534",
    "rs10497191", "rs772262", "rs7657799", "rs7745461", "rs192655",
    "rs3823159", "rs6990312", "rs2166624", "rs12439433", "rs798443",
    "rs1876482", "rs3737576", "rs7554936", "rs2125345", "rs200354",
    "rs4411548", "rs4471745", "rs6754311", "rs10007810", "rs4954681",
    "rs1407434", "rs7226659",
]


def post(host, path, ids):
    body = json.dumps({"ids": ids}).encode()
    req = urllib.request.Request(
        f"https://{host}{path}",
        data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            sys.stderr.write(f"  retry {attempt + 1} ({e})\n")
            time.sleep(2)
    return {}


def batched(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


def main():
    freqs = {}
    pos38 = {}
    pos37 = {}
    for chunk in batched(CANDIDATES, 20):
        sys.stderr.write(f"fetching freqs+GRCh38 for {len(chunk)} ids...\n")
        d = post("rest.ensembl.org", "/variation/homo_sapiens?pops=1", chunk)
        freqs.update(d)
        for rs, info in d.items():
            for m in info.get("mappings", []):
                if m.get("assembly_name") == "GRCh38" and m.get("coord_system") == "chromosome":
                    pos38[rs] = (m["seq_region_name"], int(m["start"]))
        sys.stderr.write(f"fetching GRCh37 positions for {len(chunk)} ids...\n")
        d37 = post("grch37.rest.ensembl.org", "/variation/homo_sapiens", chunk)
        for rs, info in d37.items():
            for m in info.get("mappings", []):
                if m.get("assembly_name") == "GRCh37" and m.get("coord_system") == "chromosome":
                    pos37[rs] = (m["seq_region_name"], int(m["start"]))

    panel = []
    for rs in CANDIDATES:
        info = freqs.get(rs)
        if not info or rs not in pos37 or rs not in pos38:
            continue
        contig, p37 = pos37[rs]
        _, p38 = pos38[rs]
        if contig not in [str(i) for i in range(1, 23)]:
            continue  # autosomal only
        # Collect per-super-population per-allele frequencies.
        per = {sp: {} for sp in SUPERPOPS}
        for p in info.get("populations", []):
            name = p.get("population", "")
            if not name.startswith("1000GENOMES:phase_3:"):
                continue
            sp = name.split(":")[-1]
            if sp in per and len(p.get("allele", "")) == 1:
                per[sp][p["allele"]] = float(p["frequency"])
        alleles = set()
        for sp in SUPERPOPS:
            alleles |= set(per[sp].keys())
        alleles = sorted(a for a in alleles if a in "ACGT")
        if len(alleles) != 2:
            continue
        eff, ref = alleles[0], alleles[1]
        freq = []
        ok = True
        for sp in SUPERPOPS:
            e = per[sp]
            if eff in e:
                f = e[eff]
            elif ref in e:
                f = 1.0 - e[ref]
            else:
                ok = False
                break
            freq.append(round(f, 4))
        if not ok:
            continue
        panel.append((rs, contig, p37, p38, eff, ref, freq))

    panel.sort(key=lambda x: (int(x[1]), x[2]))
    sys.stderr.write(f"\nKept {len(panel)} / {len(CANDIDATES)} AIMs.\n")

    lines = [
        "//! Ancestry-informative-marker (AIM) panel: 1000 Genomes phase 3",
        "//! super-population effect-allele frequencies (Ensembl REST) plus GRCh37 /",
        "//! GRCh38 positions. Generated by `scripts/gen_aim_panel.py`; do not edit by",
        "//! hand. Frequencies are build-independent; positions are per build.",
        "",
        "/// One ancestry-informative marker.",
        "pub struct Aim {",
        "    pub rsid: &'static str,",
        "    pub contig: &'static str,",
        "    pub pos37: u32,",
        "    pub pos38: u32,",
        "    /// Effect allele (the one `freq` is the frequency of), as an ASCII base.",
        "    pub effect: u8,",
        "    /// Effect-allele frequency per super-population, in SUPERPOPS order.",
        "    pub freq: [f64; 5],",
        "}",
        "",
        '/// Super-population labels, the index order of `Aim::freq` and estimates.',
        'pub const SUPERPOPS: [&str; 5] = ["AFR", "AMR", "EAS", "EUR", "SAS"];',
        "",
        "/// Human-readable super-population names (same order as SUPERPOPS).",
        'pub const SUPERPOP_NAMES: [&str; 5] = [',
        '    "African", "Admixed American", "East Asian", "European", "South Asian",',
        "];",
        "",
        "// Frequencies are data; some coincidentally approximate math constants.",
        "#[allow(clippy::approx_constant)]",
        "#[rustfmt::skip]",
        "pub const PANEL: &[Aim] = &[",
    ]
    for rs, contig, p37, p38, eff, ref, freq in panel:
        fr = ", ".join(f"{x}" for x in freq)
        lines.append(
            f'    Aim {{ rsid: "{rs}", contig: "{contig}", pos37: {p37}, '
            f"pos38: {p38}, effect: b'{eff}', freq: [{fr}] }},"
        )
    lines.append("];")
    lines.append("")
    out = "\n".join(lines)
    with open("crates/gx-ancestry/src/panel.rs", "w", newline="\n", encoding="utf-8") as f:
        f.write(out)
    sys.stderr.write("wrote crates/gx-ancestry/src/panel.rs\n")

    # EUR-typical demo rows: dosage = round(2 * EUR effect-allele freq).
    print("# Ancestry-informative markers (GRCh37) for the Ancestry composition view.")
    eur = SUPERPOPS.index("EUR")
    for rs, contig, p37, p38, eff, ref, freq in panel:
        dose = round(2 * freq[eur])
        geno = {2: eff + eff, 1: eff + ref, 0: ref + ref}[dose]
        print(f"{rs}\t{contig}\t{p37}\t{geno}")


if __name__ == "__main__":
    main()
