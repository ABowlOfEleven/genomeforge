# Contributing to GenomeForge

Thanks for your interest in GenomeForge. This is a local-first, native Rust application;
contributions of all sizes are welcome.

## Getting started

```sh
git clone https://github.com/ABowlOfEleven/genomeforge
cd genomeforge
cargo run -p gx-app                                       # launch
cargo run -p gx-app -- assets/sample_genome_23andme.txt   # open a demo genome
```

The toolchain is pinned to the version CI uses (currently Rust 1.96.0). `rustup` will pick
it up automatically.

## Before you open a PR

Run the same checks CI runs (all three must be clean):

```sh
cargo build --workspace --locked
cargo test --workspace --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

The bar is **zero build warnings and zero clippy warnings**. Formatting is not enforced by
CI, but matching the surrounding style is appreciated.

## Project layout

A Cargo workspace of focused crates:

| Crate | Responsibility |
|-------|----------------|
| `gx-core` | Domain model: assemblies, ranges, variants, features |
| `gx-io` | Importers: raw DNA, VCF, FASTA / GenBank |
| `gx-annotate` | Online + cached annotation (MyVariant, Ensembl, PGS, PubMed) |
| `gx-plasmid` | Restriction mapping, ORFs, primers, digest / cloning / assembly |
| `gx-crispr` | Guide design, scoring, edit / base / prime editing, AAV |
| `gx-pgs` | PGS Catalog parsing and scoring |
| `gx-app` | The eframe (egui) desktop application |

All network and disk I/O runs on a background worker thread (`gx-app/src/worker.rs`); never
block the paint loop.

## Scientific honesty

Scoring methods (Doench on-target, CFD off-target, base-edit windows, PGS percentiles) are
**transparent, documented decision-support estimates**, not a clinical pipeline. New
heuristics must be documented as such, with their assumptions stated in code.

## Releasing (maintainers)

```sh
scripts/bump-version.sh 0.1.1   # sets the version everywhere
# update CHANGELOG.md, commit, then:
git tag v0.1.1 && git push origin v0.1.1
```

The release workflow verifies on all three platforms and publishes the artifacts.
