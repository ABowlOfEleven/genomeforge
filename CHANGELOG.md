# Changelog

All notable changes to GenomeForge are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and this project adheres to
semantic versioning.

## [Unreleased]

## [0.1.2] - 2026-06-11

### Added
- **Ancestry & lineage workspace:** a new tab that works entirely on your machine.
  - **Haplogroups:** predicted maternal (mtDNA) and paternal (Y-DNA) major
    haplogroups from your SNP calls, sourced from PhyloTree Build 17 and the ISOGG
    Y-tree, with honest coarse-resolution framing and a confidence badge.
  - **Ancestry composition:** a continental admixture estimate against five 1000
    Genomes super-populations (supervised EM over an ancestry-informative-marker
    panel), shown as a proportion bar, ranked percentages, and a reference plot.
    Framed as genetic similarity, not an identity test.
- **Pharmacogenomics workspace:** a CPIC-based drug-response report. SNP star-allele
  / diplotype calling for 15 pharmacogenes with metabolizer phenotypes and per-drug
  notes; genes not fully callable from array data (CYP2D6, UGT1A1, HLA-B, G6PD, ...)
  are flagged. Informational only, with a repeated not-medical-advice banner.
- **Settings: cache management and provenance:** view the local cache's size, age,
  and entry counts; clear it; and see the live data sources behind each feature, with
  a per-annotation provenance caption.

### Fixed
- The app now reports its real version (the binary previously reported `0.1.0`
  regardless of the release), so the About dialog, status bar, and self-update check
  are correct.

## [0.1.1] - 2026-06-10

### Added
- **GRCh37 ⇄ GRCh38 liftover:** convert a variant's coordinates between assemblies
  (via the Ensembl assembly-map API).
- **More CRISPR nucleases:** SpCas9-NG, SpRY, and Cas12a (Cpf1) alongside SpCas9,
  with correct PAM placement and protospacer lengths.
- **Base- and prime-editing outcome simulation:** CBE (C→T) / ABE (A→G) window
  outcomes for a guide, and a transparent prime-editing installed-edit preview.
- **Plasmid assembly simulation:** Golden-Gate / restriction-ligation and Gibson
  assembly from fragments, into a circular product.
- **Variant "Learn more" additions:** gnomAD variant link and gene-level PubMed /
  ClinVar shortcuts.
- **Quality of life:** drag-and-drop file open, a recent-files list, CSV export of the
  variant table, copy buttons, an About dialog, and a keyboard-shortcuts reference.
- **Update check:** compares the running version against the latest GitHub release.
- **Project infrastructure:** `CHANGELOG.md`, `CONTRIBUTING.md`, issue templates, and a
  `scripts/bump-version.sh` that sets the version everywhere it is declared.

## [0.1.0] - 2026-06-10

First multiplatform release (Windows, Linux/Flatpak, macOS).

### Added
- Genome browser & variant explorer: 23andMe/AncestryDNA, VCF, FASTA, GenBank import;
  ClinVar/gnomAD/dbSNP annotation via MyVariant.info; common-health-SNP quick-jumps;
  filterable variant table; per-variant deep dive with dbSNP/ClinVar/GWAS/MedlinePlus
  links and a live PubMed search.
- Plasmid designer: circular/linear maps, restriction mapping, ORFs/translation, primer
  Tm/GC, digest/cloning preview.
- CRISPR studio: SpCas9 guide design with a Doench 2014 on-target model, CFD off-target
  scoring, NHEJ/HDR edit simulation, AAV cargo planner.
- Phenotype & risk: PGS Catalog polygenic scores with coverage and an ancestry caveat.
- Beginner / Intermediate / Expert experience tiers with per-section tutorials.
- Packaging: Windows MSI + portable exe, Linux tarball + Flatpak, macOS universal dmg;
  CI + release workflows across all three platforms.

[Unreleased]: https://github.com/ABowlOfEleven/genomeforge/compare/v0.1.2...HEAD
[0.1.2]: https://github.com/ABowlOfEleven/genomeforge/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/ABowlOfEleven/genomeforge/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/ABowlOfEleven/genomeforge/releases/tag/v0.1.0
