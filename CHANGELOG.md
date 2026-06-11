# Changelog

All notable changes to GenomeForge are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and this project adheres to
semantic versioning.

## [Unreleased]

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

[Unreleased]: https://github.com/ABowlOfEleven/genomeforge/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/ABowlOfEleven/genomeforge/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/ABowlOfEleven/genomeforge/releases/tag/v0.1.0
