# Analysis scripts

This directory contains the analysis scripts used in the
temperament–psychopathology genetic epidemiology project.

The scripts are already included in this repository.
No additional import step is required.

## GenomicSEM analyses

- `genomicsem/01_joint_ldsc.R` — earlier joint LDSC analysis.
- `genomicsem/02_joint_ldsc_7traits.R` — earlier seven-trait analysis.
- `genomicsem/03_test_rg_contrasts.R` — earlier genetic correlation contrasts.
- `genomicsem/04_joint_ldsc_developmental.R` — principal seven-trait joint LDSC analysis.
- `genomicsem/05_test_developmental_contrasts.R` — covariance-aware developmental contrasts.
- `genomicsem/06_validate_rg_se.R` — genetic correlation SE validation.
- `genomicsem/07_age_specific_contrasts.R` — 18 age-specific joint LDSC analyses and contrasts.

## GWAS quality control and preparation

- `gwas_qc/07_prepare_age_specific.R` — age-specific temperament GWAS preparation.
- `gwas_qc/clean_arisk.py` — risk-taking GWAS cleaning.
- `gwas_qc/finalise_arisk.py` — final risk-taking GWAS filtering.
- `gwas_qc/harmonise_shy_arisk.py` — shyness/risk-taking harmonisation.

## Historical scripts

The `source_project/` directory contains archived historical
analysis scripts. Their provenance and integrity are documented
in the repository's source-script manifest.

## Reproducibility

The repository includes reproduction scripts, environment
specifications, and validation reports.

All 18 GWAS inputs, the principal seven-trait LDSC analysis,
and all 18 age-specific LDSC analyses have been reproduced
against the original numerical results.

The four principal developmental contrasts have also been
checked using independent analytical gradients.

See `../results/reproducibility/` for validation reports.

## Data availability

GWAS summary statistics and LD-score reference files are
not included in this repository. These must be obtained
separately from their respective sources.

## Manuscript outputs

The age-specific analysis produces:

`results/genomicsem/age_specific/age_specific_contrasts.csv`

A separate `08_age_specific_manuscript_outputs.R` script was
mentioned in an earlier repository template but was never
present in the original analysis projects. It is not required
for the validated statistical analyses.
