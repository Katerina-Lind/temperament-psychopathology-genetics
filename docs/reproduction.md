# Reproduction workflow

## 1. Exact analysis implementation

On the author's Mac, the working analysis lives in `~/temperament_psychopathology` and the original data/preprocessing project is `~/temperament_project`. Those folders are **not** included here, and are not published automatically.

Run the import script from this repository directory to copy the **actual scripts that generated the results**; do not hand-reconstruct analysis code or substitute new thresholds without recording the change.

```bash
bash scripts/import_local_analysis.sh ~/temperament_psychopathology
```

The import script is non-destructive; it copies selected code into `analysis/` and does not copy GWAS inputs, RDS objects, logs, credentials, or machine-specific environment files. Audit the code for personal paths, verify the selected scripts correspond to final results, and update `config/data_sources.tsv` before committing. Original older QC scripts in `~/temperament_project` should be imported separately where used.

## 2. Original pipeline outline

1. Acquire external European GWAS summary statistics; verify hashes and record provenance.
2. Munge GWAS: aligned effect allele/Z, HapMap3 matching, duplicate rsID exclusions, allele and MAF checks, SNP-specific N, reference genome checks. Retain QC logs and counts.
3. Run univariate LDSC h2 with European-ancestry LD scores.
4. Run joint LDSC on EMO/ACT/SHY/CADHD/AINT/ASUD/ARISK, with consistent harmonised ARISK input; save `developmental_7traits.rds`.
5. Compute covariance-aware delta method contrasts using S/V, validate ordering.
6. Run six age-specific GWAS through same QC and univariate LDSC.
7. Run 18 trivariate joint LDSC fits (2Y/3Y/outcome), retain S/V and check eigenvalues, compute paired rg contrasts.
8. Generate tables/figures from the saved CSV/RDS artifacts, not manually transcribed console output.

## 3. Expected local outputs

```text
results/genomicsem/developmental_7traits.rds
results/genomicsem/developmental_rg_contrasts.csv
results/genomicsem/age_specific/age_specific_contrasts.csv
results/ldsc/age_specific_h2/*[23]Y.log
```

The files under `results/tables/` are **snapshot tables**, not files proved to have been produced by importing a raw GWAS pipeline into this archive.

## 4. Refresh age-specific tables/figures

```bash
python3 scripts/plot_snapshot.py --output . \
  --age-csv ~/temperament_psychopathology/results/genomicsem/age_specific/age_specific_contrasts.csv
python3 scripts/check_results.py
```

This imports full-precision age-specific contrasts; other visualisations remain based on embedded snapshot estimates until replaced with full-source exports. Never compare changed manuscript numbers without updating the code, data-source versions, and a release tag.
