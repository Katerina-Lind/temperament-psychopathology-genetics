# Early-life temperament and psychiatric genetic architecture

**Reproducibility repository — manuscript in preparation (working release, 23 November 2026).**

This repository documents European-ancestry GWAS analyses examining the genetic overlap of infant/toddler EAS temperament dimensions (emotionality, activity and shyness) with childhood and adult psychopathology phenotypes. Analyses include univariate and cross-trait LD score regression (LDSC), GenomicSEM joint covariance estimation, covariance-aware contrasts between genetic correlations, and second- versus third-year temperament sensitivity analyses.

> **Status:** The published `results/` files are **provisional manuscript snapshots**, assembled from investigator-supplied numerical outputs. The original individual analysis scripts and GWAS files reside on the investigator's Mac and **have not been copied into this archive**. Before publication, import the exact source code and regenerate the tables using the original RDS/CSV files; see [Reproduction guide](docs/reproduction.md). Do not describe this archive as full end-to-end reproducibility until those steps are complete.

## Scientific questions

1. Do emotionality, activity and shyness show different genetic correlations across internalising, ADHD, substance use and risk-taking phenotypes?
2. Do genetic correlations with the *same outcome* differ for temperament GWAS measured in the second versus third postnatal years?

Cross-phenotype comparisons are **not tests of within-person developmental change**; the age-specific GWAS may include shared participants and differ in cohort composition and measurement.

## Repository layout

```text
analysis/            Author's original R/Bash/Python analysis scripts (to be imported)
config/              Analysis manifest and source provenance
scripts/             Audit, import and snapshot plotting utilities
results/tables/      Provisional derived results tables (CSV)
results/figures/     Figure previews in vector PDF
references/          (not bundled) LDSC European LD-score reference and HM3 list
raw GWAS/            (not bundled) must be obtained from original sources
```

## Quick start

```bash
python3 scripts/check_results.py
python3 scripts/check_local_analysis.py --project ~/temperament_psychopathology
# On investigator's Mac, after inspecting the list of files:
bash scripts/import_local_analysis.sh ~/temperament_psychopathology
```

To rebuild plots from the saved snapshot (NumPy and Matplotlib required):

```bash
python3 -m pip install -r requirements-figures.txt
python3 scripts/plot_snapshot.py --output .
```

To rebuild figures using exact-precision age-specific contrasts, additionally supply `--age-csv /path/to/age_specific_contrasts.csv`. **Note:** `plot_snapshot.py` contains manually transcribed cross-phenotype inputs. Treat its products as figure previews until the results are regenerated from raw analysis artifacts.

## Included results

- Seven-trait joint LDSC genetic correlation estimates (21 pairwise coefficients).
- Four exploratory covariance-aware cross-phenotype contrasts.
- Six age-specific SNP heritability estimates.
- Eighteen second-year versus third-year temperament contrasts with six psychiatric outcomes each.
- Four figure previews derived from these tables.

The two FDR families (4 and 18 tests) were corrected separately, **after exploratory phenotype selection**. Do not interpret the adjusted q values as multiplicity correction for all previously explored outcomes or models. `docs/statistical_caveats.md` lists pending checks.

## Data access and privacy

GWAS summary statistics (~0.5 GiB per temperament file) are **not redistributed**. See `config/data_sources.tsv` for file identifiers and `docs/data_access.md` for instructions. LD-score reference files also are excluded. Local absolute paths, cloud credentials, personal identifiers, and bulk data must not be committed.

## Citation

A citation record with placeholder authors/title is provided in `CITATION.cff`. Replace these details and add the manuscript DOI **before** public release. Until then, reference this as code accompanying a manuscript in preparation and cite underlying source GWAS papers separately.

## License

No open-source licence is applied yet. Authors should choose an appropriate licence with all contributors' agreement before publishing.
