# Cross-age Hollowell GWAS reproduction check

This is an **isolated reproduction test**, not yet a complete GWAS-to-manuscript pipeline. It reuses the original historical scripts from `~/temperament_project/scripts` without editing them and does not change the validated GWAS files.

## Preflight

From the publication repository root, run:

```bash
python3 scripts/reproduce_cross_age_temperament.py
```

The preflight verifies original GWAS files, the HapMap3 reference, historical scripts, `python3`/`Rscript`, and uncompressed SHA256 hashes of the four validated final temperament files across the original and current analysis projects.

## Full reproduction (optional after preflight)

```bash
python3 scripts/reproduce_cross_age_temperament.py --run
```

Outputs are written to `~/temperament_cross_age_reproduction`, which must be empty/nonexistent. To rerun, provide `--workdir ~/temperament_cross_age_reproduction_2` rather than deleting or overwriting previous outputs.

The script stages **symlinks to original input data and reference** and local copies of original scripts. It then runs:

1. `13_prepare_emotionality.py`
2. `14_munge_emotionality.R`
3. `16_prepare_remaining_temperament.py`
4. `17_munge_remaining_temperament.R`

Logs are preserved under `logs/`, the generated files under `data/munged/`, and a JSON comparison under `reproduction_report.json`.

## Interpretation and limits

- Identical decompressed hashes establish exact textual identity of the regenerated file and the archived final analytical input.
- Differing hashes do not necessarily imply a scientific disagreement: GenomicSEM version, sorting, numeric formatting, reference changes and package behavior may matter. Investigate record-level differences before drawing conclusions.
- The original scripts read `data/raw/temperament/` and `references/ldsc/w_hm3.snplist`. The runner assumes the originals exist at those historical paths.
- This tests **cross-age temperament only**, not age-specific GWAS, psychiatric GWAS, LDSC or GenomicSEM contrasts.
- The repository's other historical scripts remain archived, not certified portable.
- This workflow deliberately does not commit generated GWAS files to Git.
