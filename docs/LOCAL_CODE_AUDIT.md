# Local manuscript code audit

Generated: 2026-11-23 12:14 UTC

This is a **static, read-only scan**, not an execution test. It deliberately omits matched values, which could include secrets.

- Imported analysis scripts found: **32**
- Code files scanned (including repository utilities): **38**
- Findings: **31** (security: 0, portability: 31, other: 0)
- Expected analysis scripts missing: **0**

## Missing analysis scripts

None.

## Findings

| Category | File | Line | Rule |
|---|---|---:|---|
| portability | `analysis/genomicsem/01_joint_ldsc.R` | 8 | `absolute_setwd` |
| portability | `analysis/genomicsem/02_joint_ldsc_7traits.R` | 8 | `absolute_setwd` |
| portability | `analysis/genomicsem/03_test_rg_contrasts.R` | 6 | `absolute_setwd` |
| portability | `analysis/genomicsem/04_joint_ldsc_developmental.R` | 8 | `absolute_setwd` |
| portability | `analysis/genomicsem/05_test_developmental_contrasts.R` | 6 | `absolute_setwd` |
| portability | `analysis/genomicsem/07_age_specific_contrasts.R` | 8 | `absolute_setwd` |
| portability | `analysis/source_project/core_preparation/03_munge.R` | 13 | `absolute_setwd` |
| portability | `analysis/source_project/core_preparation/03_munge.R` | 14 | `absolute_setwd` |
| portability | `analysis/source_project/core_preparation/30_munge_ARISK.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/core_preparation/33_prepare_remaining_gwas.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/core_preparation/40_make_strict_allele_inputs.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/ldsc_reference/37_ldsc_childhood.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/ldsc_reference/38_ldsc_temperament_childhood.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/ldsc_reference/41_ldsc_full_ten_traits.R` | 4 | `home_project_path` |
| portability | `analysis/source_project/ldsc_reference/42_test_developmental_contrasts.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/ldsc_reference/49_ldsc_full_individual_highN.R` | 4 | `home_project_path` |
| portability | `analysis/source_project/ldsc_reference/50_validate_rg_gradients.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/ldsc_reference/51_finalise_rg_results.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/quality_control/26_audit_psychopathology_gwas.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/quality_control/27_validate_gwas_statistics.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/quality_control/34_check_ldsc_coverage.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/quality_control/39_audit_cross_gwas_alleles.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/quality_control/47_audit_sample_size_ldscore.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/quality_control/52_gwas_provenance_audit.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/quality_control/53_annotate_gwas_provenance.R` | 3 | `home_project_path` |
| portability | `analysis/source_project/quality_control/54_inventory_gwas_build_metadata.py` | 4 | `home_project_path` |
| portability | `analysis/source_project/quality_control/55_inspect_original_gwas_metadata.py` | 5 | `home_project_path` |
| portability | `analysis/source_project/quality_control/56_validate_gwas_effect_directions.R` | 3 | `home_project_path` |
| portability | `scripts/check_local_analysis.py` | 6 | `home_project_path` |
| portability | `scripts/import_local_analysis.sh` | 3 | `home_project_path` |
| portability | `scripts/import_source_preprocessing.py` | 4 | `home_project_path` |

## What to do before public release

1. Review each flagged line locally; do not paste actual credentials into reports or chats.
2. Parameterise source and output roots, then rerun each affected script on a clean working directory. Keep the original analysis revisions in Git history.
3. Import missing **essential** source preparation and LDSC scripts from `~/temperament_project` into `analysis/`; record the original Git commit or SHA256 where possible.
4. Verify ordering and parameterisation of GenomicSEM sampling covariance `V` before marking contrast P-values final.
5. Capture `R --version`, `sessionInfo()`, the LDSC Python environment, and full source GWAS metadata; pin dependency versions.
6. Check `git status --short`, `git diff --cached --stat`, and `git diff --cached --check` before pushing; never commit bulk GWAS or private data.

## Limitations

Regex-based detection cannot prove that a repository is free of secrets or that all research results reproduce. False positives and false negatives are possible.
