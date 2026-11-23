# Local manuscript code audit

Generated: 2026-11-23 12:00 UTC

This is a **static, read-only scan**, not an execution test. It deliberately omits matched values, which could include secrets.

- Imported analysis scripts found: **11**
- Code files scanned (including repository utilities): **16**
- Findings: **18** (security: 0, portability: 18, other: 0)
- Expected analysis scripts missing: **0**

## Missing analysis scripts

None.

## Findings

| Category | File | Line | Rule |
|---|---|---:|---|
| portability | `analysis/genomicsem/01_joint_ldsc.R` | 3 | `home_project_path` |
| portability | `analysis/genomicsem/01_joint_ldsc.R` | 4 | `absolute_setwd` |
| portability | `analysis/genomicsem/02_joint_ldsc_7traits.R` | 3 | `home_project_path` |
| portability | `analysis/genomicsem/02_joint_ldsc_7traits.R` | 4 | `absolute_setwd` |
| portability | `analysis/genomicsem/03_test_rg_contrasts.R` | 1 | `home_project_path` |
| portability | `analysis/genomicsem/03_test_rg_contrasts.R` | 2 | `absolute_setwd` |
| portability | `analysis/genomicsem/04_joint_ldsc_developmental.R` | 3 | `home_project_path` |
| portability | `analysis/genomicsem/04_joint_ldsc_developmental.R` | 4 | `absolute_setwd` |
| portability | `analysis/genomicsem/05_test_developmental_contrasts.R` | 1 | `home_project_path` |
| portability | `analysis/genomicsem/05_test_developmental_contrasts.R` | 2 | `absolute_setwd` |
| portability | `analysis/genomicsem/07_age_specific_contrasts.R` | 3 | `home_project_path` |
| portability | `analysis/genomicsem/07_age_specific_contrasts.R` | 4 | `absolute_setwd` |
| portability | `analysis/gwas_qc/07_prepare_age_specific.R` | 3 | `home_project_path` |
| portability | `analysis/gwas_qc/07_prepare_age_specific.R` | 4 | `home_project_path` |
| portability | `analysis/gwas_qc/harmonise_shy_arisk.py` | 4 | `home_project_path` |
| portability | `analysis/ldsc/run_developmental_rg.sh` | 4 | `home_project_path` |
| portability | `scripts/check_local_analysis.py` | 6 | `home_project_path` |
| portability | `scripts/import_local_analysis.sh` | 3 | `home_project_path` |

## What to do before public release

1. Review each flagged line locally; do not paste actual credentials into reports or chats.
2. Parameterise source and output roots, then rerun each affected script on a clean working directory. Keep the original analysis revisions in Git history.
3. Import missing **essential** source preparation and LDSC scripts from `~/temperament_project` into `analysis/`; record the original Git commit or SHA256 where possible.
4. Verify ordering and parameterisation of GenomicSEM sampling covariance `V` before marking contrast P-values final.
5. Capture `R --version`, `sessionInfo()`, the LDSC Python environment, and full source GWAS metadata; pin dependency versions.
6. Check `git status --short`, `git diff --cached --stat`, and `git diff --cached --check` before pushing; never commit bulk GWAS or private data.

## Limitations

Regex-based detection cannot prove that a repository is free of secrets or that all research results reproduce. False positives and false negatives are possible.
