# Local manuscript code audit

Generated: 2026-11-23 12:03 UTC

This is a **static, read-only scan**, not an execution test. It deliberately omits matched values, which could include secrets.

- Imported analysis scripts found: **11**
- Code files scanned (including repository utilities): **17**
- Findings: **15** (security: 0, portability: 15, other: 0)
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
| portability | `scripts/check_local_analysis.py` | 6 | `home_project_path` |
| portability | `scripts/import_local_analysis.sh` | 3 | `home_project_path` |
| portability | `scripts/patch_project_paths.py` | 16 | `home_project_path` |
| portability | `scripts/patch_project_paths.py` | 17 | `absolute_setwd` |
| portability | `scripts/patch_project_paths.py` | 24 | `absolute_setwd` |
| portability | `scripts/patch_project_paths.py` | 26 | `home_project_path` |
| portability | `scripts/patch_project_paths.py` | 27 | `home_project_path` |
| portability | `scripts/patch_project_paths.py` | 39 | `home_project_path` |
| portability | `scripts/patch_project_paths.py` | 80 | `home_project_path` |

## What to do before public release

1. Review each flagged line locally; do not paste actual credentials into reports or chats.
2. Parameterise source and output roots, then rerun each affected script on a clean working directory. Keep the original analysis revisions in Git history.
3. Import missing **essential** source preparation and LDSC scripts from `~/temperament_project` into `analysis/`; record the original Git commit or SHA256 where possible.
4. Verify ordering and parameterisation of GenomicSEM sampling covariance `V` before marking contrast P-values final.
5. Capture `R --version`, `sessionInfo()`, the LDSC Python environment, and full source GWAS metadata; pin dependency versions.
6. Check `git status --short`, `git diff --cached --stat`, and `git diff --cached --check` before pushing; never commit bulk GWAS or private data.

## Limitations

Regex-based detection cannot prove that a repository is free of secrets or that all research results reproduce. False positives and false negatives are possible.
