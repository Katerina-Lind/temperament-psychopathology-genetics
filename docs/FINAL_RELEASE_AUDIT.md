# Final release audit

Run `python3 scripts/final_release_audit.py` from any directory after extracting the add-on into the project repository. Read-only by default. To save a structured report outside the repository use `--json /tmp/temperament_release_audit.json`.

A `FAIL` is a release blocker until reviewed; a `WARN` is a manual-review item. Historical scripts in `analysis/source_project/` may legitimately preserve nonportable original paths, but runnable publication scripts and committed validation reports should not expose absolute account paths.

The audit tests committed source manifests against SHA256, validation report presence and statuses, tracked large files, external symlinks, local paths, potential tokens, Git state, metadata, and script availability. It does not execute data analysis or establish validity of the statistical assumptions.

Scientific gates still requiring explicit sign-off: provenance of fixed/effective sample sizes; precise mapping of GenomicSEM V entries; whether rerunning contrasts with the same software is sufficient for inference; robustness of age-specific findings to overlapping participants and correlated phenotypes; clear exploratory multiplicity family; final manuscript tables and plots; original GWAS licensing and redistributability; and availability/DOI of raw summary statistics and LD reference resources.

Do not commit unredacted machine-local reproduction logs or private file paths. For release, prefer an archived tag with a complete README, software lock/environment file, citation metadata, LICENSE, source references, and a short smoke-test command.
