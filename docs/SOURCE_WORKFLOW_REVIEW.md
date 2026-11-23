# Source GWAS preprocessing: reproduction checklist

**Status: code inventory only; not an end-to-end reproduction.** The scripts imported from the older project are historical analysis records. They must be reviewed and adapted in a separate commit before being advertised as executable from a clean checkout.

## Intended dependency chain (verify against code)

1. Retrieve original GWAS summary-statistic files, recording provider, version, file hashes and terms of use. **Do not commit GWAS files.**
2. Obtain and verify European LD-score reference (`eur_w_ld_chr`) and HM3 list (`w_hm3.snplist`).
3. Review original GWAS preprocessing: `03_munge.R`, `30_munge_ARISK.R`, `33_prepare_remaining_gwas.R`. Check header mappings, signed effect directions, N/effective N, SNP identifiers, duplicates, and exclusion rules.
4. Review QC and provenance: `26_audit_psychopathology_gwas.R`, `27_validate_gwas_statistics.R`, `34_check_ldsc_coverage.R`, `39_audit_cross_gwas_alleles.R`, `40_make_strict_allele_inputs.R`, `52_gwas_provenance_audit.R`, `53_annotate_gwas_provenance.R`, `56_validate_gwas_effect_directions.R`.
5. For age-specific European-ancestry temperament inputs, use `analysis/gwas_qc/07_prepare_age_specific.R`; first check consistency with historical filtering conventions.
6. Produce LDSC univariate heritability and joint covariance outputs, then the covariance-aware contrasts; confirm parameter ordering and compare SE conventions.
7. Generate tables and plots from final full-precision CSV/RDS outputs, **not manually transcribed snapshots**.

## Known caveats to resolve

- `CINT` SNP-specific N reportedly represents repeated-measures effective sample size; independently verify source definitions.
- `CADHD` sample counts differ between metadata and unique-participant descriptions; document the effective-N rationale.
- ARISK allele filtering and the exclusion of `rs10038745` must be reproduced by the final processing script.
- GenomicSEM's `V_Stand` SE conventions differ from fully propagated delta-method SEs; independent contrast verification remains desirable.
- Original scripts may contain absolute local paths, and may assume specific directories or software installations. Inspect and parameterise *copies*, then compare outputs with archived results.
- Any manuscript FDR claims describe exploratory preselected test families; ensure the scope of testing is clearly stated.

## Verification gates before public release

- [ ] Every required raw file has stable source reference and checksum
- [ ] No credentials, private data, or bulk GWAS committed
- [ ] Original and refactored script SHA256 checksums documented
- [ ] Clean-environment workflow successfully regenerated munged files
- [ ] QC metrics match archived runs within declared tolerance
- [ ] LDSC and contrast outputs match archived results
- [ ] Manuscript figures and tables built from exact final outputs
- [ ] Dependencies and computing environment pinned
- [ ] Release tagged, archived and cited with version-specific DOI
