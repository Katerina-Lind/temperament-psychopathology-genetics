# Statistical and editorial caveats (working notes)

- Estimates here are drawn from investigator-shared console output and may be rounded. Regenerate the definitive tables from saved `.rds` and `.csv` outputs on the analyst's machine.
- Delta-method contrasts use GenomicSEM LDSC sampling covariance V. The covariance vector ordering and variance convention must be cross-checked against the installed GenomicSEM version before treating P values as publication-final. A direct `rgmodel()` call failed within GenomicSEM/lavaan; it did not invalidate LDSC, but it left one implementation check outstanding.
- Exploratory BH corrections apply separately to four selected cross-phenotype contrasts and 18 age-specific comparisons, not all analyses conducted during study development.
- The seven-trait S matrix is positive definite but has condition number ~46.7; shyness 2Y/3Y rg ~0.972 implies nearly collinear traits in trivariate covariance matrices.
- High genetic correlation between age-specific temperament measures does not demonstrate within-person stability or independence between the GWAS datasets.
- Observational genetic correlations cannot establish temporal causality. Cross-phenotype contrasts of ADHD and adult substance use/risk-taking do not isolate age effects.
- Adult internalising and substance use are PGC latent factors, not interchangeable with diagnosis-specific disorders. Childhood ADHD includes symptom/diagnosis measurements, which could overlap conceptually with EAS activity.
- CINT file uses SNP-specific `n` values from a repeated-measures analysis; the source documentation gives 64,641 participants and SNP-specific N up to ~132,261. Confirm theoretical suitability of N for LDSC before making it a primary inference.
- Need explicit participant overlap, phenotype definitions, standard error checks, and source provenance table before public manuscript submission.
