# Figure legends and manuscript table notes

**Dataset freeze:** 23 November 2026. Source: GWAS LDSC and GenomicSEM console results supplied in this conversation. Figures and tables are **draft manuscript assets** until the final exported R output and all parameter-ordering checks are verified.

## Figure 1. Genetic correlation structure of early temperament and psychopathology
Lower triangle of the genetic correlation matrix estimated by joint multivariate LD score regression (LDSC) for emotionality, activity, shyness, childhood ADHD symptoms, adult internalising factor, adult substance-use factor, and risk-taking. Cell values are genetic correlations (rg); diagonal is one by definition. Estimates use European-ancestry GWAS summary statistics and a European reference LD panel. The display is descriptive; cell colours are not tests of statistical significance.

## Figure 2. Covariance-aware cross-phenotype genetic correlation contrasts
Points show estimated differences in genetic correlations between the pairs of outcomes for a common temperament trait. Horizontal intervals indicate 95% confidence intervals using the joint LDSC sampling covariance matrix and a delta-method transformation. Values of zero indicate no difference. Benjamini–Hochberg q values were computed over **four exploratory contrasts**, not over the entire set of genetic correlation comparisons considered during the project. These tests compare different phenotype definitions and should not be interpreted as longitudinal changes.

## Figure 3. Age-specific genetic correlation differences
Difference between genetic correlation of the psychiatric phenotype with temperament assessed in the second versus third year (second minus third year). Each comparison fits a trivariate LDSC model including the two age-specific temperament GWAS and a common psychiatric outcome. Points and horizontal bars show estimated differences and corresponding 95% confidence intervals, calculated using the joint sampling covariance matrix. Highlighted points denote Benjamini–Hochberg q < 0.05 among 18 exploratory contrasts. The analysis does not establish changes within individuals; samples and instruments may differ by age.

## Figure 4. Genetic correlations by temperament assessment period
Paired-point plots of genetic correlation estimates for the second-year (circles) and third-year (squares) temperament GWAS with the same psychiatric outcome. Connecting lines aid visual comparison and **are not longitudinal trajectories**. No individual-correlation confidence intervals are shown, since corresponding standard errors were not included in the supplied age-specific contrast output. Formal uncertainty for their differences is shown in Figure 3.

## Tables
- **Table 1:** Working GWAS source/provenance inventory. Blank fields are intentionally unverified; do not present them as zero or as missing study participants.
- **Table 2:** 21 distinct correlations from the seven-trait joint LDSC matrix. SE values copied from the joint LDSC result text.
- **Table 3:** Four covariance-aware exploratory cross-phenotype contrasts. Q values apply only within this set of four.
- **Table 4:** All 18 age-specific contrasts, including contrast standard errors, confidence intervals, and matrix diagnostics. Q values apply only within the 18 tests.
- **Table 5:** Six age-specific SNP-heritability estimates, standard errors, z ratios, and LDSC intercepts.

## Before submission
1. Run `scripts/build_figures.py --age-csv ~/temperament_psychopathology/results/genomicsem/age_specific/age_specific_contrasts.csv` to regenerate **all age figures and Table 4 from saved source**, replacing numbers transcribed from console excerpts.
2. Independently validate the parameter ordering of GenomicSEM's `V` and the covariance-aware contrast standard errors. Positive eigenvalues alone do not validate this computation.
3. Validate allelic matching, sources, sample overlap, effective sample sizes, and phenotype operationalisation; in particular clarify childhood internalising repeated-measures N and adult factor GWAS N.
4. Define a multiplicity strategy covering all scientific comparisons or clearly label subset contrasts exploratory.
5. Compare the main claims with prior published Hollowell and PGC analyses before deciding novelty, primary-versus-supplementary figures, and journal-specific figure sizes.
