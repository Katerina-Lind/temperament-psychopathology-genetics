library(GenomicSEM)

root <- path.expand("~/temperament_project")
temperament <- c("EMO", "ACT", "SHY", "SOC")

ld <- paste0(root, "/references/ldsc/eur_w_ld_chr/")

outdir <- file.path(root, "results/ldsc/temperament_childhood")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (outcome in c("CINT", "CADHD")) {

  traits <- c(temperament, outcome)

  files <- c(
    file.path(
      root,
      "gwas/processed/sensitivity/common_highN90",
      paste0(temperament, ".sumstats.gz")
    ),
    file.path(
      root,
      "gwas/processed/munged/childhood",
      paste0(outcome, ".sumstats.gz")
    )
  )

  stopifnot(all(file.exists(files)))

  cat("\n=== ANALYSING", outcome, "===\n")

  res <- GenomicSEM::ldsc(
    traits = files,
    sample.prev = rep(NA_real_, 5),
    population.prev = rep(NA_real_, 5),
    ld = ld,
    wld = ld,
    trait.names = traits,
    ldsc.log = file.path(outdir, paste0("temperament_", outcome))
  )

  saveRDS(
    res,
    file.path(outdir, paste0("temperament_", outcome, ".rds"))
  )

  S <- res$S
  dimnames(S) <- list(traits, traits)

  rg <- cov2cor(S)

  cat("\nGenetic correlations with", outcome, ":\n")
  print(round(rg[temperament, outcome], 3))

  cat("\nS eigenvalues:\n")
  print(eigen(S, symmetric = TRUE)$values)

  cat("\nV eigenvalues:\n")
  print(eigen(res$V, symmetric = TRUE)$values)
}
