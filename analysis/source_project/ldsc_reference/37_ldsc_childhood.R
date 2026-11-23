library(GenomicSEM)

root <- path.expand("~/temperament_project")

traits <- c("CINT", "CADHD", "CAGG")

files <- file.path(
  root,
  "gwas/processed/munged/childhood",
  paste0(traits, ".sumstats.gz")
)

stopifnot(all(file.exists(files)))

ld <- paste0(
  root, "/references/ldsc/eur_w_ld_chr/"
)

outdir <- file.path(root, "results/ldsc/childhood")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n=== CHILDHOOD LDSC ANALYSIS ===\n")

res <- GenomicSEM::ldsc(
  traits = files,
  sample.prev = rep(NA_real_, length(traits)),
  population.prev = rep(NA_real_, length(traits)),
  ld = ld,
  wld = ld,
  trait.names = traits,
  ldsc.log = file.path(outdir, "childhood_joint")
)

saveRDS(
  res,
  file.path(outdir, "childhood_joint_ldsc.rds")
)

S <- res$S
dimnames(S) <- list(traits, traits)

cat("\n=== GENETIC COVARIANCE MATRIX ===\n")
print(S)

cat("\n=== SNP HERITABILITIES ===\n")

# For three traits, GenomicSEM V follows the
# lower-triangular column-wise ordering.
positions <- c(1, 4, 6)

h2 <- diag(S)
se <- sqrt(diag(res$V)[positions])

summary <- data.frame(
  trait = traits,
  h2 = h2,
  SE = se,
  Z = h2 / se
)

print(summary)

write.table(
  summary,
  file.path(outdir, "childhood_heritability.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n=== GENETIC CORRELATIONS ===\n")

if (all(is.finite(diag(S))) && all(diag(S) > 0)) {
  print(round(cov2cor(S), 3))
} else {
  cat("At least one heritability is non-positive;",
      "genetic correlations may not be interpretable.\n")
}

cat("\n=== EIGENVALUES OF S ===\n")
print(eigen(S, symmetric = TRUE)$values)

cat("\n=== EIGENVALUES OF V ===\n")
print(eigen(res$V, symmetric = TRUE)$values)

cat("\n=== RESULT COMPONENTS ===\n")
print(names(res))

cat("\nResults saved to:", outdir, "\n")
