suppressPackageStartupMessages(library(GenomicSEM))

root <- Sys.getenv("TEMPERAMENT_ANALYSIS_ROOT")
if (!nzchar(root)) {
  stop("Set TEMPERAMENT_ANALYSIS_ROOT")
}
root <- normalizePath(root, mustWork=TRUE)
setwd(root)

trait.names <- c(
  "EMO", "ACT", "SHY",
  "CADHD", "AINT", "ASUD", "ARISK"
)

traits <- c(
  "data/ldsc_munged/temperament/EMO.sumstats.gz",
  "data/ldsc_munged/temperament/ACT.sumstats.gz",
  "data/ldsc_munged/temperament/SHY.sumstats.gz",
  "data/ldsc_munged/childhood/CADHD.sumstats.gz",
  "data/ldsc_munged/adult/AINT.sumstats.gz",
  "data/ldsc_munged/adult/ASUD.sumstats.gz",
  "data/ldsc_munged/cleaned/ARISK_final.sumstats.gz"
)

stopifnot(
  length(traits) == length(trait.names),
  all(file.exists(traits))
)

dir.create(
  "results/genomicsem",
  recursive=TRUE,
  showWarnings=FALSE
)

ref <- "data/ldsc_reference/eur_w_ld_chr/"

cat("Developmental joint LDSC\n")
cat("Traits:", paste(trait.names, collapse=", "), "\n")

result <- GenomicSEM::ldsc(
  traits=traits,
  sample.prev=rep(NA_real_, length(traits)),
  population.prev=rep(NA_real_, length(traits)),
  ld=ref,
  wld=ref,
  trait.names=trait.names,
  n.blocks=200,
  stand=TRUE,
  ldsc.log="results/genomicsem/developmental_7traits"
)

S <- result$S
V <- result$V

dimnames(S) <- list(trait.names, trait.names)

saveRDS(
  result,
  "results/genomicsem/developmental_7traits.rds"
)

cat("\n=== Genetic correlation matrix ===\n")
print(round(cov2cor(S), 3))

cat("\n=== Eigenvalues of S ===\n")
eigS <- eigen(S, symmetric=TRUE)$values
print(eigS)

cat("\n=== Matrix diagnostics ===\n")
cat("S positive definite:", all(eigS > 0), "\n")
cat("Condition number S:", kappa(S, exact=TRUE), "\n")
cat("V dimensions:", paste(dim(V), collapse=" x "), "\n")
cat("V symmetric:", isSymmetric(V, tol=1e-8), "\n")
cat("V finite:", all(is.finite(V)), "\n")

eigV <- eigen(
  (V+t(V))/2,
  symmetric=TRUE,
  only.values=TRUE
)$values

cat("Minimum eigenvalue V:", min(eigV), "\n")

cat("\nSaved developmental LDSC object.\n")
