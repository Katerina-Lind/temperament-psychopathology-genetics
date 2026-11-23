library(GenomicSEM)
library(data.table)

root <- path.expand("~/temperament_project")

temperament <- c("EMO", "ACT", "SHY", "SOC")
childhood <- c("CINT", "CADHD", "CAGG")
adult <- c("AINT", "ASUD", "ARISK")

traits <- c(temperament, childhood, adult)

inputdir <- file.path(
  root, "gwas/processed/sensitivity/strict_alleles"
)

files <- file.path(
  inputdir, paste0(traits, ".sumstats.gz")
)

stopifnot(all(file.exists(files)))

ld <- file.path(root, "references/ldsc/eur_w_ld_chr/")

outdir <- file.path(root, "results/ldsc/full_ten_traits")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cat("\n=== TEN-TRAIT JOINT LDSC ===\n")
cat("Traits:", paste(traits, collapse=", "), "\n")

res <- GenomicSEM::ldsc(
  traits=files,
  sample.prev=rep(NA_real_, length(traits)),
  population.prev=rep(NA_real_, length(traits)),
  ld=ld,
  wld=ld,
  trait.names=traits,
  ldsc.log=file.path(outdir, "full_ten_traits")
)

saveRDS(
  res,
  file.path(outdir, "full_ten_traits_ldsc.rds")
)

S <- res$S
V <- res$V

dimnames(S) <- list(traits, traits)

cat("\n=== HERITABILITIES ===\n")
print(round(diag(S), 5))

cat("\n=== MATRIX DIAGNOSTICS ===\n")

eigS <- eigen(S, symmetric=TRUE)$values
eigV <- eigen(V, symmetric=TRUE)$values

cat("Minimum S eigenvalue:", min(eigS), "\n")
cat("Minimum V eigenvalue:", min(eigV), "\n")
cat("S positive definite:", all(eigS > 0), "\n")
cat("V positive definite:", all(eigV > 0), "\n")

# Genetic correlations require positive
# diagonal genetic variances.
if (any(!is.finite(diag(S))) || any(diag(S) <= 0)) {
  stop("Non-positive or invalid heritability detected. ",
       "Inspect S before calculating correlations.")
}

rg <- cov2cor(S)

cat("\n=== TEMPERAMENT × PSYCHOPATHOLOGY RG ===\n")
print(round(
  rg[temperament, c(childhood,adult),drop=FALSE],
  3
))

# GenomicSEM V uses the lower triangle of S,
# including the diagonal, in column-major order.
idx <- which(lower.tri(S, diag=TRUE), arr.ind=TRUE)

stopifnot(
  nrow(idx)==nrow(V),
  nrow(V)==ncol(V)
)

position <- function(i,j) {
  row <- max(i,j)
  col <- min(i,j)
  match_pos <- which(
    idx[,1]==row & idx[,2]==col
  )
  stopifnot(length(match_pos)==1L)
  match_pos
}

# Delta-method SE for rg = cov_ij/sqrt(h2_i*h2_j)
rg_se <- function(i,j) {

  hi <- S[i,i]
  hj <- S[j,j]
  rij <- rg[i,j]

  gradient <- numeric(nrow(V))

  gradient[position(i,j)] <- 1/sqrt(hi*hj)
  gradient[position(i,i)] <- -rij/(2*hi)
  gradient[position(j,j)] <- -rij/(2*hj)

  variance <- as.numeric(
    crossprod(gradient, V %*% gradient)
  )

  if (!is.finite(variance) || variance < 0) {
    return(NA_real_)
  }

  sqrt(variance)
}

outcomes <- c(childhood,adult)

results <- rbindlist(lapply(temperament, function(t) {
  i <- match(t,traits)

  rbindlist(lapply(outcomes, function(o) {
    j <- match(o,traits)
    se <- rg_se(i,j)
    estimate <- rg[i,j]

    data.table(
      temperament=t,
      outcome=o,
      stage=ifelse(o %in% childhood,
                   "childhood","adult"),
      rg=estimate,
      SE=se,
      Z=estimate/se,
      p=2*pnorm(-abs(estimate/se))
    )
  }))
}))

results[, FDR_24 := p.adjust(p, method="BH")]

setcolorder(
  results,
  c("temperament","outcome","stage",
    "rg","SE","Z","p","FDR_24")
)

cat("\n=== MAIN RESULTS: 24 CORRELATIONS ===\n")
print(results, digits=4)

fwrite(
  results,
  file.path(outdir, "temperament_psychopathology_24_rg.tsv"),
  sep="\t"
)

fwrite(
  data.table(trait=traits, h2=diag(S)),
  file.path(outdir, "heritabilities.tsv"),
  sep="\t"
)

fwrite(
  data.table(trait=rownames(rg), rg),
  file.path(outdir, "genetic_correlation_matrix.tsv"),
  sep="\t"
)

cat("\n=== SIGNIFICANCE SUMMARY ===\n")
cat("Nominal p < 0.05:",
    sum(results$p < 0.05, na.rm=TRUE), "\n")
cat("FDR < 0.05:",
    sum(results$FDR_24 < 0.05, na.rm=TRUE), "\n")

cat("\n=== COMPLETE ===\n")
cat("Results:", outdir, "\n")
