library(data.table)

root <- path.expand("~/temperament_project")

infile <- file.path(
  root,
  "results/ldsc/full_ten_traits/full_ten_traits_ldsc.rds"
)

res <- readRDS(infile)

traits <- c(
  "EMO", "ACT", "SHY", "SOC",
  "CINT", "CADHD", "CAGG",
  "AINT", "ASUD", "ARISK"
)

S <- res$S
V <- res$V

rownames(S) <- traits
colnames(S) <- traits

stopifnot(
  identical(dim(S), c(10L,10L)),
  identical(dim(V), c(55L,55L)),
  isSymmetric(S),
  isSymmetric(V),
  all(diag(S) > 0)
)

# GenomicSEM ordering:
# lower triangle, column-wise, including diagonal
idx <- which(lower.tri(S, diag=TRUE), arr.ind=TRUE)

position <- function(i,j) {
  row <- max(i,j)
  col <- min(i,j)
  p <- which(idx[,1] == row & idx[,2] == col)
  stopifnot(length(p) == 1L)
  p
}

# Gradient of rg(i,j) with respect to vech(S)
rg_gradient <- function(a,b) {
  i <- match(a, traits)
  j <- match(b, traits)

  hi <- S[i,i]
  hj <- S[j,j]
  covij <- S[i,j]

  rg <- covij / sqrt(hi*hj)

  g <- numeric(nrow(V))

  g[position(i,j)] <- 1/sqrt(hi*hj)
  g[position(i,i)] <- -rg/(2*hi)
  g[position(j,j)] <- -rg/(2*hj)

  list(rg=rg, gradient=g)
}

# Compare two correlations sharing one EAS trait
contrast <- function(temperament, childhood, adult) {
  x <- rg_gradient(temperament, childhood)
  y <- rg_gradient(temperament, adult)

  delta <- x$rg - y$rg
  gdiff <- x$gradient - y$gradient

  variance <- as.numeric(
    crossprod(gdiff, V %*% gdiff)
  )

  stopifnot(is.finite(variance), variance > 0)

  se <- sqrt(variance)
  z <- delta/se
  p <- 2*pnorm(-abs(z))

  data.table(
    temperament=temperament,
    childhood=childhood,
    adult=adult,
    rg_child=x$rg,
    rg_adult=y$rg,
    difference=delta,
    SE_difference=se,
    Z=z,
    p=p,
    ci_lower=delta-1.96*se,
    ci_upper=delta+1.96*se
  )
}

# Four prespecified exploratory contrasts
comparisons <- data.table(
  temperament=c("EMO","SHY","ACT","SHY"),
  childhood=c("CINT","CINT","CADHD","CADHD"),
  adult=c("AINT","AINT","ASUD","ARISK")
)

results <- rbindlist(lapply(
  seq_len(nrow(comparisons)),
  function(k) {
    contrast(
      comparisons$temperament[k],
      comparisons$childhood[k],
      comparisons$adult[k]
    )
  }
))

results[, FDR_4 := p.adjust(p, method="BH")]

cat("\n=== CHILDHOOD–ADULT GENETIC CORRELATION CONTRASTS ===\n")
print(results, digits=4)

outdir <- file.path(
  root, "results/ldsc/developmental_contrasts"
)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

fwrite(
  results,
  file.path(outdir, "developmental_contrasts.tsv"),
  sep="\t"
)

cat("\n=== SUMMARY ===\n")
cat("Contrasts tested:", nrow(results), "\n")
cat("Nominal p < 0.05:", sum(results$p < 0.05), "\n")
cat("FDR < 0.05:", sum(results$FDR_4 < 0.05), "\n")
cat("\nSaved to:", outdir, "\n")
