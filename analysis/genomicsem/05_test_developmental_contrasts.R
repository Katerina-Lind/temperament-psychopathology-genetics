root <- Sys.getenv("TEMPERAMENT_ANALYSIS_ROOT")
if (!nzchar(root)) {
  stop("Set TEMPERAMENT_ANALYSIS_ROOT")
}
root <- normalizePath(root, mustWork=TRUE)
setwd(root)

x <- readRDS(
  "results/genomicsem/developmental_7traits.rds"
)

traits <- c(
  "EMO", "ACT", "SHY",
  "CADHD", "AINT", "ASUD", "ARISK"
)

S <- x$S
V <- x$V

dimnames(S) <- list(traits, traits)

stopifnot(
  all(dim(S) == c(7, 7)),
  all(dim(V) == c(28, 28)),
  all(is.finite(S)),
  all(is.finite(V)),
  isSymmetric(S, tol=1e-8),
  isSymmetric(V, tol=1e-8)
)

eigS <- eigen(S, symmetric=TRUE)$values
eigV <- eigen(V, symmetric=TRUE)$values

cat("Minimum eigenvalue S:", min(eigS), "\n")
cat("Minimum eigenvalue V:", min(eigV), "\n")

stopifnot(min(eigS) > 0, min(eigV) > 0)

# GenomicSEM V follows the column-major lower triangle
# of S, including its diagonal.
lower <- lower.tri(S, diag=TRUE)
theta <- S[lower]

rebuild <- function(p) {
  M <- matrix(0, 7, 7)
  M[lower] <- p
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  dimnames(M) <- list(traits, traits)
  M
}

stopifnot(isTRUE(all.equal(rebuild(theta), S)))

rg <- function(M, a, b) {
  M[a,b] / sqrt(M[a,a] * M[b,b])
}

gradient <- function(fun, p) {
  g <- numeric(length(p))

  for (k in seq_along(p)) {
    h <- 1e-6 * max(1, abs(p[k]))

    plus <- minus <- p
    plus[k] <- plus[k] + h
    minus[k] <- minus[k] - h

    g[k] <- (
      fun(rebuild(plus)) - fun(rebuild(minus))
    ) / (2*h)
  }
  g
}

contrasts <- list(
  ACT_ADHD_vs_ASUD = list(
    c("ACT","CADHD"), c("ACT","ASUD")
  ),
  ACT_ADHD_vs_ARISK = list(
    c("ACT","CADHD"), c("ACT","ARISK")
  ),
  SHY_ADHD_vs_ARISK = list(
    c("SHY","CADHD"), c("SHY","ARISK")
  ),
  EMO_AINT_vs_ASUD = list(
    c("EMO","AINT"), c("EMO","ASUD")
  )
)

results <- lapply(names(contrasts), function(name) {
  pairs <- contrasts[[name]]

  f1 <- function(M) rg(M, pairs[[1]][1], pairs[[1]][2])
  f2 <- function(M) rg(M, pairs[[2]][1], pairs[[2]][2])

  g1 <- gradient(f1, theta)
  g2 <- gradient(f2, theta)

  # Check covariance-aware SE for both individual rg values
  se1 <- sqrt(drop(t(g1) %*% V %*% g1))
  se2 <- sqrt(drop(t(g2) %*% V %*% g2))

  gd <- g1 - g2
  variance <- drop(t(gd) %*% V %*% gd)

  stopifnot(is.finite(variance), variance > 0)

  d <- f1(S) - f2(S)
  se <- sqrt(variance)
  z <- d/se

  data.frame(
    contrast=name,
    rg1=f1(S),
    SE_rg1=se1,
    rg2=f2(S),
    SE_rg2=se2,
    difference=d,
    SE_difference=se,
    z=z,
    p=2*pnorm(-abs(z)),
    CI_lower=d-1.96*se,
    CI_upper=d+1.96*se
  )
})

results <- do.call(rbind, results)
results$FDR_BH <- p.adjust(results$p, method="BH")

cat("\n=== DEVELOPMENTAL CONTRAST RESULTS ===\n")
print(results, row.names=FALSE, digits=5)

write.csv(
  results,
  "results/genomicsem/developmental_rg_contrasts.csv",
  row.names=FALSE
)

cat("\nSaved developmental_rg_contrasts.csv\n")
