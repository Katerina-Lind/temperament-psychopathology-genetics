root <- path.expand("~/temperament_psychopathology")
setwd(root)

x <- readRDS(
  "results/genomicsem/joint_ldsc_7traits.rds"
)

traits <- c("EMO", "ACT", "SHY", "WOR", "AFF", "ASB", "ALC")

S <- x$S
V <- x$V

dimnames(S) <- list(traits, traits)

stopifnot(
  identical(dim(S), c(7L, 7L)),
  identical(dim(V), c(28L, 28L)),
  all(is.finite(S)),
  all(is.finite(V))
)

cat("\n=== Matrix diagnostics ===\n")

eig_S <- eigen(S, symmetric=TRUE)$values
eig_V <- eigen((V + t(V))/2, symmetric=TRUE)$values

cat("Minimum eigenvalue S:", min(eig_S), "\n")
cat("Minimum eigenvalue V:", min(eig_V), "\n")
cat("S positive definite:", min(eig_S) > 0, "\n")
cat("V symmetric:", isSymmetric(V, tol=1e-8), "\n")

stopifnot(
  isSymmetric(S, tol=1e-8),
  isSymmetric(V, tol=1e-8),
  min(eig_S) > 0,
  min(eig_V) > -1e-8
)

# GenomicSEM stores V in lower-triangle,
# column-major order, including diagonal.
lower <- lower.tri(S, diag=TRUE)
theta <- S[lower]

rebuild_S <- function(p) {
  M <- matrix(0, 7, 7)
  M[lower] <- p
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  dimnames(M) <- list(traits, traits)
  M
}

stopifnot(isTRUE(all.equal(rebuild_S(theta), S)))

# Return genetic correlation from genetic covariance matrix.
rg <- function(M, a, b) {
  M[a, b] / sqrt(M[a, a] * M[b, b])
}

# Numerical delta-method gradient with respect
# to the 28 independent elements of S.
gradient <- function(fun, p) {
  g <- numeric(length(p))

  for (k in seq_along(p)) {
    h <- 1e-6 * max(1, abs(p[k]))

    p_plus <- p
    p_minus <- p

    p_plus[k] <- p_plus[k] + h
    p_minus[k] <- p_minus[k] - h

    g[k] <- (
      fun(rebuild_S(p_plus)) -
      fun(rebuild_S(p_minus))
    ) / (2*h)
  }

  g
}

# First validate the transformed individual SEs.
cat("\n=== Selected genetic correlations ===\n")

for (pair in list(
  c("SHY", "WOR"),
  c("SHY", "AFF"),
  c("ACT", "ASB"),
  c("ACT", "ALC")
)) {
  fun <- function(M) rg(M, pair[1], pair[2])
  g <- gradient(fun, theta)

  se <- sqrt(as.numeric(t(g) %*% V %*% g))

  cat(
    pair[1], pair[2],
    "rg =", round(fun(S), 4),
    "SE =", round(se, 4), "\n"
  )
}

# Pre-specified exploratory contrasts
contrasts <- list(
  SHY_WOR_vs_AFF = list(
    a=c("SHY","WOR"),
    b=c("SHY","AFF")
  ),
  ACT_ASB_vs_ALC = list(
    a=c("ACT","ASB"),
    b=c("ACT","ALC")
  ),
  EMO_WOR_vs_AFF = list(
    a=c("EMO","WOR"),
    b=c("EMO","AFF")
  ),
  ACT_ASB_vs_WOR = list(
    a=c("ACT","ASB"),
    b=c("ACT","WOR")
  )
)

results <- lapply(names(contrasts), function(name) {

  spec <- contrasts[[name]]

  f1 <- function(M) rg(M, spec$a[1], spec$a[2])
  f2 <- function(M) rg(M, spec$b[1], spec$b[2])

  g1 <- gradient(f1, theta)
  g2 <- gradient(f2, theta)

  g_delta <- g1 - g2

  difference <- f1(S) - f2(S)

  variance <- as.numeric(
    t(g_delta) %*% V %*% g_delta
  )

  if (!is.finite(variance) || variance <= 0) {
    stop(paste("Invalid variance:", name, variance))
  }

  se <- sqrt(variance)
  z <- difference / se
  p <- 2 * pnorm(-abs(z))

  data.frame(
    contrast=name,
    rg_1=f1(S),
    rg_2=f2(S),
    difference=difference,
    SE=se,
    z=z,
    p=p,
    CI_lower=difference - 1.96*se,
    CI_upper=difference + 1.96*se
  )
})

results <- do.call(rbind, results)
results$FDR_BH <- p.adjust(results$p, method="BH")

cat("\n=== Genetic correlation contrasts ===\n")
print(results, digits=4, row.names=FALSE)

dir.create(
  "results/genomicsem",
  recursive=TRUE,
  showWarnings=FALSE
)

write.csv(
  results,
  "results/genomicsem/rg_contrasts_7traits.csv",
  row.names=FALSE
)

cat("\nResults saved successfully.\n")
