x <- readRDS("results/genomicsem/developmental_7traits.rds")

traits <- c("EMO", "ACT", "SHY", "CADHD",
            "AINT", "ASUD", "ARISK")

S <- x$S
V <- x$V

dimnames(S) <- list(traits, traits)

# GenomicSEM sampling covariance ordering:
# column-major lower triangle including diagonal
idx <- which(lower.tri(S, diag=TRUE), arr.ind=TRUE)

stopifnot(nrow(idx) == 28, all(dim(V) == c(28,28)))

gradient_rg <- function(a, b) {
  i <- match(a, traits)
  j <- match(b, traits)

  h1 <- S[i,i]
  h2 <- S[j,j]
  g  <- S[i,j]

  grad <- numeric(nrow(idx))

  for (k in seq_len(nrow(idx))) {
    r <- idx[k,1]
    c <- idx[k,2]

    if (r == i && c == i)
      grad[k] <- -g / (2 * h1^(3/2) * sqrt(h2))

    if (r == j && c == j)
      grad[k] <- -g / (2 * h2^(3/2) * sqrt(h1))

    if (r != c &&
        ((r == i && c == j) || (r == j && c == i)))
      grad[k] <- 1 / sqrt(h1*h2)
  }

  grad
}

results <- list()

for (i in 1:6) {
  for (j in (i+1):7) {

    a <- traits[i]
    b <- traits[j]

    g <- gradient_rg(a, b)

    rg <- S[i,j] / sqrt(S[i,i]*S[j,j])
    se <- sqrt(drop(t(g) %*% V %*% g))

    results[[length(results)+1]] <- data.frame(
      trait1=a,
      trait2=b,
      rg=rg,
      SE=se
    )
  }
}

results <- do.call(rbind, results)

print(results, row.names=FALSE, digits=5)

write.csv(
  results,
  "results/genomicsem/developmental_rg_delta_validation.csv",
  row.names=FALSE
)
