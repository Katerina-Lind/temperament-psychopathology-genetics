library(data.table)

root <- path.expand("~/temperament_project")

res <- readRDS(file.path(
  root,
  "results/ldsc/full_individual_highN",
  "full_ten_traits_ldsc.rds"
))

traits <- c(
  "EMO","ACT","SHY","SOC",
  "CINT","CADHD","CAGG",
  "AINT","ASUD","ARISK"
)

S <- unname(res$S)
V <- res$V
n <- length(traits)

idx <- which(lower.tri(S, diag=TRUE), arr.ind=TRUE)
stopifnot(nrow(idx)==nrow(V), isSymmetric(S))

theta <- S[idx]

reconstruct <- function(v) {
  M <- matrix(0, n, n)
  M[idx] <- v
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}

position <- function(i,j) {
  p <- which(
    idx[,1]==max(i,j) &
    idx[,2]==min(i,j)
  )
  stopifnot(length(p)==1)
  p
}

rgfun <- function(v,i,j) {
  M <- reconstruct(v)
  M[i,j]/sqrt(M[i,i]*M[j,j])
}

check <- function(a,b) {
  i <- match(a,traits)
  j <- match(b,traits)

  r <- rgfun(theta,i,j)
  h1 <- S[i,i]
  h2 <- S[j,j]

  g <- numeric(length(theta))
  g[position(i,j)] <- 1/sqrt(h1*h2)
  g[position(i,i)] <- -r/(2*h1)
  g[position(j,j)] <- -r/(2*h2)

  # Only three coordinates can have nonzero derivatives.
  active <- unique(c(
    position(i,j),position(i,i),position(j,j)
  ))

  gn <- numeric(length(theta))
  for (k in active) {
    step <- 1e-7 * max(1,abs(theta[k]))
    plus <- theta
    minus <- theta
    plus[k] <- plus[k]+step
    minus[k] <- minus[k]-step

    gn[k] <- (
      rgfun(plus,i,j)-rgfun(minus,i,j)
    )/(2*step)
  }

  se <- sqrt(as.numeric(crossprod(g,V%*%g)))

  data.table(
    temperament=a,
    outcome=b,
    rg=r,
    SE_delta=se,
    max_gradient_error=max(abs(g-gn))
  )
}

eas <- traits[1:4]
outcomes <- traits[5:10]

results <- rbindlist(lapply(eas, function(a) {
  rbindlist(lapply(outcomes, function(b) {
    check(a,b)
  }))
}))

print(results,digits=5)

cat("\nMaximum gradient discrepancy:",
    max(results$max_gradient_error),"\n")

stopifnot(
  all(is.finite(results$SE_delta)),
  max(results$max_gradient_error) < 1e-5
)

out <- file.path(
  root,"results/ldsc/full_individual_highN",
  "gradient_validation.tsv"
)

fwrite(results,out,sep="\t")

cat("\nPASS: analytical and numerical gradients agree.\n")
cat("Saved:",out,"\n")
