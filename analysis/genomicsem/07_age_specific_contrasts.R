suppressPackageStartupMessages(library(GenomicSEM))

root <- Sys.getenv("TEMPERAMENT_ANALYSIS_ROOT")
if (!nzchar(root)) {
  stop("Set TEMPERAMENT_ANALYSIS_ROOT")
}
root <- normalizePath(root, mustWork=TRUE)
setwd(root)

age_dir <- "data/ldsc_munged/temperament_age_specific"

outcomes <- c(
  CADHD = "data/ldsc_munged/childhood/CADHD.sumstats.gz",
  AINT  = "data/ldsc_munged/adult/AINT.sumstats.gz",
  ASUD  = "data/ldsc_munged/adult/ASUD.sumstats.gz",
  ARISK = "data/ldsc_munged/cleaned/ARISK_final.sumstats.gz",
  WOR   = "data/ldsc_munged/worry.sumstats.gz",
  AFF   = "data/ldsc_munged/depressed_affect.sumstats.gz"
)

ref <- "data/ldsc_reference/eur_w_ld_chr/"
outdir <- "results/genomicsem/age_specific"

# Correlation between two traits in a covariance matrix
rg <- function(S, i, j) {
  S[i,j] / sqrt(S[i,i] * S[j,j])
}

# Finite-difference gradient for a trivariate covariance matrix
gradient <- function(S, fun) {
  lower <- lower.tri(S, diag=TRUE)
  p <- S[lower]

  rebuild <- function(v) {
    M <- matrix(0, 3, 3)
    M[lower] <- v
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
    M
  }

  g <- numeric(length(p))

  for (k in seq_along(p)) {
    h <- 1e-6 * max(1, abs(p[k]))
    p1 <- p2 <- p
    p1[k] <- p1[k] + h
    p2[k] <- p2[k] - h
    g[k] <- (fun(rebuild(p1)) - fun(rebuild(p2))) / (2*h)
  }

  g
}

results <- list()
failures <- list()

for (trait in c("EMO", "ACT", "SHY")) {
  for (outcome in names(outcomes)) {

    id <- paste(trait, outcome, sep="_")
    cat("\n==========", id, "==========\n")

    paths <- c(
      file.path(age_dir, paste0(trait, "_2Y.sumstats.gz")),
      file.path(age_dir, paste0(trait, "_3Y.sumstats.gz")),
      outcomes[[outcome]]
    )

    if (!all(file.exists(paths))) {
      failures[[id]] <- "Missing input file"
      next
    }

    result <- tryCatch(
      GenomicSEM::ldsc(
        traits=paths,
        sample.prev=rep(NA_real_, 3),
        population.prev=rep(NA_real_, 3),
        ld=ref,
        wld=ref,
        trait.names=c("Y2","Y3","OUT"),
        n.blocks=200,
        stand=TRUE,
        ldsc.log=file.path(outdir, id)
      ),
      error=function(e) {
        failures[[id]] <<- conditionMessage(e)
        NULL
      }
    )

    if (is.null(result)) next

    saveRDS(result, file.path(outdir, paste0(id, ".rds")))

    S <- result$S
    V <- result$V

    if (!all(dim(S) == c(3,3)) ||
        !all(dim(V) == c(6,6)) ||
        !all(is.finite(S)) ||
        !all(is.finite(V)) ||
        any(diag(S) <= 0)) {
      failures[[id]] <- "Invalid covariance matrices"
      next
    }

    # Check age-to-age genetic correlation as a diagnostic.
    rg_age <- rg(S, 1, 2)

    f <- function(M) rg(M,1,3) - rg(M,2,3)
    g <- gradient(S, f)

    difference <- f(S)
    variance <- as.numeric(t(g) %*% V %*% g)

    if (!is.finite(variance) || variance <= 0) {
      failures[[id]] <- "Invalid contrast variance"
      next
    }

    se <- sqrt(variance)
    z <- difference / se
    p <- 2 * pnorm(-abs(z))

    results[[id]] <- data.frame(
      temperament=trait,
      outcome=outcome,
      rg_2Y=rg(S,1,3),
      rg_3Y=rg(S,2,3),
      rg_2Y_3Y=rg_age,
      difference_2Y_minus_3Y=difference,
      SE_difference=se,
      Z=z,
      P=p,
      CI_lower=difference-1.96*se,
      CI_upper=difference+1.96*se,
      min_eigenvalue_S=min(eigen(S, symmetric=TRUE)$values),
      min_eigenvalue_V=min(
        eigen((V+t(V))/2, symmetric=TRUE)$values
      )
    )
  }
}

if (length(results)) {
  tab <- do.call(rbind, results)
  tab$FDR_BH <- p.adjust(tab$P, method="BH")

  write.csv(
    tab,
    file.path(outdir, "age_specific_contrasts.csv"),
    row.names=FALSE
  )

  cat("\n=== AGE-SPECIFIC CONTRAST RESULTS ===\n")
  print(
    tab[, c(
      "temperament","outcome",
      "rg_2Y","rg_3Y",
      "difference_2Y_minus_3Y",
      "SE_difference","P","FDR_BH"
    )],
    row.names=FALSE,
    digits=4
  )
}

if (length(failures)) {
  cat("\n=== FAILED ANALYSES ===\n")
  print(failures)
}

cat("\nCompleted:", length(results), "of 18 analyses\n")
