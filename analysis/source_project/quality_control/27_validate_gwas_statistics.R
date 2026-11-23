library(data.table)

root <- path.expand("~/temperament_project")
manifest <- fread(file.path(
  root, "config/temperament_psychopathology_datasets.tsv"
))
outdir <- file.path(root, "gwas/processed/qc")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

traits <- c("CINT", "CADHD", "CAGG",
            "AINT", "ASUD", "ARISK")

spec <- list(
  CINT  = c("z_score", "",               "",               "p_value"),
  CADHD = c("",        "beta",           "standard_error", "p_value"),
  CAGG  = c("Zscore",  "",               "",               "P-value"),
  AINT  = c("",        "BETA",           "SE",             "P"),
  ASUD  = c("",        "BETA",           "SE",             "P"),
  ARISK = c("",        "Beta",           "SE",             "p_bolt_lmm")
)

results <- list()

for (trait in traits) {
  cat("\n=== ", trait, " ===\n", sep = "")

  filename <- manifest[["file"]][match(trait, manifest[["trait"]])]
  stopifnot(length(filename) == 1L)

  x <- fread(
    file.path(root, filename),
    na.strings = c("NA", "#NA", "NaN", "")
  )

  s <- spec[[trait]]

  z <- if (nzchar(s[1])) {
    as.numeric(x[[s[1]]])
  } else {
    as.numeric(x[[s[2]]]) /
      as.numeric(x[[s[3]]])
  }

  p <- as.numeric(x[[s[4]]])

  valid <- is.finite(z) & is.finite(p) & p > 0 & p <= 1

  # Expected two-sided p-values under a normal Z statistic.
  expected_logp <- pchisq(
    z[valid]^2,
    df = 1,
    lower.tail = FALSE,
    log.p = TRUE
  ) / log(10)

  observed_logp <- log10(p[valid])

  discrepancy <- abs(expected_logp - observed_logp)

  # Restrict the principal concordance check to
  # p < 0.05, avoiding misleading relative
  # comparisons near p = 1.
  associated <- p[valid] < 0.05

  results[[trait]] <- data.table(
    trait = trait,
    variants = nrow(x),
    valid_statistics = sum(valid),
    median_abs_log10p_diff = median(discrepancy),
    pct_diff_gt_0_1 = 100 * mean(discrepancy > 0.1),
    pct_diff_gt_0_1_p_lt_0_05 =
      if (any(associated)) {
        100 * mean(discrepancy[associated] > 0.1)
      } else NA_real_
  )

  print(results[[trait]])

  if (trait == "CINT") {
    n <- as.numeric(x$n)
    cat("\nCINT: fractional N values:",
        sum(is.finite(n) & abs(n - round(n)) > 1e-6),
        "of", length(n), "\n")

    cat("CINT: N exceeding 64641:",
        sum(n > 64641, na.rm = TRUE), "\n")
  }

  if (trait == "ARISK") {
    cat("\n=== Risk-taking duplicate IDs ===\n")

    d <- x[
      duplicated(rsID) |
      duplicated(rsID, fromLast = TRUE)
    ]

    cat("Rows with duplicated rsIDs:", nrow(d), "\n")
    cat("Distinct duplicated rsIDs:", uniqueN(d$rsID), "\n")

    # Count identifiers with conflicting alleles.
    allele_pairs <- d[
      , .(
        allele_pairs = uniqueN(
          paste(Effect_allele, Other_allele, sep = "/")
        ),
        distinct_beta = uniqueN(Beta),
        distinct_se = uniqueN(SE)
      ),
      by = rsID
    ]

    cat("IDs with multiple allele pairs:",
        sum(allele_pairs$allele_pairs > 1), "\n")

    cat("IDs with differing beta:",
        sum(allele_pairs$distinct_beta > 1), "\n")

    cat("IDs with differing SE:",
        sum(allele_pairs$distinct_se > 1), "\n")

    fwrite(
      allele_pairs,
      file.path(outdir, "ARISK_duplicate_summary.tsv"),
      sep = "\t"
    )
  }

  rm(x, z, p, valid, expected_logp,
     observed_logp, discrepancy)
  gc(verbose = FALSE)
}

report <- rbindlist(results)

fwrite(
  report,
  file.path(outdir, "gwas_statistical_concordance.tsv"),
  sep = "\t"
)

cat("\n=== FINAL STATISTICAL QC ===\n")
print(report)
