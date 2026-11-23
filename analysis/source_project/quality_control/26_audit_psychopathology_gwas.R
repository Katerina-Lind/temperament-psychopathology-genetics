library(data.table)

root <- path.expand("~/temperament_project")

manifest <- fread(file.path(
  root, "config/temperament_psychopathology_datasets.tsv"
))

traits <- c("CINT", "CADHD", "CAGG",
            "AINT", "ASUD", "ARISK")

manifest <- manifest[match(traits, trait)]
stopifnot(identical(manifest$trait, traits))

outdir <- file.path(root, "gwas/processed/qc")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pick <- function(cols, alternatives) {
  ix <- match(tolower(alternatives), tolower(cols))
  ix <- ix[!is.na(ix)]
  if (length(ix)) cols[ix[1]] else NA_character_
}

get_numeric <- function(x, col) {
  if (is.na(col)) return(NULL)
  suppressWarnings(as.numeric(x[[col]]))
}

describe_numeric <- function(v) {
  if (is.null(v)) return("COLUMN ABSENT")
  good <- v[is.finite(v)]
  if (!length(good)) return("NO FINITE VALUES")
  paste(
    "missing/nonfinite =", length(v) - length(good),
    "| min =", signif(min(good), 5),
    "| median =", signif(median(good), 5),
    "| max =", signif(max(good), 5)
  )
}

results <- list()
columns_report <- list()

for (i in seq_len(nrow(manifest))) {

  trait <- manifest$trait[i]
  path <- file.path(root, manifest$file[i])

  cat("\n====================================\n")
  cat("DATASET:", trait, "\n")
  cat("FILE:", path, "\n")
  cat("====================================\n")

  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }

  x <- fread(
    path,
    na.strings = c("NA", "#NA", "NaN", ""),
    showProgress = TRUE
  )

  cols <- names(x)

  cat("\nDimensions:", nrow(x), "rows x",
      ncol(x), "columns\n")

  cat("\nCOLUMN NAMES:\n")
  print(cols)

  idcol <- pick(cols, c(
    "rsid", "SNP", "MarkerName", "variant_id"
  ))

  chrcol <- pick(cols, c(
    "chromosome", "CHR", "chr"
  ))

  poscol <- pick(cols, c(
    "base_pair_location", "BP", "position", "pos"
  ))

  a1col <- pick(cols, c(
    "effect_allele", "A1", "Allele1", "EA"
  ))

  a2col <- pick(cols, c(
    "other_allele", "A2", "Allele2", "NEA"
  ))

  ncolname <- pick(cols, c(
    "n", "N", "SAMPLESIZE", "N_total"
  ))

  betacol <- pick(cols, c("beta", "BETA", "effect"))
  secol <- pick(cols, c(
    "standard_error", "SE", "STDERR", "se"
  ))

  zcol <- pick(cols, c(
    "z_score", "Zscore", "Z", "z"
  ))

  pcol <- pick(cols, c(
    "p_value", "P", "P-value", "p_bolt_lmm",
    "pvalue"
  ))

  freqcol <- pick(cols, c(
    "effect_allele_frequency",
    "MAF", "Freq1", "freq", "EAF"
  ))

  fields <- c(
    SNP = idcol,
    CHR = chrcol,
    BP = poscol,
    A1 = a1col,
    A2 = a2col,
    N = ncolname,
    BETA = betacol,
    SE = secol,
    Z = zcol,
    P = pcol,
    FREQ = freqcol
  )

  cat("\nDETECTED FIELDS:\n")
  print(fields)

  n <- get_numeric(x, ncolname)
  beta <- get_numeric(x, betacol)
  se <- get_numeric(x, secol)
  z <- get_numeric(x, zcol)
  p <- get_numeric(x, pcol)
  freq <- get_numeric(x, freqcol)

  cat("\nSAMPLE SIZE:\n", describe_numeric(n), "\n")
  cat("BETA:\n", describe_numeric(beta), "\n")
  cat("SE:\n", describe_numeric(se), "\n")
  cat("Z:\n", describe_numeric(z), "\n")
  cat("P:\n", describe_numeric(p), "\n")
  cat("FREQUENCY:\n", describe_numeric(freq), "\n")

  if (!is.null(n)) {
    cat("\nN quantiles:\n")
    if (any(is.finite(n))) {
      print(quantile(
        n[is.finite(n)],
        c(0, .01, .05, .5, .95, .99, 1)
      ))
    }
  }

  # Preview reported Z or derived beta/SE Z.
  zcheck <- z

  if (is.null(zcheck) &&
      !is.null(beta) && !is.null(se)) {
    zcheck <- beta / se
    cat("\nUsing beta/SE for Z diagnostic.\n")
  }

  valid_z <- if (is.null(zcheck)) {
    NA_integer_
  } else {
    sum(is.finite(zcheck))
  }

  mean_chisq <- if (is.null(zcheck) ||
                    !any(is.finite(zcheck))) {
    NA_real_
  } else {
    mean(zcheck[is.finite(zcheck)]^2)
  }

  duplicate_ids <- if (is.na(idcol)) {
    NA_integer_
  } else {
    sum(duplicated(x[[idcol]]))
  }

  missing_ids <- if (is.na(idcol)) {
    NA_integer_
  } else {
    sum(is.na(x[[idcol]]) |
        trimws(as.character(x[[idcol]])) == "")
  }

  invalid_p <- if (is.null(p)) {
    NA_integer_
  } else {
    sum(!is.na(p) &
        (!is.finite(p) | p < 0 | p > 1))
  }

  valid_alleles <- if (
    is.na(a1col) || is.na(a2col)
  ) {
    NA_integer_
  } else {
    a1 <- toupper(as.character(x[[a1col]]))
    a2 <- toupper(as.character(x[[a2col]]))

    sum(
      a1 %in% c("A", "C", "G", "T") &
      a2 %in% c("A", "C", "G", "T") &
      a1 != a2
    )
  }

  cat("\nFIRST TWO ROWS:\n")
  print(head(x[, seq_len(min(ncol(x), 12)),
               with = FALSE], 2))

  results[[i]] <- data.table(
    trait = trait,
    variants = nrow(x),
    columns = ncol(x),
    snp_id_column = idcol,
    duplicate_ids = duplicate_ids,
    missing_ids = missing_ids,
    valid_allele_pairs = valid_alleles,
    n_column = ncolname,
    median_N = if (!is.null(n) &&
                   any(is.finite(n))) {
      median(n[is.finite(n)])
    } else NA_real_,
    z_column = zcol,
    valid_Z = valid_z,
    mean_chisq = mean_chisq,
    p_column = pcol,
    invalid_p = invalid_p
  )

  columns_report[[i]] <- data.table(
    trait = trait,
    field = names(fields),
    source_column = unname(fields)
  )

  rm(x, n, beta, se, z, p, freq, zcheck)
  gc(verbose = FALSE)
}

report <- rbindlist(results, fill = TRUE)
mapping <- rbindlist(columns_report)

fwrite(
  report,
  file.path(outdir, "psychopathology_gwas_audit.tsv"),
  sep = "\t",
  na = "NA"
)

fwrite(
  mapping,
  file.path(outdir, "psychopathology_column_mapping.tsv"),
  sep = "\t",
  na = "NA"
)

cat("\n====================================\n")
cat("FINAL SIX-DATASET QC SUMMARY\n")
cat("====================================\n")

print(report)

cat("\nCOLUMN MAPPING:\n")
print(mapping)

cat("\nReports saved in:", outdir, "\n")
