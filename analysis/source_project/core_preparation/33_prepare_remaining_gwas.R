library(data.table)

root <- path.expand("~/temperament_project")
manifest <- fread(file.path(
  root, "config/temperament_psychopathology_datasets.tsv"
))

hm3 <- unique(fread(
  file.path(root, "references/ldsc/w_hm3.snplist"),
  header = FALSE
)[[1]])

outbase <- file.path(root, "gwas/processed/munged")
qcdir <- file.path(root, "gwas/processed/qc")

dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)

spec <- list(
  CINT = list(
    id="variant_id", a1="effect_allele",
    a2="other_allele", z="z_score",
    n="n", freq="minor_allele_frequency",
    fixedN=NA_real_, domain="childhood"
  ),
  CADHD = list(
    id="rsid", a1="effect_allele",
    a2="other_allele", beta="beta",
    se="standard_error",
    freq="effect_allele_frequency",
    fixedN=73483, domain="childhood"
  ),
  CAGG = list(
    id="MarkerName", a1="Allele1",
    a2="Allele2", z="Zscore",
    n="SAMPLESIZE", freq="Freq1",
    fixedN=NA_real_, domain="childhood"
  ),
  AINT = list(
    id="SNP", a1="A1", a2="A2",
    beta="BETA", se="SE", freq="MAF",
    fixedN=1637337, domain="adult"
  ),
  ASUD = list(
    id="SNP", a1="A1", a2="A2",
    beta="BETA", se="SE", freq="MAF",
    fixedN=313395, domain="adult"
  )
)

summary <- list()

for (trait in names(spec)) {

  s <- spec[[trait]]
  index <- match(trait, manifest[["trait"]])
  stopifnot(!is.na(index))

  input <- file.path(root, manifest[["file"]][index])

  cat("\n=== Processing", trait, "===\n")

  needed <- unique(unlist(
    s[c("id","a1","a2","z","beta","se","n","freq")]
  ))

  x <- fread(
    input,
    select=needed,
    na.strings=c("NA","#NA","NaN","")
  )

  original_n <- nrow(x)

  y <- data.table(
    SNP = as.character(x[[s$id]]),
    A1 = toupper(as.character(x[[s$a1]])),
    A2 = toupper(as.character(x[[s$a2]]))
  )

  if (!is.null(s$z)) {
    y[, Z := as.numeric(x[[s$z]])]
  } else {
    y[, Z := as.numeric(x[[s$beta]]) /
                 as.numeric(x[[s$se]])]
  }

  if (!is.null(s$n)) {
    y[, N := as.numeric(x[[s$n]])]
  } else {
    y[, N := s$fixedN]
  }

  f <- as.numeric(x[[s$freq]])
  maf <- pmin(f, 1-f)

  y[, MAF := maf]

  # Exclude every identifier occurring more than once.
  duplicated_ids <- y[
    duplicated(SNP) | duplicated(SNP, fromLast=TRUE),
    unique(SNP)
  ]

  y <- y[!SNP %chin% duplicated_ids]

  y <- y[
    SNP %chin% hm3 &
    A1 %chin% c("A","C","G","T") &
    A2 %chin% c("A","C","G","T") &
    A1 != A2 &
    !(paste0(A1,A2) %chin%
        c("AT","TA","CG","GC")) &
    is.finite(Z) &
    is.finite(N) & N > 0 &
    is.finite(MAF) & MAF >= 0.01
  ]

  stopifnot(anyDuplicated(y$SNP) == 0L)

  outputdir <- file.path(outbase, s$domain)
  dir.create(
    outputdir, recursive=TRUE,
    showWarnings=FALSE
  )

  output <- file.path(
    outputdir, paste0(trait, ".sumstats.gz")
  )

  if (file.exists(output)) {
    stop("Refusing to overwrite existing file: ", output)
  }

  # Standard LDSC column order
  y <- y[, .(SNP,N,Z,A1,A2)]

  fwrite(
    y, output,
    sep="\t", compress="gzip"
  )

  summary[[trait]] <- data.table(
    trait=trait,
    original_variants=original_n,
    retained_variants=nrow(y),
    retention_pct=round(
      100*nrow(y)/original_n, 2
    ),
    duplicate_id_groups_excluded=
      length(duplicated_ids),
    median_N=median(y$N),
    mean_chisq=mean(y$Z^2),
    output=output
  )

  print(summary[[trait]])

  rm(x,y,f,maf)
  gc(verbose=FALSE)
}

report <- rbindlist(summary)

fwrite(
  report,
  file.path(qcdir, "remaining_gwas_preparation.tsv"),
  sep="\t"
)

cat("\n=== FINAL PREPARATION SUMMARY ===\n")
print(report)

cat("\nAll outputs are provisional pending")
cat(" source and allele-reference verification.\n")
