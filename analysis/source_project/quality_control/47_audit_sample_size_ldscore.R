library(data.table)

root <- path.expand("~/temperament_project")

traits <- c("EMO", "ACT", "SHY", "SOC")

refdir <- file.path(
  root, "references/ldsc/eur_w_ld_chr"
)

outdir <- file.path(
  root, "results/ldsc/sample_size_ldscore_qc"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cat("=== READING REFERENCE LD SCORES ===\n")

ldfiles <- file.path(
  refdir,
  paste0(1:22, ".l2.ldscore.gz")
)

stopifnot(all(file.exists(ldfiles)))

ld <- rbindlist(lapply(ldfiles, function(f) {
  fread(f, select=c("CHR", "SNP", "L2"))
}))

stopifnot(anyDuplicated(ld$SNP) == 0L)

cat("Reference variants:", nrow(ld), "\n")

summaries <- list()
deciles <- list()
chromosomes <- list()

for (trait in traits) {

  cat("\n=== ", trait, " ===\n", sep="")

  file <- file.path(
    root, "gwas/processed/munged/temperament",
    paste0(trait, ".sumstats.gz")
  )

  x <- fread(file, select=c("SNP", "N", "Z"))

  stopifnot(
    anyDuplicated(x$SNP) == 0L,
    all(is.finite(x$N)),
    all(x$N > 0),
    all(is.finite(x$Z))
  )

  maxN <- max(x$N)
  cutoff <- 0.90 * maxN

  x[, retained := N >= cutoff]

  n_before <- nrow(x)

  x <- merge(x, ld, by="SNP", sort=FALSE)

  cat("Maximum N:", maxN, "\n")
  cat("90% cutoff:", cutoff, "\n")
  cat("Original variants:", n_before, "\n")
  cat("Variants with LD scores:", nrow(x), "\n")

  # Distribution of LD scores by retention group.
  summary <- x[
    ,
    .(
      n_SNP=.N,
      median_N=median(as.numeric(N)),
      mean_N=mean(N),
      median_L2=median(L2),
      mean_L2=mean(L2),
      q25_L2=as.numeric(quantile(L2, .25)),
      q75_L2=as.numeric(quantile(L2, .75)),
      mean_chisq=mean(Z^2)
    ),
    by=retained
  ]

  summary[, trait := trait]

  setcolorder(
    summary,
    c("trait", "retained", "n_SNP",
      "median_N", "mean_N",
      "median_L2", "mean_L2",
      "q25_L2", "q75_L2", "mean_chisq")
  )

  print(summary)
  summaries[[trait]] <- summary

  # LD-score deciles within each temperament GWAS.
  breaks <- unique(as.numeric(
    quantile(x$L2, probs=seq(0,1,.1), na.rm=TRUE)
  ))

  stopifnot(length(breaks) >= 3L)

  x[, LD_decile := cut(
    L2,
    breaks=breaks,
    include.lowest=TRUE,
    labels=FALSE
  )]

  dec <- x[
    ,
    .(
      n_SNP=.N,
      retained_pct=100*mean(retained),
      median_L2=median(L2),
      median_N=median(as.numeric(N))
    ),
    by=LD_decile
  ][order(LD_decile)]

  dec[, trait := trait]

  deciles[[trait]] <- dec

  # Chromosome-specific retention.
  chr <- x[
    ,
    .(
      n_SNP=.N,
      retained_pct=100*mean(retained),
      median_N=median(as.numeric(N))
    ),
    by=CHR
  ][order(CHR)]

  chr[, trait := trait]
  chromosomes[[trait]] <- chr

  cat("\nCorrelation between N and LD score:\n")
  cat("Pearson:",
      cor(x$N, x$L2), "\n")
  cat("Spearman:",
      cor(x$N, x$L2, method="spearman"), "\n")

  cat("\nRetention across LD-score deciles:\n")
  print(dec)
}

summary <- rbindlist(summaries)
decile_table <- rbindlist(deciles)
chr_table <- rbindlist(chromosomes)

fwrite(
  summary,
  file.path(outdir, "retained_vs_excluded.tsv"),
  sep="\t"
)

fwrite(
  decile_table,
  file.path(outdir, "ldscore_decile_retention.tsv"),
  sep="\t"
)

fwrite(
  chr_table,
  file.path(outdir, "chromosome_retention.tsv"),
  sep="\t"
)

cat("\n=== DIAGNOSTIC COMPLETE ===\n")
cat("Results:", outdir, "\n")
