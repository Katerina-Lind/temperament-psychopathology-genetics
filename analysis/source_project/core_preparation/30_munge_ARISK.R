library(data.table)

root <- path.expand("~/temperament_project")

input <- file.path(
  root,
  "gwas/processed/qc/ARISK_unambiguous_hm3.tsv.gz"
)

outdir <- file.path(
  root, "gwas/processed/munged/adult"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

output <- file.path(
  outdir, "ARISK.sumstats.gz"
)

stopifnot(file.exists(input))

x <- fread(input)

# Reported GWAS sample size: Clifton et al. 2018
N <- 436236L

x[, Z := BETA / SE]

x <- x[
  is.finite(Z) &
  is.finite(BETA) &
  is.finite(SE) &
  SE > 0 &
  A1 %chin% c("A", "C", "G", "T") &
  A2 %chin% c("A", "C", "G", "T") &
  A1 != A2
]

# Create LDSC-compatible summary statistics.
y <- x[, .(
  SNP,
  N = N,
  Z,
  A1,
  A2
)]

stopifnot(
  anyDuplicated(y$SNP) == 0,
  all(is.finite(y$Z)),
  all(y$N == 436236L)
)

fwrite(
  y,
  output,
  sep="\t",
  compress="gzip"
)

cat("\n=== ARISK LDSC INPUT ===\n")
cat("Variants:", nrow(y), "\n")
cat("Sample size:", N, "\n")
cat("Mean Chi-square:", mean(y$Z^2), "\n")
cat("Duplicate SNPs:", anyDuplicated(y$SNP), "\n")
cat("Output:", output, "\n")
