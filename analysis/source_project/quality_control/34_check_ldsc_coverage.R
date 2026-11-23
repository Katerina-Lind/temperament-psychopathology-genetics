library(data.table)

root <- path.expand("~/temperament_project")

ref <- file.path(root, "references/ldsc/w_hm3.snplist")
hm3 <- unique(as.character(fread(ref, header=FALSE)[[1]]))

paths <- c(
  EMO = "gwas/processed/sensitivity/common_highN90/EMO.sumstats.gz",
  ACT = "gwas/processed/sensitivity/common_highN90/ACT.sumstats.gz",
  SHY = "gwas/processed/sensitivity/common_highN90/SHY.sumstats.gz",
  SOC = "gwas/processed/sensitivity/common_highN90/SOC.sumstats.gz",
  CINT = "gwas/processed/munged/childhood/CINT.sumstats.gz",
  CADHD = "gwas/processed/munged/childhood/CADHD.sumstats.gz",
  CAGG = "gwas/processed/munged/childhood/CAGG.sumstats.gz",
  AINT = "gwas/processed/munged/adult/AINT.sumstats.gz",
  ASUD = "gwas/processed/munged/adult/ASUD.sumstats.gz",
  ARISK = "gwas/processed/munged/adult/ARISK.sumstats.gz"
)

ids <- list()
report <- list()

for (trait in names(paths)) {
  file <- file.path(root, paths[[trait]])
  stopifnot(file.exists(file))

  snps <- fread(file, select="SNP")$SNP
  ids[[trait]] <- unique(snps)

  report[[trait]] <- data.table(
    trait=trait,
    variants=length(snps),
    unique_snps=length(ids[[trait]]),
    hm3_overlap=sum(ids[[trait]] %chin% hm3),
    hm3_coverage_pct=round(
      100 * sum(ids[[trait]] %chin% hm3) / length(hm3),
      2
    )
  )
}

summary <- rbindlist(report)

cat("\n=== REFERENCE COVERAGE ===\n")
print(summary)

cat("\n=== COMMON SNPS ===\n")

common_temperament <- Reduce(
  intersect, ids[c("EMO","ACT","SHY","SOC")]
)

common_childhood <- Reduce(
  intersect, ids[c("CINT","CADHD","CAGG")]
)

common_adult <- Reduce(
  intersect, ids[c("AINT","ASUD","ARISK")]
)

common_all <- Reduce(intersect, ids)

cat("Temperament:", length(common_temperament), "\n")
cat("Childhood:", length(common_childhood), "\n")
cat("Adult:", length(common_adult), "\n")
cat("All ten:", length(common_all), "\n")

fwrite(
  summary,
  file.path(
    root,
    "gwas/processed/qc/ldsc_reference_coverage.tsv"
  ),
  sep="\t"
)
