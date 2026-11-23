library(data.table)

root <- path.expand("~/temperament_project")
dir <- file.path(root, "results/qc/gwas_provenance")

x <- fread(file.path(dir, "gwas_provenance_master.tsv"))

x[, `:=`(
  N_type = "Requires source-specific interpretation",
  verification = "PARTIAL"
)]

x[trait %in% c("EMO","ACT","SHY","SOC"),
  `:=`(
    N_type = "SNP-specific contributing sample size",
    verification = "SOURCE_MATCH"
  )]

x[trait == "CINT",
  `:=`(
    N_type = "SNP-specific effective N from repeated-measures meta-analysis",
    verification = "SOURCE_MATCH"
  )]

x[trait == "CADHD",
  `:=`(
    published_unique_N = 70953,
    N_type = "Constant N=73483; catalog sample count, verify against SSF metadata",
    verification = "PARTIAL"
  )]

x[trait == "CAGG",
  `:=`(
    N_type = "SNP-specific meta-analysis N",
    verification = "SOURCE_MATCH"
  )]

x[trait == "AINT",
  `:=`(
    N_type = "Implied N for multivariate genomic factor",
    verification = "SOURCE_MATCH"
  )]

x[trait == "ASUD",
  `:=`(
    N_type = "Implied N for multivariate genomic factor",
    verification = "SOURCE_MATCH"
  )]

x[trait == "ARISK",
  `:=`(
    N_type = "Fixed UK Biobank sample count",
    verification = "SOURCE_MATCH"
  )]

stopifnot(
  identical(nrow(x), 10L),
  all(x$genome_build_verified == FALSE),
  all(x$reference_allele_verified == FALSE)
)

fwrite(
  x,
  file.path(dir, "gwas_provenance_annotated.tsv"),
  sep="\t"
)

cat("=== SAMPLE SIZE PROVENANCE ===\n")
print(x[, .(
  trait,
  median_N,
  max_N,
  N_type,
  verification
)], width=120)

cat("\nSaved annotated provenance table.\n")
