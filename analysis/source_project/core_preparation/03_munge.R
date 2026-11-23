#!/usr/bin/env Rscript
# Requires reviewed per-trait input format, column names, N convention and reference.
if (!requireNamespace("GenomicSEM", quietly=TRUE)) stop("Install GenomicSEM first")
source("config/analysis_config.R")
if (any(!file.exists(raw_files))) stop("Missing GWAS input(s); inspect config/analysis_config.R")
if (!file.exists(hm3_file)) stop("Missing HapMap3 reference")
if (any(is.na(munge_N))) {
  message("NOTE: N=NA requires valid SNP-specific sample size fields in the corresponding raw GWAS.")
  message("CHECK each file's field names and effect-type convention before proceeding.")
}
dir.create("data/munged", recursive=TRUE, showWarnings=FALSE)
# GenomicSEM::munge writes output into the working directory.
old <- getwd(); on.exit(setwd(old),add=TRUE)
setwd("data/munged")
GenomicSEM::munge(
  files=normalizePath(file.path(old,raw_files), mustWork=TRUE),
  hm3=normalizePath(file.path(old,hm3_file), mustWork=TRUE),
  trait.names=traits,
  N=unname(munge_N[traits]),
  parallel=FALSE
)
message("Inspect ALL munging logs: allele convention, OR/LOG_OR/BETA, SNP-specific N, INFO/MAF filters.")
