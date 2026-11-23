library(data.table)

root <- path.expand("~/temperament_project")

paths <- c(
 EMO="gwas/processed/sensitivity/common_highN90/EMO.sumstats.gz",
 ACT="gwas/processed/sensitivity/common_highN90/ACT.sumstats.gz",
 SHY="gwas/processed/sensitivity/common_highN90/SHY.sumstats.gz",
 SOC="gwas/processed/sensitivity/common_highN90/SOC.sumstats.gz",
 CINT="gwas/processed/munged/childhood/CINT.sumstats.gz",
 CADHD="gwas/processed/munged/childhood/CADHD.sumstats.gz",
 CAGG="gwas/processed/munged/childhood/CAGG.sumstats.gz",
 AINT="gwas/processed/munged/adult/AINT.sumstats.gz",
 ASUD="gwas/processed/munged/adult/ASUD.sumstats.gz",
 ARISK="gwas/processed/munged/adult/ARISK.sumstats.gz"
)

outdir <- file.path(
 root, "gwas/processed/sensitivity/strict_alleles"
)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Build one allele-pair reference across all ten files.
reference <- data.table(
 SNP=character(), A1=character(), A2=character()
)

bad <- character()

for (trait in names(paths)) {
 x <- fread(file.path(root, paths[[trait]]),
            select=c("SNP","A1","A2"))

 x[, `:=`(
   A1=toupper(A1),
   A2=toupper(A2)
 )]

 stopifnot(anyDuplicated(x$SNP)==0L)

 if (nrow(reference)>0L) {
   m <- merge(
     x, reference, by="SNP",
     suffixes=c(".new",".ref")
   )

   exact <- m$A1.new==m$A1.ref &
            m$A2.new==m$A2.ref

   swapped <- m$A1.new==m$A2.ref &
              m$A2.new==m$A1.ref

   bad <- union(bad, m$SNP[!(exact|swapped)])
 }

 new <- x[!SNP %chin% reference$SNP]
 reference <- rbindlist(list(reference,new))
}

cat("\n=== DISCORDANT SNP IDS ===\n")
cat("Total:", length(bad), "\n")

fwrite(
 data.table(SNP=bad),
 file.path(outdir, "excluded_allele_discordant_snps.tsv"),
 sep="\t"
)

report <- list()

for (trait in names(paths)) {
 input <- file.path(root, paths[[trait]])
 output <- file.path(outdir, paste0(trait,".sumstats.gz"))

 x <- fread(input)
 before <- nrow(x)

 x <- x[!SNP %chin% bad]

 stopifnot(anyDuplicated(x$SNP)==0L)

 fwrite(x, output, sep="\t", compress="gzip")

 report[[trait]] <- data.table(
  trait=trait,
  original=before,
  excluded=before-nrow(x),
  retained=nrow(x)
 )
}

report <- rbindlist(report)

cat("\n=== STRICT ALLELE INPUTS ===\n")
print(report)

fwrite(
 report,
 file.path(outdir,"strict_allele_qc.tsv"),
 sep="\t"
)

cat("\nOutputs:",outdir,"\n")
