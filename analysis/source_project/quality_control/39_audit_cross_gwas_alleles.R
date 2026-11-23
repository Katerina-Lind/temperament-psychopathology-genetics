library(data.table)

root <- path.expand("~/temperament_project")

paths <- c(
  EMO   = "gwas/processed/sensitivity/common_highN90/EMO.sumstats.gz",
  ACT   = "gwas/processed/sensitivity/common_highN90/ACT.sumstats.gz",
  SHY   = "gwas/processed/sensitivity/common_highN90/SHY.sumstats.gz",
  SOC   = "gwas/processed/sensitivity/common_highN90/SOC.sumstats.gz",
  CINT  = "gwas/processed/munged/childhood/CINT.sumstats.gz",
  CADHD = "gwas/processed/munged/childhood/CADHD.sumstats.gz",
  CAGG  = "gwas/processed/munged/childhood/CAGG.sumstats.gz",
  AINT  = "gwas/processed/munged/adult/AINT.sumstats.gz",
  ASUD  = "gwas/processed/munged/adult/ASUD.sumstats.gz",
  ARISK = "gwas/processed/munged/adult/ARISK.sumstats.gz"
)

qcdir <- file.path(root, "gwas/processed/qc")
dir.create(qcdir, recursive=TRUE, showWarnings=FALSE)

cat("=== READING GWAS ALLELES ===\n")

data <- lapply(names(paths), function(trait) {
  file <- file.path(root, paths[[trait]])
  stopifnot(file.exists(file))

  x <- fread(file, select=c("SNP","A1","A2"))

  x[, `:=`(
    SNP=as.character(SNP),
    A1=toupper(A1),
    A2=toupper(A2)
  )]

  stopifnot(anyDuplicated(x$SNP) == 0L)

  cat(trait, ":", nrow(x), "variants\n")
  x
})

names(data) <- names(paths)

comp <- c(A="T", T="A", C="G", G="C")

compare_pair <- function(a, b) {

  x <- merge(
    data[[a]],
    data[[b]],
    by="SNP",
    suffixes=c(".x",".y"),
    sort=FALSE
  )

  exact <- x$A1.x == x$A1.y &
           x$A2.x == x$A2.y

  swapped <- x$A1.x == x$A2.y &
             x$A2.x == x$A1.y

  y1c <- unname(comp[x$A1.y])
  y2c <- unname(comp[x$A2.y])

  strand <- x$A1.x == y1c &
            x$A2.x == y2c

  strand_swapped <- x$A1.x == y2c &
                    x$A2.x == y1c

  category <- fifelse(
    exact, "exact",
    fifelse(
      swapped, "swapped",
      fifelse(
        strand, "strand",
        fifelse(
          strand_swapped,
          "strand_swapped",
          "incompatible"
        )
      )
    )
  )

  count <- function(label) sum(category == label)

  result <- data.table(
    pair=paste(a,b,sep="-"),
    overlap=nrow(x),
    exact=count("exact"),
    swapped=count("swapped"),
    strand=count("strand"),
    strand_swapped=count("strand_swapped"),
    incompatible=count("incompatible")
  )

  result[, incompatible_pct :=
    round(100*incompatible/overlap, 4)]

  result[, swapped_pct :=
    round(100*swapped/overlap, 2)]

  result
}

pairs <- combn(names(paths), 2, simplify=FALSE)

results <- rbindlist(lapply(pairs, function(p) {
  compare_pair(p[1], p[2])
}))

cat("\n=== CROSS-GWAS ALLELE AUDIT ===\n")

print(results)

cat("\n=== PAIRS WITH INCOMPATIBLE ALLELES ===\n")
print(results[incompatible > 0])

cat("\n=== SUMMARY ===\n")
cat("Pairs assessed:", nrow(results), "\n")
cat("Pairs with incompatible alleles:",
    sum(results$incompatible > 0), "\n")

output <- file.path(
  qcdir,
  "cross_gwas_allele_audit.tsv"
)

fwrite(results, output, sep="\t")

cat("\nSaved:", output, "\n")
