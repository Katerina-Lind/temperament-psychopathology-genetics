suppressPackageStartupMessages(library(data.table))

old_root <- Sys.getenv("TEMPERAMENT_SOURCE_ROOT")
new_root <- Sys.getenv("TEMPERAMENT_ANALYSIS_ROOT")

if (!nzchar(old_root) || !nzchar(new_root)) {
  stop("Set TEMPERAMENT_SOURCE_ROOT and TEMPERAMENT_ANALYSIS_ROOT")
}

old_root <- normalizePath(old_root, mustWork=TRUE)
new_root <- normalizePath(new_root, mustWork=TRUE)

input_dir <- file.path(
  old_root, "gwas/temperament/age_specific/raw"
)

output_dir <- file.path(
  new_root, "data/ldsc_munged/temperament_age_specific"
)

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

hm3 <- unique(fread(
  file.path(new_root, "data/ldsc_reference/w_hm3.snplist"),
  header=FALSE
)[[1]])

traits <- c(
  EMO="Emotionality",
  ACT="Activity",
  SHY="Shyness"
)

summary <- list()

for (trait in names(traits)) {
  for (age in c("2Y","3Y")) {

    name <- paste0(trait, "_", age)

    input <- file.path(
      input_dir,
      paste0(
        "EASgwas_sumstats_",
        traits[[trait]], "_", age,
        "_Hollowell2025.txt"
      )
    )

    if (!file.exists(input)) {
      stop("Missing input: ", input)
    }

    cat("\nProcessing:", name, "\n")

    x <- fread(
      input,
      select=c(
        "SNP","A1","A2","Freq1",
        "Beta","SE","N"
      )
    )

    original <- nrow(x)

    x[, `:=`(
      SNP=as.character(SNP),
      A1=toupper(A1),
      A2=toupper(A2),
      Freq1=as.numeric(Freq1),
      Beta=as.numeric(Beta),
      SE=as.numeric(SE),
      N=as.numeric(N)
    )]

    # Exclude ALL records associated with duplicated SNP IDs.
    duplicates <- x[
      duplicated(SNP) | duplicated(SNP, fromLast=TRUE),
      unique(SNP)
    ]

    x <- x[!SNP %chin% duplicates]

    # Restrict to HM3 SNPs and valid biallelic nucleotides.
    x <- x[
      SNP %chin% hm3 &
      A1 %chin% c("A","C","G","T") &
      A2 %chin% c("A","C","G","T") &
      A1 != A2
    ]

    # Exclude strand-ambiguous variants.
    x <- x[!paste0(A1,A2) %chin% c(
      "AT","TA","CG","GC"
    )]

    # Frequency and summary-statistic QC.
    x[, MAF := pmin(Freq1, 1-Freq1)]

    x <- x[
      is.finite(MAF) & MAF > 0.01 &
      is.finite(Beta) &
      is.finite(SE) & SE > 0 &
      is.finite(N) & N > 0
    ]

    x[, Z := Beta / SE]

    x <- x[is.finite(Z)]

    output <- file.path(
      output_dir,
      paste0(name, ".sumstats.gz")
    )

    fwrite(
      x[, .(SNP,A1,A2,Z,N)],
      output,
      sep="\t",
      compress="gzip"
    )

    summary[[name]] <- data.table(
      trait=name,
      original_SNPs=original,
      munged_SNPs=nrow(x),
      median_N=median(x$N),
      min_N=min(x$N),
      max_N=max(x$N)
    )

    cat(
      "Original:", original,
      "Retained:", nrow(x),
      "Median N:", median(x$N), "\n"
    )

    rm(x)
    gc(verbose=FALSE)
  }
}

result <- rbindlist(summary)

fwrite(
  result,
  file.path(
    new_root,
    "results/ldsc/age_specific_h2/preparation_qc.tsv"
  ),
  sep="\t"
)

cat("\n=== SUMMARY ===\n")
print(result)
