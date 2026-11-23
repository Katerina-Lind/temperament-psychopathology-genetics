library(data.table)

root <- path.expand("~/temperament_project")
outdir <- file.path(root,"results/qc/gwas_provenance")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

traits <- c(
  "EMO","ACT","SHY","SOC",
  "CINT","CADHD","CAGG",
  "AINT","ASUD","ARISK"
)

raw_files <- c(
  "temperament/raw/EASgwas_sumstats_Emotionality_Hollowell2025.txt",
  "temperament/raw/EASgwas_sumstats_Activity_Hollowell2025.txt",
  "temperament/raw/EASgwas_sumstats_Shyness_Hollowell2025.txt",
  "temperament/raw/EASgwas_sumstats_Sociability_Hollowell2025.txt",
  "childhood/internalising_GCST90054778/GCST90054778_buildGRCh37.tsv.gz",
  "childhood/ADHD_symptoms_GCST90568440/GCST90568440.tsv.gz",
  "childhood/aggression_Pappa2015_GCST002986/pappa 150810 EAGLE MA aggression_METAL RESULTS.txt",
  "adult/internalising_PGC2025/F4_Internalizing_2025.tsv.gz",
  "adult/substance_use_PGC2025/F5_SubstanceUse_2025.tsv.gz",
  "adult/risk_taking_Clifton2018/Risk_SumStats_Clifton_2018.csv.gz"
)

names(raw_files) <- traits

get_processed <- function(t) {
  if(t %in% traits[1:4]) {
    file.path(
      root,
      "results/ldsc/internalising_individual_highN/inputs",
      paste0(t,".sumstats.gz")
    )
  } else {
    file.path(
      root,
      "gwas/processed/sensitivity/strict_alleles",
      paste0(t,".sumstats.gz")
    )
  }
}

get_raw <- function(t) {

  f <- file.path(root,"gwas",raw_files[[t]])

  if(t %in% traits[1:4]) {
    x <- fread(f,select=c(
      "SNP","A1","A2","Beta","SE"
    ))
    setnames(x,c("Beta","SE"),c("beta","se"))

  } else if(t=="CINT") {
    x <- fread(f,select=c(
      "variant_id","effect_allele",
      "other_allele","z_score"
    ))
    setnames(
      x,
      c("variant_id","effect_allele",
        "other_allele","z_score"),
      c("SNP","A1","A2","source_Z")
    )

  } else if(t=="CADHD") {
    x <- fread(f,select=c(
      "rsid","effect_allele","other_allele",
      "beta","standard_error"
    ))
    setnames(
      x,
      c("rsid","effect_allele","other_allele",
        "standard_error"),
      c("SNP","A1","A2","se")
    )

  } else if(t=="CAGG") {
    x <- fread(f,select=c(
      "MarkerName","Allele1","Allele2","Zscore"
    ))
    setnames(
      x,
      c("MarkerName","Allele1","Allele2","Zscore"),
      c("SNP","A1","A2","source_Z")
    )

  } else if(t %in% c("AINT","ASUD")) {
    x <- fread(f,select=c(
      "SNP","A1","A2","BETA","SE"
    ))
    setnames(
      x,c("BETA","SE"),c("beta","se")
    )

  } else if(t=="ARISK") {
    x <- fread(f,select=c(
      "rsID","Effect_allele",
      "Other_allele","Beta","SE"
    ))
    setnames(
      x,
      c("rsID","Effect_allele",
        "Other_allele","Beta","SE"),
      c("SNP","A1","A2","beta","se")
    )
  }

  if(!"source_Z" %in% names(x)) {
    x[,source_Z := beta/se]
  }

  x[,`:=`(
    A1=toupper(A1),
    A2=toupper(A2)
  )]

  x[,.(SNP,A1,A2,source_Z)]
}

complement <- function(x) {
  chartr("ACGT","TGCA",x)
}

results <- list()
examples <- list()

for(k in seq_along(traits)) {

  t <- traits[k]

  cat("\n=== ",t," ===\n",sep="")

  p <- fread(
    get_processed(t),
    select=c("SNP","A1","A2","Z")
  )

  stopifnot(anyDuplicated(p$SNP)==0L)

  set.seed(5600+k)

  if(nrow(p)>10000L) {
    p <- p[sample.int(nrow(p),10000L)]
  }

  cat("Reading original GWAS...\n")

  r <- get_raw(t)

  r <- r[SNP %chin% p$SNP]

  # Protect against ambiguous source SNP identifiers.
  duplicate_ids <- unique(
    r$SNP[duplicated(r$SNP)]
  )

  r <- r[!SNP %chin% duplicate_ids]

  x <- merge(
    p,r,by="SNP",
    suffixes=c("_processed","_raw")
  )

  x[,`:=`(
    A1_processed=toupper(A1_processed),
    A2_processed=toupper(A2_processed)
  )]

  x[,orientation := fifelse(
    A1_processed==A1_raw &
      A2_processed==A2_raw,
    "direct",
    fifelse(
      A1_processed==A2_raw &
        A2_processed==A1_raw,
      "swapped",
      fifelse(
        A1_processed==complement(A1_raw) &
          A2_processed==complement(A2_raw),
        "complement_direct",
        fifelse(
          A1_processed==complement(A2_raw) &
            A2_processed==complement(A1_raw),
          "complement_swapped",
          "incompatible"
        )
      )
    )
  )]

  x[,expected_Z := fifelse(
    orientation %in% c(
      "direct","complement_direct"
    ),
    source_Z,
    fifelse(
      orientation %in% c(
        "swapped","complement_swapped"
      ),
      -source_Z,
      NA_real_
    )
  )]

  x[,abs_error := abs(Z-expected_Z)]

  # Avoid unstable sign comparisons around Z=0.
  x[,sign_mismatch := (
    is.finite(expected_Z) &
    is.finite(Z) &
    abs(expected_Z)>.05 &
    abs(Z)>.05 &
    sign(expected_Z)!=sign(Z)
  )]

  valid <- x[
    is.finite(expected_Z) & is.finite(Z)
  ]

  summary <- data.table(
    trait=t,
    sampled_processed=nrow(p),
    matched=nrow(x),
    duplicate_source_ids=length(duplicate_ids),
    direct=sum(x$orientation=="direct"),
    swapped=sum(x$orientation=="swapped"),
    complement_direct=sum(
      x$orientation=="complement_direct"
    ),
    complement_swapped=sum(
      x$orientation=="complement_swapped"
    ),
    incompatible=sum(
      x$orientation=="incompatible"
    ),
    sign_mismatches=sum(
      x$sign_mismatch,na.rm=TRUE
    ),
    correlation_Z=if(nrow(valid)>2)
      cor(valid$Z,valid$expected_Z)
      else NA_real_,
    median_abs_Z_error=if(nrow(valid)>0)
      median(valid$abs_error)
      else NA_real_,
    max_abs_Z_error=if(nrow(valid)>0)
      max(valid$abs_error)
      else NA_real_
  )

  results[[t]] <- summary

  examples[[t]] <- x[
    orientation=="incompatible" |
    sign_mismatch==TRUE,
    .(
      SNP,A1_processed,A2_processed,
      A1_raw,A2_raw,
      Z,source_Z,expected_Z,
      orientation,sign_mismatch
    )
  ][1:min(.N,30)]

  print(summary)
}

report <- rbindlist(results)
issues <- rbindlist(
  examples,
  idcol="trait",
  fill=TRUE
)

fwrite(
  report,
  file.path(outdir,"effect_direction_audit.tsv"),
  sep="\t"
)

fwrite(
  issues,
  file.path(outdir,"effect_direction_issues.tsv"),
  sep="\t"
)

cat("\n=== ALL TRAITS ===\n")
print(report,digits=6)

cat("\n=== POTENTIAL ISSUES ===\n")
print(issues)

cat("\n=== COMPLETE ===\n")
