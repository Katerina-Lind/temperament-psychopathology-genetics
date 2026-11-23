library(data.table)

root <- path.expand("~/temperament_project")

outdir <- file.path(root, "results/qc/gwas_provenance")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

traits <- c(
  "EMO","ACT","SHY","SOC",
  "CINT","CADHD","CAGG",
  "AINT","ASUD","ARISK"
)

temperament <- traits[1:4]

source <- data.table(
  trait = traits,
  phenotype = c(
    "Cross-age emotionality",
    "Cross-age activity",
    "Cross-age shyness",
    "Cross-age sociability",
    "Childhood/adolescent internalising symptoms",
    "Childhood ADHD symptom measures",
    "Childhood/adolescent aggressive behaviour",
    "Adult psychiatric Internalizing factor F4",
    "Adult psychiatric Substance Use factor F5",
    "Adult self-reported risk-taking"
  ),
  source = c(
    rep("Hollowell et al., 2026",4),
    "Jami et al., 2022",
    "van der Laan et al., 2025",
    "Pappa et al., 2015",
    "Grotzinger et al., 2026",
    "Grotzinger et al., 2026",
    "Clifton et al., 2018"
  ),
  published_unique_N = c(
    50426,50453,45656,43963,
    64641,70953,18988,
    NA,NA,436236
  ),
  published_effective_N = c(
    NA,NA,NA,NA,
    132260,NA,NA,
    NA,NA,NA
  ),
  source_url = c(
    rep("https://doi.org/10.1038/s41562-026-02486-5",4),
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC10859168/",
    "https://doi.org/10.1038/s41588-025-02295-y",
    "https://doi.org/10.1002/ajmg.b.32333",
    rep("https://doi.org/10.1038/s41586-025-09820-3",2),
    "https://doi.org/10.1038/s42003-018-0042-6"
  ),
  provenance_note = c(
    rep("European cross-age EAS meta-analysis",4),
    "Repeated assessments; N is an effective sample size",
    "Multiple raters/ages; N interpretation to verify",
    "Weak SNP heritability; exploratory",
    "Latent genomic factor; N interpretation unresolved",
    "Latent genomic factor; N interpretation unresolved",
    "UK Biobank BOLT-LMM risk-taking propensity"
  )
)

# Read estimates from our saved full LDSC model.
fit_path <- file.path(
  root,
  "results/ldsc/full_individual_highN",
  "full_ten_traits_ldsc.rds"
)

stopifnot(file.exists(fit_path))

fit <- readRDS(fit_path)
S <- unname(fit$S)
V <- fit$V

stopifnot(
  identical(dim(S), c(10L,10L)),
  identical(dim(V), c(55L,55L)),
  isSymmetric(S),
  isSymmetric(V)
)

# GenomicSEM V uses the column-major lower
# triangle of the genetic covariance matrix.
idx <- which(lower.tri(S, diag=TRUE), arr.ind=TRUE)

h2_rows <- which(idx[,1] == idx[,2])

heritability <- data.table(
  trait = traits,
  h2 = diag(S),
  h2_SE = sqrt(diag(V)[h2_rows])
)

heritability[, h2_Z := h2/h2_SE]

audit_one <- function(t) {

  if (t %in% temperament) {
    f <- file.path(
      root,
      "results/ldsc/internalising_individual_highN",
      "inputs",
      paste0(t, ".sumstats.gz")
    )
  } else {
    f <- file.path(
      root,
      "gwas/processed/sensitivity/strict_alleles",
      paste0(t, ".sumstats.gz")
    )
  }

  stopifnot(file.exists(f))

  hdr <- names(fread(f, nrows=0))

  required <- c("SNP","N","Z","A1","A2")

  stopifnot(all(required %in% hdr))

  x <- fread(f, select=required)

  stopifnot(
    anyDuplicated(x$SNP)==0L,
    all(is.finite(x$N)),
    all(x$N>0),
    all(is.finite(x$Z))
  )

  x[, `:=`(
    A1=toupper(A1),
    A2=toupper(A2)
  )]

  valid_alleles <- x$A1 %in% c("A","C","G","T") &
    x$A2 %in% c("A","C","G","T") &
    x$A1 != x$A2

  palindromic <- (
    x$A1=="A" & x$A2=="T"
  ) | (
    x$A1=="T" & x$A2=="A"
  ) | (
    x$A1=="C" & x$A2=="G"
  ) | (
    x$A1=="G" & x$A2=="C"
  )

  data.table(
    trait=t,
    input_file=f,
    n_SNP=nrow(x),
    n_unique_N=uniqueN(x$N),
    min_N=min(x$N),
    median_N=median(as.numeric(x$N)),
    mean_N=mean(x$N),
    max_N=max(x$N),
    mean_chisq=mean(x$Z^2),
    max_abs_Z=max(abs(x$Z)),
    n_invalid_alleles=sum(!valid_alleles),
    n_palindromic=sum(palindromic),
    fraction_palindromic=mean(palindromic),
    fraction_N_below_half_max=mean(x$N < .5*max(x$N)),
    genome_build_verified=FALSE,
    reference_allele_verified=FALSE
  )
}

cat("=== AUDITING TEN GWAS INPUTS ===\n")

qc <- rbindlist(lapply(traits, function(t) {
  cat("Reading:",t,"\n")
  audit_one(t)
}))

result <- merge(source,qc,by="trait")
result <- merge(result,heritability,by="trait")

result[, caution := fifelse(
  h2_Z < 4,
  "LOW_H2_SIGNAL",
  "REVIEW"
)]

result[trait %in% c("AINT","ASUD"),
       caution := "VERIFY_FACTOR_N"]

result[trait=="CINT",
       caution := "VERIFY_EFFECTIVE_N"]

result[, trait_order := match(trait, traits)]
setorder(result, trait_order)
result[, trait_order := NULL]

fwrite(
  result,
  file.path(outdir,"gwas_provenance_master.tsv"),
  sep="\t"
)

cat("\n=== SOURCE VS PROCESSED SAMPLE SIZE ===\n")
print(result[, .(
  trait,
  published_unique_N,
  published_effective_N,
  min_N,
  median_N,
  max_N,
  n_unique_N
)],digits=6)

cat("\n=== SNP AND ALLELE QC ===\n")
print(result[, .(
  trait,
  n_SNP,
  n_invalid_alleles,
  fraction_palindromic,
  mean_chisq
)],digits=5)

cat("\n=== HERITABILITY AND REVIEW FLAGS ===\n")
print(result[, .(
  trait,h2,h2_SE,h2_Z,caution
)],digits=5)

cat("\n=== COMPLETE ===\n")
cat(file.path(
  outdir,"gwas_provenance_master.tsv"
),"\n")
