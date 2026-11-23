library(GenomicSEM)

root <- path.expand("~/temperament_psychopathology")
setwd(root)

traits <- c(
  "data/ldsc_munged/temperament/EMO.sumstats.gz",
  "data/ldsc_munged/temperament/ACT.sumstats.gz",
  "data/ldsc_munged/temperament/SHY.sumstats.gz",
  "data/ldsc_munged/worry.sumstats.gz",
  "data/ldsc_munged/depressed_affect.sumstats.gz",
  "data/ldsc_munged/antisocial_2022.sumstats.gz",
  "data/ldsc_munged/alcohol_problems.sumstats.gz"
)

trait.names <- c(
  "EMO", "ACT", "SHY",
  "WOR", "AFF", "ASB", "ALC"
)

stopifnot(all(file.exists(traits)))

ld <- "data/ldsc_reference/eur_w_ld_chr/"
wld <- ld

stopifnot(file.exists(
  file.path(ld, "1.l2.ldscore.gz")
))

# Treat all eight GWAS as quantitative/continuous
# or quantitative latent-factor phenotypes.
sample.prev <- rep(NA_real_, length(traits))
population.prev <- rep(NA_real_, length(traits))

cat("Starting joint LDSC analysis\n")
cat("Traits:", paste(trait.names, collapse=", "), "\n")

LDSC7 <- GenomicSEM::ldsc(
  traits = traits,
  sample.prev = sample.prev,
  population.prev = population.prev,
  ld = ld,
  wld = wld,
  trait.names = trait.names,
  n.blocks = 200,
  stand = TRUE,
  ldsc.log = "results/genomicsem/joint_ldsc_7traits"
)

saveRDS(
  LDSC7,
  "results/genomicsem/joint_ldsc_7traits.rds"
)

cat("\nAvailable output components:\n")
print(names(LDSC7))

cat("\nGenetic covariance matrix S:\n")
print(round(LDSC7$S, 4))

cat("\nEigenvalues of S:\n")
print(eigen(LDSC7$S, symmetric=TRUE)$values)

cat("\nSampling covariance matrix V dimensions:\n")
print(dim(LDSC7$V))

cat("\nCompleted successfully.\n")
