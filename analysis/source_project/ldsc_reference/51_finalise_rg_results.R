library(data.table)

root <- path.expand("~/temperament_project")

outdir <- file.path(
  root, "results/ldsc/full_individual_highN"
)

x <- fread(file.path(
  outdir, "gradient_validation.tsv"
))

stopifnot(
  nrow(x) == 24L,
  all(is.finite(x$rg)),
  all(is.finite(x$SE_delta)),
  all(x$SE_delta > 0),
  max(x$max_gradient_error) < 1e-5,
  anyDuplicated(x[,.(temperament,outcome)]) == 0L
)

x[, `:=`(
  Z = rg / SE_delta,
  p = 2 * pnorm(-abs(rg / SE_delta)),
  CI_lower = rg - 1.96 * SE_delta,
  CI_upper = rg + 1.96 * SE_delta
)]

x[, FDR_24 := p.adjust(p, method="BH")]

x[, stage := fifelse(
  outcome %in% c("CINT","CADHD","CAGG"),
  "childhood",
  "adult"
)]

x[, significant_FDR := FDR_24 < 0.05]

setcolorder(x, c(
  "temperament","outcome","stage",
  "rg","SE_delta","CI_lower","CI_upper",
  "Z","p","FDR_24","significant_FDR",
  "max_gradient_error"
))

fwrite(
  x,
  file.path(outdir, "MASTER_24_genetic_correlations.tsv"),
  sep="\t"
)

cat("\n=== MASTER RESULTS ===\n")
print(x[, .(
  temperament, outcome, rg,
  SE=SE_delta, p, FDR_24
)], digits=5)

cat("\n=== FDR SIGNIFICANT ASSOCIATIONS ===\n")
print(
  x[significant_FDR == TRUE,
    .(temperament,outcome,rg,SE=SE_delta,FDR_24)],
  digits=5
)

cat("\n=== SUMMARY ===\n")
cat("Associations tested:", nrow(x), "\n")
cat("FDR < 0.05:",
    sum(x$significant_FDR), "\n")

cat("\nSaved master table to:\n")
cat(file.path(
  outdir, "MASTER_24_genetic_correlations.tsv"
), "\n")
