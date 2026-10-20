# Cheap synthetic checks for allele orientation and proxy selection.
# Run from the repository root with: Rscript --vanilla tests/test_analysis_extensions.R

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table before running these checks.")
}
source("scripts/run_proxy_sensitivity.R")
source("scripts/run_colocalisation.R")

exposure <- data.table::data.table(
  taxon = "Example taxon", SNP = "rsTarget",
  beta.exposure = 0.20, se.exposure = 0.05,
  effect_allele.exposure = "A", other_allele.exposure = "G",
  eaf.exposure = 0.30, pval.exposure = 0.001
)
exact_empty <- data.table::data.table(SNP = character())
proxy_outcome <- data.table::data.table(
  SNP = "rsProxy", beta.outcome = 0.10, se.outcome = 0.03,
  effect_allele.outcome = "C", other_allele.outcome = "T",
  eaf.outcome = 0.25, pval.outcome = 0.01
)
proxy_map <- data.table::data.table(
  target_SNP = "rsTarget", proxy_SNP = "rsProxy",
  r2 = 0.90, signed_ld = -sqrt(0.90),
  target_allele = "A", proxy_allele = "C", maf_eligible = TRUE
)
selected <- choose_proxy_rows(exposure, exact_empty, proxy_outcome, proxy_map)
stopifnot(nrow(selected$selected) == 1L)
aligned <- align_proxy_to_target_allele(selected$selected, exposure)
stopifnot(isTRUE(all.equal(aligned$beta.outcome, -0.10)))

exact_present <- data.table::data.table(SNP = "rsTarget")
preferred <- choose_proxy_rows(exposure, exact_present, proxy_outcome, proxy_map)
stopifnot(nrow(preferred$selected) == 0L)
stopifnot(preferred$audit$status[1] == "exact_match_preferred")

trait1 <- data.table::data.table(
  SNP = c("rs1", "rs2", "rs3"), chr = "1", pos = c(100, 200, 300),
  beta = c(0.1, 0.2, 0.3), se = c(0.02, 0.03, 0.04),
  effect_allele = c("A", "C", "G"), other_allele = c("G", "T", "A"),
  eaf = c(0.2, 0.3, 0.4), genome_build = "GRCh37"
)
trait2 <- data.table::data.table(
  SNP = c("rs1", "rs2", "rs3"), chr = "1", pos = c(100, 200, 300),
  beta = c(0.5, 0.6, 0.7), se = c(0.05, 0.06, 0.07),
  effect_allele = c("G", "C", "T"), other_allele = c("A", "T", "C"),
  eaf = c(0.3, 0.25, 0.4), genome_build = "GRCh37"
)
harmonised <- harmonise_coloc_region(trait1, trait2, "test_region")
stopifnot(nrow(harmonised$data) == 2L)
stopifnot(isTRUE(all.equal(harmonised$data$beta2[harmonised$data$SNP == "rs1"], -0.5)))
stopifnot(isTRUE(all.equal(harmonised$data$eaf2[harmonised$data$SNP == "rs1"], 0.7)))
stopifnot(harmonised$audit[included == FALSE, SNP][1] == "rs3")

message("Synthetic proxy and colocalisation helper checks passed.")
