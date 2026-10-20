# Exact-versus-proxy CRP MR sensitivity analysis.
# Run from the repository root: Rscript --vanilla scripts/run_proxy_sensitivity.R
# This script reads local, pre-downloaded data only; it does not call a GWAS API.

source("scripts/analysis_extensions_common.R")

proxy_required_exposure <- c(
  "taxon", "SNP", "beta.exposure", "se.exposure",
  "effect_allele.exposure", "other_allele.exposure",
  "eaf.exposure", "pval.exposure"
)
proxy_required_outcome <- c(
  "SNP", "beta.outcome", "se.outcome",
  "effect_allele.outcome", "other_allele.outcome",
  "eaf.outcome", "pval.outcome"
)
proxy_required_map <- c(
  "target_SNP", "proxy_SNP", "r2", "signed_ld",
  "target_allele", "proxy_allele", "maf_eligible"
)

normalize_proxy_inputs <- function(exposure, exact, proxy, pairs) {
  data.table::setDT(exposure)
  data.table::setDT(exact)
  data.table::setDT(proxy)
  data.table::setDT(pairs)

  # Accept the historical exposure column name used by the main pipeline.
  if (!"taxon" %in% names(exposure) && "taxa" %in% names(exposure)) {
    data.table::setnames(exposure, "taxa", "taxon")
  }
  assert_columns(exposure, proxy_required_exposure, "exposure_instruments.csv")
  assert_columns(exact, proxy_required_outcome, "crp_exact.csv")
  assert_columns(proxy, proxy_required_outcome, "crp_proxy.csv")
  assert_columns(pairs, proxy_required_map, "proxy_pairs.csv")

  exposure[, `:=`(
    taxon = trimws(as.character(taxon)),
    SNP = trimws(as.character(SNP)),
    effect_allele.exposure = toupper(trimws(as.character(effect_allele.exposure))),
    other_allele.exposure = toupper(trimws(as.character(other_allele.exposure)))
  )]
  exact[, `:=`(
    SNP = trimws(as.character(SNP)),
    effect_allele.outcome = toupper(trimws(as.character(effect_allele.outcome))),
    other_allele.outcome = toupper(trimws(as.character(other_allele.outcome)))
  )]
  proxy[, `:=`(
    SNP = trimws(as.character(SNP)),
    effect_allele.outcome = toupper(trimws(as.character(effect_allele.outcome))),
    other_allele.outcome = toupper(trimws(as.character(other_allele.outcome)))
  )]
  pairs[, `:=`(
    target_SNP = trimws(as.character(target_SNP)),
    proxy_SNP = trimws(as.character(proxy_SNP)),
    target_allele = toupper(trimws(as.character(target_allele))),
    proxy_allele = toupper(trimws(as.character(proxy_allele)))
  )]

  numeric_columns <- list(
    list(data = exposure, columns = c("beta.exposure", "se.exposure", "eaf.exposure", "pval.exposure"), label = "exposure"),
    list(data = exact, columns = c("beta.outcome", "se.outcome", "eaf.outcome", "pval.outcome"), label = "exact outcome"),
    list(data = proxy, columns = c("beta.outcome", "se.outcome", "eaf.outcome", "pval.outcome"), label = "proxy outcome"),
    list(data = pairs, columns = c("r2", "signed_ld"), label = "proxy map")
  )
  for (item in numeric_columns) {
    for (column in item$columns) {
      item$data[, (column) := suppressWarnings(as.numeric(get(column)))]
    }
  }

  maf_flag <- toupper(trimws(as.character(pairs$maf_eligible)))
  allowed_flags <- c("TRUE", "FALSE", "T", "F", "1", "0", "YES", "NO")
  if (anyNA(maf_flag) || any(!maf_flag %in% allowed_flags)) {
    stop("proxy_pairs.csv maf_eligible must contain explicit TRUE/FALSE values.")
  }
  pairs[, maf_eligible := maf_flag %in% c("TRUE", "T", "1", "YES")]

  if (anyNA(exposure$taxon) || any(!nzchar(exposure$taxon)) || anyNA(exposure$SNP) ||
      anyNA(exact$SNP) || anyNA(proxy$SNP) || anyNA(pairs$target_SNP) ||
      anyNA(pairs$proxy_SNP) || anyNA(pairs$r2) || anyNA(pairs$signed_ld) ||
      anyNA(pairs$target_allele) || anyNA(pairs$proxy_allele)) {
    stop("Required IDs, LD values, or alleles contain missing values.")
  }

  assert_unique(exposure, c("taxon", "SNP"), "exposure instruments")
  assert_unique(exact, "SNP", "exact CRP outcome")
  assert_unique(proxy, "SNP", "proxy CRP outcome")
  assert_unique(pairs, c("target_SNP", "proxy_SNP"), "proxy map")

  if (any(!is.finite(exposure$beta.exposure)) || any(!is.finite(exposure$se.exposure)) || any(exposure$se.exposure <= 0)) {
    stop("Exposure beta/SE values must be finite and SE must be positive.")
  }
  for (out in list(exact, proxy)) {
    if (any(!is.finite(out$beta.outcome)) || any(!is.finite(out$se.outcome)) || any(out$se.outcome <= 0)) {
      stop("Outcome beta/SE values must be finite and SE must be positive.")
    }
  }
  if (any(!is.finite(pairs$r2)) || any(pairs$r2 < 0 | pairs$r2 > 1) ||
      any(!is.finite(pairs$signed_ld)) || any(abs(pairs$signed_ld) > 1)) {
    stop("Proxy r2 and signed_ld must be finite values within [0,1] and [-1,1].")
  }
  if (any(abs(pairs$r2 - pairs$signed_ld^2) > 0.05)) {
    stop("Some proxy rows have r2 inconsistent with signed_ld^2 by more than 0.05. Check LD source and allele coding.")
  }
  if (any(!exposure$effect_allele.exposure %in% c("A", "C", "G", "T")) ||
      any(!exposure$other_allele.exposure %in% c("A", "C", "G", "T"))) {
    stop("Exposure alleles must be A/C/G/T; ambiguous or non-SNP alleles need upstream resolution.")
  }
  if (any(!pairs$target_allele %in% c("A", "C", "G", "T")) ||
      any(!pairs$proxy_allele %in% c("A", "C", "G", "T"))) {
    stop("Proxy-map alleles must be A/C/G/T.")
  }
  invisible(list(exposure = exposure, exact = exact, proxy = proxy, pairs = pairs))
}

format_proxy_exposure <- function(exposure_rows) {
  arguments <- list(
    dat = as.data.frame(exposure_rows), type = "exposure",
    snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure", pval_col = "pval.exposure",
    phenotype_col = "taxon"
  )
  TwoSampleMR::format_data(arguments$dat, type = arguments$type,
    snp_col = arguments$snp_col, beta_col = arguments$beta_col,
    se_col = arguments$se_col, effect_allele_col = arguments$effect_allele_col,
    other_allele_col = arguments$other_allele_col, eaf_col = arguments$eaf_col,
    pval_col = arguments$pval_col, phenotype_col = arguments$phenotype_col)
}

format_proxy_outcome <- function(outcome_rows, label) {
  outcome_rows <- as.data.frame(outcome_rows)
  outcome_rows$outcome <- label
  TwoSampleMR::format_data(
    dat = outcome_rows, type = "outcome", snp_col = "SNP",
    beta_col = "beta.outcome", se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    eaf_col = "eaf.outcome", pval_col = "pval.outcome",
    phenotype_col = "outcome"
  )
}

choose_proxy_rows <- function(exposure_rows, exact_rows, proxy_rows, pair_rows, r2_min = 0.80) {
  targets <- data.table::data.table(target_SNP = unique(exposure_rows$SNP))
  audit <- merge(targets, pair_rows, by = "target_SNP", all.x = TRUE, sort = FALSE)
  audit[, exact_present := target_SNP %in% exact_rows$SNP]
  audit[, status := "candidate"]
  audit[is.na(proxy_SNP), status := "no_proxy_mapping"]
  audit[!is.na(proxy_SNP) & exact_present, status := "exact_match_preferred"]
  audit[!is.na(proxy_SNP) & !exact_present & !is.na(r2) & r2 < r2_min, status := "below_r2_threshold"]
  audit[!is.na(proxy_SNP) & !exact_present & !is.na(r2) & r2 >= r2_min & !is.na(maf_eligible) & !maf_eligible, status := "maf_ineligible"]

  eligible <- audit[status == "candidate"]
  if (nrow(eligible) == 0L) {
    audit[, proxy_reused_within_taxon := FALSE]
    return(list(selected = eligible, audit = audit))
  }
  eligible <- merge(eligible, proxy_rows, by.x = "proxy_SNP", by.y = "SNP", all.x = TRUE, sort = FALSE)
  has_outcome <- !is.na(eligible$beta.outcome) & !is.na(eligible$se.outcome)
  missing_keys <- eligible[!has_outcome, .(target_SNP, proxy_SNP)]
  for (i in seq_len(nrow(missing_keys))) {
    audit[target_SNP == missing_keys$target_SNP[i] & proxy_SNP == missing_keys$proxy_SNP[i], status := "proxy_outcome_missing"]
  }
  available <- eligible[has_outcome]
  if (nrow(available) == 0L) {
    audit[, proxy_reused_within_taxon := FALSE]
    return(list(selected = available, audit = audit))
  }
  data.table::setorder(available, target_SNP, -r2, proxy_SNP)
  available[, target_rank := seq_len(.N), by = target_SNP]
  selected <- available[target_rank == 1L]
  lower <- available[target_rank > 1L, .(target_SNP, proxy_SNP)]
  for (i in seq_len(nrow(selected))) {
    audit[target_SNP == selected$target_SNP[i] & proxy_SNP == selected$proxy_SNP[i], status := "selected"]
  }
  for (i in seq_len(nrow(lower))) {
    audit[target_SNP == lower$target_SNP[i] & proxy_SNP == lower$proxy_SNP[i], status := "lower_r2_for_target"]
  }
  selected[, proxy_reused_within_taxon := .N > 1L, by = proxy_SNP]
  audit[, proxy_reused_within_taxon := FALSE]
  reused <- selected[proxy_reused_within_taxon == TRUE, unique(proxy_SNP)]
  if (length(reused) > 0L) audit[proxy_SNP %in% reused, proxy_reused_within_taxon := TRUE]
  list(selected = selected, audit = audit)
}

align_proxy_to_target_allele <- function(selected, exposure_rows) {
  if (nrow(selected) == 0L) return(selected)
  data.table::setDT(selected)
  aligned <- merge(selected, exposure_rows, by.x = "target_SNP", by.y = "SNP", all = FALSE, sort = FALSE)

  target_sign <- ifelse(aligned$target_allele == aligned$effect_allele.exposure, 1,
    ifelse(aligned$target_allele == aligned$other_allele.exposure, -1, NA_real_))
  proxy_sign <- ifelse(aligned$proxy_allele == aligned$effect_allele.outcome, 1,
    ifelse(aligned$proxy_allele == aligned$other_allele.outcome, -1, NA_real_))
  if (anyNA(target_sign) || anyNA(proxy_sign)) {
    stop("A proxy-map allele did not match the target or proxy GWAS alleles. Check the signed LD allele definition.")
  }
  ld_sign <- sign(aligned$signed_ld)
  if (any(ld_sign == 0)) stop("signed_ld cannot be zero for an eligible proxy.")

  aligned[, `:=`(
    SNP = target_SNP,
    beta.outcome = beta.outcome * ld_sign * target_sign * proxy_sign,
    effect_allele.outcome = effect_allele.exposure,
    other_allele.outcome = other_allele.exposure,
    eaf.outcome = NA_real_
  )]
  aligned
}

estimate_proxy_branch <- function(harmonised, taxon, branch, n_exact, n_proxy) {
  data.table::setDT(harmonised)
  if (nrow(harmonised) == 0L || !"mr_keep" %in% names(harmonised) || !any(harmonised$mr_keep)) {
    primary <- data.table::data.table(
      taxon = taxon, branch = branch, status = "no_harmonised_instruments",
      method = NA_character_, b = NA_real_, se = NA_real_, pval = NA_real_,
      nsnp = 0L, n_exact = n_exact, n_proxy = n_proxy
    )
    return(list(primary = primary, estimators = data.table::data.table()))
  }
  kept <- harmonised[mr_keep == TRUE]
  if (anyDuplicated(kept$SNP)) stop("Duplicate target SNPs remain after proxy selection for ", taxon, ".")
  all_mr <- data.table::as.data.table(TwoSampleMR::mr(as.data.frame(kept)))
  if (nrow(all_mr) == 0L) {
    primary <- data.table::data.table(
      taxon = taxon, branch = branch, status = "no_mr_estimator",
      method = NA_character_, b = NA_real_, se = NA_real_, pval = NA_real_,
      nsnp = uniqueN(kept$SNP), n_exact = n_exact, n_proxy = n_proxy
    )
    return(list(primary = primary, estimators = data.table::data.table()))
  }
  all_mr[, `:=`(taxon = taxon, branch = branch)]
  method_order <- c("Inverse variance weighted", "Wald ratio", "Weighted median",
                    "MR Egger", "Weighted mode", "Simple mode")
  all_mr[, method_rank := match(as.character(method), method_order)]
  all_mr[is.na(method_rank), method_rank := 999L]
  data.table::setorder(all_mr, method_rank)
  selected <- all_mr[1L]
  primary <- data.table::data.table(
    taxon = taxon, branch = branch, status = "ok",
    method = as.character(selected$method), b = as.numeric(selected$b),
    se = as.numeric(selected$se), pval = as.numeric(selected$pval),
    nsnp = as.integer(selected$nsnp), n_exact = as.integer(n_exact),
    n_proxy = as.integer(n_proxy)
  )
  all_mr[, method_rank := NULL]
  list(primary = primary, estimators = all_mr)
}

analyse_one_taxon_proxy <- function(taxon, exposure, exact, proxy, pairs, r2_min = 0.80) {
  exposure_rows <- exposure[taxon == ..taxon]
  exposure_formatted <- format_proxy_exposure(exposure_rows)

  exact_rows <- exact[SNP %in% exposure_rows$SNP]
  exact_h <- data.table::data.table()
  if (nrow(exact_rows) > 0L) {
    exact_formatted <- format_proxy_outcome(exact_rows, "CRP exact-primary")
    exact_h <- data.table::as.data.table(TwoSampleMR::harmonise_data(
      exposure_dat = exposure_formatted,
      outcome_dat = exact_formatted,
      action = 2
    ))
  }
  exact_n <- if (nrow(exact_h) > 0L) sum(exact_h$mr_keep, na.rm = TRUE) else 0L
  exact_result <- estimate_proxy_branch(exact_h, taxon, "exact_primary", exact_n, 0L)

  selected <- choose_proxy_rows(exposure_rows, exact, proxy, pairs, r2_min)
  audit <- selected$audit
  proxy_h <- data.table::data.table()
  if (nrow(selected$selected) > 0L) {
    proxy_aligned <- align_proxy_to_target_allele(selected$selected, exposure_rows)
    proxy_for_format <- proxy_aligned[, .(
      SNP, beta.outcome, se.outcome,
      effect_allele.outcome, other_allele.outcome,
      eaf.outcome, pval.outcome
    )]
    proxy_formatted <- format_proxy_outcome(proxy_for_format, "CRP proxy-sensitivity")
    # Allele phase has already been resolved from signed LD; action=1 preserves
    # the constructed target-allele orientation, including known palindromic pairs.
    proxy_h <- data.table::as.data.table(TwoSampleMR::harmonise_data(
      exposure_dat = exposure_formatted,
      outcome_dat = proxy_formatted,
      action = 1
    ))
  }
  if (nrow(exact_h) > 0L) {
    exact_for_proxy <- data.table::copy(exact_h)
    exact_for_proxy[, `:=`(id.outcome = "CRP_proxy_sensitivity", outcome = "CRP proxy-sensitivity")]
  } else {
    exact_for_proxy <- data.table::data.table()
  }
  if (nrow(proxy_h) > 0L) {
    proxy_h[, `:=`(id.outcome = "CRP_proxy_sensitivity", outcome = "CRP proxy-sensitivity")]
  }
  combined_proxy <- data.table::rbindlist(list(exact_for_proxy, proxy_h), fill = TRUE)
  proxy_n <- if (nrow(proxy_h) > 0L) sum(proxy_h$mr_keep, na.rm = TRUE) else 0L
  proxy_result <- estimate_proxy_branch(
    combined_proxy, taxon, "proxy_sensitivity", exact_n, proxy_n
  )

  primary <- data.table::rbindlist(list(exact_result$primary, proxy_result$primary), fill = TRUE)
  estimators <- data.table::rbindlist(list(
    data.table::copy(exact_result$estimators),
    data.table::copy(proxy_result$estimators)
  ), fill = TRUE)
  list(primary = primary, estimators = estimators, audit = audit)
}

run_proxy_sensitivity <- function(project_root = getwd()) {
  assert_packages(c("data.table", "TwoSampleMR"))
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(project_root, "pipeline.R"))) {
    stop("Run this command from the repository root, or set project_root to that folder.")
  }
  input_dir <- Sys.getenv("AIXBIO_PROXY_INPUT_DIR", unset = "data/private/proxy_sensitivity")
  output_root <- Sys.getenv("AIXBIO_PROXY_OUTPUT_DIR", unset = "results/proxy_sensitivity/runs")
  input_dir <- file.path(project_root, input_dir)
  output_root <- file.path(project_root, output_root)
  input_files <- file.path(input_dir, c(
    "exposure_instruments.csv", "crp_exact.csv", "crp_proxy.csv", "proxy_pairs.csv"
  ))
  missing <- input_files[!file.exists(input_files)]
  if (length(missing) > 0L) {
    stop("Proxy analysis inputs are missing. See scripts/ANALYSIS_EXTENSIONS.md. Missing: ",
         paste(missing, collapse = ", "))
  }
  code_files <- file.path(project_root, c(
    "scripts/run_proxy_sensitivity.R", "scripts/analysis_extensions_common.R"
  ))
  config <- c("r2_min=0.80", "maf_filter=explicit_proxy_pairs.maf_eligible", "exact_action=2", "proxy_action=1")
  run_id <- fingerprint_run(input_files, code_files, config, c("data.table", "TwoSampleMR"))
  run_dir <- file.path(output_root, paste0("proxy_", run_id))
  checkpoint_dir <- file.path(run_dir, "checkpoints")
  failure_dir <- file.path(run_dir, "failures")
  final_dir <- file.path(run_dir, "final")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)

  if (dir.exists(final_dir)) {
    validate_completed_manifest(final_dir, run_id)
    message("Validated completed proxy run; reusing outputs at ", final_dir)
    return(invisible(final_dir))
  }
  unlock <- acquire_run_lock(file.path(run_dir, ".lock"))
  on.exit(unlock(), add = TRUE)

  exposure <- data.table::fread(input_files[1], na.strings = c("", "NA"))
  exact <- data.table::fread(input_files[2], na.strings = c("", "NA"))
  proxy <- data.table::fread(input_files[3], na.strings = c("", "NA"))
  pairs <- data.table::fread(input_files[4], na.strings = c("", "NA"))
  inputs <- normalize_proxy_inputs(exposure, exact, proxy, pairs)
  exposure <- inputs$exposure; exact <- inputs$exact
  proxy <- inputs$proxy; pairs <- inputs$pairs
  taxa <- sort(unique(exposure$taxon))
  if (length(taxa) == 0L) stop("No taxa were found in exposure_instruments.csv.")

  unit_results <- vector("list", length(taxa))
  for (i in seq_along(taxa)) {
    unit_path <- file.path(checkpoint_dir, sprintf("taxon_%04d.rds", i))
    if (file.exists(unit_path)) {
      cached <- tryCatch(readRDS(unit_path), error = function(e) {
        stop("Unreadable checkpoint preserved at ", unit_path, ": ", conditionMessage(e))
      })
      if (!identical(cached$taxon, taxa[i]) || is.null(cached$primary) || is.null(cached$audit)) {
        stop("Checkpoint does not match current taxon ordering; preserved at ", unit_path)
      }
      unit_results[[i]] <- cached
      message("Reusing taxon checkpoint ", i, "/", length(taxa), ": ", taxa[i])
      next
    }
    message("Analysing taxon ", i, "/", length(taxa), ": ", taxa[i])
    result <- tryCatch(
      analyse_one_taxon_proxy(taxa[i], exposure, exact, proxy, pairs, r2_min = 0.80),
      error = function(e) {
        failure_path <- file.path(failure_dir, sprintf("taxon_%04d_%s.txt", i, format(Sys.time(), "%Y%m%dT%H%M%S")))
        atomic_write_lines(c(paste0("taxon=", taxa[i]), paste0("error=", conditionMessage(e))), failure_path)
        stop("Taxon failed and was not checkpointed: ", taxa[i], " — ", conditionMessage(e))
      }
    )
    result$taxon <- taxa[i]
    result$unit_id <- i
    atomic_save_rds(result, unit_path)
    unit_results[[i]] <- result
  }

  primary <- data.table::rbindlist(lapply(unit_results, `[[`, "primary"), fill = TRUE)
  estimators <- data.table::rbindlist(lapply(unit_results, `[[`, "estimators"), fill = TRUE)
  if (ncol(estimators) == 0L) {
    estimators <- data.table::data.table(
      taxon = character(), branch = character(), id.exposure = character(),
      id.outcome = character(), exposure = character(), outcome = character(),
      method = character(), nsnp = integer(), b = numeric(), se = numeric(), pval = numeric()
    )
  }
  audit <- data.table::rbindlist(lapply(unit_results, `[[`, "audit"), fill = TRUE)
  primary[, fdr := NA_real_]
  for (branch_name in unique(primary$branch)) {
    index <- which(primary$branch == branch_name & is.finite(primary$pval))
    if (length(index) > 0L) primary[index, fdr := stats::p.adjust(pval, method = "BH")]
  }
  branch_summary <- primary[, .(
    taxa_total = .N,
    taxa_with_estimate = sum(status == "ok"),
    nominal_p_lt_0_05 = sum(is.finite(pval) & pval < 0.05),
    fdr_lt_0_05 = sum(is.finite(fdr) & fdr < 0.05),
    median_instruments = if (any(nsnp > 0L)) stats::median(nsnp[nsnp > 0L]) else NA_real_
  ), by = branch]
  exact_primary <- primary[branch == "exact_primary", .(taxon, b_exact = b)]
  proxy_primary <- primary[branch == "proxy_sensitivity", .(taxon, b_proxy = b)]
  paired <- merge(exact_primary, proxy_primary, by = "taxon", all = FALSE)
  paired <- paired[is.finite(b_exact) & is.finite(b_proxy)]
  direction_concordance <- if (nrow(paired) > 0L) mean(sign(paired$b_exact) == sign(paired$b_proxy)) else NA_real_
  branch_summary[, paired_taxa_with_estimates := nrow(paired)]
  branch_summary[, paired_direction_concordance := direction_concordance]

  stage_dir <- tempfile(pattern = ".publish_proxy_", tmpdir = run_dir)
  dir.create(stage_dir, recursive = FALSE, showWarnings = FALSE)
  if (!dir.exists(stage_dir)) stop("Could not create final staging directory.")
  stage_files <- c("proxy_sensitivity_results.csv", "all_estimators.csv", "proxy_pair_audit.csv", "branch_summary.csv")
  atomic_write_csv(primary, file.path(stage_dir, stage_files[1]))
  atomic_write_csv(estimators, file.path(stage_dir, stage_files[2]))
  atomic_write_csv(audit, file.path(stage_dir, stage_files[3]))
  atomic_write_csv(branch_summary, file.path(stage_dir, stage_files[4]))
  output_hashes <- tools::md5sum(file.path(stage_dir, stage_files))
  manifest <- data.table::data.table(
    run_id = run_id, file = stage_files, md5 = unname(output_hashes),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    config = paste(config, collapse = ";"),
    TwoSampleMR = as.character(utils::packageVersion("TwoSampleMR")),
    data.table = as.character(utils::packageVersion("data.table")),
    R_version = as.character(getRversion())
  )
  atomic_write_csv(manifest, file.path(stage_dir, "COMPLETE_MANIFEST.csv"))
  publish_final_directory(stage_dir, final_dir)
  message("Proxy sensitivity analysis complete. Outputs: ", final_dir)
  invisible(final_dir)
}

proxy_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(proxy_script_arg) == 1L && grepl("run_proxy_sensitivity\\.R$", proxy_script_arg)) {
  run_proxy_sensitivity()
}
