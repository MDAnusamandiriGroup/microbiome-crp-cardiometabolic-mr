#!/usr/bin/env Rscript

# =============================================================================
# PAPER 1 — CMPB HYBRID FINAL PIPELINE
# Microbiome two-sample MR + replication-aware explainable evidence ranking
# =============================================================================
#
# This self-contained script merges the user's previously executed CRP,
# benchmarking, and T2D/SBP generalisability scripts with a strengthened
# post-MR computational layer for Computer Methods and Programs in Biomedicine.
#
# PRESERVED FROM THE ORIGINAL PIPELINE
# - the existing preselected microbiome IV file and outcome-file paths;
# - TwoSampleMR harmonisation, MR estimators, Cochran Q and Egger intercept;
# - CRP primary/replication, T2D primary/replication, SBP primary/replication,
#   and optional within-family SBP sensitivity;
# - volcano, evidence heatmap, forest, SHAP and cross-dataset figures;
# - top-decile + direction-concordant consensus shortlists.
#
# STRENGTHENED IN THIS HYBRID VERSION
# - corrected rank-normalised evidence grade;
# - explicit missingness indicators with native XGBoost NA handling;
# - pairwise XGBoost ranking with repeated 5-fold CV used only for evaluation;
# - final candidate ordering from a full-data fitted ranking model;
# - p-value, |beta|, composite, Borda and full Pareto-front baselines;
# - external replication Precision@K and NDCG@K;
# - full-model ablations;
# - 300 true XGBoost bootstrap refits with ranking and SHAP stability;
# - separate fixed-score taxon-bootstrap comparison across all six methods;
# - complete manifests, parameter audits and ready-to-paste Methods text.
#
# IMPORTANT INTERPRETATION
# The target evidence grade and model predictors share MR-derived inputs.
# Therefore, XGBoost is an auditable nonlinear reweighting of evidence, not an
# independent causal-truth learner. External replication is the main evaluation.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(xgboost)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================

BASE_DIR <- normalizePath(Sys.getenv("AIXBIO_DATA_DIR", unset = file.path("data", "private")), mustWork = FALSE)
IV_PATH  <- file.path(BASE_DIR, "qin_IV_adaptive_combined.csv")
OUT_ROOT <- normalizePath(Sys.getenv("AIXBIO_OUTPUT_DIR", unset = "results"), mustWork = FALSE)

# Preserve the preselected/clumped instrument set used in the original study.
# The thresholds below are audited, but are NOT applied unless explicitly set
# to TRUE. This prevents silent changes to the manuscript's MR estimands.
APPLY_IV_P_FILTER <- FALSE
IV_P_THRESHOLD    <- 5e-6
APPLY_F_FILTER    <- FALSE
MIN_F_STATISTIC   <- 10
HARMONISE_ACTION  <- 2

TOP_K              <- 20L
CONSENSUS_TOP_PROP <- 0.10
CV_FOLDS            <- 5L
CV_REPEATS          <- 10L
MODEL_BOOTSTRAP_REFITS <- 300L
METHOD_BOOTSTRAPS      <- 300L
MASTER_SEED         <- 20260802L

# Rank-normalised weak-supervision target.
W_STAT       <- 0.50
W_SNR        <- 0.30
W_INSTRUMENT <- 0.20
LAMBDA_HET   <- 0.10
LAMBDA_EGGER <- 0.10

DETECTED_CORES <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(DETECTED_CORES) || DETECTED_CORES < 1L) DETECTED_CORES <- 2L
XGB_NROUNDS_MAX <- 2000L
XGB_EARLY_STOP  <- 50L
XGB_PARAMS <- list(
  objective        = "rank:pairwise",
  eval_metric      = paste0("ndcg@", TOP_K),
  eta              = 0.03,
  max_depth        = 3L,
  min_child_weight = 2,
  subsample        = 0.80,
  colsample_bytree = 0.80,
  lambda           = 1,
  alpha            = 0,
  nthread          = max(1L, as.integer(DETECTED_CORES) - 1L),
  seed             = MASTER_SEED
)

OUTCOMES <- list(
  list(
    trait = "CRP", role = "primary", required = TRUE,
    tag = "CRP_primary_ebi-a-GCST90025959",
    candidates = c(
      file.path(BASE_DIR, "Paper1_ebi-a-GCST90025959", "CRP_ebi-a-GCST90025959_selected_snps.csv"),
      file.path(BASE_DIR, "CRP_ebi-a-GCST90025959_selected_snps.csv")
    )
  ),
  list(
    trait = "CRP", role = "replication", required = TRUE,
    tag = "CRP_replication_ieu-b-35",
    candidates = c(file.path(BASE_DIR, "ieu-b-35_selected_snps.csv"))
  ),
  list(
    trait = "T2D", role = "primary", required = TRUE,
    tag = "T2D_primary_ebi-a-GCST90018926",
    candidates = c(
      file.path(BASE_DIR, "Paper3_T2D__20251226_170937", "T2D_main_ebi-a-GCST90018926", "ebi-a-GCST90018926_selected_snps.csv"),
      file.path(BASE_DIR, "ebi-a-GCST90018926_selected_snps.csv")
    )
  ),
  list(
    trait = "T2D", role = "replication", required = TRUE,
    tag = "T2D_replication_finn-b-T2D",
    candidates = c(
      file.path(BASE_DIR, "Paper3_T2D__20251226_170937", "T2D_sens_finn-b-T2D", "finn-b-T2D_selected_snps.csv"),
      file.path(BASE_DIR, "finn-b-T2D_selected_snps.csv")
    )
  ),
  list(
    trait = "SBP", role = "primary", required = TRUE,
    tag = "SBP_primary_ieu-b-38",
    candidates = c(
      file.path(BASE_DIR, "ClusterB_LipidsBP", "SBP_ICBP_ieu-b-38", "ieu-b-38_selected_snps.csv"),
      file.path(BASE_DIR, "ieu-b-38_selected_snps.csv")
    )
  ),
  list(
    trait = "SBP", role = "replication", required = TRUE,
    tag = "SBP_replication_ukb-b-20175",
    candidates = c(
      file.path(BASE_DIR, "SBP_UKB_ukb-b-20175", "ukb-b-20175_selected_snps.csv"),
      file.path(BASE_DIR, "ukb-b-20175_selected_snps.csv")
    )
  ),
  list(
    trait = "SBP", role = "within_family", required = FALSE,
    tag = "SBP_within_family_ieu-b-4817",
    candidates = c(
      file.path(BASE_DIR, "SBP_WithinFamily_ieu-b-4817", "ieu-b-4817_selected_snps.csv"),
      file.path(BASE_DIR, "ieu-b-4817_selected_snps.csv")
    )
  )
)

RUN_ABLATIONS_FOR <- c("CRP_primary_ebi-a-GCST90025959")
RUN_MODEL_BOOTSTRAP_FOR <- c("CRP_primary_ebi-a-GCST90025959")
RUN_METHOD_BOOTSTRAP_FOR <- c(
  "CRP_primary_ebi-a-GCST90025959",
  "T2D_primary_ebi-a-GCST90018926",
  "SBP_primary_ieu-b-38"
)
MAKE_FIGURES <- TRUE
MAKE_DESIGN_FIGURES <- TRUE
FOREST_TOP_N <- 20L

# =============================================================================
# 2. GENERAL HELPERS
# =============================================================================

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[1])) y else x
clean_snp <- function(x) trimws(as.character(x))
clean_allele <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N", ".")] <- NA_character_
  x
}
pick1 <- function(candidates, nms) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0L) NA_character_ else hit[1]
}
first_existing <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) NA_character_ else normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}
safe_tag <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

rank01 <- function(x, higher_is_better = TRUE, missing_value = 0.5) {
  out <- rep(missing_value, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  z <- x[ok]
  if (!higher_is_better) z <- -z
  if (length(unique(z)) == 1L) out[ok] <- 0.5
  else out[ok] <- (rank(z, ties.method = "average") - 1) / (length(z) - 1)
  out
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

jaccard <- function(a, b) {
  a <- unique(as.character(a)); b <- unique(as.character(b))
  den <- length(union(a, b))
  if (den == 0L) return(NA_real_)
  length(intersect(a, b)) / den
}

precision_at_k <- function(relevance, order_index, k = TOP_K) {
  k <- min(as.integer(k), length(order_index))
  if (k < 1L) return(NA_real_)
  mean(as.integer(relevance[order_index[seq_len(k)]] > 0), na.rm = TRUE)
}

ndcg_at_k <- function(relevance, order_index, k = TOP_K) {
  k <- min(as.integer(k), length(order_index))
  if (k < 1L) return(NA_real_)
  rel <- as.numeric(relevance[order_index[seq_len(k)]])
  dcg <- sum((2^rel - 1) / log2(seq_len(k) + 1), na.rm = TRUE)
  ideal <- head(sort(as.numeric(relevance), decreasing = TRUE, na.last = TRUE), k)
  idcg <- sum((2^ideal - 1) / log2(seq_len(k) + 1), na.rm = TRUE)
  if (!is.finite(idcg) || idcg <= 0) return(NA_real_)
  dcg / idcg
}

pairwise_median_spearman <- function(mat) {
  if (is.null(dim(mat)) || nrow(mat) < 2L) return(NA_real_)
  cm <- suppressWarnings(cor(t(mat), method = "spearman", use = "pairwise.complete.obs"))
  vals <- cm[upper.tri(cm)]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) NA_real_ else median(vals)
}

safe_format_outcome <- function(dat, snp, beta, se, ea, oa, p = NA_character_,
                                eaf = NA_character_, outcome_name = "outcome") {
  dat[[snp]] <- clean_snp(dat[[snp]])
  dat[[ea]]  <- clean_allele(dat[[ea]])
  dat[[oa]]  <- clean_allele(dat[[oa]])
  dat$outcome <- outcome_name
  args <- list(
    dat = dat, type = "outcome", snp_col = snp, beta_col = beta,
    se_col = se, effect_allele_col = ea, other_allele_col = oa,
    phenotype_col = "outcome"
  )
  if (!is.na(p)) args$pval_col <- p
  if (!is.na(eaf)) args$eaf_col <- eaf
  do.call(TwoSampleMR::format_data, args)
}

make_stratified_fold_plan <- function(y, folds, repeats, seed) {
  n <- length(y)
  plan <- matrix(NA_integer_, nrow = n, ncol = repeats)
  strata <- cut(rank(y, ties.method = "average"),
                breaks = unique(quantile(seq_len(n), probs = seq(0, 1, length.out = 11), type = 1)),
                include.lowest = TRUE, labels = FALSE)
  if (length(strata) != n || anyNA(strata)) strata <- rep(1L, n)
  for (r in seq_len(repeats)) {
    set.seed(seed + r)
    fold_r <- integer(n)
    for (s in sort(unique(strata))) {
      idx <- which(strata == s)
      fold_r[idx] <- sample(rep(seq_len(folds), length.out = length(idx)))
    }
    plan[, r] <- fold_r
  }
  plan
}

make_rank_dmatrix <- function(X, y = NULL) {
  d <- xgb.DMatrix(data = X, label = y, missing = NA_real_)
  setinfo(d, "group", as.integer(nrow(X)))
  d
}

extract_best_iteration <- function(model, fallback) {
  bi <- model$best_iteration
  if (is.null(bi) || !is.finite(bi) || bi < 1L) as.integer(fallback) else as.integer(bi)
}

compute_evidence_components <- function(dt) {
  x <- copy(as.data.table(dt))
  x[, `:=`(
    pen_het = as.integer(!is.na(het_Q_p) & het_Q_p < 0.05),
    pen_egger = as.integer(!is.na(egger_p) & egger_p < 0.05),
    miss_het_Q = as.integer(is.na(het_Q)),
    miss_het_Q_p = as.integer(is.na(het_Q_p)),
    miss_egger_int = as.integer(is.na(egger_int)),
    miss_egger_p = as.integer(is.na(egger_p))
  )]
  x[, `:=`(
    stat_component = rank01(neglog10_p),
    snr_component = rank01(signal_to_noise),
    instrument_component = rank01(log1p(nsnp))
  )]
  x[, evidence_grade_raw :=
      W_STAT * stat_component +
      W_SNR * snr_component +
      W_INSTRUMENT * instrument_component -
      LAMBDA_HET * pen_het -
      LAMBDA_EGGER * pen_egger]
  x[, evidence_grade := rank01(evidence_grade_raw)]
  # XGBoost's default exponential NDCG gain requires integer labels <=31.
  x[, relevance_label := as.integer(pmin(31L, floor(32 * evidence_grade)))]
  x
}

# =============================================================================
# 3. TRANSPARENT BASELINES
# =============================================================================

pareto_front_levels <- function(objective_matrix) {
  M <- as.matrix(objective_matrix)
  M[!is.finite(M)] <- 0.5
  n <- nrow(M)
  dominates <- vector("list", n)
  domination_count <- integer(n)
  first_front <- integer(0)

  # NSGA-II-style nondominated sorting: O(n^2) dominance checks.
  for (i in seq_len(n)) {
    Si <- integer(0)
    ni <- 0L
    for (j in seq_len(n)) {
      if (i == j) next
      i_dom_j <- all(M[i, ] >= M[j, ]) && any(M[i, ] > M[j, ])
      j_dom_i <- all(M[j, ] >= M[i, ]) && any(M[j, ] > M[i, ])
      if (i_dom_j) Si <- c(Si, j)
      else if (j_dom_i) ni <- ni + 1L
    }
    dominates[[i]] <- Si
    domination_count[i] <- ni
    if (ni == 0L) first_front <- c(first_front, i)
  }

  level <- integer(n)
  front <- unique(first_front)
  current <- 1L
  while (length(front) > 0L) {
    level[front] <- current
    next_front <- integer(0)
    for (i in front) {
      for (j in dominates[[i]]) {
        domination_count[j] <- domination_count[j] - 1L
        if (domination_count[j] == 0L) next_front <- c(next_front, j)
      }
    }
    front <- unique(next_front)
    current <- current + 1L
  }
  if (any(level == 0L)) stop("Pareto sorting failed to assign all observations.")
  level
}

add_baseline_scores_and_ranks <- function(dt) {
  x <- copy(as.data.table(dt))
  x[, score_pval := neglog10_p]
  x[, score_abeta := abs_b]
  x[, score_compE := evidence_grade]
  x[, qc_risk := pen_het + pen_egger]
  x[, qc_coverage := 1 - 0.5 * (miss_het_Q_p + miss_egger_p)]
  x[, score_borda := rowMeans(cbind(
    rank01(stat_component),
    rank01(snr_component),
    rank01(instrument_component),
    rank01(qc_risk, higher_is_better = FALSE),
    rank01(qc_coverage)
  ))]
  objectives <- cbind(
    stat = stat_component,
    snr = snr_component,
    instrument = instrument_component,
    qc_low_risk = rank01(qc_risk, higher_is_better = FALSE),
    qc_coverage = rank01(qc_coverage)
  )
  x[, pareto_level := pareto_front_levels(objectives)]
  x[, score_pareto := -pareto_level + 1e-3 * evidence_grade]
  score_cols <- c(
    xgb = "score_xgb", pval = "score_pval", abeta = "score_abeta",
    compE = "score_compE", borda = "score_borda", pareto = "score_pareto"
  )
  for (nm in names(score_cols)) {
    sc <- score_cols[[nm]]
    x[, (paste0("rank_", nm)) := frank(-get(sc), ties.method = "average")]
  }
  x
}

# =============================================================================
# 4. XGBOOST RANKING, CV, ABLATION AND BOOTSTRAP
# =============================================================================

run_repeated_cv_evaluation <- function(X, y_label, y_continuous, taxa, fold_plan,
                                       params, nrounds_max, early_stop,
                                       seed_offset = 0L) {
  n <- nrow(X)
  repeats <- ncol(fold_plan)
  folds <- max(fold_plan)
  oof_percentile <- matrix(NA_real_, nrow = n, ncol = repeats)
  fold_stats <- list()
  counter <- 1L
  for (r in seq_len(repeats)) {
    for (f in seq_len(folds)) {
      valid_idx <- which(fold_plan[, r] == f)
      train_idx <- setdiff(seq_len(n), valid_idx)
      if (length(valid_idx) < 5L || length(train_idx) < 20L) next
      dtrain <- make_rank_dmatrix(X[train_idx, , drop = FALSE], y_label[train_idx])
      dvalid <- make_rank_dmatrix(X[valid_idx, , drop = FALSE], y_label[valid_idx])
      p <- params
      p$seed <- as.integer(params$seed + seed_offset + 1000L * r + f)
      model <- xgb.train(
        params = p, data = dtrain, nrounds = nrounds_max,
        watchlist = list(train = dtrain, valid = dvalid),
        early_stopping_rounds = early_stop, maximize = TRUE, verbose = 0L
      )
      pred <- predict(model, dvalid)
      # Fold-normalised percentiles are used only for comparable OOF diagnostics.
      oof_percentile[valid_idx, r] <- rank01(pred)
      ord <- order(pred, decreasing = TRUE)
      fold_stats[[counter]] <- data.table(
        repeat = r, fold = f,
        n_train = length(train_idx), n_valid = length(valid_idx),
        best_iteration = extract_best_iteration(model, nrounds_max),
        target_spearman = safe_cor(pred, y_continuous[valid_idx], "spearman"),
        target_ndcg = ndcg_at_k(y_label[valid_idx], ord, min(TOP_K, length(valid_idx)))
      )
      counter <- counter + 1L
    }
  }
  if (anyNA(oof_percentile)) stop("Repeated CV left missing OOF percentile predictions.")
  repeat_rank <- apply(oof_percentile, 2, function(z) frank(-z, ties.method = "average"))
  mean_pct <- rowMeans(oof_percentile)
  mean_rank <- frank(-mean_pct, ties.method = "average")
  repeat_stability <- rbindlist(lapply(seq_len(repeats), function(r) {
    data.table(
      repeat = r,
      rank_spearman_vs_mean = safe_cor(repeat_rank[, r], mean_rank, "spearman")
    )
  }))
  list(
    oof_percentile = oof_percentile,
    oof_percentile_mean = mean_pct,
    oof_percentile_sd = apply(oof_percentile, 1, sd),
    oof_rank_mean = mean_rank,
    oof_rank_sd = apply(repeat_rank, 1, sd),
    fold_stats = rbindlist(fold_stats, fill = TRUE),
    repeat_stability = repeat_stability
  )
}

fit_final_ranker <- function(X, y_label, params, nrounds, seed) {
  dtrain <- make_rank_dmatrix(X, y_label)
  p <- params
  p$seed <- as.integer(seed)
  xgb.train(params = p, data = dtrain,
            nrounds = as.integer(max(1L, nrounds)), verbose = 0L)
}

run_model_bootstrap_refits <- function(X, y_label, taxa, reference_score,
                                       params, nrounds, n_bootstrap, seed) {
  n <- nrow(X)
  reference_top <- taxa[order(reference_score, decreasing = TRUE)[seq_len(min(TOP_K, n))]]
  stability <- vector("list", n_bootstrap)
  shap_rows <- vector("list", n_bootstrap)
  dall <- xgb.DMatrix(X, missing = NA_real_)
  for (b in seq_len(n_bootstrap)) {
    set.seed(seed + b)
    idx <- sample.int(n, n, replace = TRUE)
    dboot <- make_rank_dmatrix(X[idx, , drop = FALSE], y_label[idx])
    p <- params
    p$seed <- as.integer(seed + 10000L + b)
    model <- xgb.train(params = p, data = dboot, nrounds = nrounds, verbose = 0L)
    pred <- predict(model, dall)
    top_b <- taxa[order(pred, decreasing = TRUE)[seq_len(min(TOP_K, n))]]
    stability[[b]] <- data.table(
      bootstrap = b,
      topK_jaccard = jaccard(top_b, reference_top),
      rank_spearman = safe_cor(pred, reference_score, "spearman")
    )
    sh <- as.matrix(predict(model, dall, predcontrib = TRUE))
    bias <- which(colnames(sh) %in% c("BIAS", "bias", "Bias"))
    if (length(bias) == 1L) sh <- sh[, -bias, drop = FALSE]
    imp <- colMeans(abs(sh), na.rm = TRUE)
    shap_rows[[b]] <- data.table(
      bootstrap = b, feature = names(imp), mean_abs_shap = as.numeric(imp)
    )
    if (b %% 25L == 0L) message("  model bootstrap refits: ", b, "/", n_bootstrap)
  }
  stability_dt <- rbindlist(stability)
  shap_long <- rbindlist(shap_rows)
  shap_summary <- shap_long[, .(
    mean_abs_shap = mean(mean_abs_shap, na.rm = TRUE),
    sd_abs_shap = sd(mean_abs_shap, na.rm = TRUE),
    median_abs_shap = median(mean_abs_shap, na.rm = TRUE),
    q025 = quantile(mean_abs_shap, 0.025, na.rm = TRUE),
    q975 = quantile(mean_abs_shap, 0.975, na.rm = TRUE)
  ), by = feature][order(-mean_abs_shap)]
  shap_wide <- dcast(shap_long, bootstrap ~ feature, value.var = "mean_abs_shap")
  shap_matrix <- as.matrix(shap_wide[, -"bootstrap"])
  list(
    stability = stability_dt,
    shap_long = shap_long,
    shap_summary = shap_summary,
    shap_pairwise_spearman = pairwise_median_spearman(shap_matrix)
  )
}

run_method_taxon_bootstrap <- function(feature_dt, n_bootstrap, seed) {
  full <- add_baseline_scores_and_ranks(compute_evidence_components(feature_dt))
  methods <- c("xgb", "pval", "abeta", "compE", "borda", "pareto")
  full_top <- lapply(methods, function(m) {
    head(full[order(get(paste0("rank_", m))), exposure], TOP_K)
  })
  names(full_top) <- methods
  out <- vector("list", n_bootstrap * length(methods))
  counter <- 1L
  for (b in seq_len(n_bootstrap)) {
    set.seed(seed + b)
    idx <- sample.int(nrow(feature_dt), nrow(feature_dt), replace = TRUE)
    boot <- copy(feature_dt[idx])
    boot <- compute_evidence_components(boot)
    boot <- add_baseline_scores_and_ranks(boot)
    score_cols <- paste0("score_", methods)
    unique_boot <- boot[, lapply(.SD, mean, na.rm = TRUE), by = exposure, .SDcols = score_cols]
    for (m in methods) {
      sc <- paste0("score_", m)
      top_b <- head(unique_boot[order(-get(sc)), exposure], TOP_K)
      out[[counter]] <- data.table(
        bootstrap = b, method = m,
        topK_jaccard = jaccard(top_b, full_top[[m]])
      )
      counter <- counter + 1L
    }
  }
  long <- rbindlist(out)
  summary <- long[, .(
    jaccard_mean = mean(topK_jaccard, na.rm = TRUE),
    jaccard_sd = sd(topK_jaccard, na.rm = TRUE),
    jaccard_median = median(topK_jaccard, na.rm = TRUE),
    q025 = quantile(topK_jaccard, 0.025, na.rm = TRUE),
    q975 = quantile(topK_jaccard, 0.975, na.rm = TRUE)
  ), by = method][order(-jaccard_mean)]
  list(long = long, summary = summary)
}

within_dataset_metrics <- function(ranking_dt) {
  dt <- copy(as.data.table(ranking_dt))
  methods <- c("pval", "abeta", "compE", "borda", "pareto")
  n <- nrow(dt)
  top_n <- max(1L, ceiling(CONSENSUS_TOP_PROP * n))
  top_xgb_k <- dt[order(rank_xgb)][seq_len(min(TOP_K, n)), exposure]
  top_xgb_dec <- dt[order(rank_xgb)][seq_len(min(top_n, n)), exposure]
  rbindlist(lapply(methods, function(m) {
    rc <- paste0("rank_", m)
    top_k <- dt[order(get(rc))][seq_len(min(TOP_K, n)), exposure]
    top_dec <- dt[order(get(rc))][seq_len(min(top_n, n)), exposure]
    data.table(
      baseline = m,
      spearman_rho = safe_cor(dt$rank_xgb, dt[[rc]], "spearman"),
      kendall_tau = safe_cor(dt$rank_xgb, dt[[rc]], "kendall"),
      topK_overlap = length(intersect(top_xgb_k, top_k)),
      topK_jaccard = jaccard(top_xgb_k, top_k),
      topdecile_overlap = length(intersect(top_xgb_dec, top_dec)),
      topdecile_jaccard = jaccard(top_xgb_dec, top_dec),
      n_taxa = n, topdecile_n = top_n
    )
  }))
}

make_model_diagnostic_plots <- function(ranked_dt, outdir, tag) {
  dt <- as.data.table(ranked_dt)
  dt <- dt[is.finite(evidence_grade) & is.finite(score_xgb)]
  if (nrow(dt) < 10L) return(invisible(NULL))
  rho <- safe_cor(dt$evidence_grade, dt$score_xgb, "spearman")
  p1 <- ggplot(dt, aes(x = evidence_grade, y = score_xgb)) +
    geom_point(alpha = 0.7) +
    labs(
      title = paste0("Ranking-model diagnostic — ", tag),
      subtitle = paste0("Spearman rho with prespecified evidence grade = ", round(rho, 3)),
      x = "Prespecified evidence grade", y = "Full-data XGBoost ranking score"
    ) + theme_bw(base_size = 12)
  ggsave(file.path(outdir, "Fig_ModelDiagnostic_evidence_grade_vs_XGBoost.png"),
         p1, width = 8, height = 6, dpi = 320)
  p2 <- ggplot(dt, aes(x = score_xgb)) + geom_histogram(bins = 40) +
    labs(title = paste0("Distribution of XGBoost ranking scores — ", tag),
         x = "XGBoost ranking score", y = "Taxa") + theme_bw(base_size = 12)
  ggsave(file.path(outdir, "Fig_ModelDiagnostic_XGBoost_score_distribution.png"),
         p2, width = 8, height = 5.5, dpi = 320)
}


# =============================================================================
# 5. FIGURE HELPERS RETAINED FROM THE USER'S PROVEN PIPELINE
# =============================================================================

pick_primary_mr <- function(mr_dt) {
  # use slice_head (not slice) to avoid masking issues
  mr_dt %>%
    dplyr::mutate(
      method = as.character(method),
      exposure = as.character(exposure),
      nsnp = as.integer(nsnp)
    ) %>%
    dplyr::group_by(exposure) %>%
    dplyr::arrange(
      dplyr::desc(method == "Inverse variance weighted"),
      dplyr::desc(method == "Wald ratio")
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
}

clip_q <- function(x, q = c(0.05, 0.95)) {
  if (all(is.na(x))) return(x)
  lo <- quantile(x, q[1], na.rm = TRUE)
  hi <- quantile(x, q[2], na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

make_volcano <- function(mr_res_dt, outdir, tag) {
  mr_tbl <- as_tibble(mr_res_dt)
  
  primary <- pick_primary_mr(mr_tbl) %>%
    dplyr::mutate(
      neglog10p = -log10(pmax(pval, 1e-300)),
      fdr = p.adjust(pval, method = "BH"),
      sig = dplyr::case_when(
        fdr < 0.05 ~ "FDR<0.05",
        fdr < 0.10 ~ "FDR<0.10",
        pval < 0.05 ~ "P<0.05",
        TRUE ~ "NS"
      ),
      nsnp_cat = ifelse(nsnp >= 2, ">=2 SNP", "1 SNP")
    )
  
  label_dt <- primary %>%
    dplyr::arrange(fdr, pval) %>%
    dplyr::slice_head(n = 10)
  
  p <- ggplot(primary, aes(x = b, y = neglog10p)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(shape = nsnp_cat, colour = sig), alpha = 0.85, size = 2.1) +
    geom_text(data = label_dt, aes(label = exposure), size = 3, vjust = -0.7, check_overlap = TRUE) +
    labs(
      title = paste0("MR screen (primary method per taxon): ", tag),
      x = "MR effect estimate (b)",
      y = expression(-log[10](p))
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(outdir, "Fig_MR_volcano_primary.png"), p, width = 9, height = 6, dpi = 320)
  ggsave(file.path(outdir, "Fig_MR_volcano_primary.pdf"), p, width = 9, height = 6)
}

make_heatmap_top20 <- function(top20_dt, outdir, tag) {
  top20 <- as_tibble(top20_dt)
  
  need_cols <- c("exposure","b","se","pval","fdr","nsnp","het_Q_p","egger_p","score_pred")
  for (cc in need_cols) if (!cc %in% names(top20)) top20[[cc]] <- NA
  
  hm <- top20 %>%
    dplyr::mutate(
      abs_b = abs(b),
      b_c = clip_q(b),
      se_c = clip_q(se),
      abs_b_c = clip_q(abs_b),
      nsnp_c = clip_q(nsnp),
      het_Q_p_c = clip_q(het_Q_p),
      egger_p_c = clip_q(egger_p),
      logp = -log10(pmax(pval, 1e-300)),
      logfdr = -log10(pmax(fdr, 1e-300)),
      logp_c = clip_q(logp),
      logfdr_c = clip_q(logfdr)
    ) %>%
    dplyr::select(
      exposure,
      `b` = b_c,
      `SE` = se_c,
      `|b|` = abs_b_c,
      `nsnp` = nsnp_c,
      `-log10(p)` = logp_c,
      `-log10(FDR)` = logfdr_c,
      `het_Q_p` = het_Q_p_c,
      `egger_p` = egger_p_c
    ) %>%
    tidyr::pivot_longer(-exposure, names_to = "metric", values_to = "value") %>%
    dplyr::group_by(metric) %>%
    dplyr::mutate(value_z = as.numeric(scale(value))) %>%
    dplyr::ungroup()
  
  taxa_order <- if ("rank_xgb" %in% names(top20)) top20 %>% dplyr::arrange(rank_xgb) %>% dplyr::pull(exposure) else top20 %>% dplyr::arrange(fdr, pval) %>% dplyr::pull(exposure)
  hm$exposure <- factor(hm$exposure, levels = rev(taxa_order))
  
  p <- ggplot(hm, aes(x = metric, y = exposure, fill = value_z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    labs(
      title = paste0("Evidence heatmap: top 20 XGBoost-prioritised taxa (nsnp ≥ 2) — ", tag),
      x = NULL, y = NULL, fill = "z-score"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
  
  ggsave(file.path(outdir, "Fig_Top20_evidence_heatmap.png"), p, width = 10, height = 7.5, dpi = 320)
  ggsave(file.path(outdir, "Fig_Top20_evidence_heatmap.pdf"), p, width = 10, height = 7.5)
}

make_forest_methods <- function(mr_all_dt, taxa_vec, outdir, tag) {
  keep_methods <- c("Inverse variance weighted", "Wald ratio", "Weighted median", "MR Egger")
  
  dt <- as.data.table(mr_all_dt)
  dt <- dt[exposure %in% taxa_vec & method %in% keep_methods]
  if (nrow(dt) == 0) return(invisible(NULL))
  
  taxa_vec <- taxa_vec[taxa_vec %in% unique(dt$exposure)]
  dt[, exposure := factor(exposure, levels = rev(taxa_vec))]
  dt[, method := factor(as.character(method), levels = keep_methods)]
  dt[, ci_lo := b - 1.96 * se]
  dt[, ci_hi := b + 1.96 * se]
  
  dodge <- position_dodge(width = 0.7)
  
  p <- ggplot(dt, aes(x = exposure, y = b, colour = method, shape = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, position = dodge, alpha = 0.85) +
    geom_point(position = dodge, size = 2.2, alpha = 0.9) +
    coord_flip() +
    labs(
      title = paste0("Forest plot (multi-method MR): shortlist taxa — ", tag),
      x = NULL,
      y = "MR effect estimate (b) with 95% CI"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  h <- max(6, 0.28 * length(taxa_vec) + 3)
  ggsave(file.path(outdir, "Fig_Forest_shortlist_multimethod.png"), p, width = 10, height = h, dpi = 320)
  ggsave(file.path(outdir, "Fig_Forest_shortlist_multimethod.pdf"), p, width = 10, height = h)
}


make_crossdataset_plots <- function(cmp_dt, out_root, spearman_val = NA_real_, prefix = "") {
  cmp <- as_tibble(cmp_dt)
  title_suffix <- if (is.finite(spearman_val)) paste0(" (Spearman = ", round(spearman_val, 3), ")") else ""
  
  p1 <- ggplot(cmp, aes(x = b_primary, y = b_rep)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.4) +
    geom_point(aes(colour = dir_concordant), alpha = 0.75, size = 2.0) +
    labs(
      title = paste0("Cross-dataset effect concordance (b): primary vs replication", title_suffix),
      x = "b (primary)",
      y = "b (replication)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(out_root, paste0(prefix, "Fig_CrossDataset_b_scatter.png")), p1, width = 7.5, height = 6, dpi = 320)
  ggsave(file.path(out_root, paste0(prefix, "Fig_CrossDataset_b_scatter.pdf")), p1, width = 7.5, height = 6)
  
  p2 <- ggplot(cmp, aes(x = rank_primary, y = rank_rep)) +
    geom_point(aes(colour = dir_concordant), alpha = 0.75, size = 2.0) +
    scale_x_reverse() + scale_y_reverse() +
    labs(
      title = paste0("Cross-dataset rank agreement: primary vs replication", title_suffix),
      x = "Rank (primary; smaller = higher priority)",
      y = "Rank (replication; smaller = higher priority)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(out_root, paste0(prefix, "Fig_CrossDataset_rank_scatter.png")), p2, width = 7.5, height = 6, dpi = 320)
  ggsave(file.path(out_root, paste0(prefix, "Fig_CrossDataset_rank_scatter.pdf")), p2, width = 7.5, height = 6)
  
  p3 <- ggplot(cmp, aes(x = abs_rank_diff)) +
    geom_histogram(bins = 40, alpha = 0.85) +
    labs(
      title = "Distribution of absolute rank differences (primary vs replication)",
      x = "|Rank(primary) - Rank(replication)|",
      y = "Count"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(out_root, paste0(prefix, "Fig_CrossDataset_abs_rank_diff_hist.png")), p3, width = 7.5, height = 5.5, dpi = 320)
  ggsave(file.path(out_root, paste0(prefix, "Fig_CrossDataset_abs_rank_diff_hist.pdf")), p3, width = 7.5, height = 5.5)
}

make_replication_focus_panels <- function(cmp_dt, top_primary_vec, shortlist_vec, out_root, label_n = 12, prefix = "") {
  cmp <- as.data.table(cmp_dt)
  
  cmp[, group := fifelse(exposure %in% shortlist_vec, "Consensus (top decile + concordant)",
                         fifelse(exposure %in% top_primary_vec, "Top decile primary only", "Other"))]
  cmp[, group := factor(group, levels = c("Consensus (top decile + concordant)", "Top decile primary only", "Other"))]
  
  lab <- head(cmp[exposure %in% shortlist_vec][order(rank_primary)], label_n)
  
  pA <- ggplot(cmp, aes(x = b_primary, y = b_rep)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.4) +
    geom_point(aes(colour = group), alpha = 0.8, size = 2.2) +
    geom_text(data = lab, aes(label = exposure), size = 3, vjust = -0.7, check_overlap = TRUE) +
    labs(
      title = "Replication focus: effect concordance (highlight top decile + consensus)",
      x = "b (primary)",
      y = "b (replication)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  pB <- ggplot(cmp, aes(x = rank_primary, y = rank_rep)) +
    geom_point(aes(colour = group), alpha = 0.8, size = 2.2) +
    scale_x_reverse() + scale_y_reverse() +
    geom_text(data = lab, aes(label = exposure), size = 3, vjust = -0.7, check_overlap = TRUE) +
    labs(
      title = "Replication focus: rank agreement (highlight top decile + consensus)",
      x = "Rank (primary; smaller = higher priority)",
      y = "Rank (replication; smaller = higher priority)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    g <- gridExtra::arrangeGrob(pA, pB, ncol = 2)
    ggsave(file.path(out_root, paste0(prefix, "Fig_ReplicationFocus_panels.png")), g, width = 14, height = 6, dpi = 320)
    ggsave(file.path(out_root, paste0(prefix, "Fig_ReplicationFocus_panels.pdf")), g, width = 14, height = 6)
  } else {
    ggsave(file.path(out_root, paste0(prefix, "Fig_ReplicationFocus_effect_scatter.png")), pA, width = 7.5, height = 6, dpi = 320)
    ggsave(file.path(out_root, paste0(prefix, "Fig_ReplicationFocus_effect_scatter.pdf")), pA, width = 7.5, height = 6)
    ggsave(file.path(out_root, paste0(prefix, "Fig_ReplicationFocus_rank_scatter.png")), pB, width = 7.5, height = 6, dpi = 320)
    ggsave(file.path(out_root, paste0(prefix, "Fig_ReplicationFocus_rank_scatter.pdf")), pB, width = 7.5, height = 6)
  }
}

make_consensus_forest_primarymethod <- function(rank_primary_dt, rank_rep_dt, shortlist_vec, out_root, topn = 20, prefix = "") {
  rp <- as.data.table(rank_primary_dt)
  rr <- as.data.table(rank_rep_dt)
  
  ss <- shortlist_vec
  ss <- ss[ss %in% intersect(rp$exposure, rr$exposure)]
  ss <- head(ss, topn)
  if (length(ss) == 0) return(invisible(NULL))
  
  rp2 <- rp[exposure %in% ss, .(exposure, dataset = "Primary", b, se)]
  rr2 <- rr[exposure %in% ss, .(exposure, dataset = "Replication", b, se)]
  dt <- rbind(rp2, rr2, fill = TRUE)
  
  dt[, ci_lo := b - 1.96 * se]
  dt[, ci_hi := b + 1.96 * se]
  dt[, exposure := factor(exposure, levels = rev(ss))]
  
  dodge <- position_dodge(width = 0.6)
  
  p <- ggplot(dt, aes(x = exposure, y = b, colour = dataset, shape = dataset)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, position = dodge, alpha = 0.85) +
    geom_point(position = dodge, size = 2.2, alpha = 0.9) +
    coord_flip() +
    labs(
      title = "Consensus taxa: primary-method estimates in primary vs replication",
      x = NULL,
      y = "MR effect estimate (b) with 95% CI"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  h <- max(6, 0.28 * length(ss) + 3)
  ggsave(file.path(out_root, paste0(prefix, "Fig_Consensus_forest_primarymethod_primary_vs_replication.png")),
         p, width = 10, height = h, dpi = 320)
  ggsave(file.path(out_root, paste0(prefix, "Fig_Consensus_forest_primarymethod_primary_vs_replication.pdf")),
         p, width = 10, height = h)
}

# ============================================================
# SHAP plotting (base R): summary + importance + dependence
# ============================================================
plot_shap_suite <- function(shap_mat, X_use, outdir) {
  if (is.null(shap_mat) || ncol(shap_mat) == 0) return(invisible(NULL))
  
  pal <- colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))
  
  imp_dt <- data.table(
    feature = colnames(shap_mat),
    mean_abs_shap = colMeans(abs(shap_mat))
  )[order(-mean_abs_shap)]
  
  topn <- min(15, ncol(shap_mat))
  top_feats <- imp_dt$feature[1:topn]
  
  sh_df <- as.data.frame(shap_mat[, top_feats, drop = FALSE])
  x_df  <- as.data.frame(X_use[, top_feats, drop = FALSE])
  
  # grayscale summary
  png(file.path(outdir, "Fig_SHAP_summary_ranking.png"), width = 1800, height = 1200, res = 200)
  par(mar = c(5, 18, 3, 2))
  ypos <- rev(seq_len(topn))
  plot(NA, xlim = range(sh_df, finite = TRUE), ylim = c(0.5, topn + 0.5),
       yaxt = "n", xlab = "SHAP value (impact on model output)", ylab = "",
       main = paste0("SHAP summary (top ", topn, " features)"))
  axis(2, at = ypos, labels = rev(top_feats), las = 1)
  abline(v = 0, lty = 2)
  for (i in seq_len(topn)) {
    f <- rev(top_feats)[i]
    shap_v <- sh_df[[f]]
    x_v <- x_df[[f]]
    if (any(!is.finite(x_v))) {
      med_x <- suppressWarnings(median(x_v[is.finite(x_v)], na.rm = TRUE))
      if (!is.finite(med_x)) med_x <- 0
      x_v[!is.finite(x_v)] <- med_x
    }
    r <- rank(x_v, ties.method = "average") / length(x_v)
    yj <- rep(i, length(shap_v)) + runif(length(shap_v), -0.25, 0.25)
    points(shap_v, yj, pch = 16, cex = 0.5, col = gray(1 - r, alpha = 0.6))
  }
  dev.off()
  
  # grayscale importance
  png(file.path(outdir, "Fig_SHAP_importance_ranking.png"), width = 1600, height = 1000, res = 200)
  par(mar = c(6, 18, 3, 1))
  topm <- min(30, nrow(imp_dt))
  barplot(rev(imp_dt$mean_abs_shap[1:topm]),
          names.arg = rev(imp_dt$feature[1:topm]),
          horiz = TRUE, las = 1,
          xlab = "Mean(|SHAP|)", main = "SHAP feature importance (ranking model)")
  dev.off()
  
  # colored summary
  png(file.path(outdir, "Fig_SHAP_summary_ranking_COLOR.png"), width = 1800, height = 1200, res = 200)
  par(mar = c(5, 18, 3, 2))
  plot(NA, xlim = range(sh_df, finite = TRUE), ylim = c(0.5, topn + 0.5),
       yaxt = "n", xlab = "SHAP value (impact on model output)", ylab = "",
       main = paste0("SHAP summary (top ", topn, " features; colored by feature value)"))
  axis(2, at = ypos, labels = rev(top_feats), las = 1)
  abline(v = 0, lty = 2)
  for (i in seq_len(topn)) {
    f <- rev(top_feats)[i]
    shap_v <- sh_df[[f]]
    x_v <- x_df[[f]]
    if (any(!is.finite(x_v))) {
      med_x <- suppressWarnings(median(x_v[is.finite(x_v)], na.rm = TRUE))
      if (!is.finite(med_x)) med_x <- 0
      x_v[!is.finite(x_v)] <- med_x
    }
    rng <- range(x_v, finite = TRUE)
    z <- if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] != rng[2]) (x_v - rng[1]) / (rng[2] - rng[1]) else rep(0.5, length(x_v))
    col_idx <- pmax(1, pmin(100, as.integer(z * 99) + 1))
    cols <- pal(100)[col_idx]
    yj <- rep(i, length(shap_v)) + runif(length(shap_v), -0.25, 0.25)
    points(shap_v, yj, pch = 16, cex = 0.55, col = adjustcolor(cols, alpha.f = 0.65))
  }
  dev.off()
  
  # colored importance
  imp2 <- imp_dt[!is.na(feature) & !is.na(mean_abs_shap)][order(-mean_abs_shap)]
  topm <- min(30, nrow(imp2))
  cols_imp <- pal(topm)
  
  png(file.path(outdir, "Fig_SHAP_importance_ranking_COLOR.png"), width = 1600, height = 1000, res = 200)
  par(mar = c(6, 18, 3, 1))
  barplot(rev(imp2$mean_abs_shap[1:topm]),
          names.arg = rev(imp2$feature[1:topm]),
          horiz = TRUE, las = 1,
          col = rev(cols_imp),
          xlab = "Mean(|SHAP|)", main = "SHAP feature importance (ranking model)")
  dev.off()
  
  # dependence plots (up to 4)
  dep_feats <- intersect(c("abs_b", "se", "het_Q_p", "nsnp"), imp2$feature)
  if (length(dep_feats) < 4) dep_feats <- imp2$feature[1:min(4, nrow(imp2))]
  
  X_plot <- as.data.frame(X_use)
  for (f in dep_feats[1:min(4, length(dep_feats))]) {
    if (!f %in% colnames(shap_mat)) next
    x_v <- X_plot[[f]]
    s_v <- shap_mat[, f]
    
    color_by <- if ("abs_b" %in% names(X_plot) && f != "abs_b") "abs_b" else if ("se" %in% names(X_plot) && f != "se") "se" else f
    c_v <- X_plot[[color_by]]
    if (any(!is.finite(c_v))) {
      med_c <- suppressWarnings(median(c_v[is.finite(c_v)], na.rm = TRUE))
      if (!is.finite(med_c)) med_c <- 0
      c_v[!is.finite(c_v)] <- med_c
    }
    
    rng <- range(c_v, finite = TRUE)
    z <- if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] != rng[2]) (c_v - rng[1]) / (rng[2] - rng[1]) else rep(0.5, length(c_v))
    col_idx <- pmax(1, pmin(100, as.integer(z * 99) + 1))
    cols <- pal(100)[col_idx]
    
    png(file.path(outdir, paste0("Fig_SHAP_dependence_", f, "_COLOR.png")), width = 1600, height = 1100, res = 200)
    par(mar = c(5, 5, 3, 2))
    plot(x_v, s_v, pch = 16, cex = 0.7,
         col = adjustcolor(cols, alpha.f = 0.65),
         xlab = f, ylab = paste0("SHAP(", f, ")"),
         main = paste0("SHAP dependence: ", f, " (colored by ", color_by, ")"))
    abline(h = 0, lty = 2)
    dev.off()
  }
  
  invisible(imp_dt)
}



# =============================================================================
# 6. LOAD THE PRESELECTED EXPOSURE INSTRUMENT FILE
# =============================================================================

if (!file.exists(IV_PATH)) stop("Exposure IV file not found: ", IV_PATH)
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

iv <- as.data.frame(fread(IV_PATH))
if (!all(c("SNP", "taxa") %in% names(iv))) stop("IV file must contain SNP and taxa columns.")
beta_col <- pick1(c("beta.exposure", "beta", "b", "Beta"), names(iv))
se_col   <- pick1(c("se.exposure", "se", "SE", "standard_error"), names(iv))
ea_col   <- pick1(c("effect_allele.exposure", "effect_allele", "ea", "A1"), names(iv))
oa_col   <- pick1(c("other_allele.exposure", "other_allele", "nea", "A2"), names(iv))
p_col    <- pick1(c("pval.exposure", "pval", "p", "p_value"), names(iv))
eaf_col  <- pick1(c("eaf.exposure", "eaf", "af", "effect_allele_frequency"), names(iv))
if (any(is.na(c(beta_col, se_col, ea_col, oa_col, p_col)))) {
  stop("Could not detect required exposure beta, SE, allele and P-value columns.")
}
iv$SNP <- clean_snp(iv$SNP)
iv[[ea_col]] <- clean_allele(iv[[ea_col]])
iv[[oa_col]] <- clean_allele(iv[[oa_col]])
exp_args <- list(
  dat = iv, type = "exposure", snp_col = "SNP", beta_col = beta_col,
  se_col = se_col, effect_allele_col = ea_col, other_allele_col = oa_col,
  pval_col = p_col, phenotype_col = "taxa"
)
if (!is.na(eaf_col)) exp_args$eaf_col <- eaf_col
exp_dat <- do.call(TwoSampleMR::format_data, exp_args)

exp_audit <- data.table(
  metric = c(
    "input_rows", "unique_taxa", "unique_SNPs", "rows_p_le_threshold",
    "rows_p_gt_threshold", "APPLY_IV_P_FILTER", "IV_P_THRESHOLD"
  ),
  value = as.character(c(
    nrow(exp_dat), uniqueN(exp_dat$exposure), uniqueN(exp_dat$SNP),
    sum(exp_dat$pval.exposure <= IV_P_THRESHOLD, na.rm = TRUE),
    sum(exp_dat$pval.exposure > IV_P_THRESHOLD, na.rm = TRUE),
    APPLY_IV_P_FILTER, IV_P_THRESHOLD
  ))
)
fwrite(exp_audit, file.path(OUT_ROOT, "AUDIT_exposure_instrument_input.csv"))
if (APPLY_IV_P_FILTER) {
  before <- nrow(exp_dat)
  exp_dat <- exp_dat[is.finite(exp_dat$pval.exposure) & exp_dat$pval.exposure <= IV_P_THRESHOLD, ]
  message("Applied IV P filter: ", nrow(exp_dat), "/", before, " rows retained.")
} else {
  message("Preserving the preselected IV file; no new IV P-value filter applied.")
}
if (nrow(exp_dat) == 0L) stop("No exposure instruments available.")

dup_audit <- as.data.table(exp_dat)[, .N, by = .(exposure, SNP)][N > 1L]
if (nrow(dup_audit) > 0L) {
  fwrite(dup_audit, file.path(OUT_ROOT, "AUDIT_duplicate_exposure_taxon_SNP.csv"))
  warning("Duplicate taxon-SNP rows detected; see audit output.")
}

# =============================================================================
# 7. ONE-OUTCOME END-TO-END PIPELINE
# =============================================================================

run_one_outcome <- function(config, out_path) {
  trait <- config$trait; role <- config$role; tag <- config$tag
  outdir <- file.path(OUT_ROOT, trait, tag)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("\n============================================================")
  message("RUNNING: ", tag)
  message("Outcome file: ", out_path)

  out <- as.data.frame(fread(out_path))
  snp_o  <- pick1(c("SNP", "rsid", "variant", "rs_id", "rsID"), names(out))
  beta_o <- pick1(c("beta.outcome", "beta", "b", "beta_EUR", "beta_EAS", "beta_AFR"), names(out))
  se_o   <- pick1(c("se.outcome", "se", "stderr", "standard_error", "se_EUR", "se_EAS", "se_AFR"), names(out))
  ea_o   <- pick1(c("effect_allele.outcome", "effect_allele", "ea", "allele1", "A1"), names(out))
  oa_o   <- pick1(c("other_allele.outcome", "other_allele", "nea", "allele0", "A2"), names(out))
  p_o    <- pick1(c("pval.outcome", "pval", "p", "p_value"), names(out))
  eaf_o  <- pick1(c("eaf.outcome", "eaf", "af", "effect_allele_frequency"), names(out))
  if (any(is.na(c(snp_o, beta_o, se_o, ea_o, oa_o)))) {
    stop("Outcome columns missing for ", tag, ". Found: ", paste(names(out), collapse = ", "))
  }
  out_dat <- safe_format_outcome(out, snp_o, beta_o, se_o, ea_o, oa_o,
                                 p = p_o, eaf = eaf_o, outcome_name = tag)
  dat_h <- TwoSampleMR::harmonise_data(exp_dat, out_dat, action = HARMONISE_ACTION)
  if (nrow(dat_h) == 0L) stop("No harmonised variants for ", tag)
  dh <- as.data.table(dat_h)
  dh[, F_snp := (beta.exposure / se.exposure)^2]

  strength_audit <- dh[, .(
    nsnp_before_optional_filter = uniqueN(SNP),
    min_F_before = suppressWarnings(min(F_snp, na.rm = TRUE)),
    mean_F_before = suppressWarnings(mean(F_snp, na.rm = TRUE)),
    median_F_before = suppressWarnings(median(F_snp, na.rm = TRUE)),
    n_F_lt10_before = sum(F_snp < MIN_F_STATISTIC, na.rm = TRUE)
  ), by = exposure]
  for (cc in c("min_F_before", "mean_F_before", "median_F_before")) {
    strength_audit[!is.finite(get(cc)), (cc) := NA_real_]
  }

  n_weak_removed <- 0L
  if (APPLY_F_FILTER) {
    keep <- !is.finite(dh$F_snp) | dh$F_snp >= MIN_F_STATISTIC
    n_weak_removed <- sum(is.finite(dh$F_snp) & dh$F_snp < MIN_F_STATISTIC)
    dh <- dh[keep]
    if (nrow(dh) == 0L) stop("No variants remain after optional F filtering for ", tag)
  }

  strength <- dh[, .(
    nsnp = uniqueN(SNP),
    min_F = suppressWarnings(min(F_snp, na.rm = TRUE)),
    mean_F = suppressWarnings(mean(F_snp, na.rm = TRUE)),
    median_F = suppressWarnings(median(F_snp, na.rm = TRUE)),
    n_weak_F_lt10 = sum(F_snp < MIN_F_STATISTIC, na.rm = TRUE)
  ), by = exposure]
  for (cc in c("min_F", "mean_F", "median_F")) strength[!is.finite(get(cc)), (cc) := NA_real_]
  fwrite(strength_audit, file.path(outdir, "00_instrument_strength_audit_before_optional_filter.csv"))
  fwrite(strength, file.path(outdir, "00_instrument_strength_used.csv"))

  mr_res <- TwoSampleMR::mr(as.data.frame(dh))
  mr_dt <- as.data.table(mr_res)
  fwrite(mr_dt, file.path(outdir, "01_MR_results_all_estimators.csv"))
  fwrite(mr_dt, file.path(outdir, "1MR_results.csv"))

  mr_dt[, `:=`(method = as.character(method), exposure = as.character(exposure))]
  pref_order <- c("Inverse variance weighted", "Wald ratio", "Weighted median",
                  "MR Egger", "Weighted mode", "Simple mode")
  mr_dt[, method_rank := match(method, pref_order)]
  mr_dt[is.na(method_rank), method_rank := 999L]
  primary <- mr_dt[order(method_rank), .SD[1], by = exposure]
  if ("nsnp" %in% names(primary)) primary[, nsnp := NULL]
  primary <- merge(primary, strength, by = "exposure", all.x = TRUE)
  primary[, `:=`(
    fdr = p.adjust(pval, method = "BH"),
    neglog10_p = -log10(pmax(pval, 1e-300)),
    abs_b = abs(b),
    signal_to_noise = fifelse(is.finite(se) & se > 0, abs(b) / se, NA_real_),
    sign_b = fifelse(b >= 0, 1L, -1L)
  )]

  exp_ge2 <- strength[nsnp >= 2L, exposure]
  het <- tryCatch(
    suppressMessages(as.data.table(TwoSampleMR::mr_heterogeneity(as.data.frame(dh[exposure %in% exp_ge2])))),
    error = function(e) data.table()
  )
  het_one <- data.table(exposure = primary$exposure, het_Q = NA_real_, het_Q_p = NA_real_)
  if (nrow(het) > 0L) {
    het[, method := as.character(method)]
    pick <- het[method == "Inverse variance weighted"]
    if (nrow(pick) == 0L) pick <- het
    pick <- pick[, .SD[1], by = exposure]
    setnames(pick, c("Q", "Q_pval"), c("het_Q", "het_Q_p"), skip_absent = TRUE)
    het_one <- merge(het_one[, .(exposure)], pick[, .(exposure, het_Q, het_Q_p)], by = "exposure", all.x = TRUE)
  }

  exp_ge3 <- strength[nsnp >= 3L, exposure]
  plei <- tryCatch(
    suppressMessages(as.data.table(TwoSampleMR::mr_pleiotropy_test(as.data.frame(dh[exposure %in% exp_ge3])))),
    error = function(e) data.table()
  )
  plei_one <- data.table(exposure = primary$exposure, egger_int = NA_real_, egger_p = NA_real_)
  if (nrow(plei) > 0L) {
    setnames(plei, c("egger_intercept", "pval"), c("egger_int", "egger_p"), skip_absent = TRUE)
    plei_one <- merge(plei_one[, .(exposure)], plei[, .(exposure, egger_int, egger_p)], by = "exposure", all.x = TRUE)
  }

  steiger_one <- data.table(exposure = primary$exposure, steiger_dir = NA)
  need_steiger <- all(c("pval.exposure", "pval.outcome", "samplesize.exposure", "samplesize.outcome") %in% names(dh)) ||
    all(c("r.exposure", "r.outcome") %in% names(dh))
  if (need_steiger) {
    steiger_one <- tryCatch({
      st <- suppressMessages(as.data.table(TwoSampleMR::directionality_test(as.data.frame(dh))))
      cand <- intersect(c("correct_causal_direction", "steiger_dir", "correct_direction"), names(st))
      if (length(cand) == 0L) st[, steiger_dir := NA] else setnames(st, cand[1], "steiger_dir")
      st[, .(exposure, steiger_dir)]
    }, error = function(e) data.table(exposure = primary$exposure, steiger_dir = NA))
  }

  feat <- Reduce(function(x, y) merge(x, y, by = "exposure", all.x = TRUE),
                 list(primary, het_one, plei_one, steiger_one))
  feat <- compute_evidence_components(feat)
  fwrite(feat, file.path(outdir, "02_evidence_feature_matrix.csv"))
  fwrite(feat, file.path(outdir, "2MR_feature_table.csv"))

  feature_names <- c(
    "neglog10_p", "b", "se", "abs_b", "signal_to_noise", "sign_b",
    "nsnp", "min_F", "mean_F", "median_F", "n_weak_F_lt10",
    "het_Q", "het_Q_p", "egger_int", "egger_p",
    "miss_het_Q", "miss_het_Q_p", "miss_egger_int", "miss_egger_p"
  )
  feature_names <- intersect(feature_names, names(feat))
  X <- as.matrix(feat[, ..feature_names])
  storage.mode(X) <- "double"
  colnames(X) <- feature_names
  y_label <- feat$relevance_label
  y_cont <- feat$evidence_grade
  taxa <- feat$exposure

  fold_plan <- make_stratified_fold_plan(
    y_cont, CV_FOLDS, CV_REPEATS, MASTER_SEED + sum(utf8ToInt(tag))
  )
  cv <- run_repeated_cv_evaluation(
    X, y_label, y_cont, taxa, fold_plan, XGB_PARAMS,
    XGB_NROUNDS_MAX, XGB_EARLY_STOP, seed_offset = sum(utf8ToInt(tag))
  )
  best_rounds <- as.integer(round(median(cv$fold_stats$best_iteration, na.rm = TRUE)))
  if (!is.finite(best_rounds) || best_rounds < 1L) best_rounds <- 300L

  final_model <- fit_final_ranker(
    X, y_label, XGB_PARAMS, best_rounds,
    MASTER_SEED + 50000L + sum(utf8ToInt(tag))
  )
  d_all <- xgb.DMatrix(X, missing = NA_real_)
  final_score <- predict(final_model, d_all)
  shap <- as.matrix(predict(final_model, d_all, predcontrib = TRUE))
  bias <- which(colnames(shap) %in% c("BIAS", "bias", "Bias"))
  if (length(bias) == 1L) shap <- shap[, -bias, drop = FALSE]
  shap_imp <- data.table(feature = colnames(shap), mean_abs_shap = colMeans(abs(shap), na.rm = TRUE))[order(-mean_abs_shap)]

  feat[, `:=`(
    score_xgb = final_score,
    score_xgb_full_model = final_score,
    cv_oof_percentile_mean = cv$oof_percentile_mean,
    cv_oof_percentile_sd = cv$oof_percentile_sd,
    cv_oof_rank_mean = cv$oof_rank_mean,
    cv_oof_rank_sd = cv$oof_rank_sd
  )]
  ranking <- add_baseline_scores_and_ranks(feat)
  ranking[, `:=`(
    score_true = evidence_grade,
    score_pred = score_xgb,
    method_primary = method
  )]
  setorder(ranking, rank_xgb)

  fwrite(ranking, file.path(outdir, "03_rankings_all_methods.csv"))
  fwrite(ranking, file.path(outdir, "3AI_ranked_taxa_regression.csv"))
  fwrite(shap_imp, file.path(outdir, "04_SHAP_final_model_importance.csv"))
  fwrite(shap_imp, file.path(outdir, "4SHAP_feature_importance_regression.csv"))
  fwrite(cv$fold_stats, file.path(outdir, "03_CV_fold_performance.csv"))
  fwrite(cv$repeat_stability, file.path(outdir, "03_CV_repeat_rank_stability.csv"))
  evidence_alias <- ranking[, .(
    taxon = exposure, beta = b, se = se, p = pval, nsnp = nsnp,
    phet = het_Q_p, pegg = egger_p, score_xgb = score_xgb
  )]
  fwrite(evidence_alias, file.path(outdir, paste0(toupper(role), "_evidence_table.csv")))
  within_metrics <- within_dataset_metrics(ranking)
  fwrite(within_metrics, file.path(outdir, "Table_AI_vs_baselines_metrics.csv"))
  if (MAKE_FIGURES) {
    corr_long <- melt(within_metrics,
                      id.vars = "baseline",
                      measure.vars = c("spearman_rho", "kendall_tau"),
                      variable.name = "metric", value.name = "value")
    pcorr <- ggplot(corr_long, aes(x = value, y = baseline)) +
      geom_vline(xintercept = 0, linetype = "dashed") + geom_point(size = 2.4) +
      facet_wrap(~metric, ncol = 1) + theme_bw() +
      labs(title = paste0("XGBoost ranking vs transparent baselines — ", tag),
           x = "Rank correlation", y = NULL)
    ggsave(file.path(outdir, "Fig_AI_vs_baselines_rank_correlations.png"),
           pcorr, width = 7, height = 6, dpi = 320)
  }
  fwrite(data.table(
    exposure = taxa,
    oof_percentile_mean = cv$oof_percentile_mean,
    oof_percentile_sd = cv$oof_percentile_sd,
    oof_rank_mean = cv$oof_rank_mean,
    oof_rank_sd = cv$oof_rank_sd
  ), file.path(outdir, "03_CV_oof_diagnostics.csv"))

  # Full-data ablations. Final ranking remains the full model.
  ablation_summary <- data.table()
  ablation_scores <- list()
  if (tag %in% RUN_ABLATIONS_FOR) {
    ablation_sets <- list(
      no_QC = setdiff(feature_names, c(
        "het_Q", "het_Q_p", "egger_int", "egger_p",
        "miss_het_Q", "miss_het_Q_p", "miss_egger_int", "miss_egger_p"
      )),
      no_instrument_support = setdiff(feature_names, c(
        "nsnp", "min_F", "mean_F", "median_F", "n_weak_F_lt10"
      ))
    )
    full_top <- taxa[order(final_score, decreasing = TRUE)[seq_len(min(TOP_K, length(taxa)))]]
    for (nm in names(ablation_sets)) {
      cols <- ablation_sets[[nm]]
      Xab <- X[, cols, drop = FALSE]
      mod <- fit_final_ranker(
        Xab, y_label, XGB_PARAMS, best_rounds,
        MASTER_SEED + 70000L + match(nm, names(ablation_sets))
      )
      pred <- predict(mod, xgb.DMatrix(Xab, missing = NA_real_))
      ablation_scores[[nm]] <- pred
      top_ab <- taxa[order(pred, decreasing = TRUE)[seq_len(min(TOP_K, length(taxa)))]]
      ablation_summary <- rbind(
        ablation_summary,
        data.table(
          model = nm, n_features = length(cols),
          rank_spearman_vs_full = safe_cor(pred, final_score, "spearman"),
          topK_jaccard_vs_full = jaccard(top_ab, full_top)
        )
      )
      ranking[, (paste0("score_xgb_", nm)) := pred]
      ranking[, (paste0("rank_xgb_", nm)) := frank(-pred, ties.method = "average")]
    }
    fwrite(ablation_summary, file.path(outdir, "05_ablation_summary.csv"))
    fwrite(ranking, file.path(outdir, "03_rankings_all_methods.csv"))
    fwrite(ranking, file.path(outdir, "3AI_ranked_taxa_regression.csv"))
  }

  model_bootstrap <- NULL
  if (tag %in% RUN_MODEL_BOOTSTRAP_FOR) {
    message("Running ", MODEL_BOOTSTRAP_REFITS, " true XGBoost bootstrap refits for ", tag)
    model_bootstrap <- run_model_bootstrap_refits(
      X, y_label, taxa, final_score, XGB_PARAMS, best_rounds,
      MODEL_BOOTSTRAP_REFITS, MASTER_SEED + 90000L
    )
    fwrite(model_bootstrap$stability, file.path(outdir, "06_model_refit_bootstrap_ranking_stability.csv"))
    fwrite(model_bootstrap$shap_long, file.path(outdir, "06_model_refit_bootstrap_SHAP_long.csv"))
    fwrite(model_bootstrap$shap_summary, file.path(outdir, "06_model_refit_bootstrap_SHAP_summary.csv"))
  }

  method_bootstrap <- NULL
  if (tag %in% RUN_METHOD_BOOTSTRAP_FOR) {
    method_bootstrap <- run_method_taxon_bootstrap(
      ranking, METHOD_BOOTSTRAPS, MASTER_SEED + 110000L + sum(utf8ToInt(tag))
    )
    fwrite(method_bootstrap$long, file.path(outdir, "07_fixed_score_taxon_bootstrap_all_methods_long.csv"))
    fwrite(method_bootstrap$summary, file.path(outdir, "07_fixed_score_taxon_bootstrap_all_methods_summary.csv"))
  }

  robust20 <- head(ranking[nsnp >= 2][order(rank_xgb)], 20L)
  single20 <- head(ranking[nsnp == 1][order(rank_xgb)], 20L)
  fwrite(robust20, file.path(outdir, "Table_top20_nsnp_ge2.csv"))
  fwrite(single20, file.path(outdir, "Table_top20_nsnp_eq1.csv"))

  if (MAKE_FIGURES) {
    make_volcano(mr_dt, outdir, tag)
    if (nrow(robust20) > 0L) make_heatmap_top20(robust20, outdir, tag)
    shortlist_for_forest <- head(ranking[order(rank_xgb), exposure], FOREST_TOP_N)
    make_forest_methods(mr_dt, shortlist_for_forest, outdir, tag)
    make_model_diagnostic_plots(ranking, outdir, tag)
    plot_shap_suite(shap, X, outdir)
    if (!is.null(model_bootstrap)) {
      pboot <- ggplot(model_bootstrap$stability, aes(y = topK_jaccard, x = factor(1))) +
        geom_boxplot() + theme_bw() +
        labs(title = paste0("XGBoost model-refit bootstrap stability — ", tag),
             x = NULL, y = paste0("Jaccard@", TOP_K)) +
        theme(axis.text.x = element_blank())
      ggsave(file.path(outdir, "Fig_model_refit_bootstrap_Jaccard.png"), pboot,
             width = 5.5, height = 5.5, dpi = 320)
    }
    if (!is.null(method_bootstrap)) {
      pmethod <- ggplot(method_bootstrap$long, aes(x = method, y = topK_jaccard)) +
        geom_boxplot() + coord_flip() + theme_bw() +
        labs(title = paste0("Taxon-resampling stability across methods — ", tag),
             x = NULL, y = paste0("Jaccard@", TOP_K))
      ggsave(file.path(outdir, "Fig_fixed_score_taxon_bootstrap_all_methods.png"), pmethod,
             width = 7, height = 5, dpi = 320)
    }
  }

  param_dt <- data.table(
    parameter = c(
      "trait", "role", "tag", "outcome_path", "n_harmonised_rows_used",
      "n_taxa", "APPLY_IV_P_FILTER", "IV_P_THRESHOLD", "APPLY_F_FILTER",
      "MIN_F_STATISTIC", "n_weak_rows_removed", "CV_FOLDS", "CV_REPEATS",
      "TOP_K", "best_rounds", "median_CV_target_spearman",
      "median_CV_target_NDCG", "median_CV_repeat_rank_spearman",
      "model_bootstrap_refits", "model_bootstrap_median_topK_jaccard",
      "model_bootstrap_median_rank_spearman", "model_bootstrap_SHAP_pairwise_spearman",
      "fixed_score_method_bootstraps"
    ),
    value = as.character(c(
      trait, role, tag, out_path, nrow(dh), nrow(feat), APPLY_IV_P_FILTER,
      IV_P_THRESHOLD, APPLY_F_FILTER, MIN_F_STATISTIC, n_weak_removed,
      CV_FOLDS, CV_REPEATS, TOP_K, best_rounds,
      median(cv$fold_stats$target_spearman, na.rm = TRUE),
      median(cv$fold_stats$target_ndcg, na.rm = TRUE),
      median(cv$repeat_stability$rank_spearman_vs_mean, na.rm = TRUE),
      if (is.null(model_bootstrap)) 0L else MODEL_BOOTSTRAP_REFITS,
      if (is.null(model_bootstrap)) NA_real_ else median(model_bootstrap$stability$topK_jaccard, na.rm = TRUE),
      if (is.null(model_bootstrap)) NA_real_ else median(model_bootstrap$stability$rank_spearman, na.rm = TRUE),
      if (is.null(model_bootstrap)) NA_real_ else model_bootstrap$shap_pairwise_spearman,
      if (is.null(method_bootstrap)) 0L else METHOD_BOOTSTRAPS
    ))
  )
  fwrite(param_dt, file.path(outdir, "00_analysis_parameters_and_audit.csv"))

  list(
    trait = trait, role = role, tag = tag, out_path = out_path, outdir = outdir,
    feature = feat, ranking = ranking, mr = mr_dt, cv = cv, shap = shap_imp,
    shap_matrix = shap, model_bootstrap = model_bootstrap,
    method_bootstrap = method_bootstrap, ablation = ablation_summary,
    within_metrics = within_metrics, best_rounds = best_rounds
  )
}

# =============================================================================
# 8. RUN ALL CONFIGURED OUTCOMES
# =============================================================================

resolved <- lapply(OUTCOMES, function(cfg) {
  path <- first_existing(cfg$candidates)
  if (is.na(path) && isTRUE(cfg$required)) {
    stop("Required outcome file not found for ", cfg$tag, ". Tried:\n",
         paste(cfg$candidates, collapse = "\n"))
  }
  list(config = cfg, path = path)
})

runs <- list(); manifest <- list()
for (item in resolved) {
  cfg <- item$config
  if (is.na(item$path)) {
    manifest[[length(manifest) + 1L]] <- data.table(
      trait = cfg$trait, role = cfg$role, tag = cfg$tag,
      status = "SKIPPED_FILE_NOT_FOUND", path = NA_character_
    )
    next
  }
  result <- tryCatch(
    run_one_outcome(cfg, item$path),
    error = function(e) {
      message("FAILED: ", cfg$tag, " | ", conditionMessage(e))
      structure(list(error = conditionMessage(e)), class = "failed_run")
    }
  )
  if (inherits(result, "failed_run")) {
    manifest[[length(manifest) + 1L]] <- data.table(
      trait = cfg$trait, role = cfg$role, tag = cfg$tag,
      status = paste0("FAILED: ", result$error), path = item$path
    )
  } else {
    runs[[cfg$tag]] <- result
    manifest[[length(manifest) + 1L]] <- data.table(
      trait = cfg$trait, role = cfg$role, tag = cfg$tag,
      status = "COMPLETED", path = item$path
    )
  }
}
manifest_dt <- rbindlist(manifest, fill = TRUE)
fwrite(manifest_dt, file.path(OUT_ROOT, "RUN_MANIFEST.csv"))
required_tags <- vapply(Filter(function(x) isTRUE(x$required), OUTCOMES), function(x) x$tag, character(1))
failed_required <- setdiff(required_tags, names(runs))
if (length(failed_required) > 0L) stop("Required run(s) did not complete: ", paste(failed_required, collapse = ", "))

# =============================================================================
# 9. PRIMARY–REPLICATION BENCHMARKING AND CONSENSUS SHORTLISTS
# =============================================================================

CORE_METHODS <- c("xgb", "pval", "abeta", "compE", "borda", "pareto")

compare_primary_replication <- function(primary_run, replication_run, trait,
                                        prefix = paste0(trait, "_")) {
  outdir <- file.path(OUT_ROOT, trait, "CROSS_DATASET")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  p <- copy(primary_run$ranking); r <- copy(replication_run$ranking)
  pkeep <- intersect(c("exposure", "b", "se", "pval", "fdr", "nsnp",
                       paste0("rank_", CORE_METHODS), paste0("score_", CORE_METHODS),
                       grep("^rank_xgb_no_", names(p), value = TRUE),
                       grep("^score_xgb_no_", names(p), value = TRUE)), names(p))
  rkeep <- intersect(c("exposure", "b", "se", "pval", "fdr", "nsnp",
                       paste0("rank_", CORE_METHODS), paste0("score_", CORE_METHODS)), names(r))
  p <- p[, ..pkeep]; r <- r[, ..rkeep]
  setnames(p, setdiff(names(p), "exposure"), paste0(setdiff(names(p), "exposure"), "_primary"))
  setnames(r, setdiff(names(r), "exposure"), paste0(setdiff(names(r), "exposure"), "_rep"))
  cmp <- merge(p, r, by = "exposure", all = FALSE)
  cmp[, `:=`(
    dir_concordant = sign(b_primary) == sign(b_rep),
    replication_nominal = is.finite(pval_rep) & pval_rep < 0.05,
    replication_supported = as.integer(
      sign(b_primary) == sign(b_rep) & is.finite(pval_rep) & pval_rep < 0.05
    ),
    abs_rank_diff = abs(rank_xgb_primary - rank_xgb_rep)
  )]

  benchmark <- rbindlist(lapply(CORE_METHODS, function(m) {
    rp <- paste0("rank_", m, "_primary"); rr <- paste0("rank_", m, "_rep")
    ord_p <- order(cmp[[rp]], decreasing = FALSE)
    ord_r <- order(cmp[[rr]], decreasing = FALSE)
    top_p <- cmp$exposure[ord_p[seq_len(min(TOP_K, nrow(cmp)))]]
    top_r <- cmp$exposure[ord_r[seq_len(min(TOP_K, nrow(cmp)))]]
    data.table(
      trait = trait, method = m, n_shared = nrow(cmp),
      n_replication_supported = sum(cmp$replication_supported, na.rm = TRUE),
      precision_at_K = precision_at_k(cmp$replication_supported, ord_p, TOP_K),
      ndcg_at_K = ndcg_at_k(cmp$replication_supported, ord_p, TOP_K),
      cross_dataset_rank_spearman = safe_cor(cmp[[rp]], cmp[[rr]], "spearman"),
      topK_overlap = length(intersect(top_p, top_r)),
      topK_jaccard = jaccard(top_p, top_r)
    )
  }))

  ablation_cols <- grep("^rank_xgb_no_.*_primary$", names(cmp), value = TRUE)
  ablation_rep <- rbindlist(lapply(ablation_cols, function(rc) {
    ord <- order(cmp[[rc]], decreasing = FALSE)
    data.table(
      trait = trait,
      ablation = sub("^rank_xgb_", "", sub("_primary$", "", rc)),
      precision_at_K = precision_at_k(cmp$replication_supported, ord, TOP_K),
      ndcg_at_K = ndcg_at_k(cmp$replication_supported, ord, TOP_K),
      rank_spearman_vs_full = safe_cor(cmp[[rc]], cmp$rank_xgb_primary, "spearman")
    )
  }), fill = TRUE)

  top_n <- ceiling(CONSENSUS_TOP_PROP * nrow(primary_run$ranking))
  top_primary <- head(primary_run$ranking[order(rank_xgb), exposure], top_n)
  consensus <- cmp[exposure %in% top_primary & dir_concordant == TRUE][order(rank_xgb_primary)]

  fwrite(cmp, file.path(outdir, paste0(trait, "_cross_dataset_comparison.csv")))
  fwrite(cmp, file.path(OUT_ROOT, paste0("5CrossDataset_compare_", trait, ".csv")))
  fwrite(benchmark, file.path(outdir, paste0(trait, "_benchmarking_metrics.csv")))
  fwrite(consensus, file.path(outdir, paste0(trait, "_consensus_top10pct_direction_concordant.csv")))
  fwrite(consensus, file.path(OUT_ROOT, paste0("6Consensus_shortlist_", trait, "_top10pct_primary_dirconcordant.csv")))
  if (nrow(ablation_rep) > 0L) fwrite(ablation_rep, file.path(outdir, paste0(trait, "_ablation_replication_metrics.csv")))
  method_stability <- primary_run$method_bootstrap$summary %||% data.table()
  benchmark_main <- copy(benchmark)
  if (nrow(method_stability) > 0L) {
    benchmark_main <- merge(benchmark_main, method_stability[, .(
      method, bootstrap_jaccard_mean = jaccard_mean,
      bootstrap_jaccard_sd = jaccard_sd
    )], by = "method", all.x = TRUE)
  }
  fwrite(benchmark_main, file.path(outdir, "benchmark_table_main.csv"))

  if (MAKE_FIGURES) {
    rho <- safe_cor(cmp$rank_xgb_primary, cmp$rank_xgb_rep, "spearman")
    plot_cmp <- copy(cmp)
    setnames(plot_cmp,
             c("rank_xgb_primary", "rank_xgb_rep"),
             c("rank_primary", "rank_rep"))
    make_crossdataset_plots(plot_cmp, outdir, spearman_val = rho, prefix = prefix)
    make_replication_focus_panels(
      plot_cmp, top_primary, consensus$exposure, outdir, prefix = prefix
    )
    rp_plot <- copy(primary_run$ranking); rr_plot <- copy(replication_run$ranking)
    make_consensus_forest_primarymethod(
      rp_plot, rr_plot, consensus$exposure, outdir, topn = FOREST_TOP_N, prefix = prefix
    )
    long <- melt(benchmark, id.vars = c("trait", "method"),
                 measure.vars = c("precision_at_K", "ndcg_at_K"),
                 variable.name = "metric", value.name = "value")
    pb <- ggplot(long, aes(x = method, y = value, fill = metric)) +
      geom_col(position = "dodge") + theme_bw() +
      labs(title = paste0(trait, " external replication benchmarking at K=", TOP_K),
           x = NULL, y = "Metric")
    ggsave(file.path(outdir, paste0("Fig_", trait, "_external_replication_benchmark.png")),
           pb, width = 8, height = 5.5, dpi = 320)
  }
  list(compare = cmp, benchmark = benchmark_main, consensus = consensus, ablation = ablation_rep)
}

cross_results <- list(); all_benchmark <- list()
for (trait in c("CRP", "T2D", "SBP")) {
  ptag <- Filter(function(x) x$trait == trait && x$role == "primary", OUTCOMES)[[1]]$tag
  rtag <- Filter(function(x) x$trait == trait && x$role == "replication", OUTCOMES)[[1]]$tag
  if (all(c(ptag, rtag) %in% names(runs))) {
    res <- compare_primary_replication(runs[[ptag]], runs[[rtag]], trait)
    cross_results[[trait]] <- res
    all_benchmark[[trait]] <- res$benchmark
  }
}
if (length(all_benchmark) > 0L) {
  fwrite(rbindlist(all_benchmark, fill = TRUE),
         file.path(OUT_ROOT, "ALL_TRAITS_benchmarking_summary.csv"))
}

# Within-family SBP is a sensitivity comparison, not independent replication.
if (all(c("SBP_primary_ieu-b-38", "SBP_within_family_ieu-b-4817") %in% names(runs))) {
  sens <- compare_primary_replication(
    runs[["SBP_primary_ieu-b-38"]], runs[["SBP_within_family_ieu-b-4817"]],
    "SBP_WITHIN_FAMILY_SENSITIVITY", prefix = "SBP_withinfamily_"
  )
  cross_results[["SBP_within_family"]] <- sens
}

# =============================================================================
# 10. WORKFLOW AND GRAPHICAL-ABSTRACT FIGURES
# =============================================================================

if (MAKE_DESIGN_FIGURES) {
  figdir_design <- file.path(OUT_ROOT, "FIGS_design")
  dir.create(figdir_design, recursive = TRUE, showWarnings = FALSE)
  gp_box <- gpar(col = "grey35", fill = "grey98", lwd = 0.8)
  draw_box <- function(x, y, w, h, label) {
    grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(w, "npc"), unit(h, "npc"),
                   r = unit(0.05, "snpc"), gp = gp_box)
    grid.text(label, unit(x, "npc"), unit(y, "npc"), gp = gpar(fontsize = 9))
  }
  draw_arrow2 <- function(x0, y0, x1, y1) {
    grid.segments(unit(x0, "npc"), unit(y0, "npc"), unit(x1, "npc"), unit(y1, "npc"),
                  arrow = arrow(type = "closed", length = unit(0.12, "inches")),
                  gp = gpar(col = "grey35"))
  }
  draw_workflow <- function() {
    grid.text("Replication-aware microbiome MR evidence-ranking workflow",
              unit(0.04, "npc"), unit(0.94, "npc"), just = c("left", "top"),
              gp = gpar(fontsize = 13, fontface = "bold"))
    xs <- c(0.10, 0.28, 0.46, 0.64, 0.82)
    labs <- c(
      "Preselected\ntaxon instruments",
      "Harmonisation\n+ two-sample MR",
      "Evidence matrix\n+ missingness",
      "Pairwise XGBoost\n+ SHAP",
      "Replication benchmark\n+ consensus shortlist"
    )
    for (i in seq_along(xs)) draw_box(xs[i], 0.52, 0.15, 0.30, labs[i])
    for (i in 1:4) draw_arrow2(xs[i] + 0.075, 0.52, xs[i + 1] - 0.075, 0.52)
    grid.text("CRP primary/replication; portability to SBP and T2D; within-family SBP sensitivity",
              unit(0.04, "npc"), unit(0.12, "npc"), just = c("left", "bottom"),
              gp = gpar(fontsize = 9, col = "grey35"))
  }
  png(file.path(figdir_design, "Fig1_WorkflowDiagram.png"), width = 11, height = 5.5, units = "in", res = 320)
  grid.newpage(); draw_workflow(); dev.off()
  cairo_pdf(file.path(figdir_design, "Fig1_WorkflowDiagram.pdf"), width = 11, height = 5.5)
  grid.newpage(); draw_workflow(); dev.off()

  png(file.path(figdir_design, "Graphical_Abstract.png"), width = 1328, height = 531, res = 300)
  grid.newpage(); draw_workflow(); dev.off()
  cairo_pdf(file.path(figdir_design, "Graphical_Abstract.pdf"), width = 1328/300, height = 531/300)
  grid.newpage(); draw_workflow(); dev.off()
}

# =============================================================================
# 11. READY-TO-PASTE METHODS TEXT AND REPRODUCIBILITY OUTPUTS
# =============================================================================

methods_text <- paste0(
"METHODS INSERT — USE ONLY AFTER RUN_MANIFEST SHOWS COMPLETED\n",
"===========================================================\n\n",
"The analysis used a previously selected and linkage-disequilibrium-clumped taxon instrument file. The final pipeline audited exposure-association P values and per-SNP F statistics but did not alter the prespecified instrument set (APPLY_IV_P_FILTER=", APPLY_IV_P_FILTER, "; APPLY_F_FILTER=", APPLY_F_FILTER, "). Exposure and outcome alleles were harmonised using TwoSampleMR with action=", HARMONISE_ACTION, ". The Wald ratio was used for single-instrument taxa and inverse-variance weighted MR for taxa with at least two instruments; complementary estimators and sensitivity diagnostics were retained when estimable.\n\n",
"For each outcome dataset, MR effect estimates, uncertainty, instrument count, per-taxon F-statistic summaries, Cochran Q, MR-Egger intercepts and explicit indicators of unavailable diagnostics were assembled into an evidence matrix. A prespecified weak-supervision grade combined rank-normalised statistical strength, signal-to-noise and instrument support with weights ", W_STAT, ", ", W_SNR, " and ", W_INSTRUMENT, ", respectively. Penalties of ", LAMBDA_HET, " were applied when P_heterogeneity<0.05 or P_Egger<0.05.\n\n",
"XGBoost used a pairwise ranking objective with one query group per outcome dataset. Model behaviour was evaluated using ", CV_REPEATS, " repeated ", CV_FOLDS, "-fold cross-validations; fold-level predictions were used only for target-reconstruction diagnostics and were converted to within-fold percentiles for cross-fold summaries. The final taxon ordering was obtained from a model refitted to the complete evidence matrix using the median early-stopped iteration selected across CV folds. SHAP values were calculated from this final fitted model.\n\n",
"The XGBoost ranking was compared with rankings based on P value, absolute effect magnitude, the deterministic composite grade, Borda aggregation and iterative Pareto-front sorting. External utility was evaluated among taxa shared by primary and replication datasets; replication-supported taxa had concordant effect directions and nominal replication evidence (P_rep<0.05). Precision@", TOP_K, " and NDCG@", TOP_K, " were computed against this set. Cross-dataset robustness was additionally quantified using Spearman rank correlation, top-", TOP_K, " overlap and Jaccard similarity. The consensus shortlist required top-decile priority in the primary dataset and concordant effect direction in replication.\n\n",
"For the primary CRP dataset, models were refitted after removing QC/sensitivity features and after removing instrument-support features. Ranking stability was assessed through ", MODEL_BOOTSTRAP_REFITS, " bootstrap XGBoost refits, with top-", TOP_K, " Jaccard similarity, global Spearman correlation and SHAP-importance concordance calculated relative to the final full-data model. A separate ", METHOD_BOOTSTRAPS, "-resample fixed-score taxon bootstrap compared shortlist-composition stability across the six ranking approaches.\n\n",
"Because the weak-supervision grade and predictors share MR-derived inputs, model performance was interpreted as evidence-integration fidelity and auditability rather than independent confirmation of causality.\n"
)
writeLines(methods_text, file.path(OUT_ROOT, "METHODS_READY_TO_PASTE_AFTER_SUCCESSFUL_RUN.txt"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(OUT_ROOT, "sessionInfo.txt"))
writeLines(capture.output(XGB_PARAMS), file.path(OUT_ROOT, "XGBoost_parameters.txt"))

message("\n============================================================")
message("HYBRID FINAL PIPELINE COMPLETED")
message("Outputs: ", OUT_ROOT)
message("Open RUN_MANIFEST.csv before updating any manuscript number.")
