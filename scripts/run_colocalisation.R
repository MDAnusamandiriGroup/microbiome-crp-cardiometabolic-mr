# Targeted ABF and optional SuSiE colocalisation.
# Run from the repository root: Rscript --vanilla scripts/run_colocalisation.R
# All summary statistics and LD matrices must be prepared locally beforehand.

source("scripts/analysis_extensions_common.R")

coloc_manifest_required <- c(
  "region_id", "taxon", "chromosome", "lead_position", "window_bp",
  "trait1_file", "trait2_file", "n1", "n2", "sdY1", "sdY2",
  "genome_build", "ld_rds"
)
coloc_summary_required <- c(
  "SNP", "chr", "pos", "beta", "se", "effect_allele",
  "other_allele", "eaf", "genome_build"
)

read_coloc_sumstats <- function(path, label, chromosome, lead_position, window_bp, genome_build) {
  data <- data.table::fread(path, na.strings = c("", "NA"))
  assert_columns(data, coloc_summary_required, label)
  data[, `:=`(
    SNP = trimws(as.character(SNP)),
    chr = as.character(chr),
    pos = suppressWarnings(as.numeric(pos)),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    effect_allele = toupper(trimws(as.character(effect_allele))),
    other_allele = toupper(trimws(as.character(other_allele))),
    eaf = suppressWarnings(as.numeric(eaf)),
    genome_build = trimws(as.character(genome_build))
  )]
  assert_unique(data, "SNP", label)
  if (anyNA(data$SNP) || any(!nzchar(data$SNP)) || anyNA(data$chr) ||
      anyNA(data$effect_allele) || anyNA(data$other_allele) ||
      anyNA(data$genome_build) || any(!nzchar(data$genome_build)) ||
      any(data$genome_build != genome_build)) {
    stop(label, " has missing variant identifiers, chromosome, alleles, or mismatched genome build.")
  }
  if (anyNA(data$pos) || any(!is.finite(data$pos)) || any(!is.finite(data$beta)) ||
      any(!is.finite(data$se)) || any(data$se <= 0)) {
    stop(label, " contains missing positions, non-finite beta, or invalid SE values.")
  }
  if (any(!is.na(data$eaf) & (data$eaf < 0 | data$eaf > 1))) {
    stop(label, " EAF values must be proportions between 0 and 1.")
  }
  chr_text <- as.character(chromosome)
  data <- data[chr == chr_text & abs(pos - lead_position) <= window_bp]
  data
}

harmonise_coloc_region <- function(trait1, trait2, region_id) {
  rename_suffix <- function(x, suffix) {
    data.table::setnames(x, setdiff(names(x), "SNP"),
      paste0(setdiff(names(x), "SNP"), suffix))
    x
  }
  first <- rename_suffix(data.table::copy(trait1), "1")
  second <- rename_suffix(data.table::copy(trait2), "2")
  joined <- merge(first, second, by = "SNP", all = FALSE, sort = TRUE)
  if (nrow(joined) == 0L) stop("No overlapping variants in region ", region_id, ".")
  if (any(joined$chr1 != joined$chr2) || any(joined$pos1 != joined$pos2)) {
    stop("Shared variant coordinates disagree between traits in region ", region_id,
         ". Check genome build and variant identifiers.")
  }

  same <- joined$effect_allele1 == joined$effect_allele2 &
    joined$other_allele1 == joined$other_allele2
  reversed <- joined$effect_allele1 == joined$other_allele2 &
    joined$other_allele1 == joined$effect_allele2
  same[is.na(same)] <- FALSE
  reversed[is.na(reversed)] <- FALSE
  aligned <- same | reversed
  audit <- data.table::data.table(
    region_id = region_id, SNP = joined$SNP,
    position = joined$pos1,
    allele_status = ifelse(same, "same", ifelse(reversed, "reversed", "unmatched")),
    included = aligned
  )
  kept <- joined[aligned]
  if (nrow(kept) < 2L) {
    stop("Fewer than two allele-aligned shared variants in region ", region_id,
         ". No coloc analysis was run.")
  }
  reverse_rows <- reversed[aligned]
  beta2 <- kept$beta2
  eaf2 <- kept$eaf2
  beta2[reverse_rows] <- -beta2[reverse_rows]
  eaf2[reverse_rows] <- 1 - eaf2[reverse_rows]
  list(
    data = data.table::data.table(
      SNP = kept$SNP, position = kept$pos1,
      beta1 = kept$beta1, se1 = kept$se1, eaf1 = kept$eaf1,
      beta2 = beta2, se2 = kept$se2, eaf2 = eaf2
    ),
    audit = audit
  )
}

make_coloc_dataset <- function(data, beta_col, se_col, eaf_col, sample_size, sd_y, label) {
  beta <- data[[beta_col]]
  se <- data[[se_col]]
  eaf <- data[[eaf_col]]
  if (!is.finite(sample_size) || sample_size <= 0) stop(label, " sample size must be positive.")
  dataset <- list(
    beta = beta,
    varbeta = se^2,
    snp = data$SNP,
    position = data$position,
    type = "quant",
    N = as.integer(sample_size)
  )
  if (is.finite(sd_y) && sd_y > 0) {
    dataset$sdY <- as.numeric(sd_y)
  } else {
    maf <- pmin(eaf, 1 - eaf)
    if (any(!is.finite(maf))) {
      stop(label, " requires finite EAF values when sdY is not supplied.")
    }
    dataset$MAF <- maf
  }
  dataset
}

extract_coloc_summary <- function(result, method, region_id, taxon, p12_value) {
  summary_object <- result$summary
  if (is.null(summary_object)) stop("Coloc returned no summary for region ", region_id, ".")
  if (is.data.frame(summary_object)) {
    summary_table <- data.table::as.data.table(summary_object)
  } else {
    summary_table <- data.table::as.data.table(as.list(unlist(summary_object, use.names = TRUE)))
  }
  summary_table[, `:=`(
    region_id = region_id, taxon = taxon, method = method,
    p1 = 1e-4, p2 = 1e-4, p12 = p12_value
  )]
  data.table::setcolorder(summary_table, c(
    "region_id", "taxon", "method", "p1", "p2", "p12",
    setdiff(names(summary_table), c("region_id", "taxon", "method", "p1", "p2", "p12"))
  ))
  summary_table
}

analyse_one_coloc_region <- function(region, input_dir, p12_values = c(1e-6, 1e-5, 1e-4)) {
  region_id <- as.character(region$region_id)
  trait1_path <- file.path(input_dir, as.character(region$trait1_file))
  trait2_path <- file.path(input_dir, as.character(region$trait2_file))
  if (!file.exists(trait1_path) || !file.exists(trait2_path)) {
    stop("Missing summary-statistic file(s) for region ", region_id, ".")
  }
  trait1 <- read_coloc_sumstats(
    trait1_path, "trait 1 summary statistics", region$chromosome,
    as.numeric(region$lead_position), as.numeric(region$window_bp),
    as.character(region$genome_build)
  )
  trait2 <- read_coloc_sumstats(
    trait2_path, "trait 2 summary statistics", region$chromosome,
    as.numeric(region$lead_position), as.numeric(region$window_bp),
    as.character(region$genome_build)
  )
  harmonised <- harmonise_coloc_region(trait1, trait2, region_id)
  coloc_data <- harmonised$data
  n1 <- as.numeric(region$n1); n2 <- as.numeric(region$n2)
  sd_y1 <- suppressWarnings(as.numeric(region$sdY1))
  sd_y2 <- suppressWarnings(as.numeric(region$sdY2))
  dataset1 <- make_coloc_dataset(coloc_data, "beta1", "se1", "eaf1", n1, sd_y1, "trait 1")
  dataset2 <- make_coloc_dataset(coloc_data, "beta2", "se2", "eaf2", n2, sd_y2, "trait 2")
  coloc::check_dataset(dataset1)
  coloc::check_dataset(dataset2)

  abf_rows <- vector("list", length(p12_values))
  for (i in seq_along(p12_values)) {
    p12_value <- p12_values[i]
    fitted <- coloc::coloc.abf(
      dataset1 = dataset1, dataset2 = dataset2,
      p1 = 1e-4, p2 = 1e-4, p12 = p12_value
    )
    abf_rows[[i]] <- extract_coloc_summary(fitted, "ABF", region_id,
      as.character(region$taxon), p12_value)
  }
  abf_summary <- data.table::rbindlist(abf_rows, fill = TRUE)

  susie_summary <- data.table::data.table()
  ld_file <- as.character(region$ld_rds)
  if (!is.na(ld_file) && nzchar(trimws(ld_file))) {
    if (!requireNamespace("susieR", quietly = TRUE)) {
      stop("Region ", region_id, " requests SuSiE; install the susieR package first.")
    }
    ld_path <- file.path(input_dir, ld_file)
    if (!file.exists(ld_path)) stop("Missing LD RDS for region ", region_id, ": ", ld_path)
    ld <- readRDS(ld_path)
    if (!is.matrix(ld) || !is.numeric(ld) || is.null(rownames(ld)) || is.null(colnames(ld))) {
      stop("LD RDS for region ", region_id, " must be a named numeric matrix.")
    }
    if (anyDuplicated(rownames(ld)) || anyDuplicated(colnames(ld)) ||
        !setequal(rownames(ld), colnames(ld))) {
      stop("LD matrix row/column names must be unique and contain the same SNP IDs.")
    }
    snps <- coloc_data$SNP
    if (!all(snps %in% rownames(ld)) || !all(snps %in% colnames(ld))) {
      stop("LD matrix does not contain every allele-aligned SNP for region ", region_id, ".")
    }
    ld <- ld[snps, snps, drop = FALSE]
    if (any(!is.finite(ld)) || max(abs(ld - t(ld))) > 1e-6 ||
        any(abs(diag(ld) - 1) > 1e-4) || any(abs(ld) > 1.0001)) {
      stop("LD matrix failed symmetry, diagonal, or range checks for region ", region_id, ".")
    }
    # LD alleles must match the forward alleles in the summary-statistic files.
    # If the LD source codes the opposite allele, align its signed correlations before use.
    dataset1$LD <- ld
    dataset2$LD <- ld
    coloc::check_dataset(dataset1, req = "LD")
    coloc::check_dataset(dataset2, req = "LD")
    susie1 <- coloc::runsusie(dataset1, suffix = paste0(region_id, "_trait1"))
    susie2 <- coloc::runsusie(dataset2, suffix = paste0(region_id, "_trait2"))
    susie_rows <- vector("list", length(p12_values))
    for (i in seq_along(p12_values)) {
      p12_value <- p12_values[i]
      fitted <- coloc::coloc.susie(
        dataset1 = susie1, dataset2 = susie2,
        p1 = 1e-4, p2 = 1e-4, p12 = p12_value
      )
      susie_rows[[i]] <- extract_coloc_summary(fitted, "SuSiE", region_id,
        as.character(region$taxon), p12_value)
    }
    susie_summary <- data.table::rbindlist(susie_rows, fill = TRUE)
  }

  list(
    region_id = region_id,
    taxon = as.character(region$taxon),
    abf = abf_summary,
    susie = susie_summary,
    audit = harmonised$audit,
    region_summary = data.table::data.table(
      region_id = region_id, taxon = as.character(region$taxon),
      chromosome = as.character(region$chromosome),
      lead_position = as.numeric(region$lead_position),
      window_bp = as.numeric(region$window_bp),
      n_trait1_in_window = nrow(trait1), n_trait2_in_window = nrow(trait2),
      n_shared_allele_aligned = nrow(coloc_data),
      n_shared_allele_excluded = sum(!harmonised$audit$included),
      susie_requested = !is.na(ld_file) && nzchar(trimws(ld_file)),
      status = "completed"
    )
  )
}

run_colocalisation <- function(project_root = getwd()) {
  assert_packages(c("data.table", "coloc"))
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(project_root, "pipeline.R"))) {
    stop("Run this command from the repository root, or set project_root to that folder.")
  }
  input_dir <- Sys.getenv("AIXBIO_COLOC_INPUT_DIR", unset = "data/private/colocalisation")
  output_root <- Sys.getenv("AIXBIO_COLOC_OUTPUT_DIR", unset = "results/colocalisation/runs")
  input_dir <- file.path(project_root, input_dir)
  output_root <- file.path(project_root, output_root)
  manifest_path <- file.path(input_dir, "region_manifest.csv")
  if (!file.exists(manifest_path)) {
    stop("Missing region_manifest.csv. See scripts/ANALYSIS_EXTENSIONS.md for the required dense regional inputs.")
  }
  regions <- data.table::fread(manifest_path, na.strings = c("", "NA"))
  assert_columns(regions, coloc_manifest_required, "region_manifest.csv")
  assert_unique(regions, "region_id", "colocalisation region manifest")
  regions[, `:=`(
    region_id = trimws(as.character(region_id)),
    taxon = trimws(as.character(taxon)),
    chromosome = as.character(chromosome),
    lead_position = suppressWarnings(as.numeric(lead_position)),
    window_bp = suppressWarnings(as.numeric(window_bp)),
    n1 = suppressWarnings(as.numeric(n1)), n2 = suppressWarnings(as.numeric(n2)),
    sdY1 = suppressWarnings(as.numeric(sdY1)), sdY2 = suppressWarnings(as.numeric(sdY2)),
    genome_build = trimws(as.character(genome_build)),
    trait1_file = trimws(as.character(trait1_file)),
    trait2_file = trimws(as.character(trait2_file)),
    ld_rds = trimws(as.character(ld_rds))
  )]
  if (nrow(regions) == 0L || anyNA(regions$region_id) || any(!nzchar(regions$region_id))) {
    stop("The region manifest must contain at least one region with a non-empty region_id.")
  }
  if (anyNA(regions$taxon) || any(!nzchar(regions$taxon))) {
    stop("Every region manifest row must name a taxon.")
  }
  if (any(!is.finite(regions$lead_position)) || any(!is.finite(regions$window_bp)) ||
      any(regions$window_bp <= 0) || any(!is.finite(regions$n1)) || any(!is.finite(regions$n2)) ||
      any(regions$n1 <= 0) || any(regions$n2 <= 0)) {
    stop("Manifest positions, windows, and sample sizes must be finite; windows and sample sizes must be positive.")
  }
  if (anyNA(regions$chromosome) || anyNA(regions$genome_build) ||
      any(!nzchar(regions$genome_build)) || anyNA(regions$trait1_file) ||
      anyNA(regions$trait2_file) || any(!nzchar(regions$trait1_file)) ||
      any(!nzchar(regions$trait2_file))) {
    stop("Region manifest needs chromosome, genome_build, and both summary-statistic file paths for every row.")
  }
  if (any(!is.na(regions$sdY1) & regions$sdY1 <= 0) || any(!is.na(regions$sdY2) & regions$sdY2 <= 0)) {
    stop("Optional sdY1/sdY2 values must be positive when supplied.")
  }
  p12_values <- c(1e-6, 1e-5, 1e-4)
  summary_files <- unique(c(
    manifest_path,
    file.path(input_dir, regions$trait1_file),
    file.path(input_dir, regions$trait2_file),
    file.path(input_dir, regions$ld_rds[!is.na(regions$ld_rds) & nzchar(regions$ld_rds)])
  ))
  missing <- summary_files[!file.exists(summary_files)]
  if (length(missing) > 0L) stop("Missing colocalisation input file(s): ", paste(missing, collapse = ", "))
  code_files <- file.path(project_root, c(
    "scripts/run_colocalisation.R", "scripts/analysis_extensions_common.R"
  ))
  config <- c("p1=1e-4", "p2=1e-4", "p12=1e-6,1e-5,1e-4", "quantitative_traits_only=TRUE")
  packages <- if (any(!is.na(regions$ld_rds) & nzchar(regions$ld_rds))) c("data.table", "coloc", "susieR") else c("data.table", "coloc")
  run_id <- fingerprint_run(summary_files, code_files, config, packages)
  run_dir <- file.path(output_root, paste0("coloc_", run_id))
  checkpoint_dir <- file.path(run_dir, "checkpoints")
  failure_dir <- file.path(run_dir, "failures")
  final_dir <- file.path(run_dir, "final")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(final_dir)) {
    validate_completed_manifest(final_dir, run_id)
    message("Validated completed colocalisation run; reusing outputs at ", final_dir)
    return(invisible(final_dir))
  }
  unlock <- acquire_run_lock(file.path(run_dir, ".lock"))
  on.exit(unlock(), add = TRUE)

  unit_results <- vector("list", nrow(regions))
  for (i in seq_len(nrow(regions))) {
    region <- regions[i]
    unit_path <- file.path(checkpoint_dir, sprintf("region_%04d.rds", i))
    if (file.exists(unit_path)) {
      cached <- tryCatch(readRDS(unit_path), error = function(e) {
        stop("Unreadable region checkpoint preserved at ", unit_path, ": ", conditionMessage(e))
      })
      if (!identical(cached$region_id, as.character(region$region_id))) {
        stop("Checkpoint region ID does not match the manifest; preserved at ", unit_path)
      }
      unit_results[[i]] <- cached
      message("Reusing region checkpoint ", i, "/", nrow(regions), ": ", region$region_id)
      next
    }
    message("Colocalising region ", i, "/", nrow(regions), ": ", region$region_id)
    result <- tryCatch(
      analyse_one_coloc_region(region, input_dir, p12_values),
      error = function(e) {
        failure_path <- file.path(failure_dir, sprintf("region_%04d_%s.txt", i, format(Sys.time(), "%Y%m%dT%H%M%S")))
        atomic_write_lines(c(paste0("region_id=", region$region_id), paste0("error=", conditionMessage(e))), failure_path)
        stop("Region failed and was not checkpointed: ", region$region_id, " — ", conditionMessage(e))
      }
    )
    atomic_save_rds(result, unit_path)
    unit_results[[i]] <- result
  }

  abf <- data.table::rbindlist(lapply(unit_results, `[[`, "abf"), fill = TRUE)
  susie <- data.table::rbindlist(lapply(unit_results, `[[`, "susie"), fill = TRUE)
  if (ncol(susie) == 0L) {
    susie <- data.table::data.table(
      region_id = character(), taxon = character(), method = character(),
      p1 = numeric(), p2 = numeric(), p12 = numeric(), nsnps = integer(),
      hit1 = character(), hit2 = character(),
      PP.H0.abf = numeric(), PP.H1.abf = numeric(), PP.H2.abf = numeric(),
      PP.H3.abf = numeric(), PP.H4.abf = numeric()
    )
  }
  audit <- data.table::rbindlist(lapply(unit_results, `[[`, "audit"), fill = TRUE)
  region_summary <- data.table::rbindlist(lapply(unit_results, `[[`, "region_summary"), fill = TRUE)
  stage_dir <- tempfile(pattern = ".publish_coloc_", tmpdir = run_dir)
  dir.create(stage_dir, recursive = FALSE, showWarnings = FALSE)
  if (!dir.exists(stage_dir)) stop("Could not create final staging directory.")
  stage_files <- c("coloc_abf_summary.csv", "coloc_susie_summary.csv", "variant_alignment_audit.csv", "region_summary.csv")
  atomic_write_csv(abf, file.path(stage_dir, stage_files[1]))
  atomic_write_csv(susie, file.path(stage_dir, stage_files[2]))
  atomic_write_csv(audit, file.path(stage_dir, stage_files[3]))
  atomic_write_csv(region_summary, file.path(stage_dir, stage_files[4]))
  output_hashes <- tools::md5sum(file.path(stage_dir, stage_files))
  manifest <- data.table::data.table(
    run_id = run_id, file = stage_files, md5 = unname(output_hashes),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    config = paste(config, collapse = ";"),
    coloc = as.character(utils::packageVersion("coloc")),
    data.table = as.character(utils::packageVersion("data.table")),
    susieR = if (requireNamespace("susieR", quietly = TRUE)) as.character(utils::packageVersion("susieR")) else NA_character_,
    R_version = as.character(getRversion())
  )
  atomic_write_csv(manifest, file.path(stage_dir, "COMPLETE_MANIFEST.csv"))
  publish_final_directory(stage_dir, final_dir)
  message("Colocalisation complete. Outputs: ", final_dir)
  invisible(final_dir)
}

coloc_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(coloc_script_arg) == 1L && grepl("run_colocalisation\\.R$", coloc_script_arg)) {
  run_colocalisation()
}
