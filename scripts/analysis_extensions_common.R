# Shared helpers for the proxy-sensitivity and colocalisation extensions.
# This file has no analysis side effects when sourced.

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Install the required R package(s) first: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

assert_columns <- function(data, required, label) {
  absent <- setdiff(required, names(data))
  if (length(absent) > 0L) {
    stop(label, " is missing required column(s): ", paste(absent, collapse = ", "))
  }
  invisible(TRUE)
}

assert_unique <- function(data, keys, label) {
  duplicates <- data[, .N, by = keys][N > 1L]
  if (nrow(duplicates) > 0L) {
    stop(label, " contains duplicate key(s) for: ", paste(keys, collapse = ", "))
  }
  invisible(TRUE)
}

md5_text <- function(lines) {
  temp_file <- tempfile(pattern = "aixbio_fingerprint_")
  on.exit(unlink(temp_file), add = TRUE)
  writeLines(as.character(lines), temp_file, useBytes = TRUE)
  unname(tools::md5sum(temp_file))
}

fingerprint_run <- function(input_files, code_files, config_lines, package_names) {
  all_files <- unique(c(input_files, code_files))
  missing <- all_files[!file.exists(all_files)]
  if (length(missing) > 0L) {
    stop("Cannot fingerprint missing file(s): ", paste(missing, collapse = ", "))
  }
  file_hashes <- tools::md5sum(all_files)
  package_versions <- vapply(package_names, function(pkg) {
    as.character(utils::packageVersion(pkg))
  }, character(1))
  fingerprint_lines <- c(
    paste(names(file_hashes), unname(file_hashes), sep = "="),
    paste0("config=", config_lines),
    paste0("package=", names(package_versions), "@", package_versions),
    paste0("R=", getRversion())
  )
  substr(md5_text(fingerprint_lines), 1L, 16L)
}

acquire_run_lock <- function(lock_dir) {
  if (!dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not acquire the single-writer lock at ", lock_dir,
         ". If no other run is active, inspect and remove the stale lock manually.")
  }
  function() unlink(lock_dir, recursive = TRUE, force = TRUE)
}

atomic_write_csv <- function(data, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing output: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = ".csv")
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  data.table::fwrite(data, temp_path, na = "NA")
  check <- data.table::fread(temp_path, na.strings = "NA")
  if (nrow(check) != nrow(data) || !identical(names(check), names(data))) {
    stop("Read-back validation failed for temporary CSV: ", path)
  }
  if (!file.rename(temp_path, path)) stop("Could not publish output file: ", path)
  invisible(path)
}

atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing checkpoint: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  saveRDS(object, temp_path, version = 3)
  check <- readRDS(temp_path)
  if (!identical(names(check), names(object))) {
    stop("Read-back validation failed for temporary RDS checkpoint: ", path)
  }
  if (!file.rename(temp_path, path)) stop("Could not publish checkpoint: ", path)
  invisible(path)
}

atomic_write_lines <- function(lines, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = ".txt")
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  writeLines(lines, temp_path, useBytes = TRUE)
  check <- readLines(temp_path, warn = FALSE)
  if (!identical(check, as.character(lines))) stop("Read-back validation failed: ", path)
  if (!file.rename(temp_path, path)) stop("Could not publish file: ", path)
  invisible(path)
}

validate_completed_manifest <- function(final_dir, expected_run_id) {
  manifest_path <- file.path(final_dir, "COMPLETE_MANIFEST.csv")
  if (!file.exists(manifest_path)) {
    stop("A final-output folder exists without a completion manifest. It has been preserved at ",
         final_dir, ". Inspect it and rerun with a new input/config fingerprint.")
  }
  manifest <- data.table::fread(manifest_path)
  assert_columns(manifest, c("run_id", "file", "md5"), "completion manifest")
  if (nrow(manifest) < 1L || any(manifest$run_id != expected_run_id)) {
    stop("Completion manifest fingerprint does not match this run. Existing outputs were preserved.")
  }
  files <- file.path(final_dir, manifest$file)
  if (any(!file.exists(files))) stop("A manifested output is missing. Existing outputs were preserved.")
  observed <- unname(tools::md5sum(files))
  if (!identical(observed, as.character(manifest$md5))) {
    stop("Output checksum mismatch in completed run. Existing outputs were preserved.")
  }
  TRUE
}

publish_final_directory <- function(staging_dir, final_dir) {
  if (dir.exists(final_dir)) stop("Refusing to overwrite existing final directory: ", final_dir)
  if (!file.rename(staging_dir, final_dir)) {
    stop("Could not publish final directory. Staging files were preserved at ", staging_dir)
  }
  invisible(final_dir)
}
