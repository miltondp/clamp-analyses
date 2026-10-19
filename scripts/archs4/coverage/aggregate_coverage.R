#!/usr/bin/env Rscript
# Summarize pathway recovery across coverage models.
suppressPackageStartupMessages({
  library(data.table)
})

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (i in seq(1L, length(args), by = 2L)) {
    if (i == length(args)) stop("Missing value for ", args[[i]])
    out[[gsub("-", "_", substring(args[[i]], 3L), fixed = TRUE)]] <- args[[i + 1L]]
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", gsub("_", "-", name))
  value
}

read_summary <- function(path) {
  row <- fread(path)
  index <- regmatches(path, regexpr("/seed[0-9]+/", path))
  row[, summary_seed := seed]
  row[, seed := as.integer(gsub("[^0-9]", "", index))]
  row
}

relabel_databases <- function(dt, labels) {
  if (is.null(labels)) return(dt)
  pairs <- strsplit(strsplit(labels, ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
  map <- setNames(vapply(pairs, `[[`, character(1), 2L), vapply(pairs, `[[`, character(1), 1L))
  dt[database %chin% names(map), database_label := map[database]]
  dt
}

aggregate_coverage <- function(args) {
  ora_root <- required_arg(args, "ora_root")
  coverage_out <- required_arg(args, "coverage_out")
  cross_out <- required_arg(args, "cross_out")
  panel_out <- required_arg(args, "panel_out")
  databases <- args$databases

  summary_files <- list.files(ora_root, pattern = "^summary\\.csv$", recursive = TRUE, full.names = TRUE)
  canonical <- "/[^/]+/rs[0-9]+/seed[0-9]+/CLAMPbase/[^./]+/summary\\.csv$"
  summary_files <- grep(canonical, summary_files, value = TRUE)
  if (!length(summary_files)) stop("No published ORA summaries found below ", ora_root)
  complete_files <- file.path(dirname(summary_files), "complete")
  if (any(!file.exists(complete_files))) {
    stop("Incomplete ORA directories: ", paste(dirname(summary_files)[!file.exists(complete_files)], collapse = ", "))
  }

  coverage <- rbindlist(lapply(summary_files, read_summary), fill = TRUE)
  coverage <- relabel_databases(coverage, args$database_labels)
  if (!is.null(databases)) {
    databases <- trimws(strsplit(as.character(databases), ",", fixed = TRUE)[[1L]])
    coverage <- coverage[database %chin% databases]
  }
  required <- c(
    "dataset", "fraction", "seed", "model", "database", "database_label",
    "n_samples", "eligible_pathways", "recovered_pathways", "recovered_percent", "fdr"
  )
  if (!all(required %in% names(coverage))) {
    stop("ORA summaries are missing: ", paste(setdiff(required, names(coverage)), collapse = ", "))
  }
  if (!nrow(coverage)) stop("No ORA summaries remain after database filtering")
  if (anyDuplicated(coverage[, .(dataset, fraction, seed, model, database)])) {
    stop("Duplicate model/database ORA summaries found")
  }
  setorder(coverage, dataset, database, fraction, model, seed)

  cross <- coverage[fraction == 100]
  panel <- coverage[, .(
    mean_recovered = mean(recovered_pathways),
    sd_recovered = if (.N > 1L) sd(recovered_pathways) else NA_real_,
    mean_percent = mean(recovered_percent),
    eligible_pathways = unique(eligible_pathways),
    n_observations = .N,
    median_samples = as.integer(median(n_samples))
  ), by = .(dataset, fraction, model, database, database_label, fdr)]

  dir.create(dirname(coverage_out), recursive = TRUE, showWarnings = FALSE)
  fwrite(coverage, coverage_out)
  fwrite(cross, cross_out)
  fwrite(panel, panel_out)
}

aggregate_coverage(parse_cli())
