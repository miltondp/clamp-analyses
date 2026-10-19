#!/usr/bin/env Rscript
# Summarize trait recovery at a fixed model rank.
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

args <- parse_cli()
clampbase_dir <- required_arg(args, "clampbase_dir")
full_clampbase_dir <- required_arg(args, "full_clampbase_dir")
fdr_cutoff <- as.numeric(required_arg(args, "fdr"))
traits_out <- required_arg(args, "traits_out")


scan_dir <- function(dir, model_label) {
  paths <- list.files(
    dir, pattern = "^sat_rs[0-9]+_k1728_seed[0-9]+(_[A-Za-z0-9]+)?\\.tsv\\.gz$", full.names = TRUE
  )
  if (!length(paths)) stop("No K=1728 trait GLS summaries found under ", dir)
  rbindlist(lapply(paths, function(path) {
    m <- regmatches(
      basename(path),
      regexec("^sat_rs([0-9]+)_k1728_seed([0-9]+)(?:_[A-Za-z0-9]+)?\\.tsv\\.gz$", basename(path))
    )[[1]]
    d <- fread(path, select = c("phenotype", "fdr"))
    data.table(
      fraction = as.integer(m[2]),
      seed = as.integer(m[3]),
      model = model_label,
      k = 1728L,
      eligible_traits = uniqueN(d$phenotype),
      recovered_traits = uniqueN(d[fdr < fdr_cutoff, phenotype]),
      fdr = fdr_cutoff
    )
  }))
}

scan_full_dir <- function(dir, model_label) {
  paths <- list.files(dir, pattern = "^cov_rs100_seed[0-9]+\\.tsv\\.gz$", full.names = TRUE)
  if (!length(paths)) stop("No 100% coverage trait GLS summaries found under ", dir)
  rbindlist(lapply(paths, function(path) {
    seed <- as.integer(sub("^cov_rs100_seed([0-9]+)\\.tsv\\.gz$", "\\1", basename(path)))
    d <- fread(path, select = c("phenotype", "lv", "fdr"))
    if (uniqueN(d$lv) != 1728L) stop("Expected 1728 LVs in ", path, ", found ", uniqueN(d$lv))
    data.table(
      fraction = 100L,
      seed = seed,
      model = model_label,
      k = 1728L,
      eligible_traits = uniqueN(d$phenotype),
      recovered_traits = uniqueN(d[fdr < fdr_cutoff, phenotype]),
      fdr = fdr_cutoff
    )
  }))
}

traits <- rbindlist(list(
  scan_dir(clampbase_dir, "CLAMPbase"),
  scan_full_dir(full_clampbase_dir, "CLAMPbase")
))
if (anyDuplicated(traits[, .(fraction, seed, model)])) {
  stop("Duplicate fraction/seed/model saturation trait GLS summaries found")
}
setorder(traits, model, fraction, seed)

dir.create(dirname(traits_out), recursive = TRUE, showWarnings = FALSE)
fwrite(traits, traits_out)
