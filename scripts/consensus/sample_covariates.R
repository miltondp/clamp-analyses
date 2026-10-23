#!/usr/bin/env Rscript
# Build per-sample technical covariates for the GTEx and recount2 technical-noise screens.
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "..", "archs4", "common.R"))
args <- parse_cli()
dataset <- required_arg(args, "dataset")
samples <- as.character(readRDS(required_arg(args, "samples_rds")))
out_path <- required_arg(args, "out")

if (dataset == "gtex") {
  attr <- fread(required_arg(args, "attributes"), sep = "\t", quote = "")
  keep <- c("SAMPID", "SMCENTER", "SMNABTCHT", "SMGEBTCHT", "SMGEBTCH", "SMRIN", "SMTSISCH",
            "SMMAPRT", "SMEXNCRT", "SMRRNART", "SMRDTTL", "SMGNSDTC")
  attr <- attr[, ..keep]
  setnames(attr, "SAMPID", "sample")
  out <- attr[match(samples, sample)]
  out[, sample := samples]
  out[, log10_SMRDTTL := log10(as.numeric(SMRDTTL) + 1)]
  out[, SMRDTTL := NULL]
} else if (dataset == "recount2") {
  files <- list.files(required_arg(args, "rse_dir"), pattern = "^rse_gene\\.Rdata$", recursive = TRUE, full.names = TRUE)
  fields <- c("project", "run", "paired_end", "biosample_submission_date", "mapped_read_count",
              "avg_read_length", "proportion_of_reads_reported_by_sra_downloaded", "auc")
  read_one <- function(f) {
    suppressPackageStartupMessages(library(SummarizedExperiment))
    e <- new.env()
    load(f, envir = e)
    cd <- as.data.frame(colData(get(ls(e)[1L], envir = e)))
    as.data.table(cd[, intersect(fields, names(cd)), drop = FALSE])
  }
  cd <- rbindlist(mclapply(files, read_one, mc.cores = as.integer(args$n_cores %||% 8L)), fill = TRUE)
  cd[, sample := paste0(project, ".", run)]
  cd <- unique(cd, by = "sample")
  out <- cd[match(samples, sample)]
  out[, sample := samples]
  out[, `:=`(paired_end = as.character(paired_end),
             submission_year = substr(as.character(biosample_submission_date), 1L, 4L),
             log10_mapped_read_count = log10(as.numeric(mapped_read_count) + 1),
             log10_auc = log10(as.numeric(auc) + 1),
             avg_read_length = as.numeric(avg_read_length),
             proportion_downloaded = as.numeric(proportion_of_reads_reported_by_sra_downloaded))]
  out <- out[, .(sample, project, paired_end, submission_year, log10_mapped_read_count,
                 avg_read_length, proportion_downloaded, log10_auc)]
} else {
  stop("Unknown dataset: ", dataset)
}

matched <- sum(complete.cases(out[, !"sample"]))
message(dataset, ": ", matched, " of ", length(samples), " samples with complete covariates")
stopifnot(matched > 0.5 * length(samples))
fwrite(out, out_path, sep = "\t", na = "")
