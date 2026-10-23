#!/usr/bin/env Rscript
# Project a compendium's samples onto a consensus model.
suppressPackageStartupMessages({
  library(bigstatsr)
  library(data.table)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "..", "archs4", "common.R"))
args <- parse_cli()
model <- list(Z = as.matrix(read.csv(required_arg(args, "z"), row.names = 1, check.names = FALSE)),
              L2 = jsonlite::read_json(required_arg(args, "manifest"))$L2)
stopifnot(is.numeric(model$L2), length(model$L2) == 1L)
n_cores <- as.integer(args$n_cores %||% 8L)

if (!is.null(args$preprocess_dir)) {
  prep <- args$preprocess_dir
  genes <- readRDS(file.path(prep, "metadata_filtered.rds"))$gene_symbols_thin
  samples <- filtered_sample_names(prep)
  fbm <- FBM(nrow = length(genes), ncol = length(samples), backingfile = file.path(prep, "fbm_filtered"),
             create_bk = FALSE, is_read_only = TRUE)
} else {
  genes <- as.character(readRDS(required_arg(args, "genes_rds")))
  samples <- as.character(readRDS(required_arg(args, "samples_rds")))
  fbm <- readRDS(required_arg(args, "fbm_rds"))
  if (!file.exists(fbm$backingfile)) fbm$backingfile <- normalizePath(required_arg(args, "fbm_backing"))
}
stopifnot(identical(dim(fbm), c(length(genes), length(samples))), all(genes %in% rownames(model$Z)))

B <- projectCLAMP(list(Z = model$Z[genes, , drop = FALSE], L2 = model$L2), newdata = fbm, ncores = n_cores)
stopifnot(all(is.finite(B)))
rownames(B) <- colnames(model$Z)
colnames(B) <- samples
fwrite(cbind(data.table(LV = rownames(B)), as.data.table(B)), required_arg(args, "out"))
message("Projected ", ncol(B), " samples onto ", nrow(B), " meta gene modules")
