#!/usr/bin/env Rscript
# Compute the recount2 SVD and estimate model rank.
suppressPackageStartupMessages({
  library(bigstatsr)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
fbm_rds <- required_arg(args, "fbm_rds")
genes_rds <- required_arg(args, "genes_rds")
samples_rds <- required_arg(args, "samples_rds")
svd_out <- required_arg(args, "svd")
k_out <- required_arg(args, "k")
n_cores <- as.integer(args$n_cores %||% 1L)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

fbm <- readRDS(fbm_rds)
genes <- readRDS(genes_rds)
samples <- readRDS(samples_rds)
if (!identical(nrow(fbm), length(genes)) || !identical(ncol(fbm), length(samples))) {
  stop("FBM dimensions ", paste(dim(fbm), collapse = "x"),
       " do not match ", length(genes), " genes and ", length(samples), " samples")
}

svd_k <- as.integer(args$svd_k %||% floor((min(nrow(fbm), ncol(fbm)) - 1L) / 4))
message("Randomized SVD: ", nrow(fbm), " x ", ncol(fbm), ", k = ", svd_k, ", cores = ", n_cores)

svd_res <- with_single_threaded_blas(
  n_cores,
  big_randomSVD(fbm, k = svd_k, ncores = n_cores)
)
saveRDS(svd_res, svd_out)

eigenvalues <- sort(svd_res$d^2 / (ncol(fbm) - 1), decreasing = TRUE)
clamp_k <- PCAtools::chooseGavishDonoho(
  .dim = c(nrow(fbm), ncol(fbm)),
  var.explained = eigenvalues,
  noise = median(eigenvalues)
) * 2
message("Inferred CLAMP K = ", clamp_k)
saveRDS(clamp_k, k_out)
