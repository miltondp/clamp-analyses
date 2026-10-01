#!/usr/bin/env Rscript

# Fits CLAMPbase on a filtered ARCHS4 FBM.

suppressPackageStartupMessages({
  library(bigstatsr)
  library(dplyr)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
meta_path <- required_arg(args, "metadata")
samples_path <- required_arg(args, "samples")
fbm_backingfile <- required_arg(args, "fbm")
svd_path <- required_arg(args, "svd")
k_path <- required_arg(args, "k")
stage <- required_arg(args, "stage")
out_dir <- required_arg(args, "out_dir")
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

if (stage != "base") stop("--stage must be 'base', got: ", stage)

meta <- readRDS(meta_path)
genes <- meta$gene_symbols_thin
n_genes <- meta$n_genes_thin
n_samples <- meta$n_samples

sample_names <- readRDS(samples_path)
validate_sample_names(sample_names, n_samples)

fbm <- attach_fbm(fbm_backingfile, nrow = n_genes, ncol = n_samples)
svd_res <- readRDS(svd_path)
clamp_k <- readRDS(k_path)
message("CLAMPbase: ", n_genes, " genes x ", n_samples, " samples, K = ", clamp_k)

ensure_dir(out_dir)

base_res <- CLAMPbase(
  Y = fbm,
  svdres = svd_res,
  trace = TRUE,
  clamp_k = clamp_k
)
write_clamp_model(
  base_res,
  model_dir = file.path(out_dir, "CLAMPbase"),
  genes = genes,
  sample_names = sample_names,
  rds_path = file.path(out_dir, "CLAMPbase.rds")
)
message("Wrote CLAMPbase to ", out_dir)
