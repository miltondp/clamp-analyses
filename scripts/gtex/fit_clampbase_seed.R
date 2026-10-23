#!/usr/bin/env Rscript
# Fit an independent CLAMPbase replicate of GTEx or another FBM compendium.
suppressPackageStartupMessages({
  library(bigstatsr)
  library(rsvd)
  library(jsonlite)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "..", "archs4", "common.R"))

args <- parse_cli()
fbm_rds <- required_arg(args, "fbm_rds")
backing <- required_arg(args, "backing")
genes <- as.character(readRDS(required_arg(args, "genes_rds")))
samples <- as.character(readRDS(required_arg(args, "samples_rds")))
out_dir <- required_arg(args, "out_dir")
seed <- as.integer(required_arg(args, "seed"))
seed_index <- as.integer(required_arg(args, "seed_index"))
expected_k <- as.integer(args$expected_k %||% NA_integer_)
dataset <- args$dataset %||% "gtex"
started <- Sys.time()

if (file.exists(file.path(out_dir, "CLAMPbase.rds"))) stop("Refusing to overwrite ", out_dir)
dir.create(file.path(out_dir, "CLAMPbase"), recursive = TRUE, showWarnings = FALSE)


template <- readRDS(fbm_rds)
Y <- FBM(nrow = nrow(template), ncol = ncol(template),
         backingfile = sub("\\.bk$", "", backing), create_bk = FALSE)
stopifnot(identical(nrow(Y), length(genes)), identical(ncol(Y), length(samples)))


set.seed(seed)
svd_k <- floor((min(nrow(Y), ncol(Y)) - 1) / 4)
message(dataset, " randomized SVD: ", nrow(Y), " x ", ncol(Y), ", k = ", svd_k, ", seed = ", seed)
svd_res <- rsvd(Y[], k = svd_k)
saveRDS(svd_res, file.path(out_dir, "svd.rds"))
eigenvalues <- sort(svd_res$d^2 / (ncol(Y) - 1), decreasing = TRUE)
clamp_k <- PCAtools::chooseGavishDonoho(
  .dim = c(nrow(Y), ncol(Y)), var.explained = eigenvalues, noise = median(eigenvalues)
) * 2
saveRDS(clamp_k, file.path(out_dir, "CLAMP_K.rds"))
message("Inferred CLAMP K = ", clamp_k)
if (!is.na(expected_k) && clamp_k != expected_k) {
  warning("Inferred K ", clamp_k, " differs from the published ", dataset, " K ", expected_k)
}

set.seed(seed)
base_res <- CLAMPbase(Y = Y, svdres = svd_res, trace = TRUE, clamp_k = clamp_k, rseed = seed)
stopifnot(identical(dim(base_res$Z), c(length(genes), as.integer(clamp_k))),
          all(is.finite(as.matrix(base_res$Z))), all(is.finite(as.matrix(base_res$B))))
rownames(base_res$Z) <- genes
colnames(base_res$B) <- samples
saveRDS(base_res, file.path(out_dir, "CLAMPbase.rds"))
write.csv(base_res$Z, file.path(out_dir, "CLAMPbase", "Z.csv"))
write.csv(base_res$B, file.path(out_dir, "CLAMPbase", "B.csv"))

finished <- Sys.time()
write_json(list(
  schema_version = 1L, status = "complete", model = "CLAMPbase", dataset = dataset,
  fraction = 100L, seed_index = seed_index, rng_seed = seed, random_start = TRUE,
  svd_k = svd_k, latent_variables = as.integer(clamp_k), genes = length(genes), samples = length(samples),
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  finished_utc = format(finished, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
  package_versions = list(R = as.character(getRversion()), CLAMP = as.character(packageVersion("CLAMP")))
), file.path(out_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
file.create(file.path(out_dir, "complete"))
message("Wrote ", dataset, " CLAMPbase replicate to ", out_dir)
