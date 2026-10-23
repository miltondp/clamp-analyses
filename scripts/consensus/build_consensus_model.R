#!/usr/bin/env Rscript
# Describe a consensus CLAMP model built from meta gene modules clustered across seed models.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "..", "archs4", "common.R"))
args <- parse_cli()
out_dir <- required_arg(args, "out_dir")
reference_rds <- required_arg(args, "reference_model")
reference_manifest <- read_json(required_arg(args, "reference_manifest"))

mz <- fread(required_arg(args, "metaz"), nrows = 0L)
genes <- nrow(fread(required_arg(args, "metaz"), select = 1L))
summ <- fread(required_arg(args, "summary"))
stopifnot(identical(names(mz)[-1L], summ$metaLV))
ref_col <- grep("^S[0-9]+_lvs$", names(summ), value = TRUE)
stopifnot(length(ref_col) == 1L)
mapping <- data.table(metaLV = summ$metaLV, cluster = summ$cluster, n_models = summ$n_models,
                      n_lvs = summ$n_lvs, min_cos = summ$min_cos, mean_cos = summ$mean_cos,
                      reference_lvs = summ[[ref_col]], members = summ$members)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(mapping, file.path(out_dir, "mgm_mapping.csv"))
write_json(list(
  schema_version = 1L, status = "complete", model = required_arg(args, "model_name"),
  dataset = required_arg(args, "dataset"), reference_model = normalizePath(reference_rds),
  L2 = readRDS(reference_rds)$L2, models = strsplit(required_arg(args, "models"), ",")[[1L]],
  min_cos = as.numeric(required_arg(args, "min_cos")), min_models = as.integer(required_arg(args, "min_models")),
  genes = genes, samples = reference_manifest$samples, latent_variables = nrow(mapping)
), file.path(out_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
message("Consensus model: ", genes, " genes x ", nrow(mapping), " meta gene modules")
