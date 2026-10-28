#!/usr/bin/env Rscript
# Historical recipe from commit 6935e139dc9ef6c4d17ea4ff17aa3d7f603051fe,
# nbs/01_model_building/04_archs4/01_archs4_preprocess.ipynb and config.R.
# Default analyses use the bundled snapshot; this writes only to a NEW path.
suppressPackageStartupMessages({library(hdf5r); library(biomaRt); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || args[[1L]] != "--h5" || args[[3L]] != "--out") {
    stop("Usage: Rscript scripts/inputs/build_archs4_gene_lengths.R --h5 <ARCHS4 v2.5 HDF5> --out <new.rds>")
}
if (file.exists(args[[4L]])) stop("Refusing to overwrite ", args[[4L]])
h5 <- H5File$new(args[[2L]], mode = "r")
gene_symbols <- h5[["/meta/genes/symbol"]]$read()
h5$close_all()
ensembl <- biomaRt::useEnsembl(biomart = "ensembl", version = 107)
mart <- biomaRt::useDataset("hsapiens_gene_ensembl", ensembl)
gene_info <- biomaRt::getBM(
    filters = "hgnc_symbol",
    attributes = c("ensembl_gene_id", "hgnc_symbol", "transcript_length"),
    values = gene_symbols, mart = mart
)
gene_lengths <- gene_info %>%
    group_by(hgnc_symbol, ensembl_gene_id) %>%
    summarize(gene_length_bp = max(transcript_length, na.rm = TRUE), .groups = "drop") %>%
    pull(gene_length_bp, name = hgnc_symbol)
stopifnot(length(gene_lengths) > 0L, all(is.finite(gene_lengths)), all(gene_lengths > 0))
dir.create(dirname(args[[4L]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(gene_lengths, args[[4L]])
message("Wrote ", length(gene_lengths), " Ensembl-107 lengths to ", args[[4L]],
        ". Compare values with the bundled snapshot before any scientific rebuild.")
