#!/usr/bin/env Rscript
# OneK1K cell-type score of one fit: per cell type, the Pearson r of its one-to-one matched LV
# (score_per_ct() in scripts/pseudobulk/common.R, as in the Fig. 2A benchmark notebook
# nbs/03_model_biology/00_pseudobulk/00_benchmark.ipynb) with the true cell-type fractions.
# Reads one B.csv, so seed and k variants get the same score as the production fits without
# running benchmark_pseudobulk (which expands over all 6 datasets).
#
# Usage: celltype_score.R --b B.csv --truth truthFrac_v0.csv --dataset PBMC_1k1k
#          --strategy main --model CLAMPfull --out per_cell_type.csv [--exclude "A,B"]
# Output: one row per cell type (strategy, model, cell_type, cor), plus a "mean" row.
suppressPackageStartupMessages(library(data.table))
source(here::here("scripts", "pseudobulk", "common.R"))

args <- parse_cli()
dataset <- required_arg(args, "dataset")
exclude <- args$exclude %||% ""
exclude <- if (isTRUE(exclude) || !nzchar(exclude)) character() else strsplit(exclude, ",", fixed = TRUE)[[1]]

truth <- filter_analysis_cell_types(
  read_csv_matrix(required_arg(args, "truth")), dataset, setNames(list(exclude), dataset)
)
fac <- read_B(required_arg(args, "b"))
message("B: ", nrow(fac), " samples x ", ncol(fac), " LVs; truth: ", ncol(truth), " cell types")
shared <- intersect(rownames(fac), rownames(truth))
if (length(shared) != nrow(truth)) {
  stop(nrow(truth) - length(shared), " truth samples are missing from B")
}

cors <- score_per_ct(fac, truth, colnames(truth))
out <- data.table(
  strategy = required_arg(args, "strategy"), model = required_arg(args, "model"),
  cell_type = c(names(cors), "mean"), cor = c(unname(cors), mean(cors, na.rm = TRUE))
)
path <- required_arg(args, "out")
ensure_parent(path)
fwrite(out, path)
message("mean r = ", signif(mean(cors, na.rm = TRUE), 4))
