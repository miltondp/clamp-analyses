#!/usr/bin/env Rscript
# Sample linear-scale expression for the log-transform comparison.
suppressPackageStartupMessages({
  library(bigstatsr)
  library(hdf5r)
  library(data.table)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "../archs4/common.R"))

args <- parse_cli()
n_samples <- as.integer(args$n_samples %||% 2000L)
seed <- as.integer(args$seed %||% 123L)
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)

fbm_sample <- function(bk, samples_rds, genes_rds) {
  ncol <- length(readRDS(samples_rds))
  X <- FBM(file.size(bk) / 8 / ncol, ncol, backingfile = sub("\\.bk$", "", bk), create_bk = FALSE)
  s <- big_apply(X, function(X, ind) {
    b <- X[, ind, drop = FALSE]
    b[is.na(b)] <- 0
    l <- log2(b + 1)
    cbind(rowSums(l), rowSums(l^2))
  }, a.combine = "plus", block.size = 500)
  m <- s[, 1] / ncol
  keep <- m >= mean_cutoff & s[, 2] / ncol - m^2 >= var_cutoff
  stopifnot(sum(keep) == length(readRDS(genes_rds)))
  set.seed(seed)
  x <- X[which(keep), sort(sample.int(ncol, n_samples))]
  x[is.na(x)] <- 0
  rownames(x) <- readRDS(genes_rds)
  x
}

archs4_sample <- function() {
  h5 <- H5File$new(required_arg(args, "archs4_h5"), mode = "r")
  on.exit(h5$close_all())
  sym <- h5[["/meta/genes/symbol"]]$read()
  meta <- readRDS(required_arg(args, "archs4_metadata"))
  keep_s <- which(h5[["/meta/samples/singlecellprobability"]]$read() < as.numeric(required_arg(args, "archs4_sc_max")))
  stopifnot(length(keep_s) == meta$n_samples)
  dset <- h5[["/data/expression"]]
  set.seed(seed)
  cand <- sample(keep_s)
  counts <- list()
  i <- 0
  while (length(counts) < n_samples) {
    i <- i + 1
    r <- tryCatch(dset[cand[i], ], error = function(e) NULL)
    if (!is.null(r)) counts[[length(counts) + 1]] <- r
  }
  message("ARCHS4: ", i - n_samples, " unreadable samples skipped")
  counts <- do.call(rbind, counts)
  gl <- readRDS(required_arg(args, "archs4_gene_lengths"))
  idx <- which(sym %in% names(gl))
  summed <- rowsum(t(counts[, idx]), group = sym[idx])
  tpm <- t(tpm_norm(t(summed), gl[rownames(summed)]))
  tpm[meta$gene_symbols_thin, ]
}

lin <- list(
  "GTEx (TPM)" = fbm_sample(required_arg(args, "gtex_bk"), required_arg(args, "gtex_samples"), required_arg(args, "gtex_genes")),
  "recount2 (RPKM)" = fbm_sample(required_arg(args, "recount2_bk"), required_arg(args, "recount2_samples"),
                                 required_arg(args, "recount2_genes")),
  "ARCHS4 (TPM)" = archs4_sample()
)
for (f in strsplit(required_arg(args, "pseudobulk"), ",")[[1]]) {
  lin[[paste0(basename(dirname(dirname(f))), " (CPM)")]] <- as.matrix(fread(f), rownames = 1)
}
ensure_dir(dirname(required_arg(args, "out")))
saveRDS(lin, required_arg(args, "out"))
print(data.table(dataset = names(lin), genes = sapply(lin, nrow), samples = sapply(lin, ncol)))
