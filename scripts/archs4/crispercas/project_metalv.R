# Project selected metaLVs into ARCHS4 samples.
suppressPackageStartupMessages({library(bigstatsr); library(data.table); library(CLAMP); library(hdf5r)})

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (i in seq(1L, length(args), by = 2L)) out[[gsub("-", "_", substring(args[[i]], 3L), fixed = TRUE)]] <- args[[i + 1L]]
  out
}
opt <- parse_cli()
n_cores <- as.integer(if (is.null(opt$n_cores)) 8L else opt$n_cores)

lvs <- fread(opt$selected)$metaLV
mz <- as.matrix(read.csv(opt$metaz, row.names = 1, check.names = FALSE))
stopifnot(!anyDuplicated(lvs), all(lvs %in% colnames(mz)))
l2 <- readRDS(opt$model)$L2
genes <- readRDS(file.path(opt$preprocess_dir, "metadata_filtered.rds"))$gene_symbols_thin
smeta <- fread(opt$sample_metadata, select = c("sample", "h5_index", "series_id"))
stopifnot(all(genes %in% rownames(mz)))
fbm <- FBM(nrow = length(genes), ncol = nrow(smeta), backingfile = file.path(opt$preprocess_dir, "fbm_filtered"),
           create_bk = FALSE, is_read_only = TRUE)
# Compute selected rows of the full-model ridge projector. Subsetting Z before
# solving changes the coefficients by reallocating signal from omitted MGMs.
z <- mz[genes, , drop = FALSE]
gram <- crossprod(z) + diag(l2, ncol(z))
rhs <- diag(ncol(z))[, match(lvs, colnames(z)), drop = FALSE]
weights <- z %*% solve(gram, rhs)
options(bigstatsr.check.parallel.blas = FALSE)
B <- t(big_cprodMat(fbm, weights, ncores = n_cores))
rownames(B) <- lvs

h5 <- H5File$new(opt$h5, mode = "r")
read_meta <- function(k) h5[[paste0("meta/samples/", k)]]$read()[smeta$h5_index]
ann <- data.table(sample = smeta$sample, series_id = smeta$series_id, source = read_meta("source_name_ch1"),
                  title = read_meta("title"), characteristics = read_meta("characteristics_ch1"))
h5$close_all()
fwrite(cbind(ann, as.data.table(t(B))), opt$out, sep = "\t")
message("projected ", ncol(B), " samples onto all ", ncol(z), " metaLVs; retained ", nrow(B),
        " selected metaLVs (L2 = ", round(l2, 4), ")")
