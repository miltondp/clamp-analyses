#!/usr/bin/env Rscript
# Export Z.csv, summary.csv and manifest.json from one saved CLAMPfull fit (.rds).
# The files are written exactly as fit_clampfull.R writes them, so an export of a
# published model is byte-identical to its Z.csv and summary.csv. B is not exported.
#
#   Rscript export_run.R --rds <fit.rds> --model S123 --md5-file MD5SUM.txt --out-dir <dir>
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L) stop("Arguments must be --key value pairs")
opt <- setNames(as.list(args[c(FALSE, TRUE)]), sub("^--", "", args[c(TRUE, FALSE)]))
for (key in c("rds", "model", "md5-file", "out-dir")) {
  if (is.null(opt[[key]]) || !nzchar(opt[[key]])) stop("Missing --", key)
}

# md5 recorded alongside the fits (not recomputed: the files were verified when staged).
md5 <- read.table(opt[["md5-file"]], col.names = c("md5", "file"), colClasses = "character")
md5 <- md5$md5[basename(md5$file) == basename(opt$rds)]
if (length(md5) != 1L) stop("No unique md5 entry for ", basename(opt$rds), " in ", opt[["md5-file"]])

res <- readRDS(opt$rds)
Z <- res$Z
S <- res$summary
if (!is.data.frame(Z) || is.null(rownames(Z)) || anyDuplicated(rownames(Z))) {
  stop("Z must be a data.frame with unique gene rownames")
}
if (!identical(colnames(Z), paste0("LV", seq_len(ncol(Z))))) stop("Z columns are not LV1..LV", ncol(Z))
if (anyNA(Z) || any(!is.finite(as.matrix(Z)))) stop("Z contains missing or non-finite values")
if (!identical(nrow(res$B), ncol(Z))) stop("B has ", nrow(res$B), " rows but Z has ", ncol(Z), " LVs")
if (!is.data.frame(S) || !nrow(S) || !"LV" %in% colnames(S)) stop("summary is empty or has no LV column")
if (!all(S$LV %in% colnames(Z))) stop("summary has LVs that are not in Z")

out_dir <- opt[["out-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(Z, file.path(out_dir, "Z.csv"))
write.csv(S, file.path(out_dir, "summary.csv"), row.names = FALSE)

manifest <- list(
  schema_version = 1L,
  model = opt$model,
  source = list(path = opt$rds, md5 = md5),
  genes = nrow(Z),
  samples = ncol(res$B),
  latent_variables = ncol(Z),
  L2 = res$L2
)
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
message("Exported ", opt$model, ": ", nrow(Z), " genes x ", ncol(Z), " LVs, ", ncol(res$B), " samples")
