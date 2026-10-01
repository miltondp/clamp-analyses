#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(hdf5r)
  library(Biostrings)
  library(parallel)
})
script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "..", "common.R"))
args <- parse_cli()
stage <- required_arg(args, "stage")
MODEL_DIR <- required_arg(args, "model_dir")
OUT <- required_arg(args, "out_dir")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MODELS <- list(canonical = list(dir = MODEL_DIR))
H5 <- args$h5 %||% "data/archs4/human_gene_v2.5.h5"
PREP <- args$preprocess_dir %||% "output/01_model_building/02_archs4/00_preprocess"
ENS_DIR <- args$ensembl_dir %||% "data/archs4/technical_noise"
ENSEMBL <- list(
  gtf = list(file = file.path(ENS_DIR, "Homo_sapiens.GRCh38.107.gtf.gz")),
  canonical_gtf = list(file = file.path(ENS_DIR, "Homo_sapiens.GRCh38.108.gtf.gz")),
  cdna = list(file = file.path(ENS_DIR, "Homo_sapiens.GRCh38.cdna.all.fa.gz")),
  ncrna = list(file = file.path(ENS_DIR, "Homo_sapiens.GRCh38.ncrna.fa.gz"))
)
CATEGORICAL <- c("platform_id", "instrument_model", "library_selection", "library_source",
                 "library_strategy", "molecule_ch1", "submission_year")
NUMERIC_QC <- c("log10_library_size", "log10_detected_genes", "mt_fraction",
                "ribo_fraction", "protein_coding_fraction", "singlecellprobability")
STUDY <- "series_id"
MIN_LEVEL_N <- 500L
SC_MAX <- 0.5
N_CORES <- as.integer(args$n_cores %||% 8L)
BLOCK_SIZE <- 1000L
MAX_SPAN <- 2500L
TH <- list(cor = 0.2, fdr = 0.2, robust_z = 4, annot_auc = 0.7, annot_fdr = 0.05)
R2_CUT <- 0.25
if (stage == "gene_properties") {
gtf_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]*)".*'), "\\1", x)
  out[!grepl(paste0(key, ' "'), x, fixed = TRUE)] <- NA_character_
  out
}
  genes <- as.character(data.table::fread(file.path(MODELS$canonical$dir, "Z.csv"), select = 1L)[[1L]])

  h5 <- H5File$new(H5, mode = "r")
  h5_genes <- data.table(symbol = h5[["/meta/genes/symbol"]]$read(),
                         ensembl_gene_id = sub("\\..*", "", h5[["/meta/genes/ensembl_gene"]]$read()))
  h5$close_all()

  gtf <- fread(cmd = paste("zcat", shQuote(ENSEMBL$gtf$file), "| grep -v '^#'"), sep = "\t", header = FALSE,
               quote = "", select = c(1L, 3L, 4L, 5L, 9L),
               col.names = c("chromosome", "feature", "start", "end", "attr"))
  gtf <- gtf[feature %in% c("gene", "transcript")]

  gene_tab <- gtf[feature == "gene", .(ensembl_gene_id = gtf_attr(attr, "gene_id"),
                                       gtf_symbol = gtf_attr(attr, "gene_name"),
                                       gene_biotype = gtf_attr(attr, "gene_biotype"),
                                       chromosome, gene_span_bp = end - start + 1L)]
  tx_tab <- gtf[feature == "transcript", .(ensembl_gene_id = gtf_attr(attr, "gene_id"),
                                           ensembl_transcript_id = gtf_attr(attr, "transcript_id"),
                                           transcript_biotype = gtf_attr(attr, "transcript_biotype"),
                                           is_basic = grepl('tag "basic"', attr, fixed = TRUE))]
  rm(gtf)

  canon <- fread(cmd = paste("zcat", shQuote(ENSEMBL$canonical_gtf$file),
                             "| awk -F'\\t' '$3 == \"transcript\"' | grep 'tag \"Ensembl_canonical\"'"),
                 sep = "\t", header = FALSE, quote = "", select = 9L, col.names = "attr")
  canonical_ids <- unique(gtf_attr(canon$attr, "transcript_id"))
  tx_tab[, is_canonical := ensembl_transcript_id %in% canonical_ids]

  read_fasta_props <- function(path) {
    seqs <- readDNAStringSet(path)
    gc <- rowSums(letterFrequency(seqs, letters = c("G", "C")))
    acgt <- rowSums(letterFrequency(seqs, letters = c("A", "C", "G", "T")))
    data.table(ensembl_transcript_id = sub("\\..*", "", sub(" .*", "", names(seqs))),
               cdna_length_bp = width(seqs), transcript_gc_percent = 100 * gc / pmax(acgt, 1))
  }
  fasta <- unique(rbind(read_fasta_props(ENSEMBL$cdna$file), read_fasta_props(ENSEMBL$ncrna$file)),
                  by = "ensembl_transcript_id")

  tx_tab <- merge(tx_tab, fasta, by = "ensembl_transcript_id")
  setorder(tx_tab, ensembl_gene_id, -is_canonical, -is_basic, -cdna_length_bp)
  tx_best <- tx_tab[, .SD[1L], by = ensembl_gene_id]

  candidates <- unique(rbind(
    h5_genes[symbol %in% genes, .(input_gene_symbol = symbol, ensembl_gene_id, map_source = "archs4_h5")],
    gene_tab[gtf_symbol %in% genes, .(input_gene_symbol = gtf_symbol, ensembl_gene_id, map_source = "ensembl_gtf")]
  ), by = c("input_gene_symbol", "ensembl_gene_id"))
  candidates[, n_candidate_genes := .N, by = input_gene_symbol]
  candidates <- merge(candidates, gene_tab[, !"gtf_symbol"], by = "ensembl_gene_id")
  candidates <- merge(candidates, tx_best, by = "ensembl_gene_id")
  candidates[, `:=`(is_protein_coding = gene_biotype == "protein_coding",
                    source_rank = as.integer(map_source != "archs4_h5"))]
  setorder(candidates, input_gene_symbol, -is_protein_coding, source_rank, -is_canonical, -cdna_length_bp)
  chosen <- candidates[, .SD[1L], by = input_gene_symbol]

  props <- chosen[match(genes, input_gene_symbol)]
  props[, input_gene_symbol := genes]
  props[, c("is_protein_coding", "source_rank") := NULL]
  fwrite(props, file.path(OUT, "gene_properties.csv"), na = "")
} else if (stage == "sample_qc") {
expected_samples <- filtered_sample_names(PREP, SC_MAX)
  meta_fields <- c(CATEGORICAL[CATEGORICAL != "submission_year"], "series_id",
                   "submission_date", "singlecellprobability", "geo_accession")
  h5 <- H5File$new(H5, mode = "r")
  meta <- as.data.table(setNames(lapply(meta_fields, function(f) h5[[paste0("/meta/samples/", f)]]$read()), meta_fields))
  gene_symbols <- h5[["/meta/genes/symbol"]]$read()
  gene_biotype <- h5[["/meta/genes/biotype"]]$read()
  h5$close_all()

  sample_idx <- match(expected_samples, meta$geo_accession)
  stopifnot(!anyNA(sample_idx), !anyDuplicated(sample_idx))
  read_idx <- sort(sample_idx)

  mt_idx <- which(grepl("^MT-", gene_symbols))
  ribo_idx <- which(grepl("^RP[LS][0-9]+[AXY]?[0-9]*$", gene_symbols))
  pc_idx <- which(gene_biotype == "protein_coding")
  cat("MT genes:", length(mt_idx), "| ribosomal:", length(ribo_idx), "| protein-coding:", length(pc_idx), "\n")

  build_blocks <- function(idx, max_n, max_rows) {
    starts <- integer(0); start <- 1L
    for (j in seq_along(idx)) {
      if (j - start + 1L > max_n || idx[j] - idx[start] + 1L > max_rows) { starts <- c(starts, start); start <- j }
    }
    starts <- c(starts, start)
    Map(function(a, b) idx[a:b], starts, c(starts[-1L] - 1L, length(idx)))
  }
  blocks <- build_blocks(read_idx, BLOCK_SIZE, MAX_SPAN)
  cat("Blocks:", length(blocks), "| max rows per read:", max(sapply(blocks, function(b) diff(range(b)) + 1L)), "\n")

  qc_for_idx <- function(idx) {
    fh <- H5File$new(H5, mode = "r")
    on.exit(fh$close_all(), add = TRUE)
    rng <- range(idx)
    x <- fh[["/data/expression"]][rng[1L]:rng[2L], , drop = FALSE]
    x <- x[idx - rng[1L] + 1L, , drop = FALSE]
    lib <- rowSums(x); denom <- pmax(lib, 1)
    data.table(h5_index = idx, library_size = lib, detected_genes = rowSums(x > 0),
               mt_fraction = rowSums(x[, mt_idx, drop = FALSE]) / denom,
               ribo_fraction = rowSums(x[, ribo_idx, drop = FALSE]) / denom,
               protein_coding_fraction = rowSums(x[, pc_idx, drop = FALSE]) / denom)
  }
  na_rows <- function(idx) data.table(h5_index = idx, library_size = NA_real_, detected_genes = NA_real_,
                                      mt_fraction = NA_real_, ribo_fraction = NA_real_,
                                      protein_coding_fraction = NA_real_)
  qc_block <- function(i) {
    idx <- blocks[[i]]
    out <- tryCatch(qc_for_idx(idx), error = function(e) NULL)
    if (!is.null(out)) return(out)
    message("block ", i, " (rows ", min(idx), "-", max(idx), ") failed; retrying in sub-blocks")
    pieces <- split(idx, ceiling(seq_along(idx) / 50))
    rbindlist(lapply(pieces, function(p) tryCatch(qc_for_idx(p), error = function(e) na_rows(p))))
  }

  qc_list <- mclapply(seq_along(blocks), qc_block, mc.cores = N_CORES, mc.preschedule = FALSE)
  bad <- which(!vapply(qc_list, function(r) inherits(r, "data.table"), logical(1)))
  for (i in bad) qc_list[[i]] <- qc_block(i)
  qc <- rbindlist(qc_list)
  stopifnot(identical(qc$h5_index, read_idx))
  qc <- qc[match(sample_idx, h5_index)]

  out <- data.table(sample = expected_samples, h5_index = sample_idx)
  out <- cbind(out, meta[sample_idx, !"geo_accession"])
  out[, submission_year := sub(".* ", "", submission_date)]
  out <- cbind(out, qc[, !"h5_index"])
  out[, `:=`(log10_library_size = log10(library_size + 1), log10_detected_genes = log10(detected_genes + 1))]

  failed <- out[is.na(library_size), .(sample, h5_index)]
  fwrite(failed, file.path(OUT, "unreadable_samples.csv"))
  fwrite(out, file.path(OUT, "sample_metadata.tsv.gz"), sep = "\t", na = "")
} else if (stage == "gene_level") {
props <- fread(file.path(OUT, "gene_properties.csv"))
safe_cor_test <- function(x, y, method, min_n = 10L) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < min_n || sd(x) == 0 || sd(y) == 0) return(c(r = NA_real_, p = NA_real_, n = n))
  tst <- suppressWarnings(cor.test(x, y, method = method, exact = if (method == "spearman") FALSE else NULL))
  c(r = unname(tst$estimate), p = tst$p.value, n = n)
}
membership_test <- function(prop, active, min_n = 10L) {
  ok <- is.finite(prop); n1 <- sum(active & ok); n0 <- sum(!active & ok)
  if (n1 < min_n || n0 < min_n) return(c(rbc = NA_real_, p = NA_real_))
  tst <- suppressWarnings(wilcox.test(prop[active & ok], prop[!active & ok], exact = FALSE))
  c(rbc = 2 * unname(tst$statistic) / (n1 * n0) - 1, p = tst$p.value)
}

gene_level_for_model <- function(model) {
  Z <- read_matrix_csv(file.path(MODELS[[model]]$dir, "Z.csv"))
  pr <- props[match(rownames(Z), input_gene_symbol)]
  len <- as.numeric(pr$cdna_length_bp); len[!is.finite(len) | len <= 0] <- NA
  gc <- as.numeric(pr$transcript_gc_percent); gc[!is.finite(gc)] <- NA
  loglen <- log10(len)

  res <- rbindlist(lapply(seq_len(ncol(Z)), function(j) {
    z <- Z[, j]; active <- is.finite(z) & z > 0
    tests <- list(
      length_pearson_all = safe_cor_test(z, loglen, "pearson"),
      GC_pearson_all = safe_cor_test(z, gc, "pearson"),
      length_spearman_all = safe_cor_test(z, len, "spearman"),
      GC_spearman_all = safe_cor_test(z, gc, "spearman"),
      length_pearson_nonzero = safe_cor_test(z[active], loglen[active], "pearson"),
      GC_pearson_nonzero = safe_cor_test(z[active], gc[active], "pearson"),
      length_spearman_nonzero = safe_cor_test(z[active], len[active], "spearman"),
      GC_spearman_nonzero = safe_cor_test(z[active], gc[active], "spearman"))
    row <- list(LV = colnames(Z)[j], n_nonzero = sum(active), fraction_nonzero = mean(active))
    for (nm in names(tests)) { row[[nm]] <- tests[[nm]][["r"]]; row[[paste0(nm, "_p")]] <- tests[[nm]][["p"]] }
    for (pn in c("length", "GC")) {
      m <- membership_test(if (pn == "length") loglen else gc, active)
      row[[paste0(pn, "_membership_rbc")]] <- m[["rbc"]]; row[[paste0(pn, "_membership_p")]] <- m[["p"]]
    }
    as.data.table(row)
  }))

  for (pc in grep("_p$", names(res), value = TRUE)) res[[sub("_p$", "_FDR", pc)]] <- p.adjust(res[[pc]], method = "BH")
  rz <- function(x) { s <- mad(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (x - median(x, na.rm = TRUE)) / s }
  stat_cols <- c(grep("_(pearson|spearman)_(all|nonzero)$", names(res), value = TRUE),
                 grep("_membership_rbc$", names(res), value = TRUE))
  for (cc in stat_cols) res[[paste0(cc, "_robust_z")]] <- rz(res[[cc]])

  pearson_stats <- c("length_pearson_all", "GC_pearson_all", "length_pearson_nonzero", "GC_pearson_nonzero")
  hits <- as.data.table(lapply(setNames(pearson_stats, pearson_stats), function(s)
    (abs(res[[s]]) >= TH$cor & res[[paste0(s, "_FDR")]] < TH$fdr & abs(res[[paste0(s, "_robust_z")]]) >= TH$robust_z) %in% TRUE))
  res[, gene_flag := Reduce(`|`, hits)]
  res[, gene_flag_tests := apply(hits, 1, function(r) paste(pearson_stats[r], collapse = ";"))]

  summary_path <- file.path(MODELS[[model]]$dir, "summary.csv")
  if (file.exists(summary_path)) {
    pw <- fread(summary_path)[AUC > TH$annot_auc & FDR < TH$annot_fdr]
    setorder(pw, LV, FDR, -AUC)
    ann <- pw[, .(n_pathways = .N, top_pathway = pathway[1L]), by = LV]
    res <- merge(res, ann, by = "LV", all.x = TRUE, sort = FALSE)
    res[is.na(n_pathways), n_pathways := 0L]
  } else {
    res[, `:=`(n_pathways = NA_integer_, top_pathway = NA_character_)]
  }
  res[, annotated := n_pathways > 0L]
  res[, model := model]
  res
}

fwrite(gene_level_for_model("canonical"), file.path(OUT, "lv_gene_properties_canonical.csv"))

} else if (stage == "sample_level") {
samples <- fread(file.path(OUT, "sample_metadata.tsv.gz"), colClasses = list(character = c(CATEGORICAL, STUDY)))
read_wide_matrix_csv <- function(path, lines_per_read = 20L) {
  con <- file(path, open = "r"); on.exit(close(con), add = TRUE)
  header <- gsub('"', "", strsplit(readLines(con, n = 1L), ",", fixed = TRUE)[[1L]], fixed = TRUE)[-1L]
  rows <- list()
  repeat {
    lines <- readLines(con, n = lines_per_read)
    if (!length(lines)) break
    for (line in lines) {
      f <- strsplit(line, ",", fixed = TRUE)[[1L]]
      rows[[length(rows) + 1L]] <- list(name = gsub('"', "", f[1L], fixed = TRUE), values = as.numeric(f[-1L]))
    }
  }
  m <- matrix(0, length(rows), length(header), dimnames = list(vapply(rows, `[[`, "", "name"), header))
  for (i in seq_along(rows)) { m[i, ] <- rows[[i]]$values; rows[i] <- list(NULL) }
  m
}

collapse_levels <- function(x) {
  x <- trimws(as.character(x)); x[is.na(x) | !nzchar(x)] <- "unknown"
  cnt <- table(x); x[x %in% names(cnt)[cnt < MIN_LEVEL_N]] <- "other"; x
}

sample_level_for_model <- function(model, chunk_size = 20000L) {
  B <- read_wide_matrix_csv(file.path(MODELS[[model]]$dir, "B.csv"))
  meta <- samples[match(colnames(B), sample)]
  stopifnot(identical(meta$sample, colnames(B)))
  n <- ncol(B)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))

  cat_vars <- lapply(setNames(CATEGORICAL, CATEGORICAL), function(v) collapse_levels(meta[[v]]))
  cat_vars <- cat_vars[vapply(cat_vars, function(x) length(unique(x)) > 1L, logical(1))]
  study <- sub("[,;].*", "", as.character(meta[[STUDY]])); study[is.na(study) | !nzchar(study)] <- "unknown"

  mu <- rowMeans(B)
  sst <- Reduce(`+`, lapply(chunks, function(j) rowSums(B[, j, drop = FALSE]^2))) - n * mu^2

  group_sums <- function(g) {
    lev <- sort(unique(g)); gi <- match(g, lev)
    s <- matrix(0, nrow(B), length(lev), dimnames = list(rownames(B), lev))
    for (j in chunks) {
      part <- rowsum(t(B[, j, drop = FALSE]), gi[j], reorder = TRUE)
      s[, as.integer(rownames(part))] <- s[, as.integer(rownames(part))] + t(part)
    }
    list(sums = s, n = tabulate(gi, length(lev)))
  }
  eta2 <- function(gs) (rowSums(sweep(gs$sums^2, 2, gs$n, "/")) - n * mu^2) / sst

  long <- list()
  for (v in names(cat_vars)) {
    gs <- group_sums(cat_vars[[v]])
    long[[v]] <- data.table(LV = rownames(B), covariate = v, measure = "eta2", value = eta2(gs))
  }
  gs_study <- group_sums(study)
  long[[STUDY]] <- data.table(LV = rownames(B), covariate = STUDY, measure = "study_eta2", value = eta2(gs_study))
  rm(gs_study)

  num_mat <- as.matrix(meta[, ..NUMERIC_QC])
  num_mat <- apply(num_mat, 2, function(x) { x[!is.finite(x)] <- median(x[is.finite(x)], na.rm = TRUE); x })
  num_rank <- apply(num_mat, 2, rank, ties.method = "average")
  spearman <- do.call(rbind, mclapply(seq_len(nrow(B)), function(i)
    cor(rank(B[i, ], ties.method = "average"), num_rank)[1L, ], mc.cores = N_CORES))
  rownames(spearman) <- rownames(B)
  for (v in NUMERIC_QC) long[[v]] <- data.table(LV = rownames(B), covariate = v, measure = "spearman", value = spearman[, v])

  design <- do.call(cbind, c(list(intercept = rep(1, n)),
    lapply(names(cat_vars), function(v) {
      x <- cat_vars[[v]]; lev <- sort(unique(x)); ref <- names(which.max(table(x)))
      d <- vapply(setdiff(lev, ref), function(l) as.numeric(x == l), numeric(n))
      colnames(d) <- paste0(v, ":", setdiff(lev, ref)); d
    }), list(scale(num_mat))))
  qrd <- qr(design); Q <- qr.Q(qrd)[, seq_len(qrd$rank), drop = FALSE]
  rm(design, qrd)
  proj <- matrix(0, nrow(B), ncol(Q))
  for (j in chunks) proj <- proj + B[, j, drop = FALSE] %*% Q[j, , drop = FALSE]
  long[["technical_joint"]] <- data.table(LV = rownames(B), covariate = "technical_joint", measure = "r2",
                                          value = (rowSums(proj^2) - n * mu^2) / sst)
  out <- rbindlist(long)
  out[, model := model]
  out
}

fwrite(sample_level_for_model("canonical"), file.path(OUT, "lv_sample_metadata_canonical.csv"))

} else if (stage == "summary") {
gene_all <- fread(file.path(OUT, "lv_gene_properties_canonical.csv"))
gene_level <- list(canonical = gene_all)
sample_all <- fread(file.path(OUT, "lv_sample_metadata_canonical.csv"))
sample_wide <- dcast(sample_all, model + LV ~ covariate, value.var = "value")
sample_flags <- sample_wide[, .(model, LV, technical_joint, series_eta2 = get(STUDY))]
sample_flags[, sample_flag := is.finite(technical_joint) & technical_joint >= R2_CUT]
cat("Sample-level cut: joint multiple-regression technical R2 >=", R2_CUT, "\n")
flagged <- merge(gene_all[, .(model, LV, gene_flag, annotated, top_pathway, n_nonzero)],
                 sample_flags, by = c("model", "LV"))
flagged[, any_flag := gene_flag | sample_flag]

cls <- flagged[model == "canonical"]
cls[, class := dplyr::case_when(!any_flag ~ "clean",
                                annotated %in% TRUE ~ "technical_correlated_biology",
                                TRUE ~ "likely_technical")]
cls[, flag_source := dplyr::case_when(gene_flag & sample_flag ~ "gene+sample",
                                      gene_flag ~ "gene", sample_flag ~ "sample", TRUE ~ "none")]
cls <- cls[order(as.integer(sub("^LV", "", LV)))]
fwrite(cls, file.path(OUT, "lv_technical_classification.csv"))

affected <- copy(cls[any_flag == TRUE,
  .(LV, gene_flag, sample_flag, technical_joint, series_eta2, n_nonzero)])
affected[, affected_group := fifelse(gene_flag & sample_flag, "Z and B",
                                     fifelse(gene_flag, "Z only", "B only"))]
stopifnot(nrow(affected) == uniqueN(affected$LV),
          all(affected$gene_flag | affected$sample_flag),
          all(affected$sample_flag == (is.finite(affected$technical_joint) &
                                      affected$technical_joint >= R2_CUT)),
          identical(sort(affected$LV[affected$gene_flag]),
                    sort(gene_level$canonical[gene_flag == TRUE, LV])))

pathway_hits <- fread(file.path(MODELS$canonical$dir, "summary.csv"))[
  LV %in% affected$LV & AUC > TH$annot_auc & FDR < TH$annot_fdr]
setorder(pathway_hits, LV, FDR, -AUC)
pathway_hits <- unique(pathway_hits, by = c("LV", "pathway"))
pathways_per_lv <- pathway_hits[, .(n_pathways = .N,
  pathways = paste(pathway, collapse = " | ")), by = LV]

trait_file <- required_arg(args, "traits")
trait_hits <- fread(cmd = paste("gzip -dc", shQuote(trait_file)))[
  lv %in% affected$LV & is.finite(fdr) & fdr < 0.05]
setorder(trait_hits, lv, fdr, pvalue)
trait_hits <- unique(trait_hits, by = c("lv", "phenotype"))
traits_per_lv <- trait_hits[, .(n_traits = .N,
  traits = paste(paste0(phenotype, ": ", phenotype_desc), collapse = " | ")), by = lv]
setnames(traits_per_lv, "lv", "LV")

affected <- merge(affected, pathways_per_lv, by = "LV", all.x = TRUE, sort = FALSE)
affected <- merge(affected, traits_per_lv, by = "LV", all.x = TRUE, sort = FALSE)
affected[is.na(n_pathways), `:=`(n_pathways = 0L, pathways = "")]
affected[is.na(n_traits), `:=`(n_traits = 0L, traits = "")]
affected[, `:=`(has_annotation = n_pathways > 0L,
                has_trait = n_traits > 0L)]
affected[, has_both_filters := has_annotation & has_trait]
setorder(affected, affected_group, -technical_joint, LV)
fwrite(affected, file.path(OUT, "affected_lv_summary.csv"))

group_levels <- c("Z only", "B only", "Z and B")
group_counts <- affected[, .(total_flagged_no_filters = .N,
  annotated_lvs = sum(has_annotation), trait_lvs = sum(has_trait),
  both_filters_lvs = sum(has_both_filters)), by = affected_group]
group_summary <- merge(data.table(affected_group = group_levels), group_counts,
                       by = "affected_group", all.x = TRUE, sort = FALSE)
for (v in setdiff(names(group_summary), "affected_group"))
  set(group_summary, which(is.na(group_summary[[v]])), v, 0L)
group_summary <- group_summary[match(group_levels, affected_group)]
stopifnot(sum(group_summary$total_flagged_no_filters) == nrow(affected),
          all(group_summary$both_filters_lvs <= group_summary$annotated_lvs),
          all(group_summary$both_filters_lvs <= group_summary$trait_lvs))
group_summary <- rbind(group_summary, data.table(
  affected_group = "Z or B",
  total_flagged_no_filters = nrow(affected),
  annotated_lvs = sum(affected$has_annotation),
  trait_lvs = sum(affected$has_trait),
  both_filters_lvs = sum(affected$has_both_filters)))
fwrite(group_summary, file.path(OUT, "affected_lv_group_summary.csv"))

} else stop("Unknown stage: ", stage)
