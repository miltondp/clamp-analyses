# Combine and simplify enriched metaLV pathways.
suppressPackageStartupMessages({library(clusterProfiler); library(enrichplot); library(data.table)})

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (i in seq(1L, length(args), by = 2L)) out[[gsub("-", "_", substring(args[[i]], 3L), fixed = TRUE)]] <- args[[i + 1L]]
  out
}
opt <- parse_cli()
th <- as.numeric(opt$jaccard)
max_size <- as.integer(opt$max_size)
lvs <- fread(opt$selected)$metaLV

go <- fread(cmd = paste("cut -f1,2", shQuote(opt$bp_gmt)), header = FALSE, col.names = c("ID", "name"))
bp <- fread(file.path(opt$ora_dir, "bp/enrichment.csv.gz"))
bp[go, Description := i.name, on = "ID"]
e <- rbind(fread(file.path(opt$ora_dir, "canonical/enrichment.csv.gz"))[, db := "canonical"], bp[, db := "bp"], fill = TRUE)
e[, size := as.integer(sub("/.*", "", BgRatio))]
e[, ID := paste(db, ID, sep = ":")]

members <- list(); top <- list()
for (m in lvs) {
  x <- as.data.frame(e[LV == m & p.adjust < 0.05 & size <= max_size][order(p.adjust)])
  x <- x[!duplicated(x$Description), ]
  rownames(x) <- x$ID
  er <- new("enrichResult", result = x, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
            organism = "UNKNOWN", ontology = "UNKNOWN", gene = unique(unlist(strsplit(x$geneID, "/"))), keytype = "SYMBOL", readable = FALSE)
  er <- pairwise_termsim(er, method = "JC", showCategory = nrow(x))
  sim <- er@termsim[x$Description, x$Description]
  sim[is.na(sim)] <- t(sim)[is.na(sim)]; diag(sim) <- 1; sim[is.na(sim)] <- 0
  x$cluster <- cutree(hclust(as.dist(1 - sim), method = "average"), h = 1 - th)[x$Description]
  members[[m]] <- data.table(metaLV = m, x[, c("db", "Description", "p.adjust", "size", "Count", "geneID", "cluster")])
  rep <- x[!duplicated(x$cluster), ]
  top[[m]] <- data.table(metaLV = m, db = rep$db, pathway = rep$Description, padj = rep$p.adjust, genes = rep$geneID)
}
fwrite(rbindlist(members), opt$members_out)
fwrite(rbindlist(top), opt$top_out)
print(rbindlist(top)[, head(.SD, 3), by = metaLV][, .(metaLV, db, pathway, padj = signif(padj, 2))])
