#!/usr/bin/env Rscript
# Prepare canonical models and project real or permuted loadings.
suppressPackageStartupMessages({library(CLAMP);library(data.table);library(jsonlite);library(bigstatsr)})
script <- dirname(normalizePath(sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(dirname(script),'common.R'))
a <- parse_cli()
stage <- required_arg(a, 'stage')
if (stage == 'prepare') {
  src <- required_arg(a,'model'); dest <- required_arg(a,'output')
  m <- read_projection_model(src)
  stopifnot(ncol(m$Z)==1728L, !anyDuplicated(rownames(m$Z)), !anyDuplicated(colnames(m$Z)))
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  saveRDS(m,dest,compress=FALSE)
  fwrite(data.table(gene=rownames(m$Z)),a$genes,sep='\t')
  write_json(list(source=normalizePath(src),md5=unname(tools::md5sum(src)),n_genes=nrow(m$Z),n_lvs=ncol(m$Z),L2=m$L2,R=as.character(getRversion()),CLAMP=as.character(packageVersion('CLAMP'))),paste0(dest,'.json'),pretty=TRUE,auto_unbox=TRUE)
} else if (stage == 'project') {
  dir.create(dirname(required_arg(a,'output')),recursive=TRUE,showWarnings=FALSE)
  m <- readRDS(required_arg(a,'model')); input <- required_arg(a,'input')
  if (endsWith(input,'.rds')) {
   x <- readRDS(input); x$is_read_only <- TRUE
   genes <- as.character(readRDS(required_arg(a,'genes'))); samples <- as.character(readRDS(required_arg(a,'samples')))
  } else {
   x <- read_matrix_csv(input); genes <- rownames(x); samples <- colnames(x)
  }
  common <- intersect(rownames(m$Z),genes)
  stopifnot(length(common)>0,!anyDuplicated(common),!anyDuplicated(samples))
  m$Z <- m$Z[common,,drop=FALSE]

  x <- x[match(common,genes),,drop=FALSE]; rownames(x)<-common; colnames(x)<-samples
  seed <- as.integer(a$seed %||% 0L)
  if(seed>0L) {

   set.seed(seed); perm <- sample(seq_along(common))
   m$Z<-m$Z[perm,,drop=FALSE];rownames(m$Z)<-common
   fwrite(data.table(gene=common,source_gene=common[perm]),paste0(a$output,'.permutation.tsv'),sep='\t')
  }
  B <- CLAMP::projectCLAMP(m,x,scale=1)
  rownames(B)<-colnames(m$Z);colnames(B)<-samples

  saveRDS(B,required_arg(a,'output'),compress=FALSE)
  write_json(list(model=normalizePath(a$model),model_md5=unname(tools::md5sum(a$model)),input=normalizePath(input),input_md5=unname(tools::md5sum(input)),seed=seed,n_genes=length(common),n_lvs=nrow(B),n_samples=ncol(B),scale=1,R=as.character(getRversion()),CLAMP=as.character(packageVersion('CLAMP'))),paste0(a$output,'.json'),pretty=TRUE,auto_unbox=TRUE)
} else {
  stop('Unknown --stage: ', stage)
}
