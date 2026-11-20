# Final-atlas adaptation of parent_lineage_preserving_null_v2.R.
# Null generation and empirical-p helpers are vendored unchanged.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: covarnet_parent_null.R prepared_input_dir")
ROOT <- normalizePath(args[1], mustWork=TRUE)
SCRIPT <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(),value=TRUE)[1]))
source(file.path(dirname(SCRIPT),'covarnet_legacy','null_statistics.R'))
readm <- function(path) as.matrix(read.delim(path,check.names=FALSE,row.names=1))
N <- as.integer(Sys.getenv('COVARNET_PARENT_NULL_PERMUTATIONS','2000'))
if(!is.finite(N)||N<1) stop('Invalid permutation count')
set.seed(as.integer(Sys.getenv("COVARNET_PARENT_NULL_SEED","20260826")))
out <- file.path(ROOT,'parent_null')
if(dir.exists(out)) stop('Refuse overwrite')
dir.create(out)
meta <- read.delim(file.path(ROOT,'CoVarNet_sample_metadata.tsv'),check.names=FALSE)
C <- readm(file.path(ROOT,'CoVarNet_node_counts.tsv'))
Wfile <- list.files(file.path(ROOT,'stable_NMF'),pattern='^module_node_weights_K[0-9]+.tsv$',full.names=TRUE)
if(length(Wfile)!=1L) stop('Unique completed primary nsNMF basis required')
W <- readm(Wfile)
C <- C[rownames(W),,drop=FALSE]
meta <- meta[match(colnames(C),meta$covarnet_sample_uid),]
stopifnot(identical(as.character(meta$covarnet_sample_uid),colnames(C)))
parentmap <- read.delim(file.path(ROOT,'CoVarNet_node_parent_map.tsv'))
parent <- parentmap$parent_lineage[match(rownames(C),parentmap$cell_state)]
stopifnot(!anyNA(parent))
P <- readm(file.path(ROOT,'CoVarNet_full_parent_counts.tsv'))[,colnames(C),drop=FALSE]
# The omitted-state remainder includes Unknown and all excluded/constant states.
# Augment each parent block, so r2dtable fixes FULL sample-parent count rather
# than only the sum of the selected nodes. Remainders are never biological nodes.
A <- C; ap <- parent
for(lin in unique(parent)) {
  remainder <- P[lin,]-colSums(C[parent==lin,,drop=FALSE])
  if(any(remainder<0)) stop('Selected counts exceed full parent')
  A <- rbind(A,remainder); ap <- c(ap,lin)
}
A <- validate_count_matrix(A)
den <- meta$n_total_atlas_cells
stopifnot(all(is.finite(den)),all(den>0))
prop <- sweep(C,2,den,'/')
within <- C; residual <- C
Q <- list()
for(lin in unique(parent)) {
  design <- model.matrix(~parent_rank+study,data.frame(parent_rank=rank(P[lin,]/den),study=factor(meta$study_canonical)))
  Q[[lin]] <- qr(design)
  if(ncol(design)>=nrow(design)||Q[[lin]]$rank>=nrow(design)-2L) stop('Parent-adjusted residual has insufficient degrees of freedom')
}
for(i in seq_len(nrow(C))) {
  within[i,] <- ifelse(P[parent[i],]>0,C[i,]/P[parent[i],],NA_real_)
  residual[i,] <- qr.resid(Q[[parent[i]]],rank(prop[i,]))
}
stat <- function(m,idx,cross=FALSE) {
  if(length(idx)<2L)return(NA_real_)
  pairs <- combn(idx,2)
  if(cross) pairs <- pairs[,parent[pairs[1,]]!=parent[pairs[2,]],drop=FALSE]
  if(!ncol(pairs))return(NA_real_)
  vals <- apply(pairs,2,function(z) suppressWarnings(cor(m[z[1],],m[z[2],],method='spearman',use='pairwise.complete.obs')))
  vals <- vals[is.finite(vals)]
  if(length(vals))median(vals) else NA_real_
}
tops <- lapply(seq_len(ncol(W)),function(k) order(W[,k],decreasing=TRUE)[seq_len(min(5L,nrow(W)))])
observed <- do.call(rbind,lapply(seq_along(tops),function(k) data.frame(module=colnames(W)[k],metric=c('raw_all','raw_cross_parent','parent_adjusted_all','parent_adjusted_cross_parent','within_parent_all'),observed=c(stat(prop,tops[[k]]),stat(prop,tops[[k]],TRUE),stat(residual,tops[[k]]),stat(residual,tops[[k]],TRUE),stat(within,tops[[k]])))))
null <- array(NA_real_,c(N,ncol(W),4))
for(b in seq_len(N)) {
  simulated <- checked_parent_null(A,ap)[seq_len(nrow(C)),,drop=FALSE]
  p <- sweep(simulated,2,den,'/'); rr <- p
  for(i in seq_len(nrow(C))) rr[i,] <- qr.resid(Q[[parent[i]]],rank(p[i,]))
  for(k in seq_along(tops))null[b,k,] <- c(stat(p,tops[[k]]),stat(p,tops[[k]],TRUE),stat(rr,tops[[k]]),stat(rr,tops[[k]],TRUE))
}
observed$empirical_P <- NA_real_; observed$n_finite_null <- 0L
for(k in seq_along(tops))for(j in seq_len(4)) {
  ix <- (k-1L)*5L+j
  observed$empirical_P[ix] <- empirical_upper(observed$observed[ix],null[,k,j])
  observed$n_finite_null[ix] <- sum(is.finite(null[,k,j]))
}
observed$FDR <- ave(observed$empirical_P,observed$metric,FUN=function(p)p.adjust(p,'BH'))
write.table(observed,file.path(out,'parent_lineage_null_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
saveRDS(null,file.path(out,'parent_null_distribution.rds'))
writeLines(c('full sample-parent counts fixed; state totals fixed; r2dtable within each parent',
 'top five nodes evaluated for ALL newly fitted modules; missing statistic remains NA',
 'within-parent proportion uses all cells of the parent, including Unknown; parent absent means NA',
 'Null inference is conditional on selected modules; it is not held-out module discovery validation',
 'within-parent statistic is descriptive; no null P value is invented'),file.path(out,'METHOD.txt'))
