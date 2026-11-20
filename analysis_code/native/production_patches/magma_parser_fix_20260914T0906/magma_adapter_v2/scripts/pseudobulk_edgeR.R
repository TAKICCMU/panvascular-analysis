# Explicit per-contrast edgeR model. Invoked only by approved dispatcher.
suppressPackageStartupMessages({library(edgeR);library(Matrix);library(jsonlite)})
a<-commandArgs(TRUE)
if(length(a)!=5L) stop('counts.mtx.gz units.tsv genes.tsv selected_design.json NEW_OUT')
spec<-fromJSON(a[4]);out<-a[5]
if(dir.exists(out))stop('Refusing existing output');dir.create(out,recursive=TRUE)
on.exit<-NULL
con<-gzfile(a[1],'rt');counts<-readMM(con);close(con)
u<-read.delim(a[2],check.names=FALSE,stringsAsFactors=FALSE);genes<-read.delim(a[3],check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(ncol(counts)==nrow(u),nrow(counts)==nrow(genes),!anyDuplicated(u$pseudobulk_id))
m<-as.data.frame(spec$rows,stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(m$pseudobulk_id),all(m$pseudobulk_id %in% u$pseudobulk_id))
m$group<-factor(m$comparison_group,levels=c(spec$reference,spec$case))
for(k in c('donor_uid','study','sex','assay_type','vascular_bed'))if(k %in% names(m))m[[k]]<-factor(m[[k]])
if('age' %in% names(m))m$age<-as.numeric(m$age)
f<-as.formula(spec$formula);design<-model.matrix(f,m)
if(nrow(design)!=nrow(m))stop('Missing covariates silently removed rows')
if(qr(design)$rank!=ncol(design))stop('Non-identifiable design; no automatic covariate removal')
if(nrow(design)-qr(design)$rank<2L)stop('Insufficient residual degrees of freedom')
groupcols<-which(grepl('^group',colnames(design)))
if(length(groupcols)!=1L||any(grepl(':',colnames(design)[groupcols])))stop('Explicit simple case-reference group coefficient required')
x<-counts[,match(m$pseudobulk_id,u$pseudobulk_id),drop=FALSE]
if(any(x@x<0)||any(x@x!=floor(x@x)))stop('Not integer pseudobulk counts')
y<-DGEList(counts=x,genes=genes)
keep<-filterByExpr(y,design=design,min.count=10)
if(sum(keep)<spec$min_tested_genes)stop('Insufficient expressed genes')
y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method='TMM')
y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
test<-glmQLFTest(fit,coef=groupcols);tab<-topTags(test,n=Inf,sort.by='none')$table
# Historical cross-context synthesis used this approximation. Never call it exact.
tab$legacy_SE_approx<-ifelse(tab$F>0,abs(tab$logFC)/sqrt(tab$F),NA_real_)
tab$SE_kind<-'APPROXIMATE_FROM_QL_F_NOT_COEFFICIENT_SE'
tab$contrast_id<-spec$contrast_id;tab$cell_state<-spec$cell_state;tab$level<-spec$level;tab$meta_context<-spec$meta_context
write.table(tab,file.path(out,'DE.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(m,file.path(out,'included_units.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(design,file.path(out,'design_matrix.tsv'),sep='\t',quote=FALSE)
write_json(list(status='COMPLETE',formula=spec$formula,coef=colnames(design)[groupcols],n_units=nrow(m),n_donors=length(unique(m$donor_uid)),FDR_scope='all tested genes within this state and approved contrast',SE_note='Approximate QL-F-derived SE; meta-analysis gated separately'),file.path(out,'MODEL.json'),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
