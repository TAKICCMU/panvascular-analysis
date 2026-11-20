suppressPackageStartupMessages({library(edgeR);library(Matrix);library(jsonlite)})
set.seed(0)
# Keep the full nuisance column space, remove only exact redundant columns.
independent_basis <- function(z) {q<-qr(z,tol=1e-9); z[,sort(q$pivot[seq_len(q$rank)]),drop=FALSE]}
construct_design <- function(m) {
 m$group<-factor(m$comparison_group,levels=c('reference','case'))
 if(anyNA(m$group))stop('Unknown comparison group')
 if(anyDuplicated(m$pseudobulk_id))stop('Duplicate pseudobulk ID')
 build <- function(mm,terms=character()) {
  repeated<-names(which(table(mm$donor_uid)>1L))
  mm$analysis_block<-factor(ifelse(mm$donor_uid %in% repeated,paste0('DONOR::',mm$donor_uid),paste0('STUDY::',mm$study_block)))
  nuisance<-model.matrix(~0+analysis_block,mm)
  # Tissue/assay remain adjustment covariates, not separately defined disease backgrounds.
  for (k in c('tissue','assay_type',terms)) {
   vals<-mm[[k]]
   if(k=='age') {if(length(unique(vals))>1L)nuisance<-cbind(nuisance,age=as.numeric(vals))}
   else if(length(unique(vals))>1L){a<-model.matrix(~0+factor(vals));colnames(a)<-paste0(k,'::',levels(factor(vals)));nuisance<-cbind(nuisance,a)}
  }
  basis<-independent_basis(nuisance)
  g<-as.numeric(mm$group=='case');design<-cbind(basis,groupcase=g)
  if(qr(design,tol=1e-9)$rank!=ncol(design))stop('NON_IDENTIFIABLE_DISEASE: condition is confounded with mandatory study/donor/tissue/assay space; no fallback')
  if(nrow(design)-ncol(design)<2L)stop('Insufficient residual degrees of freedom')
  residual<-qr.resid(qr(basis),g)
  mm$condition_design_information<-residual^2
  list(meta=mm,design=design,nuisance_full_columns=ncol(nuisance),nuisance_basis_columns=ncol(basis))
 }
 check_arm<-function(mm) all(vapply(c('case','reference'),function(g)length(unique(mm$donor_uid[mm$comparison_group==g]))>=5L,logical(1)))
 if(!check_arm(m))stop('Fewer than 5 independent donors per arm')
 ans<-build(m);accepted<-character();covlog<-list()
 m$age<-suppressWarnings(as.numeric(m$age));m$sex[!m$sex %in% c('female','male')]<-NA
 for(k in c('sex','age')) {
  keep<-!is.na(m[[k]]);fraction<-mean(keep)
  if(fraction<0.9){covlog[[k]]<-paste0('not included: completeness ',round(fraction,3),' < 0.90');next}
  cand<-m[keep,,drop=FALSE]
  if(!check_arm(cand)){covlog[[k]]<-'not included: complete cases fail independent donor gate';next}
  if(length(unique(cand[[k]])) < ifelse(k=='age',5,2)){covlog[[k]]<-'not included: insufficient covariate variation';next}
  tryans<-tryCatch(build(cand,c(accepted,k)),error=function(e)e)
  if(inherits(tryans,'error')){covlog[[k]]<-paste0('not included: ',conditionMessage(tryans));next}
  oldans<-build(cand,accepted)
  if(ncol(tryans$design)==ncol(oldans$design)){covlog[[k]]<-'already in mandatory nuisance span: no additional coefficient required';next}
  m<-cand;accepted<-c(accepted,k);ans<-tryans;covlog[[k]]<-paste0('included; complete-case retention ',round(fraction,3))
 }
 ans$covariate_log<-covlog;ans$accepted_demographics<-accepted;ans
}
run_model<-function(input,genes_path,out) {
 if(dir.exists(out))stop('Refusing existing output');dir.create(out,recursive=TRUE)
 m<-read.delim(file.path(input,'units.tsv'),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA',''))
 m$sex[is.na(m$sex)]<-'';m$age[is.na(m$age)]<-'';m$assay_type[is.na(m$assay_type)]<-'UNREPORTED'
 genes<-read.delim(genes_path,stringsAsFactors=FALSE);con<-gzfile(file.path(input,'counts.mtx.gz'),'rt');x<-readMM(con);close(con)
 stopifnot(nrow(x)==nrow(genes),ncol(x)==nrow(m),!anyDuplicated(genes$gene_id),all(x@x>=0),all(x@x==floor(x@x)))
 write.table(m,file.path(out,'all_candidate_units.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 model<-tryCatch(construct_design(m),error=function(e)e)
 if(inherits(model,'error')){write_json(list(status='NOT_ESTIMABLE',reason=conditionMessage(model)),file.path(out,'MODEL.json'),auto_unbox=TRUE,pretty=TRUE);return(invisible(NULL))}
 mm<-model$meta;design<-model$design;colnames(x)<-m$pseudobulk_id;rownames(x)<-genes$gene_id;x<-x[,match(mm$pseudobulk_id,m$pseudobulk_id),drop=FALSE]
 rownames(design)<-mm$pseudobulk_id
 write.table(mm,file.path(out,'included_units.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(design,file.path(out,'design_matrix.tsv'),sep='\t',quote=FALSE,col.names=NA)
 deleted<-m[!m$pseudobulk_id %in% mm$pseudobulk_id,,drop=FALSE];write.table(deleted,file.path(out,'demographic_complete_case_exclusions.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 y<-DGEList(counts=x,genes=genes);keep<-filterByExpr(y,design=design,min.count=10)
 write.table(data.frame(genes,passes_filterByExpr=keep),file.path(out,'gene_filter.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 if(sum(keep)<100L){write_json(list(status='NOT_ESTIMABLE',reason='Fewer than 100 expressed genes'),file.path(out,'MODEL.json'),auto_unbox=TRUE);return(invisible(NULL))}
 y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method='TMM')
 ck<-which(genes$gene_symbol=='CKB');stopifnot(length(ck)==1L)
 cm<-cpm(as.matrix(x[ck,,drop=FALSE]),lib.size=y$samples$lib.size*y$samples$norm.factors,log=TRUE,prior.count=2)
 cktab<-cbind(mm,CKB_raw_count=as.numeric(x[ck,]),CKB_log2_CPM=as.numeric(cm),CKB_passes_filter=keep[ck],TMM_norm_factor=y$samples$norm.factors,effective_library_size=y$samples$lib.size*y$samples$norm.factors)
 write.table(cktab,file.path(out,'CKB_unit_expression.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE,legacy=TRUE);test<-glmQLFTest(fit,coef=which(colnames(design)=='groupcase'));tab<-topTags(test,n=Inf,sort.by='none')$table
 stopifnot(all(is.finite(tab$PValue)),all(tab$PValue>=0&tab$PValue<=1),!anyDuplicated(tab$gene_id))
 tab$FDR<-p.adjust(tab$PValue,method='BH')
 write.table(tab,file.path(out,'DE.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(tab[tab$gene_symbol=='CKB',,drop=FALSE],file.path(out,'CKB_result.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 info<-aggregate(condition_design_information~study,mm,sum);write.table(info,file.path(out,'study_condition_information.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
 write_json(list(status='COMPLETE',direction='positive logFC means disease higher than healthy',n_candidate_units=nrow(m),n_units=nrow(mm),n_case_donors=length(unique(mm$donor_uid[mm$comparison_group=='case'])),n_reference_donors=length(unique(mm$donor_uid[mm$comparison_group=='reference'])),n_studies=length(unique(mm$study)),n_tested_genes=nrow(tab),CKB_tested=any(tab$gene_symbol=='CKB'),covariates=model$covariate_log,accepted_demographics=model$accepted_demographics,mandatory_terms='study/repeated-donor blocks + tissue + assay; only redundant columns removed, never their column space',nuisance_full_columns=model$nuisance_full_columns,nuisance_basis_columns=model$nuisance_basis_columns,residual_df=nrow(design)-ncol(design),repeated_donor_blocks=length(unique(mm$analysis_block[grepl('^DONOR::',mm$analysis_block)])),FDR_scope='BH over all tested genes within each state/context',edgeR=as.character(packageVersion('edgeR')),seed=0),file.path(out,'MODEL.json'),auto_unbox=TRUE,pretty=TRUE)
}
if(sys.nframe()==0L){a<-commandArgs(TRUE);stopifnot(length(a)==3L);run_model(a[1],a[2],a[3])}
