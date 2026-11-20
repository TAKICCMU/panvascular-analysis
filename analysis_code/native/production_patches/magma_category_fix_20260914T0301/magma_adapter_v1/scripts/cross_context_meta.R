suppressPackageStartupMessages({library(metafor);library(data.table);library(jsonlite)})
a<-commandArgs(TRUE)
if(length(a)!=3L)stop('input.tsv output.tsv APPROVED_parameters.json required; no hard-coded threshold fallback')
stopifnot(!file.exists(a[2]))
p<-fromJSON(a[3]);stopifnot(isTRUE(p$approve_legacy_ql_se_approximation),p$meta_min_contexts>=2,p$strict_shared_min_contexts>=p$meta_min_contexts,p$direction_consistency>=0.5,p$direction_consistency<=1,p$fdr_threshold>0,p$fdr_threshold<1)
if(p$context_synthesis_strategy!='one_approved_model_per_context')stop('Unapproved within-context synthesis route')
x<-fread(a[1]);required<-c('level','cell_state','gene_id','meta_context','logFC','legacy_SE_approx')
stopifnot(all(required %in% names(x)))
if(anyDuplicated(x[,.(level,cell_state,gene_id,meta_context)]))stop('Repeated effects per context: cannot treat pooled and study effects as independent')
out<-x[,{
 d<-.SD[is.finite(logFC)&is.finite(legacy_SE_approx)&legacy_SE_approx>0]
 n<-uniqueN(d$meta_context);consistency<-if(n>0)max(mean(d$logFC>0),mean(d$logFC<0)) else NA_real_
 if(n<p$meta_min_contexts)list(status='NOT_TESTED_INSUFFICIENT_CONTEXTS',effect=NA_real_,SE=NA_real_,P=NA_real_,I2=NA_real_,n_context=n,direction_consistency=consistency,n_invalid_effects=.N-n) else {
  values<-data.frame(effect=as.numeric(d$logFC),se=as.numeric(d$legacy_SE_approx))
  fit<-tryCatch(rma.uni(yi=effect,sei=se,data=values,method='REML'),error=function(e)e)
  if(inherits(fit,'error'))list(status=paste0('FAILED:',conditionMessage(fit)),effect=NA_real_,SE=NA_real_,P=NA_real_,I2=NA_real_,n_context=n,direction_consistency=consistency,n_invalid_effects=.N-n) else list(status='COMPLETE',effect=as.numeric(fit$b),SE=as.numeric(fit$se),P=as.numeric(fit$pval),I2=as.numeric(fit$I2),n_context=n,direction_consistency=consistency,n_invalid_effects=.N-n)
 }
},by=.(level,cell_state,gene_id)]
out[,target_set_meta_FDR:=p.adjust(P,'BH'),by=.(level,cell_state)]
out[,`:=`(method='REML',SE_kind='APPROXIMATE_FROM_QL_F_NOT_COEFFICIENT_SE',effect_scale='log2_fold_change',meta_min_contexts=p$meta_min_contexts,strict_shared_min_contexts=p$strict_shared_min_contexts,required_direction_consistency=p$direction_consistency)]
out[,passes_main_definition:=status=='COMPLETE' & n_context>=p$meta_min_contexts & direction_consistency>=p$direction_consistency & !is.na(target_set_meta_FDR) & target_set_meta_FDR<p$fdr_threshold]
out[,passes_strict_definition:=passes_main_definition & n_context>=p$strict_shared_min_contexts]
fwrite(out,a[2],sep='\t')
fwrite(out[passes_main_definition==TRUE],sub('\\.tsv$','_main.tsv',a[2]),sep='\t')
fwrite(out[passes_strict_definition==TRUE],sub('\\.tsv$','_strict.tsv',a[2]),sep='\t')
write_json(list(status='COMPLETE',method='REML',SE_kind='APPROXIMATE_FROM_QL_F_NOT_COEFFICIENT_SE',parameters=p,n_input=nrow(x),n_results=nrow(out)),paste0(a[2],'.receipt.json'),pretty=TRUE,auto_unbox=TRUE)
