suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(1)
a<-commandArgs(TRUE);B<-a[1];mode<-a[2]
expr<-parse(file.path(B,'evidence/original_pathway_enrichment.R'))
for(e in expr)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('meta_group')))eval(e)
stopifnot(exists('meta_group'))
if(mode=='test'){
  for(v in list(list(b=c(-1,-1.3),se=c(.2,.3)),list(b=c(-.7,.1,1.2),se=c(.3,.2,.4)))){
    got<-meta_group(v$b,v$se);ref<-metafor::rma.uni(yi=v$b,sei=v$se,method='DL')
    stopifnot(abs(got$beta-as.numeric(ref$b))<1e-10,abs(got$se-ref$se)<1e-10,abs(got$p-ref$pval)<1e-10)
  }
  stopifnot(is.na(meta_group(1,.2)$p))
  write_json(list(status='PASS',native_DL_matches_metafor=TRUE,single_context_not_meta=TRUE,standard_error='inherited approximate abs(logFC)/sqrt(F), not exact QL SE'),file.path(B,'tests/DL_PASS.json'),auto_unbox=TRUE,pretty=TRUE)
  quit('no')
}
stopifnot(mode=='run')
d<-fread(file.path(B,'inputs/all_manuscript_scope_DE.tsv.gz'))
stopifnot(!anyDuplicated(d[,.(annotation_level,cell_state,disease_context,gene_id)]))
valid<-is.finite(d$logFC)&is.finite(d$F)&d$F>0&!is.na(d$gene_symbol)&d$gene_symbol!=''
d[,SE_approx:=abs(logFC)/sqrt(pmax(F,0))]
valid<-valid & is.finite(d$SE_approx)&d$SE_approx>0
fwrite(d[!valid,.(annotation_level,cell_state,disease_context,gene_id,gene_symbol,logFC,F,SE_approx)],file.path(B,'tables/meta_input_exclusions.tsv.gz'),sep='\t')
d<-d[valid]
z<-d[,meta_group(logFC,SE_approx),by=.(annotation_level,cell_state,gene_id,gene_symbol)]
fwrite(z[k<2],file.path(B,'tables/genes_with_only_one_context.tsv.gz'),sep='\t')
z<-z[k>=2 & is.finite(z)]
z[,FDR_all_tested_genes:=p.adjust(p,'BH'),by=.(annotation_level,cell_state)]
z[,SE_source:='original_pathway_approximation_abs_logFC_over_sqrt_F']
fwrite(z,file.path(B,'tables/full_gene_DL_meta.tsv.gz'),sep='\t')
write_json(list(status='COMPLETE_FULL_GENE_DL_META',genes=nrow(z),states=uniqueN(z[,.(annotation_level,cell_state)]),standard_error='original approximation; not exact QL uncertainty',primary_DE_source='manuscript strategy replication',target_set_REML_modified=FALSE),file.path(B,'META_COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
