suppressPackageStartupMessages({library(data.table);library(Matrix);library(edgeR);library(jsonlite)})
setDTthreads(1);set.seed(0)
a<-commandArgs(TRUE);B<-a[1];mode<-a[2]
# Evaluate only original scalar constants and function definitions. No original paths or top-level jobs execute.
original<-parse(file.path(B,'evidence/original_05b.R'))
allowed<-c('MIN_CASE_SAMPLES','MIN_CONTROL_SAMPLES','MIN_DISEASE_CONTEXTS_FOR_META','MIN_DIRECTION_CONSISTENCY','MIN_COVARIATE_COMPLETENESS','GWAS_ALIGNED_CONTEXTS','bh','normalise_condition','prepare_design','run_one_de','meta_one')
for(e in original)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]]) %in% allowed)eval(e)
if(mode=='test'){
 m<-data.table(pseudobulk_id=paste0('p',1:12),condition=rep(c('healthy','disease'),each=6),study_canonical=rep(c('A','B'),6),donor_uid_v2=rep(paste0('d',1:6),2),tissue_v2='aorta',sex_v2=NA_character_,age_v2=NA_real_)
 p<-prepare_design(m);stopifnot(p$n_paired_donor_tissues==6L,qr(p$design)$rank==ncol(p$design),grepl('analysis_block',p$design_formula))
 counts<-matrix(rnbinom(600*12,mu=80,size=10),600,12);counts[1,7:12]<-counts[1,7:12]*3;features<-data.table(gene_id=paste0('g',1:600),gene_symbol=c('CKB',paste0('G',2:600)));rownames(counts)<-features$gene_id;colnames(counts)<-m$pseudobulk_id
 got<-run_one_de(counts,features,m,'L1','test','Atherosclerosis');y<-DGEList(counts=counts,genes=data.frame(gene_symbol=features$gene_symbol,row.names=features$gene_id));keep<-filterByExpr(y,design=p$design,min.count=10);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);y<-estimateDisp(y,p$design,robust=TRUE);q<-glmQLFTest(glmQLFit(y,p$design,robust=TRUE),coef=ncol(p$design));ref<-topTags(q,n=Inf,sort.by='none')$table
 stopifnot(isTRUE(all.equal(got$PValue,ref$PValue,tolerance=1e-12)),isTRUE(all.equal(got$FDR,p.adjust(ref$PValue,'BH'),tolerance=1e-12)),got[gene_symbol=='CKB',logFC]>0)
 xx<-data.table(logFC=c(-1,-1.3),F=c(8,10),PValue=c(.02,.01),FDR=c(.2,.1),disease_context=c('Aortic_dissection','Atherosclerosis'));z<-meta_one(copy(xx));r<-metafor::rma.uni(yi=xx$logFC,sei=abs(xx$logFC)/sqrt(xx$F),method='REML');stopifnot(abs(z$meta_P-r$pval)<1e-12)
 write_json(list(status='PASS',tests=c('original function extraction','paired blocks','direct native QLF numerical identity','all-gene BH','effect direction','native REML numerical identity'),edgeR=as.character(packageVersion('edgeR')),metafor=as.character(packageVersion('metafor')),glmQLFit_default_formals=capture.output(formals(edgeR:::glmQLFit.default))),file.path(B,'tests/PASS.json'),pretty=TRUE,auto_unbox=TRUE);quit('no')
}
roster<-fread(file.path(B,'tables/model_roster.tsv'));genes<-fread(file.path(B,'genes.tsv'));gwas<-fread(file.path(B,'evidence/original_GWAS_master.tsv'));priority<-unique(toupper(gwas[as.character(primary_FDR) %chin% c('TRUE','True')|as.character(primary_Bonferroni) %chin% c('TRUE','True'),SYMBOL]));stopifnot(length(priority)==2199L,'CKB' %in% priority)
if(mode=='fit'){
 id<-a[3];r<-roster[model_id==id];stopifnot(nrow(r)==1L);out<-file.path(B,'results',id);if(dir.exists(out))stop('Existing output');dir.create(out)
 inp<-file.path(B,'inputs',id);m<-fread(file.path(inp,'units.tsv'),na.strings=c('NA',''));con<-gzfile(file.path(inp,'counts.mtx.gz'),'rt');x<-as(readMM(con),'dgCMatrix');close(con);rownames(x)<-genes$gene_id;colnames(x)<-m$pseudobulk_id
 prep<-prepare_design(m);if(is.null(prep)){write_json(list(status='NOT_ESTIMABLE',reason='original design gate'),file.path(out,'MODEL.json'),auto_unbox=TRUE);quit('no')}
 fwrite(prep$meta,file.path(out,'included_units.tsv'),sep='\t');write.table(prep$design,file.path(out,'design.tsv'),sep='\t',quote=FALSE)
 result<-tryCatch(run_one_de(x,genes,m,r$level,r$cell_state,r$context),error=function(e)e)
 if(inherits(result,'error')||is.null(result)){write_json(list(status='NOT_ESTIMABLE',reason=if(is.null(result))'original expressed gene gate' else conditionMessage(result)),file.path(out,'MODEL.json'),auto_unbox=TRUE);quit('no')}
 fwrite(result,file.path(out,'DE.tsv.gz'),sep='\t');fwrite(result[gene_symbol=='CKB'],file.path(out,'CKB.tsv'),sep='\t');fwrite(result[toupper(gene_symbol) %chin% priority],file.path(out,'priority.tsv.gz'),sep='\t')
 losaudit<-list();if(uniqueN(m$study_canonical)>=2L)for(s in unique(m$study_canonical)){
  z<-tryCatch(run_one_de(x,genes,m,r$level,r$cell_state,r$context,s,priority),error=function(e)e)
  if(is.null(z)||inherits(z,'error')){losaudit[[s]]<-list(status='NOT_ESTIMABLE',reason=if(is.null(z))'original thresholds' else conditionMessage(z));next}
  fwrite(z,file.path(out,paste0('LOS_',s,'.tsv.gz')),sep='\t');losaudit[[s]]<-list(status='COMPLETE',genes=nrow(z))
 }
 paired<-prep$n_paired_donor_tissues;nr<-uniqueN(prep$meta[condition=='healthy',donor_uid_v2]);nc<-uniqueN(prep$meta[condition=='disease',donor_uid_v2])
 write_json(list(status='COMPLETE',n_genes=nrow(result),CKB_tested=any(result$gene_symbol=='CKB'),n_case=sum(prep$meta$condition=='disease'),n_control=sum(prep$meta$condition=='healthy'),case_donors=nc,reference_donors=nr,n_studies=uniqueN(prep$meta$study_canonical),paired_donor_tissues=paired,formula=prep$design_formula,block_omitted=uniqueN(prep$meta$analysis_block)>1&&!grepl('analysis_block',prep$design_formula),LOS=losaudit),file.path(out,'MODEL.json'),pretty=TRUE,auto_unbox=TRUE)
}
if(mode=='meta'){
 state<-a[3];level<-a[4];rr<-roster[roster$cell_state==state & roster$level==level & roster$status=='READY'];xs<-list()
 for(id in rr$model_id){p<-file.path(B,'results',id,'priority.tsv.gz');if(file.exists(p))xs[[id]]<-fread(p)}
 if(!length(xs))quit('no');d<-rbindlist(xs);out<-list()
 for(g in unique(d$gene_id)){v<-d[gene_id==g];ans<-meta_one(copy(v));if(is.null(ans))next;ans[,`:=`(annotation_level=level,cell_state=state,gene_id=g,gene_symbol=v$gene_symbol[1])];out[[g]]<-ans}
 if(length(out)){z<-rbindlist(out);z[,target_set_meta_FDR:=p.adjust(meta_P,'BH')];fn<-a[5];fwrite(z,fn,sep='\t')}
}
