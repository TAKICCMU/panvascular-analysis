a<-commandArgs(TRUE);source(a[1]);out<-a[2];dir.create(out,recursive=TRUE,showWarnings=FALSE);set.seed(0)
n<-24L;m<-data.frame(pseudobulk_id=paste0('u',1:n),comparison_group=rep(c('reference','case'),each=12),donor_uid=paste0('d',1:n),study=rep(rep(c('s1','s2'),each=6),2),study_block=rep(rep(c('s1','s2'),each=6),2),tissue='aorta',assay_type='scRNA-seq',sex=rep(c('female','male'),12),age=rep(seq(30,52,2),2),stringsAsFactors=FALSE)
m$donor_uid[13:18]<-m$donor_uid[1:6]
d<-construct_design(m);stopifnot(nrow(d$design)==24,qr(d$design)$rank==ncol(d$design),tail(colnames(d$design),1)=='groupcase',all(grepl('DONOR',d$meta$analysis_block[c(1:6,13:18)])))
bad<-m;bad$donor_uid<-paste0('b',1:n);bad$study_block<-bad$study<-ifelse(bad$comparison_group=='case','case_only','control_only');z<-tryCatch(construct_design(bad),error=function(e)e);stopifnot(inherits(z,'error'),grepl('NON_IDENTIFIABLE',conditionMessage(z)))
# Controls-only study included without silently dropping study adjustment.
ext<-m[1:6,];ext$pseudobulk_id<-paste0('ext',1:6);ext$donor_uid<-paste0('ed',1:6);ext$study_block<-ext$study<-'external_controls';ext$age<-c(6,6,11,46,57,47)
p<-construct_design(rbind(m,ext));stopifnot(nrow(p$meta)==30,all(ext$pseudobulk_id %in% p$meta$pseudobulk_id))
x<-matrix(rnbinom(800*n,mu=50,size=20),nrow=800);x[1,m$comparison_group=='case']<-x[1,m$comparison_group=='case']*4L;genes<-data.frame(gene_id=paste0('testgene',1:800),gene_symbol=c('CKB',paste0('g',2:800)))
inp<-file.path(out,'input');dir.create(inp);
# Matrix::writeMM does not reliably write compressed connection on older Matrix.
plain<-file.path(inp,'plain.mtx');writeMM(as(x,'sparseMatrix'),plain);zz<-gzfile(file.path(inp,'counts.mtx.gz'),'wb');writeBin(readBin(plain,'raw',n=file.info(plain)$size),zz);close(zz)
write.table(m,file.path(inp,'units.tsv'),sep='\t',row.names=FALSE,quote=FALSE);write.table(genes,file.path(out,'genes.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
run_model(inp,file.path(out,'genes.tsv'),file.path(out,'results'))
ck<-read.delim(file.path(out,'results/CKB_result.tsv'));stopifnot(nrow(ck)==1,ck$logFC>1,ck$PValue<0.05)
tt<-read.delim(file.path(out,'results/DE.tsv'));stopifnot(isTRUE(all.equal(tt$FDR,p.adjust(tt$PValue,'BH'),tolerance=1e-12)))
# Independent direct edgeR calculation must reproduce wrapper values.
y<-DGEList(counts=x,genes=genes);kp<-filterByExpr(y,design=d$design,min.count=10);y<-calcNormFactors(y[kp,,keep.lib.sizes=FALSE]);y<-estimateDisp(y,d$design,robust=TRUE);f<-glmQLFit(y,d$design,robust=TRUE,legacy=TRUE);q<-topTags(glmQLFTest(f,coef=ncol(d$design)),n=Inf,sort.by='none')$table
stopifnot(isTRUE(all.equal(tt$logFC,q$logFC,tolerance=1e-10)),isTRUE(all.equal(tt$PValue,q$PValue,tolerance=1e-10)))
write_json(list(status='PASS',tests=c('paired donor block','full nuisance rank','fully confounded condition rejected','external controls including children retained','native edgeR direct equality','CKB direction and BH scope'),edgeR=as.character(packageVersion('edgeR'))),file.path(out,'PASS.json'),auto_unbox=TRUE,pretty=TRUE)
