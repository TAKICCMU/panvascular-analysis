suppressPackageStartupMessages(library(data.table))
a<-commandArgs(TRUE);stopifnot(length(a)==3L);root<-a[1];previous<-a[2];out<-a[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
m<-fread(file.path(previous,'aligned_latest_metadata.tsv'))
readm<-function(p){d<-fread(p);rn<-d[[1]];d[[1]]<-NULL;m<-as.matrix(d);storage.mode(m)<-'double';rownames(m)<-rn;m}
donor_means<-function(p,donors){
 stopifnot(is.matrix(p),ncol(p)==length(donors),!anyNA(donors))
 keys<-unique(donors)
 if(!length(keys))return(matrix(numeric(),nrow=0,ncol=nrow(p),dimnames=list(character(),rownames(p))))
 z<-vapply(keys,function(d)rowMeans(p[,donors==d,drop=FALSE]),numeric(nrow(p)))
 z<-matrix(z,nrow=nrow(p),ncol=length(keys),dimnames=list(rownames(p),keys))
 t(z)
}
remeta<-function(y,v){ok<-is.finite(y)&is.finite(v)&v>0;y<-y[ok];v<-v[ok];k<-length(y);if(k<2)return(list(k=k,effect=NA_real_,se=NA_real_,p=NA_real_,tau2=NA_real_,status='NOT_ESTIMABLE_FEWER_THAN_TWO_STUDIES'))
 w<-1/v;mu<-sum(w*y)/sum(w);Q<-sum(w*(y-mu)^2);tau<-max(0,(Q-k+1)/(sum(w)-sum(w^2)/sum(w)));w<-1/(v+tau);mu<-sum(w*y)/sum(w);se<-sqrt(1/sum(w));list(k=k,effect=mu,se=se,p=2*pnorm(-abs(mu/se)),tau2=tau,status='ESTIMATED_DL_APPROXIMATION')}
# Top-five-node within-study/tissue Spearman, collapsing repeat profiles by donor first.
prop<-readm(file.path(root,'CoVarNet_node_proportion_of_all_atlas_cells.tsv'));prop<-prop[,m$covarnet_sample_uid,drop=FALSE]
bas<-readm(file.path(root,'stable_NMF/module_node_weights_K14.tsv'));cr<-list();audit<-list()
for(cm in colnames(bas)){
 nodes<-rownames(bas)[order(bas[,cm],decreasing=TRUE)[1:5]]
 for(study in unique(m$study_canonical))for(tissue in unique(m[study_canonical==study,tissue_v2])){
 ix<-which(m$study_canonical==study&m$tissue_v2==tissue)
 pp<-donor_means(prop[nodes,ix,drop=FALSE],m$donor[ix]);n<-nrow(pp);
 audit[[length(audit)+1]]<-data.table(module=cm,study=study,tissue=tissue,n_profiles=length(ix),n_donors=n,status=if(n<5)'NOT_ESTIMABLE_FEWER_THAN_FIVE_DONORS' else 'ELIGIBLE');if(n<5)next
 for(pair in combn(nodes,2,simplify=FALSE)){
 rho<-suppressWarnings(cor(pp[,pair[1]],pp[,pair[2]],method='spearman'));if(!is.finite(rho))next
 rho<-max(-.999999,min(.999999,rho))
 cr[[length(cr)+1]]<-data.table(module=cm,study=study,tissue=tissue,node1=pair[1],node2=pair[2],n_donors=n,rho=rho,z=atanh(rho),variance=1/(n-3))
 }}}
fwrite(rbindlist(audit),file.path(out,'spearman_stratum_eligibility.tsv'),sep='\t')
corr<-if(length(cr))rbindlist(cr) else data.table(module=character(),study=character(),tissue=character(),node1=character(),node2=character(),n_donors=integer(),rho=double(),z=double(),variance=double());fwrite(corr,file.path(out,'within_study_tissue_donor_spearman.tsv'),sep='\t')
if(nrow(corr)){
 cm<-corr[,remeta(z,variance),by=.(module,tissue,node1,node2)];for(tis in unique(m[,.(n=uniqueN(study_canonical)),by=.(tissue_v2,donor)][n>1,tissue_v2]))cm[tissue==tis,`:=`(effect=NA_real_,se=NA_real_,p=NA_real_,status='META_BLOCKED_CROSS_STUDY_DONOR_OVERLAP')];cm[,`:=`(pooled_rho=tanh(effect),FDR=p.adjust(p,'BH'))];fwrite(cm,file.path(out,'within_tissue_spearman_meta.tsv'),sep='\t')}

writeLines(c('Disease statistics reused from content-verified prior run; no refit.','Spearman donor averaging explicitly preserves matrix dimensions for zero/one donor strata.','Fewer than five donors skipped with eligibility audit; constant nodes skipped.','Fisher-z Spearman variances are approximate; meta blocked when donors repeat across studies.'),file.path(out,'PATCH_METHODS.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines('COMPLETED',file.path(out,'ASSOCIATIONS_COMPLETED.txt'))
