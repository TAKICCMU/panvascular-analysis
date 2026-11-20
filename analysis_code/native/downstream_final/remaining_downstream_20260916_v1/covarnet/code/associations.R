.libPaths(c('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1/Rlib',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(clubSandwich)})
args<-commandArgs(TRUE);stopifnot(length(args)==3L)
root<-args[1];sidecar<-args[2];out<-args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
readm<-function(p){d<-fread(p);rn<-d[[1]];d[[1]]<-NULL;m<-as.matrix(d);storage.mode(m)<-'double';rownames(m)<-rn;m}
m<-fread(file.path(root,'CoVarNet_sample_metadata.tsv'))
s<-fread(sidecar);stopifnot(!anyDuplicated(s$analysis_profile_uid),!anyDuplicated(m$covarnet_sample_uid))
i<-match(m$covarnet_sample_uid,s$analysis_profile_uid);stopifnot(!anyNA(i));s<-s[i]
stopifnot(!any(m$study_canonical=='GSE220774'),!any(grepl('GSM7025944|GSM7025945',m$covarnet_sample_uid)))
m[,`:=`(background=s$proposed_background_id,group=s$group_binary,donor=s$donor_block)]
stopifnot(!anyNA(m$donor),all(nzchar(m$donor)),all(m$group %in% c('healthy','disease')))
m[,association_scope_eligible:=background %in% c('G1','G2','G3') & tolower(as.character(s$main_analysis_scope_candidate))=='true']
fwrite(m,file.path(out,'scope_audit.tsv'),sep='\t')
m<-m[association_scope_eligible==TRUE]
stopifnot(nrow(m)>0)
fwrite(m,file.path(out,'aligned_latest_metadata.tsv'),sep='\t')
h<-fread(file.path(root,'stable_NMF/sample_module_relative_coefficient_K14.tsv'));h<-h[match(m$covarnet_sample_uid,h$covarnet_sample_uid)]
stopifnot(identical(h$covarnet_sample_uid,m$covarnet_sample_uid))
mods<-grep('^CM',names(h),value=TRUE);hh<-as.matrix(h[,..mods]);stopifnot(all(is.finite(hh)),all(hh>=0))
# CLR of relative NMF coefficients, with explicit fixed pseudocount sensitivity.
cluster_effect<-function(y,g,donor){
 x<-cbind(1,as.numeric(g=='disease'));G<-length(unique(donor))
 if(G<6 || length(unique(donor[g=='healthy']))<3 || length(unique(donor[g=='disease']))<3 || qr(x)$rank<2) return(list(status='NOT_ESTIMABLE_DONOR_COUNTS',effect=NA_real_,variance=NA_real_,p=NA_real_,n_donors=G,df=NA_real_))
 fit<-lm(y~I(g=='disease'));ct<-tryCatch(coef_test(fit,vcov='CR2',cluster=donor,test='Satterthwaite'),error=function(e)NULL)
 if(is.null(ct))return(list(status='NOT_ESTIMABLE_CR2',effect=NA_real_,variance=NA_real_,p=NA_real_,n_donors=G,df=NA_real_))
 list(status=if(is.finite(ct$p_Satt[2])) 'ESTIMATED_CR2_SATTERTHWAITE' else 'NOT_ESTIMABLE_CR2',effect=ct$beta[2],variance=ct$SE[2]^2,p=ct$p_Satt[2],n_donors=G,df=ct$df_Satt[2])
}
# Each background x tissue x study is estimated separately; no cross-bed disease pool.
rows<-list();strata<-unique(m[,.(background,tissue_v2,study_canonical)])
raw<-cbind(m,as.data.table(hh));rawdesc<-melt(raw,id.vars=c("background","tissue_v2","study_canonical","donor","group"),measure.vars=mods,variable.name="module",value.name="relative_H")[,.(mean_relative_H=mean(relative_H),median_relative_H=median(relative_H),n_profiles=.N,n_donors=uniqueN(donor)),by=.(background,tissue_v2,study_canonical,group,module)];fwrite(rawdesc,file.path(out,"raw_H_descriptive_by_group.tsv"),sep="\t")
for(eps in c(1e-6,1e-4)){
 z<-log(hh+eps);z<-z-rowMeans(z)
 for(a in seq_len(nrow(strata))){st<-strata[a];ix<-which(m$background==st$background&m$tissue_v2==st$tissue_v2&m$study_canonical==st$study_canonical)
 for(k in seq_along(mods)) rows[[length(rows)+1]]<-cbind(st,data.table(module=mods[k],pseudocount=eps,n_profiles=length(ix)),as.data.table(cluster_effect(z[ix,k],m$group[ix],m$donor[ix])))
 }}
r<-rbindlist(rows,fill=TRUE);r[,FDR:=p.adjust(p,'BH'),by=pseudocount];fwrite(r,file.path(out,'disease_within_study_tissue_background.tsv'),sep='\t')
remeta<-function(y,v){ok<-is.finite(y)&is.finite(v)&v>0;y<-y[ok];v<-v[ok];k<-length(y);if(k<2)return(list(k=k,effect=NA_real_,se=NA_real_,p=NA_real_,tau2=NA_real_,status='NOT_ESTIMABLE_FEWER_THAN_TWO_STUDIES'))
 w<-1/v;mu<-sum(w*y)/sum(w);Q<-sum(w*(y-mu)^2);tau<-max(0,(Q-k+1)/(sum(w)-sum(w^2)/sum(w)));w<-1/(v+tau);mu<-sum(w*y)/sum(w);se<-sqrt(1/sum(w));list(k=k,effect=mu,se=se,p=2*pnorm(-abs(mu/se)),tau2=tau,status='ESTIMATED_DL_APPROXIMATION')}
donor_overlap<-m[,.(n_studies=uniqueN(study_canonical)),by=.(background,tissue_v2,donor)][n_studies>1];fwrite(donor_overlap,file.path(out,'cross_study_donor_overlap.tsv'),sep='\t')
rr<-r[,remeta(effect,variance),by=.(background,tissue_v2,module,pseudocount)];for(j in seq_len(nrow(donor_overlap))){d<-donor_overlap[j];rr[background==d$background&tissue_v2==d$tissue_v2,`:=`(effect=NA_real_,se=NA_real_,p=NA_real_,tau2=NA_real_,status='META_BLOCKED_CROSS_STUDY_DONOR_OVERLAP')]}
rr[,FDR:=p.adjust(p,'BH'),by=pseudocount];fwrite(rr,file.path(out,'disease_same_tissue_cross_study_meta.tsv'),sep='\t')
# Top-five-node within-study/tissue Spearman, collapsing repeat profiles by donor first.
prop<-readm(file.path(root,'CoVarNet_node_proportion_of_all_atlas_cells.tsv'));prop<-prop[,m$covarnet_sample_uid,drop=FALSE]
bas<-readm(file.path(root,'stable_NMF/module_node_weights_K14.tsv'));cr<-list()
for(cm in colnames(bas)){
 nodes<-rownames(bas)[order(bas[,cm],decreasing=TRUE)[1:5]]
 for(study in unique(m$study_canonical))for(tissue in unique(m[study_canonical==study,tissue_v2])){
 ix<-which(m$study_canonical==study&m$tissue_v2==tissue)
 pp<-sapply(nodes,function(node)tapply(prop[node,ix],m$donor[ix],mean));n<-nrow(pp);if(n<5)next
 for(pair in combn(nodes,2,simplify=FALSE)){
 rho<-suppressWarnings(cor(pp[,pair[1]],pp[,pair[2]],method='spearman'));if(!is.finite(rho))next
 rho<-max(-.999999,min(.999999,rho))
 cr[[length(cr)+1]]<-data.table(module=cm,study=study,tissue=tissue,node1=pair[1],node2=pair[2],n_donors=n,rho=rho,z=atanh(rho),variance=1/(n-3))
 }}}
corr<-rbindlist(cr);fwrite(corr,file.path(out,'within_study_tissue_donor_spearman.tsv'),sep='\t')
if(nrow(corr)){
 cm<-corr[,remeta(z,variance),by=.(module,tissue,node1,node2)];for(tis in unique(m[,.(n=uniqueN(study_canonical)),by=.(tissue_v2,donor)][n>1,tissue_v2]))cm[tissue==tis,`:=`(effect=NA_real_,se=NA_real_,p=NA_real_,status='META_BLOCKED_CROSS_STUDY_DONOR_OVERLAP')];cm[,`:=`(pooled_rho=tanh(effect),FDR=p.adjust(p,'BH'))];fwrite(cm,file.path(out,'within_tissue_spearman_meta.tsv'),sep='\t')}
writeLines(c('H values are relative NMF coefficients, not absolute cell abundance.','New exploratory disease adapter (not a verbatim original-paper reproduction): CLR with pseudocount 1e-6 primary and 1e-4 sensitivity. Within background x tissue x study OLS, CR2 donor-cluster variance with Satterthwaite degrees of freedom. Minimum 3 donors/group and 6 overall.','Meta-analysis combines studies only within identical background and tissue; DL normal approximation is exploratory, especially at small study counts.','Spearman meta uses donor-mean node proportions within study and tissue; Fisher-z variance 1/(n_donors-3) is an approximation for Spearman, not exact inference.','No disease coefficient pooling across vascular beds. Original diagnosis and provenance retained.'),file.path(out,'METHODS.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines('COMPLETED',file.path(out,'ASSOCIATIONS_COMPLETED.txt'))
