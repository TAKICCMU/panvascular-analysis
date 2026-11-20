# Non-production method comparison. NB Fisher SE is conditional on dispersion,
# not approved as a replacement for QL inference or as an exact marginal SE.
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
a<-commandArgs(TRUE);stopifnot(length(a)==1,!file.exists(a[1]))
set.seed(76);n<-20L;ng<-180L;group<-factor(rep(c('reference','case'),each=10),levels=c('reference','case'))
x<-matrix(rnbinom(n*ng,mu=30,size=10),ng,n);x[1:30,group=='case']<-x[1:30,group=='case']*4
D<-model.matrix(~group);y<-DGEList(x);y<-calcNormFactors(y);y<-estimateDisp(y,D,robust=TRUE);q<-glmQLFit(y,D,robust=TRUE);test<-glmQLFTest(q,coef=2)$table
f<-glmFit(y,D,prior.count=0);phi<-y$tagwise.dispersion
se<-vapply(seq_len(ng),function(i){w<-f$fitted.values[i,]/(1+phi[i]*f$fitted.values[i,]);sqrt(solve(crossprod(D,w*D))[2,2])/log(2)},numeric(1))
out<-data.frame(gene_id=seq_len(ng),legacy_log2FC=test$logFC,QL_F=test$F,QL_F_derived_SE=abs(test$logFC)/sqrt(test$F),NB_unshrunk_log2FC=f$coefficients[,2]/log(2),NB_conditional_Fisher_SE=se)
stopifnot(all(is.finite(out$NB_conditional_Fisher_SE)),median(out$legacy_log2FC[1:30])>1,any(abs(out$QL_F_derived_SE-out$NB_conditional_Fisher_SE)>1e-4))
write.table(out,a[1],sep='\t',row.names=FALSE,quote=FALSE)
write_json(list(status='PASS_SYNTHETIC_ONLY',n_genes=ng,n_donors=n,scale='log2',legacy='shrunken QL logFC divided by sqrt QL F; approximation',alternative='unshrunk NB coefficient and conditional Fisher SE using estimated dispersion; ignores QL variance uncertainty; not approved',median_absolute_SE_difference=median(abs(out$QL_F_derived_SE-out$NB_conditional_Fisher_SE)),production_approval=FALSE),paste0(a[1],'.json'),auto_unbox=TRUE,pretty=TRUE)
