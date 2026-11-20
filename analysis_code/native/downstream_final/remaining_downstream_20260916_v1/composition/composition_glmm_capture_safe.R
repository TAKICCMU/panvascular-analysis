# Reviewed beta-binomial models. Inputs are biological units, never individual cells.
suppressPackageStartupMessages({library(glmmTMB);library(jsonlite)})
a<-commandArgs(TRUE);stopifnot(length(a)==3L,!dir.exists(a[3]))
d<-read.delim(a[1],check.names=FALSE,stringsAsFactors=FALSE);spec<-fromJSON(a[2]);out<-a[3];dir.create(out,recursive=TRUE)
stopifnot(isTRUE(spec$approved))
if(any(d$n_state<0 | d$n_state>d$n_parent))stop('Invalid integer denominator')
if(any(is.na(d$analysis_profile_uid)) || any(duplicated(d$analysis_profile_uid)))stop('Duplicate or missing analysis observation UID')
if(any(is.na(d$donor_uid)|d$donor_uid==''))stop('Unknown donor not an independent person')
f<-as.formula(spec$formula)
# Do not permit a repeated donor to be treated as independent samples.
if(anyDuplicated(d$donor_uid)&&!grepl('1 *\\| *donor_uid',spec$formula))stop('Repeated donor requires donor random effect')
for(k in c('study','donor_uid','sex','comparison_group','vascular_bed','tissue','assay_type'))if(k %in% names(d))d[[k]]<-factor(d[[k]])
if(!is.null(spec$reference_levels))for(k in names(spec$reference_levels))d[[k]]<-relevel(factor(d[[k]]),ref=spec$reference_levels[[k]])
if('age' %in% names(d)){
 d$age<-as.numeric(d$age);if('age_z' %in% all.vars(f)){
 if(is.null(spec$age_center)||length(spec$age_center)!=1L||!is.numeric(spec$age_center)||!is.finite(spec$age_center)||is.null(spec$age_scale)||length(spec$age_scale)!=1L||!is.numeric(spec$age_scale)||!is.finite(spec$age_scale)||spec$age_scale<=0)stop('age_z requires explicitly approved finite age_center and positive age_scale; no scope re-estimation')
 d$age_z<-(d$age-spec$age_center)/spec$age_scale
 }}
required<-all.vars(f);original<-d;n0<-nrow(d)
if('sex' %in% required){
 support<-aggregate(donor_uid~study+sex,d,FUN=function(x)length(unique(x)));wide<-xtabs(donor_uid~study+sex,support)
 if(!all(c('female','male')%in%colnames(wide)))stop('No two-sex study support')
 eligible<-rownames(wide)[wide[,'female']>=5 & wide[,'male']>=5];d<-d[d$study %in% eligible,,drop=FALSE]
}
d<-d[complete.cases(d[,required,drop=FALSE]) & d$n_parent>=spec$min_denominator,,drop=FALSE]
write.table(original[!original$analysis_profile_uid %in% d$analysis_profile_uid,,drop=FALSE],file.path(out,'excluded_units.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
if(grepl('1 *\\| *study',spec$formula)&&length(unique(d$study))<2L)stop('Random study effect needs multiple studies')
if(nrow(d)<spec$min_units||length(unique(d$donor_uid))<spec$min_donors)stop('Insufficient independent units')
fixed<-lme4::nobars(f);mm<-model.matrix(fixed,d)
if(qr(mm)$rank<ncol(mm))stop('Confounded fixed-effect design; no automatic term deletion')
fitone<-function(control=NULL){tryCatch(if(is.null(control))glmmTMB(f,data=d,family=betabinomial('logit'))else glmmTMB(f,data=d,family=betabinomial('logit'),control=control),error=function(e)e)}
fit<-fitone();ok<-function(x)!inherits(x,'error')&&isTRUE(x$fit$convergence==0)&&isTRUE(x$sdr$pdHess)
optimizer<-'nlminb'
if(!ok(fit)){fit<-fitone(glmmTMBControl(optimizer=optim,optArgs=list(method='BFGS')));optimizer<-'BFGS'}
if(!ok(fit))stop('Both optimizers failed or Hessian not positive definite')
tab<-as.data.frame(summary(fit)$coefficients$cond);tab$term<-rownames(tab);names(tab)[1:4]<-c('estimate','SE','z','P');tab<-tab[tab$term %in% spec$tested_terms,,drop=FALSE]
if(!nrow(tab))stop('Declared tested terms absent')
tab$model_id<-spec$model_id;tab$FDR_family<-spec$FDR_family;tab$state<-spec$state;tab$level<-spec$level
write.table(tab,file.path(out,'effects.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(d,file.path(out,'included_units.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write_json(list(status='COMPLETE',formula=spec$formula,n_before=n0,n_after=nrow(d),n_donors=length(unique(d$donor_uid)),optimizer=optimizer,age_transform=if('age_z' %in% all.vars(f))list(center=spec$age_center,scale=spec$age_scale,source='approved model specification; fixed across scopes')else list(applied=FALSE),correction='BH after collecting all tested states in declared family'),file.path(out,'MODEL.json'),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
