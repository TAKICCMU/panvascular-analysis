Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1");options(mc.cores=1)
library(data.table);setDTthreads(1);suppressPackageStartupMessages(library(NMF))
w<-"/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1"
o<-file.path(w,"overnight_manuscript_integration_20260917/covarnet_final_tiers")
f<-readRDS(file.path(w,"covarnet/stable_NMF/final_nmf_K14.rds"))
h<-coef(f);h<-sweep(h,2,colSums(h),"/");h[!is.finite(h)]<-0
x<-fread(file.path(w,"covarnet/stable_NMF/sample_module_relative_coefficient_K14.tsv"))
stopifnot(identical(x[[1]],colnames(h)),max(abs(t(as.matrix(x[,-1]))-h))<1e-12)
b<-basis(f);b<-sweep(b,2,colSums(b),"/");b[!is.finite(b)]<-0
x<-fread(file.path(w,"covarnet/stable_NMF/module_node_weights_K14.tsv"));stopifnot(identical(x[[1]],rownames(b)),max(abs(as.matrix(x[,-1])-b))<1e-12)
t<-fread(file.path(o,"14_final_module_evidence_and_tiers.tsv"))
pass<-with(t,adjacent_rank_median>=.75 & leave_one_study>=.75 & leave_one_tissue>=.75 & condition_only>=.70 & hierarchical_study_bootstrap>=.65 & transformation>=.75)
stopifnot(identical(pass,t$stable_structure),sum(t$module_tier=="high_confidence_pan_vascular")==1,all(t$n_top5_pairs_meta_analyzed==10))
writeLines("PASS: exported W/H match frozen NMF RDS within 1e-12; six gates recomputed; all 14 modules and 140 top-five pairs accounted for; no NMF fit executed.",file.path(o,"VALIDATION_PASS.txt"))
