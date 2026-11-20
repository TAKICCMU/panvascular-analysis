suppressPackageStartupMessages(library(data.table))
a<-commandArgs(TRUE);patch<-a[1];fixture<-a[2];out<-a[3]
ev<-new.env(parent=globalenv());for(e in parse(patch))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e,ev)
p<-matrix(1:30,5,6,dimnames=list(paste0('n',1:5),NULL))
for(k in c(0,1,2,5,6)){d<-if(k)paste0('d',seq_len(k)) else character();z<-ev$donor_means(p[,seq_len(k),drop=FALSE],d);stopifnot(is.matrix(z),identical(dim(z),c(as.integer(k),5L)))}
z<-ev$donor_means(p,c('a','a','b','b','b','b'));stopifnot(nrow(z)==2,ncol(z)==5,all(z['a',]==rowMeans(p[,1:2,drop=FALSE])))
dir.create(out,recursive=TRUE,showWarnings=FALSE);previous<-file.path(out,'previous');dir.create(previous,showWarnings=FALSE)
m<-fread(file.path(fixture,'associations/aligned_latest_metadata.tsv'));m[1:2,study_canonical:='single_donor'];m[3:6,study_canonical:='two_donors'];fwrite(m,file.path(previous,'aligned_latest_metadata.tsv'),sep='\t')
code<-system2('/usr/bin/Rscript',c('--vanilla',patch,file.path(fixture,'fixture'),previous,file.path(out,'results')));stopifnot(code==0)
x<-fread(file.path(out,'results/spearman_stratum_eligibility.tsv'));stopifnot(all(x[study=='single_donor',n_donors]==1),all(x[study=='two_donors',n_donors]==2),all(x[study %in% c('single_donor','two_donors'),status]=='NOT_ESTIMABLE_FEWER_THAN_FIVE_DONORS'))
writeLines('PASS: donor matrix retains dimensions at 0,1,2,5,6 donors; repeated donor means exact; full sparse-stratum Spearman fixture succeeds.',file.path(out,'TEST_PASS.txt'))
