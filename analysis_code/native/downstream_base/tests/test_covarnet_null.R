a <- commandArgs(trailingOnly=TRUE)
source(file.path(a[1],'src/covarnet_legacy/null_statistics.R'))
set.seed(10)
x <- matrix(c(1L,2L,3L,4L,8L,3L,4L,1L,5L,2L,7L,1L),4,3)
p <- c('F','F','M','M')
for(i in 1:20) {
 y <- checked_parent_null(x,p)
 stopifnot(identical(rowSums(x),rowSums(y)))
 for(l in unique(p))stopifnot(identical(colSums(x[p==l,]),colSums(y[p==l,])))
}
stopifnot(is.na(empirical_upper(NA,c(1,2))),is.na(empirical_upper(1,c(NA,NA))),empirical_upper(2,c(1,2,3,NA))==0.75)
for(f in c('covarnet_nmf.R','covarnet_parent_null.R'))parse(file.path(a[1],'src',f))
cat('PASS fixed margins across20draws; unavailable statistic staysNA; empiricalplusone; Rsyntax\n')
