args <- commandArgs(trailingOnly=TRUE)

this_script<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
source(file.path(dirname(normalizePath(this_script)),'external_native_gate.R'))
stopifnot(length(args)==3L)
suppressPackageStartupMessages({library(CellChat);library(Matrix)})
d <- args[1]; min_cells <- as.integer(args[3]); set.seed(0)
x <- readMM(file.path(d,'expression.mtx'))
rownames(x)<-read.delim(file.path(d,'genes.tsv'),header=FALSE)[[1]]
colnames(x)<-read.delim(file.path(d,'cells.tsv'),header=FALSE)[[1]]
m<-read.delim(file.path(d,'obs.tsv'),row.names=1,check.names=FALSE)
stopifnot(identical(rownames(m),colnames(x)),!anyDuplicated(rownames(x)))
m$samples<-factor(m$samples)
a<-createCellChat(x,meta=m,group.by='_communication_group');a@DB<-readRDS(args[2])
a<-subsetData(a);a<-identifyOverExpressedGenes(a,do.fast=FALSE)
a<-identifyOverExpressedInteractions(a)
a<-computeCommunProb(a,type='truncatedMean',trim=0.1,population.size=FALSE)
a<-filterCommunication(a,min.cells=min_cells)
a<-computeCommunProbPathway(a);a<-aggregateNet(a)
# Retain nonsignificant edges; never condition disease contrasts on p<0.05 within samples.
z<-subsetCommunication(a,thresh=1)
write.table(z,file.path(d,'cellchat_interactions.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
saveRDS(a,file.path(d,'cellchat.rds'))
writeLines(capture.output(sessionInfo()),file.path(d,'CellChat_sessionInfo.txt'))
