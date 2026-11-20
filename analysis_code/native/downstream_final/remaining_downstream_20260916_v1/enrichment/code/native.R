# Independent explicitly approved 2026-09-16 branch; Python wrapper verifies content-bound receipt before native execution.
.libPaths(c(Sys.getenv("R_LIBS_USER"),.libPaths()))
# New evaluated gene universe and ranks; exact inherited fgsea parameters.
# No remeta from multi-df edgeR F, no online msigdbr download, no old targets.
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=3L)stop('Usage: external_pathway.R ranks.tsv frozen_genesets.tsv output_dir')
suppressPackageStartupMessages({library(data.table);library(fgsea)})
set.seed(20260822)
if(dir.exists(a[3]))stop('Non-overwrite policy: output directory exists')
dir.create(a[3],recursive=TRUE)
r<-fread(a[1]);gs<-fread(a[2])
stopifnot(all(c('analysis_id','gene','statistic','selected_up','selected_down') %in% names(r)))
stopifnot(all(c('database','pathway','gene','db_version') %in% names(gs)))
if(anyDuplicated(r[,.(analysis_id,gene)]))stop('Reviewed unique gene-symbol ranks required')
allowed<-c('GO_BP','GO_CC','GO_MF','KEGG_LEGACY','KEGG_MEDICUS','REACTOME')
if(!all(gs$database %in% allowed))stop('Unapproved pathway collection')
versions<-unique(gs$db_version)
if(length(versions)!=1L)stop('Mixed pathway resource versions')
writeLines(capture.output(sessionInfo()),file.path(a[3],'sessionInfo.txt'))
fwrite(data.table(database_version=versions),file.path(a[3],'gene_set_version.tsv'),sep='\t')
gsea<-list();ora<-list()
for(id in unique(r$analysis_id)){
  rr<-r[analysis_id==id & is.finite(statistic) & !is.na(gene) & gene!='']
  universe<-unique(rr$gene);stats<-sort(setNames(rr$statistic,rr$gene),decreasing=TRUE)
  for(db in unique(gs$database)){
    sets<-split(gs[database==db]$gene,gs[database==db]$pathway)
    z<-as.data.table(fgseaMultilevel(pathways=sets,stats=stats,minSize=10,maxSize=500,eps=1e-50,scoreType='std', BPPARAM=BiocParallel::SerialParam()))
    if(nrow(z)){z[,`:=`(analysis_id=id,database=db,n_ranked_genes=length(stats))];z[,leadingEdge:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))];gsea[[length(gsea)+1]]<-z}
    for(direction in c('up','down')){
      flag<-paste0('selected_',direction);selected<-rr[tolower(as.character(get(flag))) %in% c('true','1'),gene]
      selected<-intersect(unique(selected),universe)
      # All nonempty testable sets enter BH before displaying overlap>=3.
      q<-rbindlist(lapply(names(sets),function(path){genes<-intersect(unique(sets[[path]]),universe);k<-length(intersect(genes,selected));data.table(pathway=path,pathway_size=length(genes),overlap=k,pval=if(length(genes)>0)phyper(k-1,length(genes),length(universe)-length(genes),length(selected),lower.tail=FALSE) else NA_real_)}))
      q<-q[pathway_size>0];q[,padj:=p.adjust(pval,method='BH')]
      q[,`:=`(analysis_id=id,database=db,direction=direction,n_background_genes=length(universe),n_selected_genes=length(selected),display_overlap_ge3=overlap>=3)]
      ora[[length(ora)+1]]<-q
    }
  }
}
fwrite(rbindlist(gsea,fill=TRUE),file.path(a[3],'fgsea_all.tsv.gz'),sep='\t')
fwrite(rbindlist(ora,fill=TRUE),file.path(a[3],'ora_all_tested_before_display_filter.tsv.gz'),sep='\t')
