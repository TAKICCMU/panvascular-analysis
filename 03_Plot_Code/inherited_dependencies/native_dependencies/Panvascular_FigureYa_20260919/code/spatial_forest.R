# Adapted from FigureYa6rHRsV2 (Dekang Lv; FigureYa, Lu et al. 2025).
# Upstream commit f627917b79f28558779fb2e3ea2014cede41b89e, CC BY-NC-SA 4.0.
# Plotting only: already fitted coefficients/CI/FDR; no model or hypothesis test.
args=commandArgs(trailingOnly=TRUE); root=args[1]
.libPaths(c(file.path(root,'R_library'),.libPaths()))
library(forestplot);library(svglite);library(grid)
data=read.delim(file.path(root,'source_data/spatial/primary_shared29.tsv'),check.names=FALSE)
lineage=c(Fibroblast='Fib',Macrophage='Mac','Endothelial cell'='EC')[data$cell_state]
tabletext=cbind(c('Gene',data$gene_symbol),c('State',lineage),c('FDR',ifelse(data$FDR_BH_all_29_pairs_x_3_terms<.001,'<0.001',sprintf('%.3f',data$FDR_BH_all_29_pairs_x_3_terms))))
cols=ifelse(data$FDR_BH_all_29_pairs_x_3_terms<.05,ifelse(data$beta>0,'#278D87','#D18A32'),'#62676B')
for(ext in c('svg')) {
 out=file.path(root,'panels',paste0('Figure5c_FigureYa.',ext))
 if(ext=='svg') svglite(out,width=98/25.4,height=96/25.4) else cairo_pdf(out,width=98/25.4,height=96/25.4,family='Arial')
 plot_obj=forestplot(labeltext=tabletext,mean=c(NA,data$beta),lower=c(NA,data$CI_low),upper=c(NA,data$CI_high),
   graph.pos=3,graphwidth=unit(48,'mm'),zero=0,xlog=FALSE,
   fn.ci_norm='fpDrawCircleCI',boxsize=.15,lwd.ci=.65,ci.vertices=TRUE,ci.vertices.height=.045,
   col=fpColors(box=cols,lines=cols,zero='#90999D'),
   shapes_gp=fpShapesGp(box=lapply(c('#222222',cols),function(x)gpar(fill=x,col=x)),lines=lapply(c('#222222',cols),function(x)gpar(col=x,lwd=.65))),
   lwd.zero=.65,lwd.xaxis=.55,xticks=c(-.4,-.2,0,.2,.4,.6),clip=c(-.48,.73),
   xlab='Coefficient (95% CI)',
   hrzl_lines=list('2'=gpar(lwd=.65,col='#B5BEC4')),
   txt_gp=fpTxtGp(label=gpar(fontfamily='Arial',cex=.5,col='#222222'),ticks=gpar(fontfamily='Arial',cex=.5,col='#222222'),xlab=gpar(fontfamily='Arial',cex=.5,col='#222222')),
   lineheight=unit(2.8,'mm'),colgap=unit(1,'mm'),align=c('l','c','r'),mar=unit(c(1,1,1,1),'mm'),new_page=FALSE)
 print(plot_obj)
 dev.off()
}
writeLines(capture.output(sessionInfo()),file.path(root,'notes/spatial_R_session.txt'))
