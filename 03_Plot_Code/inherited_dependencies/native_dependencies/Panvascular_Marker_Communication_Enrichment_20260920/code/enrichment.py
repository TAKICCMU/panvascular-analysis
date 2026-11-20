from common import *
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import shutil,textwrap
f=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Content_Alignment_20260918/checks/Figure3c_seven_GSEA_ORA_verified.tsv');d=pd.read_csv(f,sep='\t');shutil.copy2(f,O/'source_data'/f.name)
names=['Lipid homeostasis','Endothelial migration (positive regulation)','Sterol transport','Sterol homeostasis','Lipoprotein assembly / remodelling / clearance','Lipoprotein clearance','LDL clearance']
assert len(d)==7
# Existing order and all seven pathways retained. No enrichment calculation.
fig=plt.figure(figsize=(183/25.4,170/25.4));y=np.arange(7)
def axes(x,top,w,h):return fig.add_axes([x/183,1-(top+h)/170,w/183,h/170])
a=axes(80,51,34,20);b=axes(130,51,30,20)
norm=colors.Normalize(0,float(max(-np.log10(d.GSEA_FDR.min()),-np.log10(d.ORA_FDR.min()))));cmap=plt.get_cmap('YlOrRd')
a.scatter(np.full(7,-.037),y,marker='s',s=5,c=['#356b8c' if v=='GO_BP' else '#a28a42' for v in d.database],clip_on=False);a.barh(y,d.NES,height=.58,color=cmap(norm(-np.log10(d.GSEA_FDR))),linewidth=0);a.set_yticks(y);a.set_yticklabels(names,fontsize=6);a.tick_params(axis='y',length=0,pad=5);a.set_ylim(6.5,-.5);a.set_xlim(0,2.05);a.set_xticks([0,1,2]);a.set_xlabel('GSEA NES',fontsize=6.5,labelpad=1);a.tick_params(axis='x',labelsize=6,length=2,pad=1)
area=38/d.overlap.max();b.scatter(d.OR,y,s=d.overlap*area,c=-np.log10(d.ORA_FDR),cmap=cmap,norm=norm,edgecolors='#555555',linewidth=.25);b.set_ylim(6.5,-.5);b.set_xlim(0,max(d.OR)*1.16);b.set_yticks([]);b.set_xlabel('ORA odds ratio',fontsize=6.5,labelpad=1);b.tick_params(axis='x',labelsize=6,length=2,pad=1);b.axvline(1,color='#bbbbbb',ls=':',lw=.4,zorder=0)
cax=axes(166,52,2,15);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax);cb.set_label('−log10(FDR)',fontsize=6,labelpad=2);cb.ax.tick_params(labelsize=6,length=2,pad=1);cb.outline.set_linewidth(.4)
# Main panel size key and category key are compact, outside the seven rows.
for j,(count,xx) in enumerate([(5,118),(10,123)]):
 k=axes(xx,73.5,4,2);k.scatter(.2,.5,s=count*area,c='#999999',edgecolors='#555555',lw=.25);k.text(.65,.5,str(count),fontsize=5.5,va='center');k.set_xlim(0,1.5);k.set_ylim(0,1);k.axis('off')
fig.text(118/183,1-73/170,'Overlap',fontsize=5.5)
r=ET.parse(B/'editable_svg/Figure3.svg').getroot();remove_ids(r,{'figure3_c_content'})
# Remove old c GO/Reactome key; categories are now source-table & standalone keys.
parents={c:p for p in r.iter() for c in p}
for g in list(r.iter(q('g'))):
 if g.get('id','').startswith('legend_') and set(t.text for t in g.iter(q('text'))) == {'GO BP','Reactome'}:parents[g].remove(g)
embed_plot(r,fig,'enrichc_');export(r,'Figure3')
# Larger readable companion, exactly same seven records.
fig,axs=plt.subplots(1,2,figsize=(183/25.4,95/25.4),gridspec_kw={'width_ratios':[1,1]});fig.subplots_adjust(left=.48,right=.88,bottom=.18,top=.91,wspace=.4)
for ax in axs:ax.set_ylim(6.5,-.5);ax.set_yticks(y);ax.tick_params(axis='y',length=0)
axs[0].set_yticklabels(names,fontsize=7);axs[1].set_yticklabels([])
axs[0].barh(y,d.NES,color=cmap(norm(-np.log10(d.GSEA_FDR))),height=.62);axs[0].set_xlabel('GSEA NES');axs[1].scatter(d.OR,y,s=d.overlap*7,c=-np.log10(d.ORA_FDR),cmap=cmap,norm=norm,edgecolors='#555555',lw=.4);axs[1].set_xlabel('ORA odds ratio');axs[1].axvline(1,lw=.5,ls=':',c='#888888')
fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.94,.35,.015,.45]),label='−log10(FDR)');axs[1].legend(handles=[Line2D([],[],marker='o',ls='',color='#777777',markersize=np.sqrt(v*7),label=str(v)) for v in [5,10]],title='Overlapping genes',fontsize=6,title_fontsize=6,frameon=False,ncol=2,bbox_to_anchor=(.5,-.16),loc='upper center');fig.savefig(O/'companion/Genetic_enrichment_seven_pathways.pdf');fig.savefig(O/'companion/Genetic_enrichment_seven_pathways.svg');plt.close(fig)
# Manuscript-primary pseudobulk: reuse stored, already selected top pathways, all states included.
f=Path('/Users/taki/Desktop/Codex任务1/pseudobulk_manuscript_strategy_enrichment_v1/completed_review_v1/primary_top_pathways_per_state_database.tsv.gz');p=pd.read_csv(f,sep='\t');shutil.copy2(f,O/'source_data'/f.name)
with PdfPages(O/'companion/Pseudobulk_primary_enrichment.pdf') as book:
 for (level,state),z in p.groupby(['level','cell_state'],sort=False):
  z=z.copy();z['leading_edge_n']=z.leadingEdge.fillna('').map(lambda s:len(s.split(';')) if s else 0);z['leading_edge_fraction']=z.leading_edge_n/z['size']
  # Keep every row of each existing top-pathway set; never rank anew by attractive direction.
  n=len(z);fig,axs=plt.subplots(1,2,figsize=(183/25.4,max(90,6*n+20)/25.4),gridspec_kw={'width_ratios':[1,1]});fig.subplots_adjust(left=.49,right=.88,top=.9,bottom=.13,wspace=.38)
  labs=[textwrap.fill(v.removeprefix('GOBP_').removeprefix('REACTOME_').removeprefix('KEGG_').replace('_',' ').lower(),50) for v in z.pathway];ys=np.arange(n);fq=-np.log10(z.padj.clip(lower=1e-300));nn=colors.Normalize(0,max(fq.max(),1))
  for ax in axs:ax.set_ylim(n-.5,-.5);ax.set_yticks(ys);ax.tick_params(axis='y',length=0);ax.axvline(0,color='#aaaaaa',lw=.4)
  axs[0].set_yticklabels(labs,fontsize=6.3);axs[1].set_yticklabels([]);axs[0].barh(ys,z.NES,color=cmap(nn(fq)),height=.62);axs[0].set_xlabel('GSEA NES')
  axs[1].scatter(z.leading_edge_fraction,ys,s=z.leading_edge_n/z.leading_edge_n.max()*72,c=fq,cmap=cmap,norm=nn,edgecolors='#666666',lw=.3);axs[1].set_xlabel('Leading-edge / gene-set size');axs[1].set_xlim(-.025,1.05)
  fig.text(.03,.96,f'{level} · {state}',fontsize=9,va='top');fig.colorbar(plt.cm.ScalarMappable(norm=nn,cmap=cmap),cax=fig.add_axes([.935,.45,.015,.37]),label='−log10(FDR)');
  handles=[Line2D([],[],marker='o',ls='',markersize=np.sqrt(v/z.leading_edge_n.max()*72),markerfacecolor='#999999',markeredgewidth=.3,color='#555555',label=str(v)) for v in sorted(set([int(z.leading_edge_n.min()),int(z.leading_edge_n.max())]))]
  axs[1].legend(handles=handles,title='Leading-edge genes',fontsize=6,title_fontsize=6,frameon=False,loc='upper left',bbox_to_anchor=(0,-.11),ncol=2)
  book.savefig(fig);plt.close(fig)
print('Genetic seven pathways preserved; pseudobulk companion pages',p.groupby(['level','cell_state']).ngroups)
