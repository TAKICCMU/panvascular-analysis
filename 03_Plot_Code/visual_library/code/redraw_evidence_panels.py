from redraw_data_panels import *
from export_native_panels import extract,native_save

def pathways():
 base=S/'Current_Display_20260920/Figure3c_current';p=base/'MAGMA_display9.tsv';q=base/'pseudobulk_display9_all7_lineages.tsv';lp=base/'pathway_display_labels.tsv';m=read(p);b=read(q);lab=read(lp).set_index('pathway_id').display_label.to_dict();paths=b.pathway.drop_duplicates().tolist();states=b.cell_state.drop_duplicates().tolist();assert len(paths)==9 and len(states)==7
 f=fig(155,143);a=axmm(f,58,12,64,52.2);a.barh(range(9),m.NES,color=C['method']['MAGMA'],height=.66);a.set(ylim=(8.5,-.5),xlim=(0,2.15),xlabel='NES');a.set_yticks(range(9),['\n'.join(textwrap.wrap(lab[p],36)) for p in m.pathway]);clean(a);a.grid(axis='x',c='#E5E8EA',lw=.45);a.set_axisbelow(True)
 for i,row in m.iterrows():f.text(130/155,1-(12+(i+.5)*5.8)/143,f'{row.padj:.3f}',va='center',fontsize=7)
 f.text(130/155,1-9/143,'Adjusted P',fontsize=7);f.text(58/155,1-7/143,'MAGMA GSEA',fontsize=7,ha='left')
 v=b.pivot(index='pathway',columns='cell_state',values='NES').reindex(index=paths,columns=states);sig=b.pivot(index='pathway',columns='cell_state',values='padj').reindex(index=paths,columns=states);vm=float(b.NES.abs().max())
 a=axmm(f,58,88,49,45);a.imshow(v,aspect='auto',cmap='RdBu_r',vmin=-vm,vmax=vm,interpolation='none');yy,xx=np.where(sig.to_numpy()<.05);a.scatter(xx,yy,c='#20252A',s=5,edgecolors='none');a.set_yticks(range(9),['\n'.join(textwrap.wrap(lab[p],36)) for p in paths]);a.set_xticks(range(7),[{'Endothelial cell':'EC','Fibroblast':'Fib','Smooth muscle cell':'SMC','Macrophage':'Mac','Monocyte':'Mo','T cell':'T','Natural killer cell':'NK'}[s] for s in states]);a.xaxis.tick_top();a.tick_params(length=0,pad=4);[sp.set_visible(False) for sp in a.spines.values()]
 f.text(58/155,1-78/143,'Pseudobulk GSEA',fontsize=7)
 for yy,pitch in [(12,5.8),(88,5.)]:
  for j,db in enumerate(['GO BP','KEGG','Reactome']):
   x=4/155;y=1-(yy+(j*3+1.5)*pitch)/143;f.text(x,y,db,rotation=90,ha='center',va='center',fontsize=6.5)
   f.add_artist(Line2D([7/155]*2,[1-(yy+(j*3)*pitch+.5)/143,1-(yy+(j*3+3)*pitch-.5)/143],transform=f.transFigure,c='#747C86',lw=.7))
 cb=f.colorbar(plt.cm.ScalarMappable(norm=Normalize(-vm,vm),cmap='RdBu_r'),cax=axmm(f,120,93,2.3,21));cb.outline.set_visible(False);cb.set_label('NES',labelpad=4);f.text(115/155,1-129/143,'● Adjusted P < 0.05',fontsize=7)
 save(f,'Figure3c',[p,q,lp],'independent genetic/expression blocks stacked for readable full pathway names; fixed9+9 pathways; all63 matrix values retained',{'MAGMA_rows':9,'PB_cells':63,'PB_color_limits':[-vm,vm],'no_new_tests':True})

def coloc():
 p=S/'Current_Display_20260920/source_data/Figure3f_same5_contexts_both_scenarios.tsv';d=read(p);pairs=d[['gene','context']].drop_duplicates();assert len(pairs)==5;f=fig(175,53)
 for j,(field,title,lim,ticks) in enumerate([('PP.H4','PP.H4',(0,1),[0,.5,1]),('PP.H3','PP.H3',(0,.4),[0,.2,.4]),('n_overlap','Shared SNPs',(0,4500),[0,2000,4000])]):
  a=axmm(f,62+j*36,6,29,28)
  for i,r in pairs.reset_index(drop=True).iterrows():
   for sc,off,mark,col in [('primary',-.12,'o',C['method']['primary']),('leave_one_CAD',.12,'s',C['method']['leave_CAD'])]:
    row=d[d.gene.eq(r.gene)&d.context.eq(r.context)&d.scenario.eq(sc)].iloc[0];strong=row.classification=='strong_colocalisation';a.scatter(row[field],i+off,s=19,marker=mark,facecolors=col if strong else 'white',edgecolors=col,lw=.75,zorder=3,clip_on=False)
  a.set(ylim=(4.5,-.5),xlim=lim);a.set_xticks(ticks,['0','2k','4k'] if j==2 else [str(x) for x in ticks]);a.set_yticks(range(5),[r.gene+' | '+r.context.replace('_',' ') for _,r in pairs.iterrows()] if j==0 else []);a.grid(axis='x',c='#E5E8EA',lw=.45);clean(a);a.set_xlabel(title)
 f.legend(handles=[Line2D([],[],marker='o',ls='',mfc=C['method']['primary'],mec=C['method']['primary'],label='Primary'),Line2D([],[],marker='s',ls='',mfc=C['method']['leave_CAD'],mec=C['method']['leave_CAD'],label='Leave-CAD'),Line2D([],[],marker='o',ls='',mfc='#747C86',mec='#747C86',label='Strong'),Line2D([],[],marker='o',ls='',mfc='white',mec='#747C86',label='Other')],ncol=4,loc='lower center',frameon=False)
 save(f,'Figure3f',[p],'5 contexts ×2 scenarios preserved; three evidence columns; separate scenario and strength symbols',{'observations':10,'strength':'stored classification, no PP.H4-only recoding'})

def tissue():
 p=S/'Figure1/L1_tissue_sample_mean_fraction.tsv';d=read(p);states=read(S/'Figure1/L1_palette_bound.tsv').sort_values('order').L1.tolist();f=fig(142,79);a=axmm(f,37,5,98,47);left=np.zeros(len(d))
 for st in states:a.barh(range(len(d)),d[st],left=left,height=.73,color=L1[st],edgecolor='none');left+=d[st].to_numpy()
 a.set_yticks(range(len(d)),d.tissue);a.set(ylim=(len(d)-.5,-.5),xlim=(0,1),xlabel='Mean sample fraction');a.set_xticks([0,.25,.5,.75,1]);clean(a)
 f.legend(handles=[Patch(fc=L1[s],label=s.replace(' cell','')) for s in states],loc='lower center',ncol=4,frameon=False,fontsize=7,columnspacing=1,handlelength=1)
 save(f,'Figure1d',[p],'fixed12tissues and original fractions; standalone16lineage key close to bars',{'cells':len(d)*len(states),'denominator':'original sample mean fractions, unchanged'})

def evidence_native_labels():
 # Read existing data-cell attributes; do not infer significance again.
 orig=P/'图件编辑资料/editable_svg/Figure3.svg';r=E.parse(str(orig)).getroot();g=r.xpath('//*[@id="FigureYa_vector_matrix_Figure3e"]')[0];cells=[e for e in g if e.get('data-row') is not None];print('evidencecell',len(cells),[dict(x.attrib) for x in cells[:1]])

if __name__=='__main__':
 try:pathways();coloc();tissue();evidence_native_labels()
 finally:
  (O/'review/evidence_redraw_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2));(O/'review/evidence_source_hashes.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2))
