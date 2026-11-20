from redraw_data_panels import *
from export_native_panels import extract,native_save,NS,Q

def evidence_matrix():
 p=P/'图件编辑资料/editable_svg/Figure3.svg';r=E.parse(str(p));g=r.xpath('//*[@id="FigureYa_vector_matrix_Figure3e"]')[0];els=[e for e in g if e.get('data-row')];assert len(els)==145
 d=read(S/'Figure3/Fig3d_e_all29.tsv');f=fig(100,122);a=axmm(f,37,17,50,93);a.set(xlim=(-.5,4.5),ylim=(28.5,-.5));a.set_xticks(range(5),['Primary\nGWAS','Primary\nPB','Leave-\nCAD','Regional\nABF','Spatial']);a.xaxis.tick_top();a.set_yticks(range(29),[r.gene_symbol+' | '+{'Fibroblast':'Fib','Endothelial cell':'EC','Macrophage':'Mac'}[r.cell_state] for _,r in d.iterrows()]);a.tick_params(length=0,pad=4);[sp.set_visible(False) for sp in a.spines.values()]
 for el in els:
  i,j=int(el.get('data-row'))-1,int(el.get('data-col'))-1
  a.add_patch(Rectangle((j-.43,i-.43),.86,.86,facecolor=el.get('fill'),edgecolor='none'))
 f.legend(handles=[Patch(fc=C['evidence'][key],label=lab) for key,lab in [('support','Support'),('inverse','Inverse'),('no_support','No support')]],ncol=3,loc='lower center',frameon=False,fontsize=7)
 save(f,'Figure3e',[p,S/'Figure3/Fig3d_e_all29.tsv'],'standalone gene/state names restored; exact145 frozen categorical cells; no evidence reclassification',{'data_cells':145,'source_values':[el.get('data-source-value') for el in els]})

def ED10():
 order=read(S/'Figure3/Fig3d_e_all29.tsv');p=S/'Figure5_final_companion/all29_original_meta.tsv';d=read(p)
 # Same two q05 sensitivity fits as original c; unit-level fits not recomputed.
 f=fig(119,124);a=axmm(f,36,9,74,101)
 for k,ss in enumerate(['all_sections','exclude_QC_fallback']):
  q=d[d.posterior.eq('q05')&d.term.eq('cell_abundance')&d.sensitivity.eq(ss)].set_index('gene').loc[order.gene_symbol].reset_index();assert len(q)==29
  col=['#3B6F91','#C1824A'][k]
  for i,row in q.iterrows():a.plot([row.CI_low,row.CI_high],[i+(k-.5)*.19]*2,c=col,lw=.7);a.scatter(row.beta,i+(k-.5)*.19,s=12,c=col,marker=['o','s'][k],edgecolors='none')
 a.set_yticks(range(29),[r.gene_symbol+' | '+{'Fibroblast':'Fib','Endothelial cell':'EC','Macrophage':'Mac'}[r.cell_state] for _,r in order.iterrows()]);a.set(ylim=(28.5,-.5),xlim=(-.48,.69),xlabel='Coefficient (95% CI)');a.axvline(0,c='#AAB1B7',lw=.55);clean(a);f.legend(handles=[Line2D([],[],c=c,marker=m,ms=3,label=l) for c,m,l in [('#3B6F91','o','All sections'),('#C1824A','s','Exclude QC-fallback')]],loc='lower center',ncol=2,frameon=False);save(f,'ExtendedData_Figure10c',[p],'29rows, complete stored sensitivity CIs, independent sensitivity legend',{'posterior':'q05'})
 p2=S/'ExtendedData10/ED10e.tsv';q=read(p2).set_index('pair');pairs=[r.gene_symbol+' | '+r.cell_state for _,r in order.iterrows()];q=q.loc[pairs];vals=-np.log10(q.to_numpy())
 f=fig(101,127);a=axmm(f,35,20,50,97);im=a.imshow(vals,aspect='auto',interpolation='none',cmap='Blues',vmin=0,vmax=10);a.set_yticks(range(29),[r.gene_symbol+' | '+{'Fibroblast':'Fib','Endothelial cell':'EC','Macrophage':'Mac'}[r.cell_state] for _,r in order.iterrows()]);a.set_xticks(range(4),['Mean\nAll','Mean\nExclude QC-\nfallback','q05\nAll','q05\nExclude QC-\nfallback']);a.xaxis.tick_top();a.tick_params(length=0,pad=4);[sp.set_visible(False) for sp in a.spines.values()]
 cb=f.colorbar(im,cax=axmm(f,90,27,2.5,28));cb.set_ticks([0,5,10],labels=['0','5','≥10']);cb.set_label('−log10 FDR',labelpad=4);cb.outline.set_visible(False);save(f,'ExtendedData_Figure10e',[p2],'standalone29 gene/state labels; explicit foursettings; original -log10 display transform and 0–10 clip',{'matrix_cells':116})
 # b only needs original shared row names. Match y positions to actual14 glyphs, not guessed page ticks.
 root=E.parse(str(P/'图件编辑资料/editable_svg/ExtendedData_Figure10.svg')).getroot();r=extract(root,['axes_2']);ag=r.xpath('//*[@id="axes_2"]')[0];uses=ag.xpath('.//*[local-name()="use" and starts-with(@*[local-name()="href"],"#C1_")]');assert len(uses)==14
 for i,u in enumerate(uses):
  el=E.SubElement(ag,Q('text'),x='119',y=str(float(u.get('y'))+2.1),style='font: 6.5px Arial;fill:#20252A;text-anchor:end');el.text=f'CM{i+1:02d}'
 rec=native_save(r,'ExtendedData_Figure10','b',['axes_2','restored module names']);rec['change']='exact original14 point coordinates retained; standalone CM01–CM14 row names restored';records.append(rec)

def program_heatmap():
 p=S/'ExtendedData2/ED2_program_matrix.tsv';d=read(p).set_index('Unnamed: 0');f=fig(127,104);a=axmm(f,34,9,74,74);im=a.imshow(d.to_numpy(),vmin=-2,vmax=2,cmap='RdBu_r',aspect='auto',interpolation='none');a.set_xticks(range(15),[x.replace(' cell','') for x in d.columns],rotation=65,ha='right');a.set_yticks(range(16),[x.replace(' cell','') for x in d.index]);a.tick_params(length=0,pad=3);[sp.set_visible(False) for sp in a.spines.values()];cb=f.colorbar(im,cax=axmm(f,113,18,2.5,30));cb.set_ticks([-2,0,2],labels=['≤−2','0','≥2']);cb.set_label('Mean gene-standardised expression',fontsize=7,labelpad=4);cb.outline.set_visible(False);save(f,'ExtendedData_Figure2a',[p],'240 originalprogram cells; explicit original statistic and ±2 clipped scale; independent labels',{'matrix_shape':list(d.shape)})

if __name__=='__main__':
 try:evidence_matrix();ED10();program_heatmap()
 finally:
  (O/'review/context_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2));(O/'review/context_source_hashes.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2))
