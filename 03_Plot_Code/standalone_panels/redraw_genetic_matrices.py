from redraw_data_panels import *
A=P/'Public_Release/Figure_Code/native_dependencies/Panvascular_FigureYa_20260919/code/genetics_assets'
traits=['AA','AAA','TAAD','CAD','Stroke','PAD','CA','CAVS'];labels=['AD','AAA','TAA','CAD','Stroke','PAD','CA','CAVS']
def matrix(name,d,lim,label,pvalues=None):
 f=fig(106,94);a=axmm(f,16,9,66,66);arr=d.loc[traits,traits].to_numpy();norm=Normalize(*lim);cm=plt.colormaps['RdBu_r'];rgba=cm(norm(arr))
 # Keep original upper/lower orientation and all values; lighter upper repeat avoids duplicate numbers.
 for i in range(8):
  for j in range(8):
   if j>i:rgba[i,j,3]=.25
 a.imshow(rgba,interpolation='none');a.set_xticks(range(8),labels,rotation=45,ha='right');a.set_yticks(range(8),labels);a.tick_params(length=0,pad=4);[sp.set_visible(False) for sp in a.spines.values()]
 for i in range(8):
  for j in range(i+1):
   a.text(j,i,f'{arr[i,j]:.2f}',ha='center',va='center',fontsize=7,color='white' if np.dot(cm(norm(arr[i,j]))[:3],[.2126,.7152,.0722])<.35 else '#20252A')
   if pvalues is not None and i>j and pvalues.loc[traits[i],traits[j]]<.05:a.scatter(j+.29,i-.30,s=4,c='#20252A',edgecolors='none')
 cb=f.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cm),cax=axmm(f,88,18,2.5,34));cb.set_label(label,labelpad=4);cb.outline.set_visible(False)
 if pvalues is not None:f.text(88/106,1-63/94,'● Nominal\n   P < 0.05',fontsize=7,va='top')
 return f
if __name__=='__main__':
 p=S/'Figure2/Figure2_Fig3a_LDSC_genetic_correlations.tsv';d=read(p).set_index('trait');q=A/'rg_with_P.tsv';ps=read(q).pivot(index='trait1',columns='trait2',values='rg_p');ps=ps.combine_first(ps.T);assert ps.loc[traits,traits].notna().all().all();assert sum(ps.loc[traits[i],traits[j]]<.05 for i in range(8) for j in range(i))==20
 f=matrix('Figure2b',d,(-1,1),'Genetic correlation ($r_g$)',ps);save(f,'Figure2b',[p,q],'8.25mm square cells; centered values and corner significance dots; independent colourbar; originaltrait order',{'matrix_cells':64,'p_threshold':'original nominal .05; unchanged'})
 p=A/'LDSC_intercept.tsv';d=read(p).set_index('trait');f=matrix('ExtendedData_Figure7f',d,(-.1,.1),'Cross-trait intercept');save(f,'ExtendedData_Figure7f',[p],'larger square cells and centered values; ±0.1 original clip; diagonal values preserved',{'matrix_cells':64,'diagonal':'original univariate intercepts, not cross-trait estimates'})
 (O/'review/genetic_matrix_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2));(O/'review/genetic_matrix_source_hashes.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2))
