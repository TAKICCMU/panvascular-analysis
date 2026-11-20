from common import *
from matplotlib.patches import Rectangle
import textwrap,shutil
alignment=json.load(open('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919/audit/alignment_panels.json'))
caps=pd.read_csv(A/'audit/atlas/independent_checks.tsv',sep='\t').dropna(subset=['cohort']).set_index('cohort').cap_Q98.to_dict()
stp=pd.read_csv(A/'panels/atlas/frozen_state_id_to_hex.tsv',sep='\t');print(stp.columns)
idx=pd.read_csv(S/'ExtendedData3_4_5/ED345_all_named_state_index.tsv',sep='\t')
# Exact fixed colors; no category re-encoding.
statecol=stp.set_index('level2_annotation')['hex'].to_dict()
stateindex=idx.set_index('level2_annotation')
records=[]

def axmm(fig,x,y,w,h):return fig.add_axes([x/180,1-(y+h)/170,w/180,h/170])
def dotpanel(fig,d,states,genes,box,legend_y,vmax,color_mode='mean'):
 x,y,w,h=box;ax=axmm(fig,x,y,w,h);nr,nc=len(states),len(genes)
 means=d.pivot(index='level2_annotation',columns='gene_id',values='mean_log1p').reindex(index=states,columns=genes);det=d.pivot(index='level2_annotation',columns='gene_id',values='detection_fraction').reindex(index=states,columns=genes)
 # Dot area (not radius) directly proportional to the recorded detection fraction.
 dia=.73*min(w*mm/nc,h*mm/nr);size=dia**2
 xx,yy=np.meshgrid(np.arange(nc),np.arange(nr));norm=colors.Normalize(0,vmax)
 positive=(det.notna() & means.notna() & det.gt(0)).to_numpy()
 ax.scatter(xx[positive],yy[positive],s=(det.to_numpy()*size)[positive],c=means.to_numpy()[positive],cmap='viridis',norm=norm,linewidths=.15,edgecolors='#555555')
 na=means.isna()|det.isna()
 if na.to_numpy().any():ax.scatter(xx[na],yy[na],marker='x',s=8,color='#666666',linewidths=.4)
 ax.set_xlim(-.5,nc-.5);ax.set_ylim(nr-.5,-.5);ax.set_xticks(range(nc));ax.set_xticklabels(genes,rotation=90,ha='center',fontsize=6.2,fontstyle='italic');ax.set_yticks(range(nr));labels=[]
 for i,st in enumerate(states):
  l1=d.loc[d.level2_annotation.eq(st),'level1_annotation'].iloc[0]
  short=st
  if st!=l1:
   for suffix in [' fibroblast',' smooth muscle cell',' macrophage',' dendritic cell',' monocyte',' neutrophil',' endothelial cell',' Schwann cell',' pericyte',' fibromyocyte',' NK cell',' T cell']:
    if short.lower().endswith(suffix.lower()):short=short[:-len(suffix)];break
  
  if st!=l1:
   if st.endswith(' NK cell'):short+=' NK'
   elif st.endswith(' T cell'):short+=' T'
   ordinal=i+1
   short=f'{ordinal} {short}'
  labels.append(short)
 ax.set_yticklabels(labels,fontsize=6.2);ax.tick_params(axis='x',length=0,pad=3);ax.tick_params(axis='y',length=0,pad=9)
 for sp in ax.spines.values():sp.set_visible(False)
 # State colors remain adjacent to matrix; parent band is separately placed beside names.
 for i,st in enumerate(states):
  ax.add_patch(Rectangle((-.80,i-.43),.15,.86,facecolor=statecol.get(st,pal.get(st,'#888888')),edgecolor='none',clip_on=False))
 # Parent-group horizontal separators, no cell grid.
 groups=d.drop_duplicates('level2_annotation').set_index('level2_annotation').reindex(states).level1_annotation
 for i in range(1,nr):
  if groups.iloc[i]!=groups.iloc[i-1] and 'display_row' in d:ax.axhline(i-.5,color='#dddddd',lw=.5)
 if 'display_row' in d:
  for parent in groups.drop_duplicates():
   ids=np.flatnonzero(groups.to_numpy()==parent);ym=(ids.min()+ids.max())/2
   # Parent text and color are distinct from adjacent individual-state color strip.
   fx=(3 if x<70 else 65)/180;fy=1-(y+h*(ym+.5)/nr)/170
   fig.text(fx,fy,parent.replace(' cell',''),rotation=90,ha='center',va='center',fontsize=6,color='#222222')
 # Dot size key in its own strip.
 lax=axmm(fig,max(3,x-40),legend_y,58,6);lax.set_xlim(0,58);lax.set_ylim(0,6);lax.axis('off')
 lax.text(0,5.5,'Detected cells (%)',fontsize=6.2,va='bottom')
 for j,f in enumerate([.25,.5,.75,1.]):
  lx=3+j*13;lax.scatter(lx,2,s=size*f,c='#78888c',edgecolors='#555555',linewidths=.15);lax.text(lx+2.5,2,str(int(100*f)),va='center',fontsize=6)
 cbax=axmm(fig,140,legend_y,34,1.8);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cbax,orientation='horizontal');cb.ax.tick_params(labelsize=6,length=2,pad=1);cb.outline.set_linewidth(.35);cb.ax.set_title('Mean log1p(CP10k)',fontsize=6.2,loc='left',pad=4)
 for i,st in enumerate(states):
  for j,g in enumerate(genes):records.append(dict(state_id=st,display_label=labels[i],gene_id=g,mean_log1p=means.iloc[i,j],detection_fraction=det.iloc[i,j],area_pt2=size*det.iloc[i,j],vmin=0,vmax=vmax,palette='viridis',zero_area_preserved=True))
 return size

roots={n:ET.parse(B/'editable_svg'/f'ExtendedData_Figure{n}.svg').getroot() for n in [2,3,4,5]}
figs={n:plt.figure(figsize=(180/25.4,170/25.4)) for n in roots}
# L1: real full cross-lineage matrix, two fixed canonical genes per known lineage.
full=pd.read_csv(S/'ExtendedData2/ED2_all_selected_markers.tsv',sep='\t');full.to_csv(O/'source_data/ED2_L1_complete_matrix.tsv',sep='\t',index=False)
genes='PECAM1 VWF DCN LUM MYH11 ACTA2 FN1 TCF21 RGS5 ABCC9 C1QA CD68 FCN1 S100A9 FCGR3B CSF3R CD1C FCER1A CD3D TRAC NKG7 GNLY MS4A1 CD79A MZB1 JCHAIN TPSAB1 KIT S100B PLP1'.split();genes=[g for g in genes if g in set(full.gene_id)]
d=full.rename(columns={'level1_annotation':'level2_annotation'});d['level1_annotation']=d.level2_annotation
states=[x for x in P.L1 if x in set(d.level2_annotation)];vmax=float(full.mean_log1p.max());l1_genes=genes.copy()
dotpanel(figs[2],d,states,genes,[46,98,127,52],86,vmax)
selections=[]
for n,panel,box,ly in [(3,'b',[55,64,119,87],53),(4,'b',[99,28,75,40],17),(4,'d',[100,108,74,43],96),(5,'b',[55,74,119,81],64)]:
 f=A/f'panels/atlas/ExtendedData_Figure{n}_{panel}_full_matrix.tsv';d=pd.read_csv(f,sep='\t');shutil.copy2(f,O/'source_data'/f.name)
 states=d.sort_values('display_row').level2_annotation.drop_duplicates().tolist();colgenes=d.sort_values('display_column').gene_id.drop_duplicates().tolist()
 chosen=[]
 for st in states:
  prefix=st.split(' ')[0];provided=str(stateindex.loc[st,'selected_genes']).split(';') if st in stateindex.index else []
  first=prefix if prefix in colgenes else next((g for g in provided if g in colgenes),None)
  if first:chosen.append(first)
 chosen+= {'mural':['DCN','ACTA2','MYH11','RGS5','S100B'],'endothelial':['PECAM1','VWF','PROX1','CA4','ACKR1'],'lymphoid':['CD3D','NKG7','GNLY','MKI67','FOXP3'],'myeloid':['C1QC','CD68','LST1','FCN1','S100A8','CD1C','LILRA4','CSF3R']}[d.cohort.iloc[0]]
 genes=[g for g in colgenes if g in set(chosen)]
 # All original states; one identity-named or original first marker/state plus lineage core markers.
 dotpanel(figs[n],d,states,genes,box,ly,caps[d.cohort.iloc[0]])
 selections.append({'figure':n,'panel':panel,'states':len(states),'original_genes':len(colgenes),'display_genes':len(genes),'genes':genes})
for n,r in roots.items():
 images=[ET.tostring(e) for e in r.iter(q('image'))]
 for rec in alignment:
  if rec['figure']==f'ExtendedData_Figure{n}' and rec['panel'] in (['b','d'] if n==4 else ['b']):remove_ids(r,set(rec['members']))
 embed_plot(r,figs[n],f'newED{n}_');export(r,f'ExtendedData_Figure{n}')
 assert all(ET.tostring(e) in images for e in r.iter(q('image')) if not e.get('id','').startswith('newED')), 'Unexpected original image change'
pd.DataFrame(records).to_csv(O/'source_data/display_marker_values.tsv',sep='\t',index=False)
(O/'notes/marker_selection.json').write_text(json.dumps({'L1_genes':l1_genes,'L2_selections':selections,'L1_encoding_change':'cross-lineage full canonical marker matrix of unchanged mean log1p; replaces within-lineage-only standardized facets, not a recomputed expression statistic','complete_matrices_retained':True,'state_palette_unchanged':True},indent=2))
print(selections)
