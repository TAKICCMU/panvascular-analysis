from f1_common import *
from matplotlib.patches import Rectangle,FancyArrowPatch
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import base64
W,H=183,170
fig=plt.figure(figsize=(W/25.4,H/25.4));fig.patch.set_facecolor('white')
def axmm(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def ftext(x,y,s,size=6.3,**kw):return fig.text(x/W,1-y/H,s,fontsize=size,va='center',**kw)
def letter(x,y,s):ftext(x,y,s,8,fontweight='bold')
a=axmm(0,0,183,20);a.set(xlim=(0,183),ylim=(20,0));a.axis('off');letter(1,2,'a')
for x,w in [(6,49),(64,64),(137,43)]:a.add_patch(Rectangle((x,2),w,16,fc='#f2f5f7',ec='#c3cfd6',lw=.5))
a.text(8,5,'Public human vascular datasets',fontsize=6.4,weight='bold',va='center');a.text(8,9,'26 studies · 247 sample IDs',fontsize=6.3,va='center');a.text(8,13,'37,336 genes',fontsize=6.3,va='center')
for x,y in [(45,9),(47,10),(49,11)]:a.add_patch(Rectangle((x,y),3,4,fc='white',ec='#356B8C',lw=.4))
a.text(67,5,'Integration and annotation',fontsize=6.4,weight='bold',va='center');a.text(67,9,'Harmony · 1,175,115 cells',fontsize=6.3,va='center');a.text(67,13,'16 lineages · 116 named states',fontsize=6.3,va='center')
a.text(140,5,'Downstream analyses',fontsize=6.4,weight='bold',va='center');a.text(140,9,'Genetics · Pseudobulk',fontsize=6.3,va='center');a.text(140,13,'CoVarNet · Spatial',fontsize=6.3,va='center')
for l,r in [(55,64),(128,137)]:a.add_patch(FancyArrowPatch((l+1,10),(r-1,10),arrowstyle='-|>',mutation_scale=7,color='#748b98',lw=.6))
D=pd.read_csv(S/'Figure1/full_cells_frozen_umap.tsv.gz',sep='\t');assert len(D)==1175115
Mpath=Path('/Users/taki/Desktop/Codex任务1/CoVarNet_pseudobulk_full_review_20260915/updated_scope_audit_v2/sample_sidecar_updated.tsv');M=pd.read_csv(Mpath,sep='\t');assert M.sample_uid.is_unique
D=D.merge(M[['sample_uid','tissue','vascular_bed']],on='sample_uid',how='left',validate='many_to_one');assert D.tissue.notna().all()
broad={'ascending aorta':'Aorta','descending aorta':'Aorta','abdominal aorta':'Aorta','aortic arch':'Aorta','aorta':'Aorta','coronary artery':'Coronary artery','carotid artery':'Carotid artery','pulmonary artery':'Pulmonary artery','aortic valve':'Aortic valve','internal thoracic artery':'Other systemic arteries','iliac artery':'Other systemic arteries','anterior tibial artery':'Other systemic arteries','brachial artery':'Other systemic arteries'}
D['Tissue']=D.vascular_bed.map(broad);assert D.Tissue.notna().all();tpal=dict(zip(['Aorta','Coronary artery','Carotid artery','Pulmonary artery','Aortic valve','Other systemic arteries'],['#4E79A7','#E69F45','#6A9C78','#9A76B5','#CB7689','#789B9E']))
D[['sample_uid','tissue','vascular_bed','Tissue']].drop_duplicates().assign(source_file=str(Mpath)).to_csv(O/'source_data/Figure1b_tissue_sample_mapping.tsv',sep='\t',index=False);pd.DataFrame(list(tpal.items()),columns=['Tissue','hex']).to_csv(O/'source_data/Figure1b_tissue_palette.tsv',sep='\t',index=False)
x0,x1=D.UMAP1.min(),D.UMAP1.max();y0,y1=D.UMAP2.min(),D.UMAP2.max();dx=x1-x0;dy=y1-y0;xlim=(x0-.04*dx,x1+.04*dx);ylim=(y0-.04*dy,y1+.04*dy)
letter(1,21,'b')
for x,y,w,h,col in [(34,27,49,33,'level1_annotation'),(125,23,51,43,'Tissue')]:
 ax=axmm(x,y,w,h);ax.set(xlim=xlim,ylim=ylim);ax.set_aspect('equal');ax.axis('off');p=pal if col=='level1_annotation' else tpal
 ax.scatter(D.UMAP1,D.UMAP2,c=D[col].map(p),s=.055,lw=0,alpha=.8,rasterized=True)
 if col=='Tissue':ftext(x+w/2,y-1.5,'Tissue',6.5,ha='center')
 else:mainax=ax
fig.canvas.draw()
# Existing actual cell anchors; only the labels are assigned to four peripheral sides.
anchors=pd.read_csv('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/reproducibility_dependencies/Panvascular_Figure1b_Complete_Labels_20260920/notes/Figure1b_16_label_mapping.tsv',sep='\t')
xy=mainax.transData.transform(anchors[['UMAP1','UMAP2']].to_numpy());xy=fig.transFigure.inverted().transform(xy);xy=np.c_[xy[:,0]*W,(1-xy[:,1])*H]
slots=[]
for side,x in [('left',31.5),('right',84.5)]:
 for y in np.linspace(29,59,6):slots.append((x,y,side))
slots +=[(44,24,'top'),(72,24,'top'),(44,64,'bottom'),(72,64,'bottom')]
pos=np.array([s[:2] for s in slots]);ri,assign=linear_sum_assignment(cdist(xy,pos));
def intersect(a,b,c,d):
 def cross(x,y,z):return np.cross(y-x,z-x)
 return cross(a,b,c)*cross(a,b,d)<-1e-8 and cross(c,d,a)*cross(c,d,b)<-1e-8
def score(perm):
 n=sum(intersect(xy[i],pos[perm[i]],xy[j],pos[perm[j]]) for i in range(16) for j in range(i))
 return n*1e5+np.sum(np.linalg.norm(xy-pos[perm],axis=1)),n
for step in range(100):
 old=score(assign)[0];best=old;bestp=None
 for i in range(16):
  for j in range(i):
   p=assign.copy();p[i],p[j]=p[j],p[i];s=score(p)[0]
   if s<best:best=s;bestp=p
 if bestp is None:break
 assign=bestp
labelax=axmm(0,0,W,H);labelax.set(xlim=(0,W),ylim=(H,0));labelax.axis('off')
for i,row in anchors.iterrows():
 x,y,side=slots[assign[i]];ha={'left':'right','right':'left','top':'center','bottom':'center'}[side];end=np.array([x,y]);labelax.plot([xy[i,0],x],[xy[i,1],y],c='#777777',lw=.3,zorder=1);labelax.text(x,y,row.L1,ha=ha,va='center',fontsize=6.1,bbox=dict(facecolor='white',edgecolor='none',pad=.35),zorder=2)
 anchors.loc[i,'new_label_x_mm']=x;anchors.loc[i,'new_label_y_mm']=y;anchors.loc[i,'side']=side
anchors.to_csv(O/'source_data/Figure1b_all16_direct_labels.tsv',sep='\t',index=False)
# c/f preserve exact cell counts and unique study counts; d uses original source table.
counts=D.level1_annotation.value_counts().reindex(P.L1);studies=D.groupby('level1_annotation').study_canonical.nunique().reindex(P.L1);labels=[x.replace(' cell','') for x in P.L1];labels=[x.replace('Natural killer','Natural killer') for x in labels]
for x,w,val,cap,panel in [(32,36,counts,'Cells per lineage','c'),(85,20,studies,'Studies per lineage','f')]:
 ax=axmm(x,69,w,30);yy=np.arange(16);ax.barh(yy,val,color=P.hex,height=.65);ax.invert_yaxis();ax.set_ylim(15.6,-.6);ax.set_yticks(yy);ax.set_yticklabels(labels if panel=='c' else ['']*16,fontsize=6);ax.tick_params(axis='y',length=0,pad=2);ax.tick_params(axis='x',labelsize=5.7,length=2);ax.set_xlabel(cap,fontsize=6,labelpad=1)
 ax.spines[['left','right','top']].set_visible(False)
 if panel=='c':ax.set_xticks([0,100000,200000],['0','100k','200k']);ax.set_xlim(0,280000)
 else:ax.set_xticks([0,13,26]);ax.set_xlim(0,28)
 for j,v in enumerate(val):ax.text((val.max()*1.06),j,f'{v:,}' if panel=='c' else str(v),va='center',fontsize=5.7)
 letter(1 if panel=='c' else 80,67,panel)
df=pd.read_csv(S/'Figure1/Fig1d.tsv',sep='\t').set_index('tissue');ax=axmm(145,69,32,30);left=np.zeros(len(df))
for k in P.L1:ax.barh(np.arange(len(df)),df[k],left=left,color=pal[k],height=.8);left+=df[k]
ax.invert_yaxis();ax.set_yticks(range(len(df)),df.index,fontsize=5.7);ax.tick_params(axis='y',length=0,pad=2);ax.tick_params(axis='x',labelsize=5.7,length=2);ax.set(xlim=(0,1));ax.set_xlabel('Mean sample fraction',fontsize=6,labelpad=1);ax.spines[['right','top','left']].set_visible(False);letter(112,67,'d')
# e frozen SVG objects added after native plot export; no recomputation of subwindows.
letter(1,106,'e')
# g major lineage canonical markers: exact complete cross-lineage values from existing source.
full=pd.read_csv(S/'ExtendedData2/ED2_all_selected_markers.tsv',sep='\t');states=[s for s in P.L1 if s not in ['Unknown','Schwann cell']];genes='PECAM1 VWF DCN LUM MYH11 ACTA2 FN1 TCF21 RGS5 ABCC9 C1QA CD68 FCN1 S100A9 FCGR3B CSF3R CD1C FCER1A CD3D TRAC NKG7 GNLY MS4A1 CD79A MZB1 JCHAIN TPSAB1 KIT'.split();d=full[full.level1_annotation.isin(states)&full.gene_id.isin(genes)].copy();assert len(d)==len(states)*len(genes)
d.to_csv(O/'source_data/Figure1g_major_lineage_marker_values.tsv',sep='\t',index=False)
v=d.pivot(index='level1_annotation',columns='gene_id',values='mean_log1p').reindex(index=states,columns=genes);f=d.pivot(index='level1_annotation',columns='gene_id',values='detection_fraction').reindex(index=states,columns=genes);vmax=float(full.mean_log1p.max());ax=axmm(37,137,125,24);X,Y=np.meshgrid(range(28),range(14));size=(.72*min(125*mm/28,24*mm/14))**2;valid=v.notna()&f.notna()&f.gt(0);ax.scatter(X[valid],Y[valid],s=(size*f.to_numpy())[valid],c=v.to_numpy()[valid],cmap='viridis',norm=colors.Normalize(0,vmax),linewidths=.12,edgecolors='#555555');ax.set(xlim=(-.5,27.5),ylim=(13.5,-.5));ax.set_xticks(range(28),genes,rotation=90,fontsize=6.0,fontstyle='italic');ax.set_yticks(range(14),[s.replace(' cell','') for s in states],fontsize=6.0);ax.tick_params(length=0,pad=2);ax.tick_params(axis='y',pad=9);ax.tick_params(axis='x',pad=-.3)
for sp in ax.spines.values():sp.set_visible(False)
for j,s in enumerate(states):ax.add_patch(Rectangle((-.95,j-.42),.14,.84,fc=pal[s],ec='none',clip_on=False))
letter(1,135,'g');cbax=axmm(165,139,2,15);cb=fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(0,vmax),cmap='viridis'),cax=cbax);cb.ax.tick_params(labelsize=5.8,length=2);ftext(164,135,'Mean log1p',5.8);ftext(164,137,'(CP10k)',5.8)
ftext(38,134,'Detected cells (%)',5.8)
for j,z in enumerate([.25,.5,.75,1]):labelax.scatter(67+j*13,134,s=size*z,c='#75898a',lw=.1);ftext(69+j*13,134,str(int(z*100)),5.8)
# Export editable text + rasterized points; append unchanged lineage view objects by translation only.
buf=io.BytesIO();fig.savefig(buf,format='svg');plt.close(fig);r=ET.fromstring(buf.getvalue());src=ET.parse('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/reproducibility_dependencies/Panvascular_Marker_Communication_Enrichment_20260920/editable_svg/Figure1.svg').getroot();members={'axes_7','axes_8','axes_9','axes_10','text_186','text_187','text_188'};allids={e.get('id') for e in src.iter() if e.get('id')};idmap={i:'f1e_'+i for i in allids}
wrap=ET.SubElement(r,q('g'),{'id':'preserved_e','transform':f'translate(0 {-34*mm})'})
for e in src.iter():
 if e.get('id') in members:
  c=copy.deepcopy(e)
  for n in c.iter():
   for k,val in list(n.attrib.items()):
    if k=='id':n.set(k,idmap[val])
    else:
     for old,new in idmap.items():val=val.replace('url(#'+old+')','url(#'+new+')')
     n.set(k,val)
  
  wrap.append(c)
# Include referenced image definitions/clips after namespacing.
for defs in src.findall(q('defs')):
 c=copy.deepcopy(defs)
 for n in c.iter():
  if n.get('id') in idmap:n.set('id',idmap[n.get('id')])
 r.insert(0,c)
export(r,'Figure1')
(O/'notes/Figure1_checks.json').write_text(json.dumps({'cells':len(D),'sample_mapped':D.sample_uid.nunique(),'tissue_groups':D.Tissue.nunique(),'l1_labels':16,'line_crossings':int(score(assign)[1]),'marker_rows':len(states),'marker_columns':len(genes),'marker_cells':len(d),'marker_vmax':vmax,'all_expression_statistics_reused':True,'e_svg_translation_only_mm':-34,'new_panel':'g','visual_review':'PENDING'},indent=2))
print('DONE',score(assign))
