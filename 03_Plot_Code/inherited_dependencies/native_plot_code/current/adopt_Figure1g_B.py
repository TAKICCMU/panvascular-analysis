from common_main import *
from matplotlib.patches import Rectangle
import shutil
PK=Path('/Users/taki/Desktop/投稿包_最新版_20260919');BK=PK/'图件编辑资料/轻量备份/Figure1g_B采用前_20260921';BK.mkdir(parents=True,exist_ok=True)
# Only affected artifacts; never duplicate the package.
for rel in ['Figures/Figure1.pdf','Figures/Figure5.pdf','图件编辑资料/editable_svg/Figure1.svg','图件编辑资料/editable_svg/Figure5.svg','全部图件_15页审阅.pdf','图件采用版本.tsv','文件来源与校验.tsv']:
 src=PK/rel;dst=BK/rel
 if not dst.exists():dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
for rel in ['FIGURE_VERSIONS.tsv','SOURCE_DATA_INDEX.tsv','本轮修改说明.md']:
 dst=BK/'work_records'/rel
 if not dst.exists():dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(O/rel,dst)
T=O/'marker_audit_20260921/tables'
M=pd.read_csv(T/'original_mean_log1p.tsv',sep='\t',index_col=0,float_precision='round_trip')
F=pd.read_csv(T/'original_detection_fraction.tsv',sep='\t',index_col=0,float_precision='round_trip')
Z=pd.read_csv(T/'per_gene_scaled_mean.tsv',sep='\t',index_col=0,float_precision='round_trip')
assert M.shape==F.shape==Z.shape==(14,28) and M.index.equals(Z.index) and M.columns.equals(Z.columns) and F.index.equals(Z.index) and F.columns.equals(Z.columns)
assert np.nanmin(Z)>=0 and np.nanmax(Z)<=1
# Read approved colours; no recomputation of expression, detection or min-max.
root=ET.parse(BK/'图件编辑资料/editable_svg/Figure1.svg').getroot();before=copy.deepcopy(root)
remove={'axes_8','axes_9',*[f'text_{i}' for i in range(154,161)],*[f'b0921_PathCollection_{i}' for i in range(1,5)]};remove_ids(root,remove)
f=plt.figure(figsize=(183/25.4,170/25.4));x,y,w,h=37,136,110,25.3
ax=f.add_axes([x/183,1-(y+h)/170,w/183,h/170]);xx,yy=np.meshgrid(np.arange(28),np.arange(14));maxdiam=.76*min(w*mm/28,h*mm/14);size=maxdiam**2
valid=np.isfinite(Z.to_numpy())&np.isfinite(F.to_numpy())&F.to_numpy().astype(float).__gt__(0)
ax.scatter(xx[valid],yy[valid],s=(size*F.to_numpy())[valid],c=Z.to_numpy()[valid],cmap='viridis',norm=colors.Normalize(0,1),linewidths=0)
missing=~np.isfinite(M.to_numpy())|~np.isfinite(F.to_numpy());ax.scatter(xx[missing],yy[missing],marker='x',s=5,c='#777',linewidths=.4)
# Neutral glyphs for zero-range genes only if identified by the previous audit.
ranges=pd.read_csv(T/'gene_range.tsv',sep='\t')
for gene in ranges.loc[ranges.zero_range.eq(True),'gene']:
 j=M.columns.get_loc(gene);good=np.isfinite(F.iloc[:,j])&F.iloc[:,j].gt(0);ax.scatter(np.full(good.sum(),j),np.flatnonzero(good),s=size*F.iloc[:,j][good],c='#bbbbbb',linewidths=0)
ax.set(xlim=(-.5,27.5),ylim=(13.5,-.5));ax.set_xticks(range(28),M.columns,rotation=90,fontsize=6,fontstyle='italic');ax.set_yticks(range(14),[s.replace(' cell','') for s in M.index],fontsize=6);ax.tick_params(length=0,pad=1.5);ax.tick_params(axis='y',pad=8);ax.tick_params(axis='x',pad=-1)
for s in ax.spines.values():s.set_visible(False)
for j,st in enumerate(M.index):ax.add_patch(Rectangle((-.98,j-.40),.15,.80,fc=pal[st],ec='none',clip_on=False))
cbax=f.add_axes([155/183,1-159/170,2/183,18/170]);cb=f.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(0,1),cmap='viridis'),cax=cbax);cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6,length=2)
f.text(154/183,1-134.8/170,'Per-gene scaled',fontsize=6,va='center');f.text(154/183,1-137.3/170,'mean expression',fontsize=6,va='center')
# Same area mapping for plotted points and legend.
handles=[ax.scatter([],[],s=size*v,c='#75898a',lw=0,label=str(int(v*100))) for v in [.25,.5,.75,1]]
f.legend(handles=handles,loc='upper left',bbox_to_anchor=(164/183,1-144/170),ncol=1,frameon=False,fontsize=5.8,handletextpad=.35,labelspacing=.55,borderaxespad=0)
f.text(165/183,1-140/170,'Detected',fontsize=5.8,va='center')
f.text(165/183,1-142.3/170,'cells (%)',fontsize=5.8,va='center')
embed_plot(root,f,'adopted_B_g_');export(root,'Figure1')
# Everything outside g is identical; axes_4 contains the old g-size key only in its last four collections.
find=lambda r,i:next(e for e in r.iter() if e.get('id')==i)
for id in ['axes_1','axes_2','axes_3','axes_5','axes_6','axes_7','preserved_e']:
 assert ET.tostring(find(root,id))==ET.tostring(find(before,id)),id
b4=copy.deepcopy(find(before,'axes_4'));remove_ids(b4,remove);assert ET.tostring(b4)==ET.tostring(find(root,'axes_4'))
# Section IDs in existing H&E headers identify the entire corresponding row; no space reallocation.
r=ET.parse(BK/'图件编辑资料/editable_svg/Figure5.svg').getroot();r0=copy.deepcopy(r);targets=[]
for g in r.iter(q('g')):
 if g.get('id') in ['r2_text_24','r2_text_28']:
  targets.append(g)
assert len(targets)==2
for g,lab in zip(targets,['H&E · GSE248608 / CA1','H&E · GSE314851 / T1077']):
 ts=list(g.iter(q('text')));assert len(ts)==1 and ts[0].text=='H&E';ts[0].text=lab
export(r,'Figure5')
# Prove all non-label SVG objects remain identical, including pixels and spot artists.
for rr in [r,r0]:
 for g in rr.iter(q('g')):
  if g.get('id') in ['r2_text_24','r2_text_28']:
   for t in g.iter(q('text')):t.text='H&E'
assert ET.tostring(r)==ET.tostring(r0)
(O/'notes/Figure1g_B_adoption_checks.json').write_text(json.dumps({'display':'B_ADOPTED','matrix_shape':[14,28],'source_tables':{str(T/f):sha(T/f) for f in ['original_mean_log1p.tsv','original_detection_fraction.tsv','per_gene_scaled_mean.tsv']},'scaled_values_read_not_recomputed':True,'area_linear_detection':True,'zero_area_at_zero':True,'max_dot_diameter_pt':maxdiam,'row_pitch_pt':h*mm/14,'column_pitch_pt':w*mm/28,'native_export_mm':[183,170],'no_174x83_page_scaling':True,'a_f_unchanged':True,'Figure5_only_two_HE_header_strings_changed':True,'annotation_validated':False,'measurement_coverage_newly_validated':False,'ORA':'PENDING'},indent=2))
print('Figure1 B and Figure5 row IDs exported')
