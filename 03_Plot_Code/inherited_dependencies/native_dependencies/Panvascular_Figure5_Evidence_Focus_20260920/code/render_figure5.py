from pathlib import Path
import os,sys,copy,json,hashlib
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import numpy as np,pandas as pd,matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
import fitz,cairosvg,xml.etree.ElementTree as E
import cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
O=Path(__file__).resolve().parents[1];B=Path('/Users/taki/Desktop/Panvascular_Marker_Spatial_Review_20260919');D=O/'source_data';P=json.loads((O/'notes/display_parameters.json').read_text());z=pd.read_csv(D/'CA1_existing_spot_values.tsv',sep='\t');U=pd.read_csv(D/'focus_all_contributing_units.tsv',sep='\t');M=pd.read_csv(D/'focus_summary.tsv',sep='\t').set_index('gene');im=Image.open(P['image']);W,H=im.size
mpl.rcParams.update({'font.family':'Arial','font.size':6,'text.color':'#222222','axes.labelcolor':'#222222','axes.linewidth':.45,'pdf.fonttype':42,'svg.fonttype':'none','xtick.major.width':.4,'ytick.major.width':.4,'xtick.major.size':2,'ytick.major.size':0})
cm=mpl.colors.LinearSegmentedColormap.from_list('retained_spatial',['#482878','#3B528B','#2196C4','#35B779','#FDE725','#F68B32','#B2182B']);ZERO='#8B8F95';COL={'FN1':'#DB8C27','SCARB1':'#278F90'};records=[]
def axes(f,x,y,w,h):return f.add_axes([x/183,1-(y+h)/170,w/183,h/170])
def tx(f,x,y,t,**kw):return f.text(x/183,1-y/170,t,va='top',**kw)
def mapplot(f,box,key,label,cbbox=None,size=.6):
 a=axes(f,*box);a.set_xlim(0,W);a.set_ylim(H,0);a.set_aspect('equal');a.set_xticks([]);a.set_yticks([]);a.spines[:].set_visible(False)
 if key=='HE':a.imshow(im,origin='upper',extent=[0,W,H,0],interpolation='nearest')
 else:
  a.set_facecolor('#050708');v=z['abundance__Fibroblast' if key=='Fib_abundance' else 'expr__'+key].to_numpy();finite=np.isfinite(v);zero=finite&(v==0);normal=finite&~zero;lo,hi=P['limits'][key];norm=mpl.colors.Normalize(lo,hi,clip=True);palette='viridis' if key=='Fib_abundance' else cm
  a.scatter(z.hires_x[zero],z.hires_y[zero],s=size,c=ZERO,lw=0)
  a.scatter(z.hires_x[normal],z.hires_y[normal],s=size,c=v[normal],cmap=palette,norm=norm,lw=0)
  if (~finite).any():a.scatter(z.hires_x[~finite],z.hires_y[~finite],marker='x',s=size*3,c='white',lw=.4)
  if cbbox:
   ca=axes(f,*cbbox);cb=f.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=palette),cax=ca,orientation='horizontal');cb.set_ticks([lo,hi]);cb.set_ticklabels([f'{lo:.2g}',f'{hi:.2g}'] if key=='Fib_abundance' else [f'≤{lo:.2g}',f'≥{hi:.2g}']);cb.ax.tick_params(labelsize=6,pad=1)
 tx(f,box[0],box[1]-3,label,fontsize=6)
 f.canvas.draw();b=a.get_position();records.append({'key':key,'allocated_mm':box,'data_mm':[b.x0*183,(1-b.y1)*170,b.width*183,b.height*170],'field':[0,W,H,0],'spot_count':1779 if key!='HE' else None,'interpolation':'nearest image only; no spot interpolation'})
 return a
# Main replacements: all untouched c-f groups remain byte-identical XML elements.
f=plt.figure(figsize=(183/25.4,170/25.4))
tx(f,7,2.4,'GSE248608 · GSM7917894_CA1',fontsize=6.2)
# 2x2 shared-context display, original full field and equal physical scale.
mapplot(f,[7,9,35,19.5],'HE','H&E')
mapplot(f,[54,9,35,19.5],'Fib_abundance','Fib abundance (q05)',[60,29.0,23,1.2])
mapplot(f,[7,37,35,19.5],'FN1','FN1 · log1p(CP10k)',[13,57,23,1.2])
mapplot(f,[54,37,35,19.5],'SCARB1','SCARB1 · log1p(CP10k)',[60,57,23,1.2])
# Data units are points, NOT per-study pooled effects. All 60 contributions shown.
a=axes(f,122,10,54,38);studies=sorted(U.GSE.unique());coord=[]
for gi,g in enumerate(['FN1','SCARB1']):
 for i,study in enumerate(studies):
  q=U[U.gene.eq(g)&U.GSE.eq(study)].sort_values('biological_unit_id')
  if q.empty:continue
  offs=np.linspace(-.15,.15,len(q)) if len(q)>1 else np.array([0.]);yy=i+(-.22 if gi==0 else .22)+offs
  a.scatter(q.beta,yy,s=4.8,c=COL[g],marker='o' if gi==0 else '^',linewidths=0,zorder=3)
  for (_,r),y in zip(q.iterrows(),yy):coord.append({'gene':g,'unit':r.biological_unit_id,'beta':r.beta,'y_display':float(y),'study':study})
 q=M.loc[g];yy=7.3+gi*.75
 a.errorbar(q.overall_beta,yy,xerr=[[q.overall_beta-q.overall_CI_low],[q.overall_CI_high-q.overall_beta]],fmt='D',ms=3,color=COL[g],lw=.8,capsize=0)
a.axvline(0,lw=.55,c='#9BA1A4',zorder=0);a.axhline(6.65,lw=.4,c='#CFD3D5');a.set_ylim(8.6,-.65);a.set_xlim(-1,.55);a.set_yticks(list(range(7))+[7.3,8.05]);a.set_yticklabels(studies+['Pooled FN1','Pooled SCARB1']);a.tick_params(axis='y',pad=3,labelsize=6);a.set_xticks([-1,-.5,0,.5]);a.set_xlabel('Conditional Fib-abundance coefficient',labelpad=2,fontsize=6);a.spines[['top','right','left']].set_visible(False)
f.legend(handles=[Line2D([],[],marker='o',ls='',ms=3,mfc=COL['FN1'],mec='none',label='FN1 (34 units)'),Line2D([],[],marker='^',ls='',ms=3.5,mfc=COL['SCARB1'],mec='none',label='SCARB1 (26 units)')],loc='upper left',bbox_to_anchor=(108/183,1-2.5/170),ncol=2,frameon=False,fontsize=6,borderaxespad=0,columnspacing=.8,handletextpad=.3)
tx(f,108,56.5,'Points: recorded units; diamonds: pooled 95% CI',fontsize=6)
# zero symbol for the expression maps, separate from the abundance encoding.
f.legend(handles=[Line2D([],[],marker='o',ls='',ms=3,mfc=ZERO,mec='none',label='Expression zero')],loc='upper left',bbox_to_anchor=(7/183,1-31/170),frameon=False,fontsize=6,borderaxespad=0,handletextpad=.3)
f.savefig(O/'review/new_ab_components.svg');plt.close(f)
NS='http://www.w3.org/2000/svg';N={'s':NS};E.register_namespace('',NS);E.register_namespace('xlink','http://www.w3.org/1999/xlink')
tree=E.parse(B/'editable_svg/Figure5.svg');root=tree.getroot();fg=root.find("s:g[@id='figure_1']",N)
# Preserve the former b in a separate, genuine editable SVG/PDF, not an overwritten or lost result.
oldb=copy.deepcopy(root);ofg=oldb.find("s:g[@id='figure_1']",N)
for el in list(ofg):
 if el.get('id') not in ['axes_5','axes_6','text_132']:ofg.remove(el)
# Crop by viewBox only; this companion is not a new ED panel.
oldb.set('viewBox',f'{100*72/25.4} 0 {83*72/25.4} {61*72/25.4}');oldb.set('width','83mm');oldb.set('height','61mm');E.ElementTree(oldb).write(O/'review/former_Figure5b_detection_matrix.svg',encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(O/'review/former_Figure5b_detection_matrix.svg'),write_to=str(O/'review/former_Figure5b_detection_matrix.pdf'))
for gid in ['Figure5a_separate_HE_spots','axes_5','axes_6']:
 el=fg.find(f"s:g[@id='{gid}']",N);assert el is not None;fg.remove(el)
new=E.parse(O/'review/new_ab_components.svg').getroot();ids={el.get('id'):'focus_'+el.get('id') for el in new.iter() if el.get('id')}
for el in new.iter():
 for k,val in list(el.attrib.items()):
  if k=='id':el.set(k,ids[val]);continue
  for old,nn in ids.items():val=val.replace('url(#'+old+')','url(#'+nn+')')
  if val.startswith('#') and val[1:] in ids:val='#'+ids[val[1:]]
  el.set(k,val)
wrap=E.SubElement(fg,'{'+NS+'}g',{'id':'Figure5_evidence_focus_ab'})
for el in new.find("s:g[@id='focus_figure_1']",N):
 if el.get('id')!='focus_patch_1':wrap.append(copy.deepcopy(el))
for defs in new.findall('s:defs',N):root.append(copy.deepcopy(defs))
svg=O/'figures/Figure5.svg';tree.write(svg,encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(svg),write_to=str(O/'figures/Figure5.pdf'))
# Detailed reviewer prototype: one complete unit forest plus exact same-section views per gene.
with PdfPages(O/'figures/Figure5_focus_prototype.pdf') as pdf:
 for g in ['FN1','SCARB1']:
  fig=plt.figure(figsize=(183/25.4,170/25.4));q=M.loc[g];d=U[U.gene.eq(g)].sort_values(['GSE','biological_unit_id']).reset_index(drop=True);n=len(d)
  tx(fig,5,3,g+' · conditional Fib-abundance association',fontsize=10,fontweight='bold');tx(fig,5,10,f'{int(q.units)} recorded units / {int(q.studies)} studies; original q05, all sections',fontsize=7)
  ax=axes(fig,63,22,47,115);yy=np.arange(n)
  ax.errorbar(d.beta,yy,xerr=d.se,fmt='o',ms=2,color='#657079',ecolor='#A4ADB3',lw=.55,capsize=0)
  labels=[r.biological_unit_id.replace('54yo_male_aortic_plaque','serial donor') for _,r in d.iterrows()]
  ax.set_yticks(yy);ax.set_yticklabels(labels,fontsize=6);ax.set_ylim(n-.3,-1);ax.axvline(0,color='#949A9D',lw=.55);ax.set_xlim(min(-1,float((d.beta-d.se).min())-.04),max(.6,float((d.beta+d.se).max())+.04));ax.spines[['top','right','left']].set_visible(False);ax.set_xlabel('Unit coefficient ± original SE',fontsize=6);ax.tick_params(axis='y',pad=3)
  for idx in range(1,n):
   if d.GSE.iloc[idx]!=d.GSE.iloc[idx-1]:ax.axhline(idx-.5,color='#E0E3E5',lw=.5,zorder=0)
  tx(fig,5,146,f'Pooled β {q.overall_beta:+.4f} (95% CI {q.overall_CI_low:+.4f} to {q.overall_CI_high:+.4f})',fontsize=7)
  tx(fig,5,151,f'Original FDR {q.overall_FDR:.3g}; I² {q.I2:.1f}%; {int(q.unit_positive)} positive / {int(q.unit_negative)} inverse units',fontsize=7)
  tx(fig,5,159,'Unit CIs were not stored; bars show SE, not 95% CI. GSE314851 IDs are libraries, not verified donors.',fontsize=6)
  for y,key,label in [(24,'HE','H&E'),(69,g,g+' · log1p(CP10k)'),(118,'Fib_abundance','Fib abundance (q05; inferred)')]:
   mapplot(fig,[127,y,48,36],key,label,None if key=='HE' else [133,y+37.5,36,1.5],size=1.4)
  tx(fig,121,13,'GSE248608 · GSM7917894_CA1',fontsize=6.3)
  # Section numbers clearly labelled; no pooled FDR adjacent to the tissue map.
  tx(fig,121,64,f'Section β {q.section_beta:+.3f}; nominal P={q.section_nominal_P:.3g}',fontsize=6)
  tx(fig,123,112,f'{int(q.n_positive):,}/{int(q.n_spots):,} positive spots; grey = zero',fontsize=6)
  fig.savefig(O/'figures'/f'Focus_{g}.svg');fig.savefig(O/'review'/f'Focus_{g}.png',dpi=180);pdf.savefig(fig);plt.close(fig)
# Scoped verification: untouched scientific panels retain exact element serialization.
old=E.parse(B/'editable_svg/Figure5.svg').getroot();checks={}
for gid in ['FigureYa6_native_forest','axes_8','axes_9','axes_10','legend_2','legend_3','legend_4']+[f'text_{i}' for i in range(136,157)]:
 a=old.find(f".//s:g[@id='{gid}']",N);b=root.find(f".//s:g[@id='{gid}']",N);checks[gid]=E.tostring(a)==E.tostring(b);assert checks[gid],gid
pd.DataFrame(coord).to_csv(D/'unit_point_display_coordinates.tsv',sep='\t',index=False)
for path in [O/'figures/Figure5.pdf']:
 doc=fitz.open(path);page=doc[0];assert abs(page.rect.width*25.4/72-183)<.01 and abs(page.rect.height*25.4/72-170)<.01;page.get_pixmap(matrix=fitz.Matrix(2.4,2.4)).save(O/'review/Figure5_after.png')
(O/'notes/render_checks.json').write_text(json.dumps({'Figure5_mm':[183,170],'all_60_unit_points':len(coord),'c_to_f_identical_SVG_elements':checks,'spatial_views':records,'new_fits':False,'new_P_FDR':False,'old_b_preserved':'review/former_Figure5b_detection_matrix.pdf','visual_review':'PENDING_ACTUAL_INSPECTION'},indent=2))
print('Figure5 + two-page quantitative/spatial prototype exported; all c-f groups unchanged.')
