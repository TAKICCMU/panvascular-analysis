from pathlib import Path
import os,sys,json,copy,hashlib
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import numpy as np,pandas as pd,matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
import fitz,cairosvg,xml.etree.ElementTree as E
import cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
O=Path(__file__).resolve().parents[1];B=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/reproducibility_dependencies/Panvascular_Figure5_AB_Finalization_20260920');F=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920/reproducibility_dependencies/Panvascular_Figure5_Evidence_Focus_20260920');S=Path('/Users/taki/Desktop/投稿包_最新版_20260919/Source_Data_Figures');R=Path('/Users/taki/Desktop/Codex任务1/spatial_visium_geo_local/extracted');D=O/'source_data/Figure5';D.mkdir(exist_ok=True)
read=lambda p:pd.read_csv(p,sep='\t',float_precision='round_trip')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
P=json.loads((F/'notes/display_parameters.json').read_text());z1=read(B/'Source_Data/Figure5/CA1_existing_spot_values.tsv')
section='GSE314851::GSM9414346_T1077';stem=section.split('::')[1];d=R/'GSE314851';cols=['spot_id','section_id','library_size_original','x','y','expr__FN1','abundance__Fibroblast']
z2=pd.concat([c[c.section_id.eq(section)] for c in pd.read_csv(S/'Figure5/spot_joint_source.tsv.gz',sep='\t',usecols=cols,chunksize=20000,float_precision='round_trip')],ignore_index=True)
pos=pd.read_csv(d/f'{stem}_tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','fullres_y','fullres_x']);z2['barcode']=z2.spot_id.str.split('::').str[-1];z2=z2.merge(pos,on='barcode',validate='one_to_one');assert len(z2)==1363 and z2.in_tissue.eq(1).all() and z2.library_size_original.ge(500).all();assert np.array_equal(z2.x,z2.fullres_x) and np.array_equal(z2.y,z2.fullres_y)
sf=json.loads((d/f'{stem}_scalefactors_json.json').read_text());z2['hires_x']=z2.x*sf['tissue_hires_scalef'];z2['hires_y']=z2.y*sf['tissue_hires_scalef'];z2.to_csv(D/'T1077_existing_spot_values.tsv',sep='\t',index=False)
# Display-only limits for the added section: same established positive 2–98% convention, no new models.
P2={'section':section,'image':str(d/f'{stem}_tissue_hires_image.png'),'limits':{'FN1':np.quantile(z2.loc[z2.expr__FN1.gt(0),'expr__FN1'],[.02,.98]).tolist(),'Fib_abundance':[float(z2.abundance__Fibroblast.min()),float(z2.abundance__Fibroblast.max())]},'scalefactors':sf}
sets=[(z1,P,'Anterior tibial artery · CA1'),(z2,P2,'Coronary artery · T1077')]
U=read(B/'Source_Data/Figure5/Figure5b_60_gene_unit_estimates.tsv');meta=read(B/'Source_Data/Figure5/all29_original_meta.tsv');M=meta[meta.gene.isin(['FN1','SCARB1'])&meta.posterior.eq('q05')&meta.term.eq('cell_abundance')&meta.sensitivity.eq('all_sections')].set_index('gene');assert len(U)==60 and U.biological_unit_id.nunique()==34
mpl.rcParams.update({'font.family':'Arial','font.size':6,'text.color':'#222222','axes.labelcolor':'#222222','axes.linewidth':.45,'pdf.fonttype':42,'svg.fonttype':'none','xtick.major.width':.4,'ytick.major.width':.4,'xtick.major.size':2,'ytick.major.size':0})
cm=mpl.colors.LinearSegmentedColormap.from_list('retained_spatial',['#482878','#3B528B','#2196C4','#35B779','#FDE725','#F68B32','#B2182B']);ZERO='#8B8F95';COL={'FN1':'#DB8C27','SCARB1':'#278F90'}
def ax(f,x,y,w,h,W=183,H=170):return f.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def tx(f,x,y,t,**kw):return f.text(x/183,1-y/170,t,va='top',**kw)
def spatial(a,z,spec,key,point=.65):
 im=Image.open(spec['image']);W,H=im.size;a.set_xlim(0,W);a.set_ylim(H,0);a.set_aspect('equal');a.set_xticks([]);a.set_yticks([]);a.spines[:].set_visible(False)
 if key=='HE':a.imshow(im,origin='upper',extent=[0,W,H,0],interpolation='nearest');return None
 a.set_facecolor('#050708');v=z['abundance__Fibroblast' if key=='Fib_abundance' else 'expr__FN1'].to_numpy();fin=np.isfinite(v);zero=fin&(v==0);positive=fin&(v!=0);lo,hi=spec['limits'][key];norm=mpl.colors.Normalize(lo,hi,clip=True);pal='viridis' if key=='Fib_abundance' else cm
 a.scatter(z.hires_x[zero],z.hires_y[zero],s=point,c=ZERO,lw=0);a.scatter(z.hires_x[positive],z.hires_y[positive],s=point,c=v[positive],cmap=pal,norm=norm,lw=0)
 if (~fin).any():a.scatter(z.hires_x[~fin],z.hires_y[~fin],s=point*3,marker='x',c='white',lw=.4)
 return mpl.cm.ScalarMappable(norm=norm,cmap=pal)
f=plt.figure(figsize=(183/25.4,170/25.4));geoms=[]
for row,(z,spec,identity) in enumerate(sets):
 top=7.5+30*row
 tx(f,53,top-7.1,identity,ha='center',fontsize=6.2)
 for x,key,label in [(5,'HE','H&E'),(37,'Fib_abundance','Inferred Fib abundance'),(73,'FN1','FN1')]:
  a=ax(f,x,top,24,21);tx(f,x+12,top-3.1,label,ha='center',fontsize=6);sm=spatial(a,z,spec,key)
  if sm is not None:
   ca=ax(f,x+25,top+2,1.3,16);cb=f.colorbar(sm,cax=ca,orientation='vertical');lo,hi=spec['limits'][key];cb.set_ticks([lo,hi]);cb.set_ticklabels([f'{lo:.2g}',f'{hi:.2g}'] if key=='Fib_abundance' else [f'≤{lo:.2g}',f'≥{hi:.2g}']);cb.ax.tick_params(labelsize=5.7,pad=1,length=1.3);cb.ax.yaxis.set_ticks_position('right')
  f.canvas.draw();bb=a.get_position();geoms.append({'section':spec['section'],'key':key,'field':[0,Image.open(spec['image']).width,Image.open(spec['image']).height,0],'data_rect_mm':[bb.x0*183,(1-bb.y1)*170,bb.width*183,bb.height*170],'spots':len(z),'label_center_mm':x+12})
# Units are common to both rows; discrete zero is separate from positive color encoding.
f.text(49/183,1-59/170,'Σ marginal q05',ha='center',va='top',fontsize=5.8);f.text(85/183,1-59/170,'log1p(CP10k)',ha='center',va='top',fontsize=5.8)
f.legend(handles=[Line2D([],[],marker='o',ls='',mfc=ZERO,mec='none',ms=3,label='Zero')],loc='upper left',bbox_to_anchor=(4/183,1-59/170),frameon=False,borderaxespad=0,handletextpad=.3,fontsize=5.8)
a=ax(f,123,10,39,38);studies=sorted(U.GSE.unique());points=[]
for gi,g in enumerate(['FN1','SCARB1']):
 for i,study in enumerate(studies):
  q=U[U.gene.eq(g)&U.GSE.eq(study)].sort_values('biological_unit_id')
  if q.empty:continue
  offs=np.linspace(-.15,.15,len(q)) if len(q)>1 else np.array([0.]);yy=i+(-.22 if gi==0 else .22)+offs
  a.scatter(q.beta,yy,s=4.8,c=COL[g],marker='o' if gi==0 else '^',linewidths=0,zorder=3)
  for (_,r),y in zip(q.iterrows(),yy):points.append({'gene':g,'GSE':study,'biological_unit_id':r.biological_unit_id,'beta':r.beta,'display_y':float(y)})
 q=M.loc[g];yy=7.3+gi*.75;a.errorbar(q.beta,yy,xerr=[[q.beta-q.CI_low],[q.CI_high-q.beta]],fmt='D',ms=3,color=COL[g],lw=.8,capsize=0)
 py=10+38*(yy+.65)/9.25;fdr=q.FDR_BH_all_29_pairs_x_3_terms;label=f'{fdr:.2e}' if fdr<.001 else f'{fdr:.4f}'
 f.text(167.2/183,1-py/170,label,ha='center',va='center',fontsize=5.8);f.text(178/183,1-py/170,f'{q.I2:.1f}%',ha='center',va='center',fontsize=5.8)
a.axvline(0,lw=.55,c='#9BA1A4',zorder=0);a.axhline(6.65,lw=.4,c='#CFD3D5');a.set_ylim(8.6,-.65);a.set_xlim(-1,.55);a.set_yticks(list(range(7))+[7.3,8.05]);a.set_yticklabels(studies+['FN1','SCARB1']);a.tick_params(axis='y',pad=3,labelsize=6);a.set_xticks([-1,-.5,0,.5]);a.set_xlabel('Conditional Fib-abundance coefficient',labelpad=2,fontsize=6);a.spines[['top','right','left']].set_visible(False)
tx(f,167.2,38.7,'FDR',ha='center',fontsize=5.8);tx(f,178,38.7,'I²',ha='center',fontsize=5.8)
zero_x=123+39/1.55
# Separate identity legends on the requested sides of the zero line.
for g,x,mark in [('FN1',133,'o'),('SCARB1',163,'^')]:
 ll=f.legend(handles=[Line2D([],[],marker=mark,ls='',mfc=COL[g],mec='none',ms=3.5,label=g)],loc='upper center',bbox_to_anchor=(x/183,1-3/170),frameon=False,borderaxespad=0,handletextpad=.4,fontsize=6)
 f.add_artist(ll)
leg=ax(f,116,56,6,2);leg.plot([.15,.85],[.5,.5],c='#555',lw=.8);leg.plot(.5,.5,'D',ms=3,c='#555');leg.set(xlim=(0,1),ylim=(0,1));leg.axis('off');tx(f,124,56,'Pooled estimate and 95% CI',fontsize=6)
f.savefig(O/'review/Figure5_new_ab.svg');plt.close(f)
NS='http://www.w3.org/2000/svg';N={'s':NS};E.register_namespace('',NS);E.register_namespace('xlink','http://www.w3.org/1999/xlink');tree=E.parse(B/'figures/Figure5.svg');root=tree.getroot();fg=root.find("s:g[@id='figure_1']",N);old=fg.find("s:g[@id='Figure5_ab_finalization']",N);assert old is not None;fg.remove(old)
new=E.parse(O/'review/Figure5_new_ab.svg').getroot();ids={el.get('id'):'r2_'+el.get('id') for el in new.iter() if el.get('id')}
for el in new.iter():
 for k,v in list(el.attrib.items()):
  if k=='id':el.set(k,ids[v]);continue
  for old,nn in ids.items():v=v.replace('url(#'+old+')','url(#'+nn+')')
  if v.startswith('#') and v[1:] in ids:v='#'+ids[v[1:]]
  el.set(k,v)
wrap=E.SubElement(fg,'{'+NS+'}g',{'id':'Figure5_two_sections_ab'})
for el in new.find("s:g[@id='r2_figure_1']",N):
 if el.get('id')!='r2_patch_1':wrap.append(copy.deepcopy(el))
for defs in new.findall('s:defs',N):root.append(copy.deepcopy(defs))
svg=O/'editable_svg/Figure5.svg';tree.write(svg,encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(svg),write_to=str(O/'figures/Figure5.pdf'))
before=fitz.open(B/'figures/Figure5.pdf');after=fitz.open(O/'figures/Figure5.pdf');clip=fitz.Rect(0,66*72/25.4,before[0].rect.width,before[0].rect.height);same=before[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=clip).samples==after[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=clip).samples;assert same
for name,doc in [('before',before),('after',after)]:doc[0].get_pixmap(matrix=fitz.Matrix(3,3)).save(O/f'review/Figure5_{name}.png')
# Full-resolution companion preserves both complete fields and exact original values.
from matplotlib.backends.backend_pdf import PdfPages
with PdfPages(O/'companion/Figure5_two_complete_sections.pdf') as pdf:
 for z,spec,identity in sets:
  f,axs=plt.subplots(1,3,figsize=(183/25.4,85/25.4));f.subplots_adjust(left=.01,right=.96,bottom=.2,top=.82,wspace=.3);f.suptitle(identity,fontsize=8,y=.97)
  for a,key,label in zip(axs,['HE','Fib_abundance','FN1'],['H&E','Inferred Fib abundance','FN1 log1p(CP10k)']):
   sm=spatial(a,z,spec,key,point=2);a.set_title(label,fontsize=7)
   if sm is not None:cb=f.colorbar(sm,ax=a,orientation='horizontal',fraction=.05,pad=.04);cb.ax.tick_params(labelsize=6)
  f.text(.5,.02,f"{spec['section']} · {len(z):,} original included spots · full field",ha='center',fontsize=6);pdf.savefig(f);plt.close(f)
sec=read(B/'Source_Data/Figure5/all29_original_section_models.tsv');q=sec[sec.section_id.eq(section)&sec.gene.eq('FN1')&sec.posterior.eq('q05')&sec.term.eq('cell_abundance')];q.to_csv(D/'T1077_original_section_coefficient.tsv',sep='\t',index=False)
pd.DataFrame(points).to_csv(D/'Figure5b_unchanged_point_mapping.tsv',sep='\t',index=False)
records={'original':P,'added':P2,'selection':'Coronary vascular bed added to anterior tibial example; T1077 original complete H&E inspected, no QC fallback, 1363 original QC spots, negative FN1 conditional coefficient consistent with overall. Selection after results review, not blinded, not minimum P. Tissue identity from frozen section manifest. Unit is a moderate-confidence library ID, not confirmed donor.','views':geoms,'gene_unit_estimates':60,'distinct_recorded_unit_IDs':34,'zero_line_page_x_mm':zero_x,'FN1_legend_center_mm':133,'SCARB1_legend_center_mm':163,'c_f_render_identical_below66mm':same,'ORA':'PENDING','histology':'Not independently confirmed','science_recomputed':False,'visual_review':'PENDING'}
(O/'notes/Figure5_checks.json').write_text(json.dumps(records,indent=2));print('Figure5 exported',same)
