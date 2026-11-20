from pathlib import Path
import os,sys,copy,json,hashlib,shutil
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
O=Path(__file__).resolve().parents[1];B=Path('/Users/taki/Desktop/Panvascular_Figure5_Evidence_Focus_20260920');S=Path('/Users/taki/Desktop/投稿包_最新版_20260919/Source_Data_Figures');D=O/'Source_Data/Figure5';P=json.loads((B/'notes/display_parameters.json').read_text());records=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')
source_map=[]
def keep(src,name,panels,meaning):
 dst=D/name;shutil.copyfile(src,dst);source_map.append(dict(panel=panels,file=str(dst.relative_to(O)),source=str(src),source_sha256=sha(src),output_sha256=sha(dst),meaning=meaning));return dst
keep(B/'source_data/CA1_existing_spot_values.tsv','CA1_existing_spot_values.tsv','5a; companion','All 1779 original spots; both FN1 and SCARB1 remain in data')
keep(B/'source_data/all29_original_meta.tsv','all29_original_meta.tsv','5b/c/d; P86 Table32','Unchanged all terms and settings, including original P/CI/FDR')
keep(B/'source_data/all29_original_unit_models.tsv','all29_original_unit_models.tsv','5b; companion','Original unit beta and SE; no stored unit CI')
keep(B/'source_data/all29_original_section_models.tsv','all29_original_section_models.tsv','companion','Original section estimates, including SCARB1 CA1')
keep(S/'Figure5/Fig5b.tsv','Former_5b_detection_fraction.tsv','former 5b supplementary data','Original nonzero-fraction matrix; no recalculation')
keep(B/'review/former_Figure5b_detection_matrix.pdf','Former_5b_detection_matrix.pdf','former 5b supplementary display','Preserved matrix PDF, no new ED number')
keep(B/'review/former_Figure5b_detection_matrix.svg','Former_5b_detection_matrix.svg','former 5b supplementary display','Editable preserved matrix')
keep(S/'Figure5/Fig5c.tsv','Fig5c_all29.tsv','5c','Frozen 29 original coefficients')
for name in ['module_coherence_meta.tsv','module_dataset_stats.tsv','module_program_direction_meta.tsv']:
 keep(S/'Figure5/module_spatial'/name,name,'5e/f','Unchanged module result source')
for name in ['Figure5_focus_prototype.pdf','Focus_FN1.svg','Focus_SCARB1.svg']:
 shutil.copyfile(B/'figures'/name,O/'companion'/name)
assert sha(B/'figures/Figure5_focus_prototype.pdf')==sha(O/'companion/Figure5_focus_prototype.pdf')
z=read(D/'CA1_existing_spot_values.tsv');full=read(D/'all29_original_unit_models.tsv');U=full[full.gene.isin(['FN1','SCARB1'])&full.posterior.eq('q05')&full.term.eq('cell_abundance')&full.sensitivity.eq('all_sections')].copy();meta=read(D/'all29_original_meta.tsv');M=meta[meta.gene.isin(['FN1','SCARB1'])&meta.posterior.eq('q05')&meta.term.eq('cell_abundance')&meta.sensitivity.eq('all_sections')].set_index('gene')
assert len(U)==60 and U.biological_unit_id.nunique()==34 and not U.duplicated(['gene','GSE','biological_unit_id']).any();assert U.groupby('gene').size().to_dict()=={'FN1':34,'SCARB1':26};U.to_csv(D/'Figure5b_60_gene_unit_estimates.tsv',sep='\t',index=False)
ids=U.pivot_table(index=['GSE','biological_unit_id'],columns='gene',values='beta',aggfunc='first').reset_index();ids.to_csv(D/'Figure5b_34_recorded_unit_IDs.tsv',sep='\t',index=False,na_rep='NOT_ESTIMABLE_OR_NOT_PROVIDED')
im=Image.open(P['image']);W,H=im.size;assert len(z)==1779;assert z.in_tissue.eq(1).all();assert z.hires_x.between(0,W).all() and z.hires_y.between(0,H).all()
mpl.rcParams.update({'font.family':'Arial','font.size':6,'text.color':'#222222','axes.labelcolor':'#222222','axes.linewidth':.45,'pdf.fonttype':42,'svg.fonttype':'none','xtick.major.width':.4,'ytick.major.width':.4,'xtick.major.size':2,'ytick.major.size':0})
cm=mpl.colors.LinearSegmentedColormap.from_list('retained_spatial',['#482878','#3B528B','#2196C4','#35B779','#FDE725','#F68B32','#B2182B']);ZERO='#8B8F95';COL={'FN1':'#DB8C27','SCARB1':'#278F90'}
def ax(f,x,y,w,h):return f.add_axes([x/183,1-(y+h)/170,w/183,h/170])
def tx(f,x,y,t,**kw):return f.text(x/183,1-y/170,t,va='top',**kw)
f=plt.figure(figsize=(183/25.4,170/25.4));tx(f,7,3,'GSE248608 · GSM7917894_CA1',fontsize=6.2)
for x,key,label in [(5,'HE','H&E'),(37,'Fib_abundance','Inferred Fib\nabundance'),(69,'FN1','FN1')]:
 a=ax(f,x,16,28,29.8);a.set_xlim(0,W);a.set_ylim(H,0);a.set_aspect('equal');a.set_xticks([]);a.set_yticks([]);a.spines[:].set_visible(False)
 tx(f,x,10,label,fontsize=6)
 if key=='HE':a.imshow(im,origin='upper',extent=[0,W,H,0],interpolation='nearest')
 else:
  a.set_facecolor('#050708');v=z['abundance__Fibroblast' if key=='Fib_abundance' else 'expr__FN1'].to_numpy();finite=np.isfinite(v);zero=finite&(v==0);positive=finite&~zero;lo,hi=P['limits'][key];norm=mpl.colors.Normalize(lo,hi,clip=True);palette='viridis' if key=='Fib_abundance' else cm
  a.scatter(z.hires_x[zero],z.hires_y[zero],s=.6,c=ZERO,lw=0);a.scatter(z.hires_x[positive],z.hires_y[positive],s=.6,c=v[positive],cmap=palette,norm=norm,lw=0)
  if (~finite).any():a.scatter(z.hires_x[~finite],z.hires_y[~finite],s=1.8,marker='x',c='white',lw=.4)
  ca=ax(f,x+1,47,26,1.2);cb=f.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=palette),cax=ca,orientation='horizontal');cb.set_ticks([lo,hi]);cb.set_ticklabels([f'{lo:.2g}',f'{hi:.2g}'] if key=='Fib_abundance' else [f'≤{lo:.2g}',f'≥{hi:.2g}']);cb.ax.tick_params(labelsize=6,pad=1)
  tx(f,x+14,53,'Sum of marginal q05' if key=='Fib_abundance' else 'log1p(CP10k)',ha='center',fontsize=6)
 f.canvas.draw();bb=a.get_position();records.append(dict(key=key,field=[0,W,H,0],allocated_mm=[x,16,28,29.8],data_mm=[bb.x0*183,(1-bb.y1)*170,bb.width*183,bb.height*170],spots=1779 if key!='HE' else None))
# Shared dedicated legend, separated from FN1 label and the color bars.
f.legend(handles=[Line2D([],[],marker='o',ls='',mfc=ZERO,mec='none',ms=3,label='Expression zero')],loc='upper left',bbox_to_anchor=(7/183,1-56/170),frameon=False,borderaxespad=0,handletextpad=.4,fontsize=6)
a=ax(f,118,10,38,38);studies=sorted(U.GSE.unique());points=[]
for gi,g in enumerate(['FN1','SCARB1']):
 for i,study in enumerate(studies):
  q=U[U.gene.eq(g)&U.GSE.eq(study)].sort_values('biological_unit_id')
  if q.empty:continue
  offs=np.linspace(-.15,.15,len(q)) if len(q)>1 else np.array([0.]);yy=i+(-.22 if gi==0 else .22)+offs
  a.scatter(q.beta,yy,s=4.8,c=COL[g],marker='o' if gi==0 else '^',linewidths=0,zorder=3)
  for (_,r),y in zip(q.iterrows(),yy):points.append({'gene':g,'GSE':study,'biological_unit_id':r.biological_unit_id,'beta':r.beta,'display_y':float(y)})
 q=M.loc[g];yy=7.3+gi*.75;a.errorbar(q.beta,yy,xerr=[[q.beta-q.CI_low],[q.CI_high-q.beta]],fmt='D',ms=3,color=COL[g],lw=.8,capsize=0)
 # Fields share exactly the plotted pooled row centers in page coordinates.
 py=10+38*(yy+.65)/9.25
 fdr=q.FDR_BH_all_29_pairs_x_3_terms;label=f'{fdr:.2e}' if fdr<.001 else f'{fdr:.4f}'
 f.text(163.6/183,1-py/170,label,ha='center',va='center',fontsize=6)
 f.text(177.6/183,1-py/170,f'{q.I2:.1f}%',ha='center',va='center',fontsize=6)
a.axvline(0,lw=.55,c='#9BA1A4',zorder=0);a.axhline(6.65,lw=.4,c='#CFD3D5');a.set_ylim(8.6,-.65);a.set_xlim(-1,.55);a.set_yticks(list(range(7))+[7.3,8.05]);a.set_yticklabels(studies+['FN1','SCARB1']);a.tick_params(axis='y',pad=3,labelsize=6);a.set_xticks([-1,-.5,0,.5]);a.set_xlabel('Conditional Fib-abundance coefficient',labelpad=2,fontsize=6);a.spines[['top','right','left']].set_visible(False)
tx(f,163.6,38.9,'FDR',ha='center',fontsize=6);tx(f,177.6,38.9,'I²',ha='center',fontsize=6)
f.legend(handles=[Line2D([],[],marker='o',ls='',mfc=COL['FN1'],mec='none',ms=3,label='FN1 (34 estimates)'),Line2D([],[],marker='^',ls='',mfc=COL['SCARB1'],mec='none',ms=3.5,label='SCARB1 (26 estimates)')],loc='upper left',bbox_to_anchor=(108/183,1-2.5/170),ncol=2,frameon=False,borderaxespad=0,handletextpad=.3,columnspacing=.7,fontsize=6)
# Symbol legend for pooled CI (not unit SE).
leg=ax(f,108,56,6,2);leg.plot([.15,.85],[.5,.5],c='#555',lw=.8);leg.plot(.5,.5,'D',ms=3,c='#555');leg.set(xlim=(0,1),ylim=(0,1));leg.axis('off');tx(f,116,56,'Pooled estimate and 95% CI',fontsize=6)
f.savefig(O/'review/new_ab_components.svg');plt.close(f)
NS='http://www.w3.org/2000/svg';N={'s':NS};E.register_namespace('',NS);E.register_namespace('xlink','http://www.w3.org/1999/xlink');tree=E.parse(B/'figures/Figure5.svg');root=tree.getroot();fg=root.find("s:g[@id='figure_1']",N);old=fg.find("s:g[@id='Figure5_evidence_focus_ab']",N);assert old is not None;fg.remove(old)
new=E.parse(O/'review/new_ab_components.svg').getroot();ids={el.get('id'):'finish_'+el.get('id') for el in new.iter() if el.get('id')}
for el in new.iter():
 for k,v in list(el.attrib.items()):
  if k=='id':el.set(k,ids[v]);continue
  for old,nn in ids.items():v=v.replace('url(#'+old+')','url(#'+nn+')')
  if v.startswith('#') and v[1:] in ids:v='#'+ids[v[1:]]
  el.set(k,v)
wrap=E.SubElement(fg,'{'+NS+'}g',{'id':'Figure5_ab_finalization'})
for el in new.find("s:g[@id='finish_figure_1']",N):
 if el.get('id')!='finish_patch_1':wrap.append(copy.deepcopy(el))
for defs in new.findall('s:defs',N):root.append(copy.deepcopy(defs))
svg=O/'figures/Figure5.svg';tree.write(svg,encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(svg),write_to=str(O/'figures/Figure5.pdf'))
orig=E.parse(B/'figures/Figure5.svg').getroot();cf={}
for gid in ['FigureYa6_native_forest','axes_8','axes_9','axes_10','legend_2','legend_3','legend_4']+[f'text_{i}' for i in range(136,157)]:
 cf[gid]=E.tostring(orig.find(f".//s:g[@id='{gid}']",N))==E.tostring(root.find(f".//s:g[@id='{gid}']",N));assert cf[gid]
# Rendering of untouched lower panels must also remain identical.
before=fitz.open(B/'figures/Figure5.pdf');after=fitz.open(O/'figures/Figure5.pdf');r=fitz.Rect(0,62*72/25.4,before[0].rect.width,before[0].rect.height);same=before[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=r).samples==after[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=r).samples;assert same
for label,d in [('before',before),('after',after)]:d[0].get_pixmap(matrix=fitz.Matrix(2.4,2.4)).save(O/'review'/f'Figure5_{label}.png')
pd.DataFrame(points).to_csv(D/'Figure5b_point_mapping.tsv',sep='\t',index=False)
pd.DataFrame(source_map).to_csv(O/'Source_Data_index.tsv',sep='\t',index=False)
inputs={str(p):sha(p) for p in [B/'figures/Figure5.pdf',B/'figures/Figure5.svg',B/'figures/Figure5_focus_prototype.pdf',B/'source_data/all29_original_meta.tsv',B/'source_data/all29_original_unit_models.tsv',Path(P['image'])]}
(O/'notes/checks.json').write_text(json.dumps({'gene_unit_estimates':60,'unique_recorded_unit_IDs':34,'gene_counts':{'FN1':34,'SCARB1':26},'overlapping_IDs':len(set(U[U.gene.eq('FN1')].biological_unit_id)&set(U[U.gene.eq('SCARB1')].biological_unit_id)),'c_f_SVG_identical':cf,'c_f_render_identical':same,'companion_PDF_byte_identical':True,'spatial_views':records,'same_field_and_physical_scale':np.allclose([r['data_mm'][2:] for r in records],records[0]['data_mm'][2:]),'source_scales_unchanged':P['limits'],'science_recomputed':False,'visual_review':'PENDING_RENDER_INSPECTION','status':'PAUSED_FOR_FIGURE5_ADOPTION_REVIEW','ORA':'PENDING','donor_identity_pending':True,'source_hashes':inputs},indent=2))
print('60 gene-unit estimates; 34 distinct IDs; c-f and detailed prototype unchanged.')
