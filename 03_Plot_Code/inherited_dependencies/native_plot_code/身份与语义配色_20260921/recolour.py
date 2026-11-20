from pathlib import Path
import os,sys,copy,re,json,io,base64,hashlib
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
import numpy as np,pandas as pd,xml.etree.ElementTree as E
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from PIL import Image
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import cairosvg,cairosvg.image,fitz
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
W=Path(__file__).resolve().parents[1];O=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920');P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');S=P/'Source_Data_Figures'
C=json.loads((W/'semantic_colors.json').read_text());L1=C['L1'];mm=72/25.4
NS='http://www.w3.org/2000/svg';XL='{http://www.w3.org/1999/xlink}href';q=lambda t:'{'+NS+'}'+t
E.register_namespace('',NS);E.register_namespace('xlink','http://www.w3.org/1999/xlink')
read=lambda p:pd.read_csv(p,sep='\t',float_precision='round_trip');find=lambda r,i:next(e for e in r.iter() if e.get('id')==i)
oldL1=read(S/'Figure1/L1_palette_bound.tsv').set_index('L1').hex.to_dict();identity={v.lower():L1[k] for k,v in oldL1.items()}
log=[];checks=[]
plt.rcParams.update({'font.family':'Arial','font.size':6,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.45,'text.color':C['style']['text'],'axes.labelcolor':C['style']['text'],'xtick.labelcolor':C['style']['text'],'ytick.labelcolor':C['style']['text']})
BL=lambda t:plt.colormaps[C['continuous']['marker']['cmap']](C['continuous']['marker']['lookup_range'][0]+np.diff(C['continuous']['marker']['lookup_range'])[0]*np.asarray(t))
MA=lambda t:plt.colormaps[C['continuous']['FN1']['cmap']](C['continuous']['FN1']['lookup_range'][0]+np.diff(C['continuous']['FN1']['lookup_range'])[0]*np.asarray(t))
def scope(r,id,mp,semantic):
 g=find(r,id);n=0
 for e in g.iter():
  if e.tag==q('text'):continue
  for k in ['fill','stroke','style']:
   v=e.get(k)
   if v:
    new=re.sub(r'#[0-9a-fA-F]{6}',lambda m:mp.get(m[0].lower(),m[0]),v)
    if new!=v:e.set(k,new);n+=1
 log.append({'figure':NAME,'object':id,'semantic':semantic,'changed_attributes':n})
def setfill(e,c):
 s=e.get('style','')
 if 'fill:' in s:s=re.sub(r'fill:\s*[^;]+','fill: '+c,s);e.set('style',s)
 else:e.set('fill',c)
def gradient(r,id,cmap):
 for im in find(r,id).iter(q('image')):
  wi,he=float(im.get('width')),float(im.get('height'));arr=np.linspace(0,1,512)[None,:] if wi>he else (np.linspace(0,1,512) if 'scale(1 -1)' in im.get('transform','') else np.linspace(1,0,512))[:,None]
  a=np.uint8(np.round(cmap(arr)*255));b=io.BytesIO();Image.fromarray(a).save(b,format='PNG');im.set(XL,'data:image/png;base64,'+base64.b64encode(b.getvalue()).decode());im.set('preserveAspectRatio','none')
 log.append({'figure':NAME,'object':id,'semantic':'continuous colourbar lookup only','changed_attributes':1})
def textstyle(r):
 # Axis/tick/guide objects are identified structurally, not by a page-wide colour search.
 for g in r.iter(q('g')):
  gid=g.get('id','')
  if 'matplotlib.axis' in gid:
   for e in g.iter():
    if e.tag==q('text'):continue
    for k in ['style','stroke','fill']:
     if e.get(k):e.set(k,re.sub(r'#[0-9a-fA-F]{6}',lambda m: {'#000000':C['style']['axis'],'#e9ecef':C['style']['guide_line'],'#e5e8eb':C['style']['guide_line']}.get(m[0].lower(),m[0]),e.get(k)))
  elif 'patch_' in gid:
   for e in g:
    st=e.get('style','')
    if 'fill: none' in st:e.set('style',re.sub(r'stroke: #(000000|222222)', 'stroke: '+C['style']['axis'],st))
 for e in r.iter(q('text')):
  s=e.get('style','');explicit=re.search(r'fill:\s*(#[0-9a-fA-F]{6}|white)',s);v=(explicit[1].lower() if explicit else e.get('fill','').lower())
  if v in ['#ffffff','white']:continue # numeric contrast on dark cells
  if 'fill:' in s:s=re.sub(r'fill:\s*[^;]+','fill: '+C['style']['text'],s)
  else:s+=';fill:'+C['style']['text']
  e.set('style',s);e.attrib.pop('fill',None)
def embed(root,fig,prefix):
 b=io.BytesIO();fig.savefig(b,format='svg',transparent=True);plt.close(fig);r=E.fromstring(b.getvalue());ids={e.get('id'):prefix+e.get('id') for e in r.iter() if e.get('id')}
 for e in r.iter():
  for k,v in list(e.attrib.items()):
   if k=='id':e.set(k,ids[v])
   else:
    for old,new in ids.items():
     v=v.replace('url(#'+old+')','url(#'+new+')');v='#'+new if v=='#'+old else v
    e.set(k,v)
 for e in r:
  if e.tag in [q('g'),q('defs')]:root.append(e)
def remove(root,id):
 for p in root.iter():
  for e in list(p):
   if e.get('id')==id:p.remove(e)
def recolor_points(g,vals,norm,cmap,oldcmap,semantic):
 els=[e for e in g if e.tag==q('path')]
 assert len(els)==len(vals),(NAME,g.get('id'),len(els),len(vals))
 for el,v in zip(els,vals):
  orig=re.search(r'fill:\s*(#[0-9a-fA-F]{6})',el.get('style',''))[1].lower();expected=colors.to_hex(oldcmap(norm(v)))
  assert orig==expected,(NAME,g.get('id'),v,orig,expected)
  setfill(el,colors.to_hex(cmap(norm(v))))
 checks.append({'figure':NAME,'object':g.get('id'),'n':len(vals),'original_values_and_geometry_unchanged':True,'encoding':semantic})
def export(r,name):
 # Bounded approved visibility finishing rules
 import importlib.util
 vp=P/'图件编辑资料/点线可见性收尾_20260921/code/visibility.py'
 spec=importlib.util.spec_from_file_location('approved_visibility',vp);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.apply(r,name,C)
 textstyle(r);p=W/'editable_svg'/f'{name}.svg';E.ElementTree(r).write(p,encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(p),write_to=str(W/'figures'/f'{name}.pdf'));d=fitz.open(W/'figures'/f'{name}.pdf');d[0].get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(W/'review'/f'{name}.png');print('exported',name,flush=True)

# Figure1: existing g circle geometry is preserved exactly; rasterized L1 map is regenerated from frozen coordinates.
NAME='Figure1';r=E.parse(W/'backup'/f'{NAME}.svg').getroot()
Z=read(O/'marker_audit_20260921/tables/per_gene_scaled_mean.tsv').set_index('cell_type') if False else pd.read_csv(O/'marker_audit_20260921/tables/per_gene_scaled_mean.tsv',sep='\t',index_col=0,float_precision='round_trip')
F=pd.read_csv(O/'marker_audit_20260921/tables/original_detection_fraction.tsv',sep='\t',index_col=0,float_precision='round_trip');good=np.isfinite(Z.to_numpy())&np.isfinite(F.to_numpy())&(F.to_numpy()>0)
recolor_points(find(r,'adopted_B_g_PathCollection_1'),Z.to_numpy()[good],colors.Normalize(0,1),BL,plt.colormaps['viridis'],'saved per-gene minmax; original circle area')
gradient(r,'adopted_B_g_axes_2',BL)
for j,st in enumerate(Z.index):scope(r,'adopted_B_g_patch_'+str(3+j),{oldL1[st].lower():L1[st]},'L1 '+st)
for id in ['axes_5','axes_6','axes_7']:scope(r,id,identity,'L1 identity')
D=read(S/'Figure1/full_cells_frozen_umap.tsv.gz');shown=D[D.level1_annotation.ne('Unknown')]
x0,x1=D.UMAP1.min(),D.UMAP1.max();y0,y1=D.UMAP2.min(),D.UMAP2.max();dx=x1-x0;dy=y1-y0
f=plt.figure(figsize=(183/25.4,170/25.4));a=f.add_axes([34/183,1-60/170,49/183,33/170]);a.set(xlim=(x0-.04*dx,x1+.04*dx),ylim=(y0-.04*dy,y1+.04*dy));a.set_aspect('equal');a.axis('off');a.scatter(shown.UMAP1,shown.UMAP2,c=shown.level1_annotation.map(L1),s=.055,lw=0,alpha=.8,rasterized=True)
b=io.BytesIO();f.savefig(b,format='svg');plt.close(f);nr=E.fromstring(b.getvalue());new=next(nr.iter(q('image')));old=next(find(r,'axes_2').iter(q('image')))
for k in ['x','y','width','height','transform']:assert old.get(k)==new.get(k),(k,old.get(k),new.get(k))
old.set(XL,new.get(XL));checks.append({'figure':NAME,'object':'axes_2','n':len(shown),'frozen_coordinates_rerendered':True,'same_image_geometry':True});del D,shown
export(r,NAME)

# Figure3: exact native records, same axes, same significance sizes/fill. Only label offsets are new.
NAME='Figure3';r=E.parse(W/'backup'/f'{NAME}.svg').getroot();remove(r,'aligned_a_axes_1')
f=plt.figure(figsize=(183/25.4,170/25.4));a=f.add_axes([15/183,1-36/170,63/183,27/170]);d=read(S/'Figure3/Fig3a_L1.tsv')
for _,v in d.iterrows():
 sig=v.FDR_BH<.05 and v['primary_r0.1__assoc_mcp_BH']<.05
 a.scatter(v.BETA_STD,v['primary_r0.1__assoc_mcz'],s=25+18*sig,facecolors=L1[v.group] if sig else 'white',edgecolors=L1[v.group],lw=.8,zorder=3)
 if sig or v.group=='Smooth muscle cell':
  lab={'Endothelial cell':'EC','Smooth muscle cell':'Smooth muscle'}.get(v.group,v.group)
  off={'Pericyte':(-5,1),'Fibromyocyte':(-6,4),'Smooth muscle':(-4,6),'EC':(-6,-8)}[lab];ha='right'
  a.annotate(lab,(v.BETA_STD,v['primary_r0.1__assoc_mcz']),xytext=off,textcoords='offset points',ha=ha,fontsize=6,arrowprops=dict(arrowstyle='-',lw=.35,color=C['style']['secondary_text']),zorder=4)
a.axhline(0,c=C['style']['guide_line'],lw=.55);a.axvline(0,c=C['style']['guide_line'],lw=.55);a.set_xlim(-.0215,.040);a.set_ylim(-3.3,6.85);a.set_xticks([-.02,0,.02,.04]);a.set_yticks([-2,0,2,4,6]);a.tick_params(labelsize=6);a.set_xlabel('MAGMA standardised beta',labelpad=2,fontsize=6);a.set_ylabel('scDRS association Z',labelpad=2,fontsize=6)
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([],[],ls='',marker='o',ms=3.5,mfc='#888',mec='#888',label='Both FDR < 0.05')],loc='upper left',bbox_to_anchor=(-.03,-.33),frameon=False,fontsize=6,borderaxespad=0,handletextpad=.4)
embed(r,f,'identity_a_')
for id in ['figure3_b_content','figure3d_compact_forest']:scope(r,id,identity,'parent L1 identity, not significance')
scope(r,'c_magma_axes_1',{'#356b8c':C['method']['MAGMA']},'MAGMA positive NES single blue')
EVID={'#278d87':C['evidence']['support'],'#d18a32':C['evidence']['inverse'],'#dce1e5':C['evidence']['no_support']}
for id in ['FigureYa_vector_matrix_Figure3e','figure3e_three_symbol_key']:scope(r,id,EVID,'evidence category')
scope(r,'figure3f_three_evidence_columns',{'#356b8c':C['method']['primary'],'#555555':C['method']['leave_CAD']},'scenario identity; original shape/fill retained')
export(r,NAME)

# Figure5: expression colours assigned from saved expression values in exact original spot order.
NAME='Figure5';r=E.parse(W/'backup'/f'{NAME}.svg').getroot();spec=json.loads((O/'notes/Figure5_checks.json').read_text());oldcm=colors.LinearSegmentedColormap.from_list('old',['#482878','#3B528B','#2196C4','#35B779','#FDE725','#F68B32','#B2182B'])
paths=[O/'reproducibility_dependencies/Panvascular_Figure5_AB_Finalization_20260920/Source_Data/Figure5/CA1_existing_spot_values.tsv',O/'source_data/Figure5/T1077_existing_spot_values.tsv']
for row,(path,axid,cbid,collid,key) in enumerate(zip(paths,['r2_axes_4','r2_axes_9'],['r2_axes_5','r2_axes_10'],['r2_PathCollection_4','r2_PathCollection_8'],['original','added'])):
 z=read(path);v=z['expr__FN1'].to_numpy();vv=v[np.isfinite(v)&(v>0)];lo,hi=spec[key]['limits']['FN1'];norm=colors.Normalize(lo,hi,clip=True);els=list(find(r,collid).iter(q('use')));assert len(els)==len(vv)
 for e,val in zip(els,vv):
  orig=re.search(r'fill:\s*(#[0-9a-fA-F]{6})',e.get('style'))[1];assert orig==colors.to_hex(oldcm(norm(val)))
  setfill(e,colors.to_hex(MA(norm(val))))
 gradient(r,cbid,MA);checks.append({'figure':NAME,'object':axid,'n_spots':len(v),'positive':len(vv),'limits':[lo,hi],'unchanged_zero_coordinates_HE':True})
MOD={'#657c92':C['module']['CM01'],'#8d7a9f':C['module']['CM03'],'#398f8b':C['module']['CM09'],'#929fa9':C['module']['Other'],'#8d999f':C['module']['Other'],'#929da4':C['module']['Other']}
for id in ['FigureYa6_native_forest','axes_8','legend_2','legend_3']:scope(r,id,{**EVID,'#62676b':C['evidence']['no_support']},'signed conditional association evidence')
for id in ['axes_9','axes_10']:scope(r,id,MOD,'module identity')
# CM02 median is originally dark grey, distinct from grey dataset points.
scope(r,'axes_10',{'#555555':C['module']['Other']},'other-module median')
for id in ['r2_axes_11','r2_legend_2','r2_legend_3']:
 if any(e.get('id')==id for e in r.iter()):scope(r,id,{'#db8c27':C['gene']['FN1'],'#278f90':C['gene']['SCARB1']},'gene identity')
# Legend duplicates can be separate named artists; target only their original gene-specific colours inside ab group.
scope(r,'Figure5_two_sections_ab',{'#db8c27':C['gene']['FN1'],'#278f90':C['gene']['SCARB1']},'FN1/SCARB1 gene-identity glyphs')
export(r,NAME)

# ED marker panels: keep every existing path and every saved statistic. No marker selection or expression re-computation.
M=read(O/'reproducibility_dependencies/Panvascular_Marker_Communication_Enrichment_20260920/source_data/display_marker_values.tsv')
for n in [2,3,4,5]:
 NAME=f'ExtendedData_Figure{n}';r=E.parse(W/'backup'/f'{NAME}.svg').getroot()
 for ai,ci in ([(1,3),(4,6)] if n==4 else [(1,3)]):
  ag=find(r,f'newED{n}_axes_{ai}');coll=next(e for e in ag if 'PathCollection' in e.get('id',''))
  xaxis=next(e for e in ag if 'matplotlib.axis' in e.get('id',''));genes=[e.text for e in xaxis.iter(q('text'))]
  ys=[e for e in ag if 'matplotlib.axis' in e.get('id','')][1];labs=[e.text for e in ys.iter(q('text'))]
  dd=M[M.gene_id.isin(genes)&M.display_label.isin(labs)].copy();dct=dd.drop_duplicates(['display_label','gene_id']).set_index(['display_label','gene_id']);assert len(dct)==len(genes)*len(labs),(NAME,len(dct),len(labs),len(genes))
  vals=[]
  for lab in labs:
   for gene in genes:
    t=dct.loc[(lab,gene)];assert t.vmin==0
    if np.isfinite(t.mean_log1p) and t.detection_fraction>0:vals.append(t.mean_log1p)
  caps=dd.vmax.unique();assert len(caps)==1
  recolor_points(coll,np.array(vals),colors.Normalize(0,caps[0]),BL,plt.colormaps['viridis'],'original mean log1p and detection area; range unchanged')
  gradient(r,f'newED{n}_axes_{ci}',BL)
  if n==2:
   # These 16 strips are L1, while ED3-5 strips are L2 and must remain the fixed state palette.
   for j,lab in enumerate(labs):scope(r,f'newED2_patch_{j+3}',{oldL1[lab].lower():L1[lab]},'L1 identity strip '+lab)
 export(r,NAME)

# Remaining pages: use only declared semantic scopes. Unrelated continuous scales and L2 raster maps are retained.
for NAME in ['Figure2','Figure4','ExtendedData_Figure1','ExtendedData_Figure6','ExtendedData_Figure7','ExtendedData_Figure8','ExtendedData_Figure9','ExtendedData_Figure10']:
 r=E.parse(W/'backup'/f'{NAME}.svg').getroot()
 if NAME=='Figure4':
  for id in ['axes_2','axes_5','axes_6','axes_7','axes_10','axes_11','axes_12','Figure4f_module_key_0921']:scope(r,id,MOD,'module identity')
  for id in ['axes_3','axes_5','legend_1']:scope(r,id,EVID,'criterion met / not met')
 if NAME=='ExtendedData_Figure6':
  for id in ['axes_1','axes_2','axes_6','legend_1']:scope(r,id,{'#278d87':C['method']['neutral_comparison']},'neutral estimate/count, not L1 or support identity')
 if NAME=='ExtendedData_Figure8':
  for id in [f'axes_{i}' for i in range(2,10)]:scope(r,id,{'#278d87':C['method']['neutral_comparison']},'scenario comparison points, identity-neutral')
  for id in ['axes_10','axes_11','legend_1']:scope(r,id,{'#354f78':C['method']['MAGMA'],'#278d87':C['method']['scDRS']},'MAGMA/scDRS method')
 if NAME=='ExtendedData_Figure9':
  for id in ['axes_1','axes_2','axes_3','legend_1','legend_2']:scope(r,id,MOD,'module identity')
 if NAME=='ExtendedData_Figure10':
  for id in ['axes_1','axes_2']:scope(r,id,MOD,'module identity')
 export(r,NAME)
(W/'tables/semantic_object_changes.json').write_text(json.dumps(log,indent=2));(W/'tables/colour_data_checks.json').write_text(json.dumps(checks,indent=2));print('ALL RENDERED')
