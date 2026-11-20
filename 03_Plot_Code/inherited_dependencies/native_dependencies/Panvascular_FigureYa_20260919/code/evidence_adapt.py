"""FigureYa plotting-grammar port to frozen, data-bearing SVG coordinates.
No statistics, ranks, clustering or selection are run. See notes/evidence_methods.md.
"""
from pathlib import Path
import xml.etree.ElementTree as E
import re,io,base64,copy,json,hashlib,os,sys
import numpy as np,pandas as pd
from PIL import Image
import matplotlib as mpl
import fitz
O=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');B=Path('/Users/taki/Desktop/投稿包_最新版_20260919');S=B/'Source_Data_Figures';A=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919')
NS='http://www.w3.org/2000/svg';XL='http://www.w3.org/1999/xlink';E.register_namespace('',NS);E.register_namespace('xlink',XL)
for d in ['main','extended_data','editable_svg','notes','review']:(O/d).mkdir(exist_ok=True)
sys.path.insert(0,str(A/'Panvascular_Local_Touchup_20260918/tools'))
import cairosvg,cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def tag(s):return '{'+NS+'}'+s
def find(r,i):return next(e for e in r.iter() if e.get('id')==i)
def axis(r,i):return find(r,'axes_'+str(i))
def numbers(s):return [float(x) for x in re.findall(r'-?\d+(?:\.\d+)?(?:e[-+]?\d+)?',s)]
def box(a):
 p=next(e for e in list(a)[0] if e.tag==tag('path'));v=np.array(numbers(p.get('d'))).reshape(-1,2);return v[:,0].min(),v[:,1].min(),v[:,0].max(),v[:,1].max()
def source(rel):
 p=S/rel;return pd.read_csv(p,sep='\t',float_precision='round_trip')
C=source('Figure3/Fig3d_e_all29.tsv');T=source('Figure4/14_final_module_evidence_and_tiers.tsv');TOP=source('Figure4/05_top_nodes_reference.tsv')
assert len(C)==29 and T.loc[T.module_tier.eq('high_confidence_pan_vascular'),'module'].tolist()==['CM09']
SP=source('Figure5/target_meta_descriptive.tsv');J=C[['gene_symbol','cell_state']].merge(SP.query("posterior=='q05' and sensitivity=='all_sections' and term=='cell_abundance'"),left_on=['gene_symbol','cell_state'],right_on=['gene','cell_context'],validate='one_to_one',how='left')
M=np.zeros((29,5),int)
for i,r in C.iterrows():
 M[i]=[3 if bool(r.primary_Bonferroni) else 1,3 if r.target_set_meta_FDR<.05 else 1,3 if r['FDR_BH__leave_one_CAD_r0.1']<.05 else 1,3 if r.primary_recomputed_n_strong>0 else 1 if r.primary_recomputed_n_contexts_tested>0 else 0,3 if J.iloc[i].FDR_BH_all_29_pairs_x_3_terms<.05 and J.iloc[i].beta>0 else 2 if J.iloc[i].FDR_BH_all_29_pairs_x_3_terms<.05 else 1]
EDP=A/'Panvascular_Figure_Redesign_20260917_v3/panels/ED10e_source.tsv';ED=pd.read_csv(EDP,sep='\t');EV=-np.log10(ED.iloc[:,1:].to_numpy(float).clip(1e-10));assert ED.iloc[:,0].tolist()==[g+' | '+s for g,s in zip(C.gene_symbol,C.cell_state)]
checks=[];coords=[];changes=[]
def diamonds(r,ai):
 """Port fpDrawDiamondCI: same CI endpoints and point centres, fixed equal area."""
 a=axis(r,ai);defs={e.get('id'):e for e in r.iter() if e.tag==tag('path') and e.get('id')};refs={u.get('{'+XL+'}href')[1:] for u in a.iter() if u.tag==tag('use')};count=0
 for ref in refs:
  p=defs.get(ref)
  if p is None or 'C' not in p.get('d',''):continue
  d=p.get('d');rad=abs(numbers(d)[1]);rr=rad*np.sqrt(np.pi/2)
  p.set('d',f'M 0 {rr} L {rr} 0 L 0 {-rr} L {-rr} 0 z');count+=sum(u.tag==tag('use') and u.get('{'+XL+'}href')=='#'+ref for u in a.iter())
 return count

def lollipops(r,ai):
 """FigureYa11 geom_point(x,value) plus a zero-baseline guide, fixed rows."""
 a=axis(r,ai);out=[]
 for g in list(a):
  if not g.get('id','').startswith('patch_'):continue
  pp=list(g)
  if len(pp)!=1 or pp[0].tag!=tag('path'):continue
  p=pp[0];st=p.get('style','');fill=re.search(r'fill: (#[0-9a-fA-F]{6})',st)
  if not fill or fill.group(1).lower() not in ['#657c92','#8d7a9f','#398f8b']:continue
  v=np.array(numbers(p.get('d'))).reshape(-1,2);x0,x1=v[:,0].min(),v[:,0].max();yy=v[:,1].mean();color=fill.group(1)
  g.remove(p);E.SubElement(g,tag('path'),{'d':f'M {x0} {yy} L {x1} {yy}','style':f'fill:none;stroke:{color};stroke-width:1.1','clip-path':p.get('clip-path')})
  E.SubElement(g,tag('circle'),{'cx':str(x1),'cy':str(yy),'r':'2.0','fill':color,'stroke':'white','stroke-width':'.35'})
  out.append((x0,x1,yy))
 return out

def inset_mesh(r,ai,factor=.94):
 """Port FigureYa83 geom_tile(width=.9,height=.9); conservative .94 for readability."""
 a=axis(r,ai);n=0
 for g in a:
  if not g.get('id','').startswith('QuadMesh_'):continue
  for p in g:
   if p.tag!=tag('path'):continue
   v=np.array(numbers(p.get('d'))).reshape(-1,2);x0,x1=v[:,0].min(),v[:,0].max();y0,y1=v[:,1].min(),v[:,1].max();dx=(x1-x0)*(1-factor)/2;dy=(y1-y0)*(1-factor)/2
   p.set('d',f'M {x0+dx} {y0+dy} L {x1-dx} {y0+dy} L {x1-dx} {y1-dy} L {x0+dx} {y1-dy} Z');n+=1
 return n

def vector_matrix(r,ai,values,cmap,vmin,vmax,name):
 """Numerical source defines colour. Pixel source is independent identity QA, not inferred data."""
 a=axis(r,ai);parent=next(e for e in a.iter() if any(ch.tag==tag('image') for ch in e));im=next(e for e in parent if e.tag==tag('image'));png=np.asarray(Image.open(io.BytesIO(base64.b64decode(im.get('{'+XL+'}href').split(',')[1]))).convert('RGBA'))
 rgb=(cmap(mpl.colors.Normalize(vmin,vmax)(values))*255).astype(np.uint8);assert png.shape==rgb.shape
 # SVG native imshow PNG uses byte truncation, including unmodified colormap.
 assert np.array_equal(png,rgb),(name,int(np.abs(png.astype(int)-rgb.astype(int)).max()))
 sx,_,_,sy,x,y=numbers(im.get('transform'));assert sx>0 and sy>0
 parent.remove(im);g=E.SubElement(parent,tag('g'),{'id':'FigureYa_vector_matrix_'+name})
 for i in range(values.shape[0]):
  for j in range(values.shape[1]):
   color='#'+''.join(f'{v:02x}' for v in rgb[i,j,:3]);E.SubElement(g,tag('rect'),{'x':str(x+j*sx+.3),'y':str(y+i*sy+.2),'width':str(sx-.6),'height':str(sy-.4),'fill':color,'data-row':str(i+1),'data-col':str(j+1),'data-source-value':str(values[i,j])})
 checks.append({'matrix':name,'shape':list(values.shape),'all_source_colours_equal_frozen_pixels':True,'cell_count':values.size})
 return x,y,sx,sy

for name in ['Figure3','Figure4','ExtendedData_Figure6','ExtendedData_Figure9','ExtendedData_Figure10']:
 inp=B/'图件编辑资料/editable_svg'/(name+'.svg');r=E.parse(inp).getroot();orig=copy.deepcopy(r);before_text=[e.text for e in r.iter() if e.tag==tag('text')];change={}
 if name=='Figure3':
  change['CI_diamond_points']=diamonds(r,5);x,y,sx,sy=vector_matrix(r,6,M,mpl.colors.ListedColormap(['white','#DCE1E5','#D18A32','#278D87']),0,3,'Figure3e')
  bx=box(axis(r,5));
  for i,g in enumerate(C.gene_symbol):coords.append(dict(figure=name,row=i+1,gene=g,forest_center_pt=bx[1]+(bx[3]-bx[1])*(i+.5)/29,matrix_center_pt=y+sy*(i+.5)))
  # Precisely retain the two frozen seven-pathway point panels and all legends.
  assert E.tostring(axis(r,3))==E.tostring(axis(orig,3)) and E.tostring(axis(r,4))==E.tostring(axis(orig,4))
 if name=='Figure4':
  change['module_weight_lollipops']=len(lollipops(r,6))+len(lollipops(r,7));change['similarity_tiles']=inset_mesh(r,8)
  for i in range(1,5):assert E.tostring(axis(r,i))==E.tostring(axis(orig,i))
 if name=='ExtendedData_Figure6':
  diamonds(r,1);diamonds(r,2);change['CI_diamond_points']=16
  # Row banding supports long L2 labels and comparison; inferred only from fixed rendered centres.
  for ai in [1,2]:
   a=axis(r,ai);x0,y0,x1,y1=box(a);g=E.Element(tag('g'),{'id':f'FigureYa_fixed_row_guides_{ai}'})
   for i in range(8):
    if i%2==0:E.SubElement(g,tag('rect'),{'x':str(x0),'y':str(y0+(y1-y0)*(i+.20)/8.3),'width':str(x1-x0),'height':str((y1-y0)*.9/8.3),'fill':'#f5f6f7'})
   a.insert(1,g)
  change['row_guides']=8
 if name=='ExtendedData_Figure9':change['module_weight_lollipops']=len(lollipops(r,1));change['similarity_tiles']=inset_mesh(r,4)
 if name=='ExtendedData_Figure10':
  change['CI_diamond_points']=diamonds(r,5)
  # Match All sections legend to the diamond; QC-fallback retains its square.
  for lg in [e for e in r.iter() if e.get('id','').startswith('legend_')]:
   for u in lg.iter():
    if u.tag==tag('use') and '#356b8c' in u.get('style',''):
     ref=u.get('{'+XL+'}href')[1:];p=find(r,ref)
     if 'C' in p.get('d',''):
      rr=abs(numbers(p.get('d'))[1])*np.sqrt(np.pi/2);p.set('d',f'M 0 {rr} L {rr} 0 L 0 {-rr} L {-rr} 0 z')
  x,y,sx,sy=vector_matrix(r,6,EV,mpl.colormaps['Blues'],0,10,'ED10e');bx=box(axis(r,5))
  for i,g in enumerate(C.gene_symbol):coords.append(dict(figure=name,row=i+1,gene=g,forest_center_pt=bx[1]+(bx[3]-bx[1])*(i+.5)/29,matrix_center_pt=y+sy*(i+.5)))
 # The referenced marker definitions include matching figure-level legends, updated together.
 for e in r.iter():
  if e.tag==tag('text'):
   st=e.get('style','');e.set('style',re.sub(r'\b5\.5px\b','6px',st))
 assert before_text==[e.text for e in r.iter() if e.tag==tag('text')]
 out=O/'editable_svg'/(name+'.svg');E.ElementTree(r).write(out,encoding='utf-8',xml_declaration=True)
 pdf=O/('main' if name.startswith('Figure') else 'extended_data')/(name+'.pdf');cairosvg.svg2pdf(url=str(out),write_to=str(pdf));d=fitz.open(pdf);pg=d[0];pg.get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review'/('evidence_'+name+'_after.png'))
 bad=[];sizes=[]
 for block in pg.get_text('dict')['blocks']:
  for line in block.get('lines',[]):
   for span in line['spans']:
    sizes.append(span['size']);rect=fitz.Rect(span['bbox'])
    if rect.x0<-.5 or rect.y0<-.5 or rect.x1>pg.rect.width+.5 or rect.y1>pg.rect.height+.5:bad.append(span['text'])
 change.update(figure=name,source_svg_sha256=sha(inp),pdf_sha256=sha(pdf),svg_sha256=sha(out),width_mm=round(pg.rect.width*25.4/72,3),height_mm=round(pg.rect.height*25.4/72,3),all_text_content_unchanged=True,outside_page_text=bad,minimum_font_pt=min(sizes));changes.append(change)
rows=pd.DataFrame(coords);rows['delta_pt']=rows.forest_center_pt-rows.matrix_center_pt;assert rows.delta_pt.abs().max()<.00002
rows.to_csv(O/'review/evidence_shared29_absolute_geometry.tsv',sep='\t',index=False)
(O/'review/evidence_checks.json').write_text(json.dumps({'figures':changes,'matrix_checks':checks,'max_shared29_displacement_pt':rows.delta_pt.abs().max(),'no_statistical_computation':True},indent=2))
print(json.dumps(changes,indent=2))
