from pathlib import Path
import os,sys,json,copy,hashlib,re,csv
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import pandas as pd,numpy as np,xml.etree.ElementTree as E
import cairosvg,fitz
O=Path(__file__).resolve().parents[1];B=Path('/Users/taki/Desktop/Panvascular_Adopted_Figure5_Merged_20260920');S=Path('/Users/taki/Desktop/投稿包_最新版_20260919/Source_Data_Figures/Figure1')
rows=list(csv.DictReader((B/'ADOPTED_FIGURE_VERSIONS.tsv').open(),delimiter='\t'));r=rows[0];assert r['figure']=='Figure1';src=Path(r['source_svg']);NS='http://www.w3.org/2000/svg';N={'s':NS};E.register_namespace('',NS);E.register_namespace('xlink','http://www.w3.org/1999/xlink')
tree=E.parse(src);root=tree.getroot();axis=root.find(".//s:g[@id='axes_1']",N);old_images=[hashlib.sha256(E.tostring(e)).hexdigest() for e in root.findall('.//s:image',N)];old_texts=[t.text for t in axis.findall('.//s:text',N)];assert len(old_texts)==6
P=pd.read_csv(S/'L1_palette_bound.tsv',sep='\t').sort_values('order');D=pd.read_csv(S/'full_cells_frozen_umap.tsv.gz',sep='\t',usecols=['_index','UMAP1','UMAP2','level1_annotation']);assert len(D)==1175115 and D._index.is_unique and set(D.level1_annotation)==set(P.L1)
# Recover the original plot's frozen affine transform; only label anchors are selected.
xmin,xmax=D.UMAP1.min(),D.UMAP1.max();ymin,ymax=D.UMAP2.min(),D.UMAP2.max();dx=xmax-xmin;dy=ymax-ymin
width=36*dx/dy;left=4+(63-width)/2;top=43;height=36
# Existing SVG's whole-axis shift remains untouched. New label group uses final page coordinates.
nums=list(map(float,re.findall(r'[-+]?\d*\.?\d+',axis.get('transform'))));sx,sy=nums[:2];mm=72/25.4
anchors=[]
for _,p in P.iterrows():
 g=D[D.level1_annotation.eq(p.L1)];med=g[['UMAP1','UMAP2']].median().to_numpy();dist=((g[['UMAP1','UMAP2']].to_numpy()-med)**2).sum(1);q=g.iloc[int(np.argmin(dist))]
 px=(left+width*(q.UMAP1-(xmin-.045*dx))/(1.09*dx))*mm+sx;py=(top+height*(1-(q.UMAP2-(ymin-.045*dy))/(1.09*dy)))*mm+sy
 anchors.append({'number':int(p.order),'L1':p.L1,'hex':p.hex,'anchor_cell_id':q['_index'],'UMAP1':float(q.UMAP1),'UMAP2':float(q.UMAP2),'anchor_page_x_pt':px,'anchor_page_y_pt':py,'n_cells':len(g)})
A=pd.DataFrame(anchors);rank=A.sort_values('anchor_page_x_pt').index;A['side']='right';A.loc[rank[:8],'side']='left'
for side in ['left','right']:
 idx=A[A.side.eq(side)].sort_values('anchor_page_y_pt').index
 A.loc[idx,'label_page_y_pt']=np.linspace(43.5,81.8,len(idx))*mm;A.loc[idx,'label_page_x_pt']=(5.5 if side=='left' else 57.5)*mm
# Remove the six old annotations, keeping the image object unchanged.
for child in list(axis):
 if not child.tag.endswith('image'):axis.remove(child)
fg=root.find("s:g[@id='figure_1']",N);new=E.SubElement(fg,'{'+NS+'}g',{'id':'Figure1b_all16_annotations'})
for _,a in A.iterrows():
 x,y=a.label_page_x_pt,a.label_page_y_pt;end=x+(4 if a.side=='left' else -4)
 E.SubElement(new,'{'+NS+'}path',{'d':f'M {end:.6f} {y:.6f} L {a.anchor_page_x_pt:.6f} {a.anchor_page_y_pt:.6f}','fill':'none','stroke':'#777777','stroke-width':'0.35'})
 E.SubElement(new,'{'+NS+'}circle',{'cx':f'{a.anchor_page_x_pt:.6f}','cy':f'{a.anchor_page_y_pt:.6f}','r':'.55','fill':'#444444'})
 E.SubElement(new,'{'+NS+'}rect',{'x':f'{x-4.8:.6f}','y':f'{y-4:.6f}','width':'9.6','height':'8','rx':'2','fill':'white','stroke':'none'})
 t=E.SubElement(new,'{'+NS+'}text',{'x':f'{x:.6f}','y':f'{y+2.2:.6f}','text-anchor':'middle','style':"fill:#222222;font:6.5px 'Arial'"});t.text=str(int(a.number))
# Complete names in the existing legend; neither order nor category colors change.
short={p:p.replace(' cell','') for p in P.L1};short['Smooth muscle cell']='Smooth muscle';short['Natural killer cell']='Natural killer';changed=[]
for _,p in P.iterrows():
 found=[]
 for t in root.findall('.//s:text',N):
  if t.text==short[p.L1] and float(t.get('x','-1'))>190 and 115<float(t.get('y','-1'))<237:found.append(t)
 assert len(found)==1,(p.L1,len(found));t=found[0];t.text=f'{int(p.order)}  {p.L1}';changed.append(p.L1)
svg=O/'figures/Figure1.svg';tree.write(svg,encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(svg),write_to=str(O/'figures/Figure1.pdf'))
assert old_images==[hashlib.sha256(E.tostring(e)).hexdigest() for e in root.findall('.//s:image',N)]
A.sort_values('number').to_csv(O/'notes/Figure1b_16_label_mapping.tsv',sep='\t',index=False)
for side,path in [('before',Path(r['source_pdf'])),('after',O/'figures/Figure1.pdf')]:
 d=fitz.open(path);d[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review'/f'Figure1_{side}.png')
(O/'notes/checks.json').write_text(json.dumps({'previous_in_map_labels':old_texts,'new_labels':16,'complete_legend_names':changed,'all_embedded_images_unchanged':True,'n_cells':len(D),'anchor_rule':'actual cell nearest within-L1 coordinate median; display anchor only, no new clustering or annotation','no_scientific_analysis':True,'source_svg':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'visual_review':'PENDING'},indent=2))
print(A[['number','L1','anchor_page_x_pt','anchor_page_y_pt','side']].to_string(index=False))
