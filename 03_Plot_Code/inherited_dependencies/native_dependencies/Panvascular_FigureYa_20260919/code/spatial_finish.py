from pathlib import Path
import xml.etree.ElementTree as E,copy,json,hashlib,subprocess,os
import pandas as pd,numpy as np,fitz
O=Path(__file__).resolve().parents[1];P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');V=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Figure_Redesign_20260917_v3')
for d in ['main','editable_svg','review','notes']:(O/d).mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
# This script composes a native-R FigureYa forest table with protected current panels.
ns='http://www.w3.org/2000/svg';E.register_namespace('',ns);E.register_namespace('xlink','http://www.w3.org/1999/xlink')
r=E.parse(P/'图件编辑资料/editable_svg/Figure5.svg').getroot();f=next(x for x in r.iter() if x.get('id')=='figure_1');base=copy.deepcopy(r)
old=next(x for x in f if x.get('id')=='axes_7');idx=list(f).index(old);f.remove(old)
q=E.parse(O/'panels/Figure5c_FigureYa.svg').getroot()
# Namespacing all native grid SVG IDs avoids collisions with matplotlib definitions.
idmap={e.get('id'):'fy5c_'+e.get('id') for e in q.iter() if e.get('id')}
for e in q.iter():
 for k,v in list(e.attrib.items()):
  if k=='id':e.set(k,idmap[v])
  else:
   for old_id,new_id in idmap.items():v=v.replace('url(#'+old_id+')','url(#'+new_id+')').replace('#'+old_id+'"','#'+new_id+'"')
   e.set(k,v)
M=72/25.4
# Native forest output is placed 1:1 (98x96 mm); no distorted text/scaled axes.
g=E.Element('{'+ns+'}g',{'id':'FigureYa6_native_forest','transform':f'translate({1*M},{73*M})'})
for e in list(q):g.append(copy.deepcopy(e))
f.insert(idx,g)
# The rest of current Figure5 is protected exactly, including actual histology/spot image panels.
check=copy.deepcopy(r);cf=next(x for x in check.iter() if x.get('id')=='figure_1');cf.remove(next(x for x in cf if x.get('id')=='FigureYa6_native_forest'))
cf0=next(x for x in base.iter() if x.get('id')=='figure_1');cf0.remove(next(x for x in cf0 if x.get('id')=='axes_7'))
assert E.tostring(check)==E.tostring(base)
svg=O/'editable_svg/Figure5.svg';E.ElementTree(r).write(svg,encoding='utf-8',xml_declaration=True)
env=os.environ.copy();env['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
subprocess.run(['/Users/taki/Desktop/Codex任务1/.venv_figure/bin/python','/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/code/export_svg_pdf.py',str(svg),str(O/'main/Figure5.pdf')],env=env,check=True)
a=fitz.open(P/'Figures/Figure5.pdf');b=fitz.open(O/'main/Figure5.pdf');assert a[0].rect==b[0].rect
b[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(O/'review/spatial_Figure5.png')
from PIL import Image,ImageDraw
ims=[]
for d in [a,b]:
 pix=d[0].get_pixmap(matrix=fitz.Matrix(1.6,1.6));im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples);ims.append(im)
can=Image.new('RGB',(ims[0].width*2,ims[0].height+24),'white');can.paste(ims[0],(0,24));can.paste(ims[1],(ims[0].width,24));dr=ImageDraw.Draw(can);dr.text((10,4),'CURRENT BASELINE',fill='black');dr.text((ims[0].width+10,4),'FIGUREYA ADAPTED: NATIVE R FOREST TABLE',fill='black');can.save(O/'review/spatial_Figure5_before_after.png')
j=pd.read_csv(O/'source_data/spatial/primary_shared29.tsv',sep='\t');assert len(j)==29
assert sha(O/'source_data/spatial/primary_shared29.tsv')==sha(V/'source_data/derived/spatial_primary_shared29.tsv')
# Record actual source bindings without numerical inference from a rendered plot.
record={'Figure5a_b_d_e_f':'baseline SVG groups unchanged; all spatial image/spot coordinates, colour limits, matrix, 14 modules, 35 points and 5 median glyphs retained byte-for-byte', 'Figure5c':{'rows':29,'source':str(V/'source_data/derived/spatial_primary_shared29.tsv'),'sha256':sha(O/'source_data/spatial/primary_shared29.tsv'),'fields':['gene_symbol','cell_state','beta','CI_low','CI_high','FDR_BH_all_29_pairs_x_3_terms'],'new_display':'FDR numeric column; separate gene and lineage columns; aligned native forestplot table','template':'FigureYa6rHRsV2/FigureYa6rHRsV2.Rmd','original_lines':'113-160','adaptations':['linear coefficient axis and zero=0 (not hazard ratio log axis)','source significance colours retained','fixed 29-row order; no sorting','6pt text; native R forestplot 3.2.0; svglite 2.2.2','original model results only; no fit or testing']},'size_mm':[183,170],'scientific_analysis_run':False,'visual_review':'PENDING_ROOT_REVIEW'}
(O/'notes/spatial_source_mapping.json').write_text(json.dumps(record,ensure_ascii=False,indent=2))
(O/'notes/spatial_legend_changes.md').write_text('Figure5c: The forest plot now separates gene and cell-state labels and reports the original BH FDR alongside each estimate and 95% confidence interval. Fib, fibroblast; EC, endothelial cell; Mac, macrophage. No estimates or multiplicity family were changed.\nOther Figure5 panels retain the current adopted display and legend meaning.\n')
print('Figure5 native R template adaptation completed')
