from pathlib import Path
import os,sys,hashlib,shutil,json,copy,re
os.environ['DYLD_FALLBACK_LIBRARY_PATH']='/Users/taki/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/lib'
sys.path.insert(0,'/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Local_Touchup_20260918/tools')
import xml.etree.ElementTree as E,pandas as pd,fitz,cairosvg,cairosvg.image
from cairosvg.surface import cairo
cairosvg.image.IMAGE_RENDERING.update({'pixelated':cairo.FILTER_NEAREST,'crisp-edges':cairo.FILTER_NEAREST})
from visibility import apply,Q,find
W=Path(__file__).resolve().parents[1];P=W.parents[1];O=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920');PRE=P/'图件编辑资料/身份与语义配色_20260921';names=['Figure5']+[f'ExtendedData_Figure{i}' for i in [2,3,4,5]];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();mm=72/25.4
for d in ['backup','figures','editable_svg','review']: (W/d).mkdir(exist_ok=True)
E.register_namespace('','http://www.w3.org/2000/svg');E.register_namespace('xlink','http://www.w3.org/1999/xlink')
v=pd.read_csv(P/'图件采用版本.tsv',sep='\t')
for _,r in v[v.figure.isin(names)].iterrows():
 for field in ['package_pdf','package_svg']:
  src=Path(r[field]);dst=W/'backup'/src.name
  if not dst.exists():shutil.copy2(src,dst)
for src in [PRE/'semantic_colors.json',PRE/'CURRENT_PENDING_CAPTIONS.md',PRE/'code/recolour.py',P/'图件采用版本.tsv',O/'FIGURE_VERSIONS.tsv']:
 dst=W/'backup'/src.name
 if not dst.exists():shutil.copy2(src,dst)
conf=json.loads((PRE/'semantic_colors.json').read_text());conf.update({'NO_SUPPORT_FILL':'#D9DDE1','NO_SUPPORT_POINT':'#747C86','NO_SUPPORT_CI':'#747C86'});(PRE/'semantic_colors.json').write_text(json.dumps(conf,indent=2,ensure_ascii=False));shutil.copy2(PRE/'semantic_colors.json',P/'图件编辑资料/semantic_colors_current.json')
# Preserve this bounded finishing rule in the current rendering entrypoint, without invoking the full renderer.
p=PRE/'code/recolour.py';s=p.read_text();marker='# Bounded approved visibility finishing rules'
if marker not in s:
 s=s.replace('def export(r,name):',"""def export(r,name):
 # Bounded approved visibility finishing rules
 import importlib.util
 vp=P/'图件编辑资料/点线可见性收尾_20260921/code/visibility.py'
 spec=importlib.util.spec_from_file_location('approved_visibility',vp);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.apply(r,name,C)""")
 p.write_text(s)
rows=[]
for name in names:
 root=E.parse(W/'backup'/f'{name}.svg').getroot();before=copy.deepcopy(root);changes=apply(root,name,conf);assert changes
 out=W/'editable_svg'/f'{name}.svg';E.ElementTree(root).write(out,encoding='utf-8',xml_declaration=True);cairosvg.svg2pdf(url=str(out),write_to=str(W/'figures'/f'{name}.pdf'))
 # Reversing these exact attributes must reproduce the complete baseline SVG, including all data geometry/images/text.
 for e,old,new,group,role in changes:e.set('style',old)
 assert E.tostring(root)==E.tostring(before),name
 row={'figure':name,'changed_attributes':len(changes),'only_targeted_style_changed':True,'main_marker_paths':0}
 if name=='Figure5':
  row['grey_points']=sum(e.tag==Q('circle') for e,*_ in changes);row['grey_CI_lines']=sum(e.tag==Q('polyline') for e,*_ in changes);assert row['grey_points']==19 and row['grey_CI_lines']==57
 else:
  row['main_marker_paths']=sum(e.tag==Q('path') and 'PathCollection' in group and 'clip-path' in e.attrib for e,old,new,group,role in changes)
 rows.append(row)
 d=fitz.open(W/'figures'/f'{name}.pdf');assert abs(d[0].rect.height/mm-170)<.001;d[0].get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(W/'review'/f'{name}.png')
# Native-scale clipped comparisons; both sides use the same clipping rectangle and 1:1 point dimensions.
comp=fitz.open();clips=[('Figure5',(0,63,97,170)),('ExtendedData_Figure2',(0,82,180,166)),('ExtendedData_Figure3',(0,49,180,164)),('ExtendedData_Figure4',(63,12,180,162)),('ExtendedData_Figure5',(0,60,180,166))]
for name,coords in clips:
 c=fitz.Rect(*(x*mm for x in coords));page=comp.new_page(width=2*c.width+8*mm,height=c.height+8*mm)
 for side,folder in enumerate(['backup','figures']):
  d=fitz.open(W/folder/f'{name}.pdf');x=side*(c.width+8*mm);page.show_pdf_page(fitz.Rect(x,8*mm,x+c.width,8*mm+c.height),d,0,clip=c,keep_proportion=True);page.insert_text((x+2*mm,5*mm),'Before' if side==0 else 'After',fontsize=7)
comp.save(W/'review/同尺寸局部前后对照.pdf',garbage=4,deflate=True)
# Pending captions only. No figure3 re-export or numerical labels added.
(W/'图注补充_待采用.md').write_text('''Figure3a待采用补充：Filled circles indicate that both the original MAGMA and scDRS FDR values are <0.05. Open circles indicate that this joint criterion is not met; they do not imply that neither analysis is significant. Fibromyocyte and Smooth muscle remain distinct L1 categories.

Figure5c待采用补充：Estimates without support under the existing criterion and their confidence intervals are shown in medium grey. The classification, estimates and intervals are unchanged.

Marker待采用补充：Marker dots have no independent outline; area and fill retain the original detection-fraction and expression encodings.

正文未写回，ORA与其他科学待办保持原状态。
''')
# No new audit system: one small change count file records scope.
(W/'changes.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2));print('READY_FOR_VISUAL_CHECK')
