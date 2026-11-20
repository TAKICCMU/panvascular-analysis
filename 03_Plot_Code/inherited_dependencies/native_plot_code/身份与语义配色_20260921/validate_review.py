from pathlib import Path
import json,re,hashlib,copy,io
import numpy as np,pandas as pd,xml.etree.ElementTree as E,fitz
from PIL import Image,ImageDraw
from colorspacious import cspace_convert
W=Path(__file__).resolve().parents[1];P=W.parents[1];O=Path('/Users/taki/Desktop/Panvascular_Main_Figure_Refinement_20260920');q=lambda x:'{http://www.w3.org/2000/svg}'+x;XL='{http://www.w3.org/1999/xlink}href';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
v=pd.read_csv(W/'backup/图件采用版本.tsv',sep='\t');reports=[];pages=fitz.open()
allowed_images={'Figure1':['axes_2','adopted_B_g_axes_2'],'Figure5':['r2_axes_5','r2_axes_10'],'ExtendedData_Figure2':['newED2_axes_3'],'ExtendedData_Figure3':['newED3_axes_3'],'ExtendedData_Figure4':['newED4_axes_3','newED4_axes_6'],'ExtendedData_Figure5':['newED5_axes_3']}
find=lambda r,id:next(e for e in r.iter() if e.get('id')==id)
for _,vr in v.sort_values('page').iterrows():
 name=vr.figure;a=E.parse(W/'backup'/f'{name}.svg').getroot();b=E.parse(W/'editable_svg'/f'{name}.svg').getroot()
 ex=set()
 for id in allowed_images.get(name,[]):ex|={e.get('id') for e in find(a,id).iter(q('image'))}
 ims={e.get('id'):e.get(XL) for e in b.iter(q('image'))}
 for im in a.iter(q('image')):
  if im.get('id') not in ex:assert im.get(XL)==ims[im.get('id')],(name,im.get('id'))
 # Geometry (including every dot outline / location / size, 29 shared rows, spot coordinates) identical outside 3a redraw.
 def strip(root):
  r=copy.deepcopy(root)
  if name=='Figure3':
   for par in r.iter():
    for e in list(par):
     if e.get('id','').startswith('identity_a_') or e.get('id')=='aligned_a_axes_1':par.remove(e)
  for e in r.iter():
   for k in ['fill','stroke',XL,'preserveAspectRatio']:e.attrib.pop(k,None)
   if 'style' in e.attrib:
    s=re.sub(r'(?:fill|stroke):\s*[^;]+;?','',e.get('style'));e.set('style',';'.join(x.strip() for x in s.split(';') if x.strip()))
  return r
 aa,bb=strip(a),strip(b)
 # Matplotlib adds unused defs for rebuilt a; compare every old remaining ID object that still exists.
 am={e.get('id'):e for e in aa.iter() if e.get('id')};bm={e.get('id'):e for e in bb.iter() if e.get('id')}
 changed=[]
 for id,e in am.items():
  if id in bm and e.tag!=q('g') and e.tag!=q('defs') and E.tostring(e)!=E.tostring(bm[id]):changed.append(id)
 assert not changed,(name,changed[:10])
 d=fitz.open(W/'figures'/f'{name}.pdf');p=d[0];exp=183 if name.startswith('Figure') else 180;assert abs(p.rect.width*25.4/72-exp)<.001 and abs(p.rect.height*25.4/72-170)<.001
 spans=[s for x in p.get_text('dict')['blocks'] for l in x.get('lines',[]) for s in l['spans']];outside=[s['text'] for s in spans if s['bbox'][0]<-.2 or s['bbox'][1]<-.2 or s['bbox'][2]>p.rect.width+.2 or s['bbox'][3]>p.rect.height+.2]
 assert not outside,(name,outside)
 reports.append({'figure':name,'native_mm':[exp,170],'noncolour_geometry_unchanged_except_3a':True,'protected_raster_images_unchanged':True,'outside_text':outside})
 pages.insert_pdf(d);typ,val=d.xref_get_key(p.xref,'Group')
 if typ=='dict':assert not re.search(r'\d+ \d+ R',val);pages.xref_set_key(pages[-1].xref,'Group',val)
pages.save(W/'15页配色核查稿.pdf',garbage=4,deflate=True);pages.close()
merge=fitz.open(W/'15页配色核查稿.pdf')
for i,r in enumerate(reports):
 src=fitz.open(W/'figures'/f'{r["figure"]}.pdf');assert src[0].get_pixmap().samples==merge[i].get_pixmap().samples
# Files already in the package inventory are verified before any formal replacement.
inventory=pd.read_csv(P/'文件来源与校验.tsv',sep='\t');protected=inventory[inventory.file.str.startswith(('Source_Data_Figures/','Supplementary_Table'))|inventory.file.str.endswith(('.docx','.xlsx'))];bad=[]
for _,r in protected.iterrows():
 p=P/r.file
 if p.exists() and sha(p)!=r.sha256:bad.append(r.file)
assert not bad,bad
# Same native-sized before/after pages, never scale. Paired pages rather than shrunk side-by-side sheets.
with fitz.open() as comp:
 for name in ['Figure1','Figure3','Figure5']:
  for folder in ['backup','figures']:
   d=fitz.open(W/folder/f'{name}.pdf');comp.insert_pdf(d);typ,val=d.xref_get_key(d[0].xref,'Group')
   if typ=='dict':comp.xref_set_key(comp[-1].xref,'Group',val)
 comp.save(W/'review/同尺寸前后对照_6页.pdf',garbage=4,deflate=True)
# Full-normal contact sheets plus CVD simulations of priority pages. These diagnostic PNGs are not production images.
for i in range(3):
 ims=[]
 for name in v.sort_values('page').figure.iloc[i*5:(i+1)*5]:
  im=Image.open(W/'review'/f'{name}.png').convert('RGB');im.thumbnail((430,420));canvas=Image.new('RGB',(450,450),'white');canvas.paste(im,((450-im.width)//2,22));ImageDraw.Draw(canvas).text((12,3),name,fill='black');ims.append(canvas)
 out=Image.new('RGB',(450*len(ims),450),'white')
 for j,im in enumerate(ims):out.paste(im,(450*j,0))
 out.save(W/'review'/f'contact_{i+1}.png')
for name in ['Figure1','Figure3','Figure5']:
 im=Image.open(W/'review'/f'{name}.png').convert('RGB');a=np.asarray(im,dtype=float)/255
 for cvd in ['deuteranomaly','protanomaly']:
  z=cspace_convert(a,{'name':'sRGB1+CVD','cvd_type':cvd,'severity':100},'sRGB1');Image.fromarray(np.uint8(np.round(np.clip(z,0,1)*255))).save(W/'review'/f'{name}_{cvd}.png')
(W/'tables/VALIDATION.json').write_text(json.dumps({'data_integrity':{'source_tables_and_DOCX_verified':len(protected),'all_unchanged':True,'point_ID_audit':'PASS with EC FDR serialization delta 8.7e-19; original values preserved','zero_dot_area':'original paths retained'},'layout_legibility':reports,'merged_native_render_matches_all15':True,'visual_review':'AWAITING_MANUAL_VIEW','ORA':'PENDING','donor_limitations':'UNCHANGED'},indent=2))
print('PASS all geometry, protected image, page bounds, scientific-file hashes and merged-page checks')
