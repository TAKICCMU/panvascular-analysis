"""Restore the native BioRender vessel/body geometry; all other content is reused."""
from pathlib import Path
from lxml import etree as E
import fitz, numpy as np, json, copy, hashlib
P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
R=P/'图件编辑资料/Figure1a_人体血管对齐'
B=R/'backup'
N='http://www.w3.org/2000/svg'
native=fitz.open(R/'BioRender_original_export.pdf')
info={x['xref']:fitz.Rect(x['bbox']) for x in native[0].get_image_info(xrefs=True)}
body,vessel=info[27],info[28]
target=fitz.Rect(14.16,24.604724,30.64,51.964724)
sx=target.width/body.width;sy=target.height/body.height
new=fitz.Rect(target.x0+(vessel.x0-body.x0)*sx,target.y0+(vessel.y0-body.y0)*sy,target.x0+(vessel.x1-body.x0)*sx,target.y0+(vessel.y1-body.y0)*sy)
old=fitz.Rect(17.28,24.764724,29.20,51.804724)
old_cmd=b'11.92 0 0 27.04 17.28 24.764724 cm\n/a0 gs /x8 Do'
new_cmd=(f'{new.width:.9f} 0 0 {new.height:.9f} {new.x0:.9f} {new.y0:.9f} cm\n/a0 gs /x8 Do').encode()
checks=[]
for rel in ['Figure/Figure1.pdf','Figures/Figure1.pdf']:
 src=B/rel;d=fitz.open(src);before=d[0].get_pixmap(matrix=fitz.Matrix(4,4),alpha=False);text=d[0].get_text()
 count=0
 for xr in d[0].get_contents():
  s=d.xref_stream(xr);count+=s.count(old_cmd)
  if old_cmd in s:d.update_stream(xr,s.replace(old_cmd,new_cmd))
 assert count==1
 out=R/'staged'/rel;out.parent.mkdir(parents=True,exist_ok=True);d.save(out,garbage=4,deflate=True)
 after=fitz.open(out);assert after[0].get_text()==text;assert after[0].rect==d[0].rect
 a=np.frombuffer(before.samples,np.uint8).reshape(before.height,before.width,3)
 px=after[0].get_pixmap(matrix=fitz.Matrix(4,4),alpha=False)
 b=np.frombuffer(px.samples,np.uint8).reshape(px.height,px.width,3)
 delta=np.any(a!=b,axis=2);allchanged=int(delta.sum());u=old|new
 delta[max(0,int(u.y0*4)-2):int(u.y1*4)+3,max(0,int(u.x0*4)-2):int(u.x1*4)+3]=False
 assert not delta.any(),int(delta.sum())
 checks.append({'file':rel,'changed_pixels_within_vessel_region':allchanged,'changed_pixels_outside_region':int(delta.sum()),'text_identical':True,'page_size_identical':True})
for rel in ['Figure/Editable_SVG/Figure1.svg','Figures/Figure1.svg']:
 root=E.parse(str(B/rel));im=root.xpath('//*[@id="a_image08f78ba088"]')[0]
 im.set('x',f'{new.x0:.9f}');im.set('y',f'{-new.y0:.9f}');im.set('width',f'{new.width:.9f}');im.set('height',f'{new.height:.9f}');im.set('transform',f'scale(1 -1) translate(0 {-new.height:.9f})')
 out=R/'staged'/rel;out.parent.mkdir(parents=True,exist_ok=True);root.write(str(out),encoding='utf-8',xml_declaration=True)
 if rel.startswith('Figure/'):
  a=copy.deepcopy(root.xpath('//*[@id="panel_a"]')[0]);s=E.Element('{'+N+'}svg',nsmap={None:N,'xlink':'http://www.w3.org/1999/xlink'},width='18.3cm',height='2.9cm',viewBox='0 0 518.740157480315 82.204724409449')
  s.append(a);E.ElementTree(s).write(str(R/'Figure1a_aligned.svg'),encoding='utf-8',xml_declaration=True)
q=fitz.open(R/'staged/Figure/Figure1.pdf');pg=q[0]
pg.get_pixmap(matrix=fitz.Matrix(16,16),clip=fitz.Rect(12,22,32,56)).save(R/'人体血管_修正后.png')
pg.get_pixmap(matrix=fitz.Matrix(3,3)).save(R/'Figure1_preview.png')
(R/'alignment_check.json').write_text(json.dumps({'native_body_bbox':list(body),'native_vessel_bbox':list(vessel),'current_body_bbox_preserved':list(target),'old_vessel_bbox':list(old),'new_vessel_bbox':list(new),'method':'Transfer relative placement from native BioRender group; body unchanged; no pixel values or other figure objects changed','checks':checks},ensure_ascii=False,indent=2))
print(json.dumps({'new_vessel':list(new),'checks':checks},ensure_ascii=False))
