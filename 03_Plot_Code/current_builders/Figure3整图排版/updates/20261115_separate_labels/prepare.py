from pathlib import Path
import fitz, json, hashlib, shutil, datetime
from lxml import etree as E
import numpy as np

A=Path('/Users/taki/Desktop/投稿包_最新版_20260919/图件编辑资料/BioRender_当前素材/Figure3整图排版')
R=Path(__file__).resolve().parent
O=R/'panels_without_labels'; U=R/'upload_ready'; B=R/'backup_before_split'
for p in [O,U,B]:p.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for p in [A/'review/biorender_status.json',A/'panel_placement.json',A/'README.md',A/'review/Figure3_BioRender.png',A/'review/Figure3_BioRender.pdf']:
 q=B/p.relative_to(A);q.parent.mkdir(parents=True,exist_ok=True)
 if not q.exists():shutil.copy2(p,q)
report={'status':'LABEL_FREE_ASSETS_READY','created_at':datetime.datetime.now().astimezone().isoformat(),'panels':{},'data_and_formal_figure_unchanged':True}
for k in 'abcdef':
 src=A/'panels'/f'Figure3{k}.pdf'; old=fitz.open(src); doc=fitz.open(src); p=doc[0]
 candidates=[s for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['text']==k and s['bbox'][0]<15 and s['bbox'][1]<15]
 assert len(candidates)==1,(k,candidates)
 s=candidates[0];rect=fitz.Rect(s['bbox']);rect.x0-=.15;rect.x1+=.15;rect.y0=max(0,rect.y0-.15);rect.y1+=.15
 p.add_redact_annot(rect,fill=(1,1,1));p.apply_redactions(images=0,graphics=0)
 pdf=O/f'Figure3{k}_nolabel.pdf';doc.save(pdf,garbage=4,deflate=True);p=doc[0]
 assert not [s for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['text']==k and s['bbox'][0]<15 and s['bbox'][1]<15]
 dpi=600; pix=p.get_pixmap(dpi=dpi,alpha=False); png=O/f'Figure3{k}_nolabel.png';pix.save(png)
 a=old[0].get_pixmap(dpi=dpi,alpha=False);b=pix
 ar=np.frombuffer(a.samples,dtype=np.uint8).reshape(a.height,a.width,3);br=np.frombuffer(b.samples,dtype=np.uint8).reshape(b.height,b.width,3)
 mask=np.any(ar!=br,axis=2);ys,xs=np.where(mask)
 rr=rect*fitz.Matrix(dpi/72,dpi/72)
 assert xs.min()>=int(rr.x0)-2 and xs.max()<=int(rr.x1)+2 and ys.min()>=int(rr.y0)-2 and ys.max()<=int(rr.y1)+2,(k,xs.min(),xs.max(),ys.min(),ys.max(),rr)
 svg=O/f'Figure3{k}_nolabel.svg'
 if k in 'def':
  root=E.parse(str(A/'panels'/f'Figure3{k}.svg'))
  found=[e for e in root.xpath('//*[local-name()="text"]') if ''.join(e.itertext()).strip()==k]
  assert len(found)==1
  found[0].getparent().remove(found[0]);root.write(str(svg),encoding='utf-8',xml_declaration=True)
 else:svg.write_text(p.get_svg_image(text_as_path=True))
 upload=svg if k in 'de' else png
 shutil.copy2(upload,U/upload.name)
 report['panels'][k]={'source_pdf':str(src),'source_sha256':sha(src),'label_rect_pt':list(rect),'font':s['font'],'font_pt':s['size'],'dimensions_pt':list(p.rect),'only_label_pixels_changed':True,'changed_pixels':int(mask.sum()),'upload':str(U/upload.name),'upload_sha256':sha(U/upload.name)}
(R/'audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'prepared':len(report['panels']),'only_label_pixels_changed':True,'upload_directory':str(U)},ensure_ascii=False))
