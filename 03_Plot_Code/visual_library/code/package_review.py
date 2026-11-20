"""Native-size review pages and field integrity checks. No scientific computation."""
from pathlib import Path
import json,hashlib,platform,importlib.metadata as md
import pandas as pd
import pymupdf as fitz
ROOT=Path(__file__).resolve().parents[1];P=ROOT.parents[1];MM=72/25.4

def pp(name):
 return ROOT/('review' if name=='Figure3de_shared_rows' else name[:-1])/f'{name}.pdf'
def size(name):
 with fitz.open(pp(name)) as d:return d[0].rect.width,d[0].rect.height

def layout(groups,columns=False):
 dims=[[size(n) for n in g] for g in groups];gap=2.5*MM;pad=3*MM;boxes=[]
 if columns:
  H=max(sum(h for w,h in dd)+gap*(len(dd)-1) for dd in dims)+2*pad;x=pad
  for row,dd in zip(groups,dims):
   cw=max(w for w,h in dd);y=pad
   for name,(w,h) in zip(row,dd):boxes.append((name,fitz.Rect(x,y,x+w,y+h)));y+=h+gap
   x+=cw+gap
  return x-gap+pad,H,boxes
 W=max(sum(w for w,h in dd)+gap*(len(dd)-1) for dd in dims)+2*pad;y=pad
 for row,dd in zip(groups,dims):
  total=sum(w for w,h in dd)+gap*(len(dd)-1);x=pad;rh=max(h for w,h in dd)
  for name,(w,h) in zip(row,dd):boxes.append((name,fitz.Rect(x,y,x+w,y+h)));x+=w+gap
  y+=rh+gap
 return W,y-gap+pad,boxes

def figure5_compact():
 # Review only: one native-size top row; forest + settings/program column +
 # module scatter below. No source panel is scaled or cropped.
 gap=2.5*MM;pad=3*MM;boxes=[];x=pad;top_h=0
 for name in ['Figure5a','Figure5b']:
  w,h=size(name);boxes.append((name,fitz.Rect(x,pad,x+w,pad+h)));x+=w+gap;top_h=max(top_h,h)
 W=x-gap+pad;yt=pad+top_h+gap;x=pad;H=yt
 for group in [['Figure5c'],['Figure5d','Figure5f'],['Figure5e']]:
  y=yt;cw=max(size(n)[0] for n in group)
  for name in group:
   w,h=size(name);boxes.append((name,fitz.Rect(x,y,x+w,y+h)));y+=h+gap
  H=max(H,y-gap);x+=cw+gap
 return max(W,x-gap+pad),H+pad,boxes

layouts={
'Figure1':layout([['Figure1a','Figure1b'],['Figure1c','Figure1f','Figure1d'],['Figure1e','Figure1g']]),
'Figure2':layout([['Figure2a'],['Figure2b','Figure2c'],['Figure2d'],['Figure2e','Figure2f']]),
'Figure3':layout([['Figure3a','Figure3b','Figure3f'],['Figure3c','Figure3de_shared_rows']],True),
'Figure4':layout([['Figure4a'],['Figure4b'],['Figure4c','Figure4d'],['Figure4e','Figure4f']]),
'Figure5':figure5_compact(),
'ExtendedData_Figure2':layout([['ExtendedData_Figure2a'],['ExtendedData_Figure2b']]),
'ExtendedData_Figure3':layout([['ExtendedData_Figure3a'],['ExtendedData_Figure3b']]),
'ExtendedData_Figure4':layout([['ExtendedData_Figure4a','ExtendedData_Figure4b'],['ExtendedData_Figure4c','ExtendedData_Figure4d']]),
'ExtendedData_Figure5':layout([['ExtendedData_Figure5a'],['ExtendedData_Figure5b']])}
allpdf=fitz.open();mainpdf=fitz.open();edpdf=fitz.open();placements=[]
for figure,(W,H,boxes) in layouts.items():
 doc=fitz.open();page=doc.new_page(width=W,height=H)
 for name,rect in boxes:
  with fitz.open(pp(name)) as src:
   assert abs(rect.width-src[0].rect.width)<1e-6 and abs(rect.height-src[0].rect.height)<1e-6
   page.show_pdf_page(rect,src,0,keep_proportion=True)
   placements.append(dict(figure=figure,panel=name,source_pdf=str(pp(name)),scale_x=1.,scale_y=1.,x_mm=rect.x0/MM,y_mm=rect.y0/MM,width_mm=rect.width/MM,height_mm=rect.height/MM))
 doc.set_metadata(dict(title=f'{figure} panel review',author='',creator='',producer=''))
 doc.save(ROOT/f'review/{figure}_整图审阅.pdf',garbage=4,deflate=True)
 page.get_pixmap(matrix=fitz.Matrix(1800/W,1800/W),alpha=False).save(ROOT/f'review/{figure}_整图审阅.png')
 allpdf.insert_pdf(doc);(edpdf if figure.startswith('Extended') else mainpdf).insert_pdf(doc)
for doc,fn in [(mainpdf,'主图1–5_整图审阅预览.pdf'),(edpdf,'ExtendedData2–5_整图审阅预览.pdf'),(allpdf,'主图及ED_整图审阅预览.pdf')]:
 doc.set_metadata(dict(title='Panvascular independent panel review',author='',creator='',producer=''));doc.save(ROOT/fn,garbage=4,deflate=True)
pd.DataFrame(placements).to_csv(ROOT/'review/preview_placements.tsv',sep='\t',index=False)
records=json.loads((ROOT/'review/panel_manifest.json').read_text());assert len(records)==41 and len({r['figure']+r['panel'] for r in records})==41
hashfile=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rows=[];checks=[];used={}
def sortkey(r):return (r['figure'].startswith('Extended'),int(r['figure'].split('Figure')[-1]),r['panel'])
for r in sorted(records,key=sortkey):
 name=r['figure']+r['panel'];pdf=pp(name);svg=pdf.with_suffix('.svg');doc=fitz.open(pdf);page=doc[0]
 spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
 outside=[s['text'] for s in spans if not (page.rect+(-.5,-.5,.5,.5)).contains(fitz.Rect(s['bbox']))]
 checks.append(dict(panel=name,pages=len(doc),outside_text=outside,pdf_text_spans=len(spans),minimum_font_pt=min((s['size'] for s in spans),default=None),width_mm=page.rect.width/MM,height_mm=page.rect.height/MM,pdf_sha256=hashfile(pdf),svg_sha256=hashfile(svg)))
 for source in r['source_files']:
  if Path(source).exists():used[source]=hashfile(source)
 rows.append(dict(panel=name,operation=r['mode'],source_files=' | '.join(r['source_files']),change=r['change'],width_mm=round(page.rect.width/MM,3),height_mm=round(page.rect.height/MM,3),pdf=str(pdf.relative_to(ROOT)),svg=str(svg.relative_to(ROOT))))
pd.DataFrame(rows).to_csv(ROOT/'逐小图_来源与修改.tsv',sep='\t',index=False)
before=json.loads((ROOT/'review/this_round_protected.json').read_text());assert all(hashfile(p)==v for p,v in before.items())
checksout=dict(panels=checks,protected_files_unchanged=len(before),source_hashes=used,scientific_analysis_run=False,preview_scaling='1:1 native PDF points; no anisotropic scaling or cropping',shared_rows=json.loads((ROOT/'review/shared_rows.json').read_text()),visual_review='See 逐图视觉检查.tsv for actual human-readable inspection findings; script checks do not substitute visual review')
(ROOT/'review/export_checks.json').write_text(json.dumps(checksout,ensure_ascii=False,indent=2))
env={k:md.version(k) for k in ['numpy','pandas','matplotlib','PyMuPDF','lxml','CairoSVG']};env['python']=platform.python_version();(ROOT/'code/environment.json').write_text(json.dumps(env,indent=2))
style='body{margin:0;background:#f0f2f4;color:#25343c;font:16px Arial,sans-serif}header{background:white;padding:24px 4vw;border-bottom:1px solid #ccd5da}main{max-width:1500px;margin:auto;padding:24px}h2{margin-top:45px}nav a{margin-right:20px}a{color:#276a92}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(390px,1fr));gap:16px}.panel{background:white;padding:16px;border-radius:8px;box-shadow:0 2px 8px #0001}.panel img{width:100%;height:430px;object-fit:contain}.links{display:flex;gap:18px}.preview{display:block;width:100%;max-height:1000px;object-fit:contain;background:white}p{line-height:1.5}.badge{font-size:13px;color:#637078}'
body=['<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>泛血管小图审阅</title><style>'+style+'</style><header><h1>泛血管小图审阅 · 当前版</h1><p>31个主图小图；同步重导出ED2–5的五张marker图，其余配套地图/热图沿用。点击图片查看原生尺寸PDF。</p><nav><a href="主图及ED_整图审阅预览.pdf">主图及ED九页预览</a><a href="主图1–5_整图审阅预览.pdf">主图五页预览</a><a href="review/重点小图_同倍率前后对照.pdf">同倍率对照</a></nav></header><main>']
for figure in layouts:
 body.append(f'<h2 id="{figure}">{figure}</h2><details><summary>打开整图审阅预览</summary><a href="review/{figure}_整图审阅.pdf"><img class="preview" src="review/{figure}_整图审阅.png"></a></details><div class="grid">')
 for r in rows:
  if r['panel'][:-1]!=figure:continue
  name=r['panel'];body.append(f'<article class="panel"><h3>{name}</h3><a href="{r["pdf"]}"><img loading="lazy" src="{figure}/{name}.png"></a><p class="badge">{r["width_mm"]} × {r["height_mm"]} mm</p><div class="links"><a href="{r["pdf"]}">PDF</a><a href="{r["svg"]}">SVG</a></div></article>')
 body.append('</div>')
body.append('</main></html>');(ROOT/'打开小图审阅.html').write_text(''.join(body))
comparison=fitz.open()
for name in ['Figure1g','Figure3c','Figure4c','Figure5a','ExtendedData_Figure3b','ExtendedData_Figure5b']:
 old=ROOT/'review/本轮修改前'/name[:-1]/f'{name}.pdf'
 if not old.exists():old=P/'独立小图_最新版'/name[:-1]/f'{name}.pdf'
 with fitz.open(old) as a,fitz.open(pp(name)) as b:
  w1,h1=a[0].rect.width,a[0].rect.height;w2,h2=b[0].rect.width,b[0].rect.height;pad=4*MM;gap=5*MM;top=10*MM
  pg=comparison.new_page(width=w1+w2+gap+2*pad,height=max(h1,h2)+top+pad)
  pg.insert_text((pad,6*MM),f'{name} | Before',fontsize=8);pg.insert_text((pad+w1+gap,6*MM),f'{name} | After',fontsize=8)
  pg.show_pdf_page(fitz.Rect(pad,top,pad+w1,top+h1),a,0);pg.show_pdf_page(fitz.Rect(pad+w1+gap,top,pad+w1+gap+w2,top+h2),b,0)
comparison.set_metadata(dict(title='Panels before and after at native physical size',author='',creator='',producer=''));comparison.save(ROOT/'review/重点小图_同倍率前后对照.pdf',garbage=4,deflate=True)
print(json.dumps(dict(panels=len(records),main_panels=sum(r['figure'].startswith('Figure') for r in records),ED_redrawn=5,protected_unchanged=len(before),outside_text={r['panel']:r['outside_text'] for r in checks if r['outside_text']},review_pages=len(allpdf)),ensure_ascii=False))
