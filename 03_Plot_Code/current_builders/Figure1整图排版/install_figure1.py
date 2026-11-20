"""Install only Figure1 and current review aggregates. No scientific analysis."""
from pathlib import Path
import csv, io, json, hashlib, shutil, zipfile
import pymupdf as fitz
import numpy as np

P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
W=Path(__file__).parent
B=P/'图件编辑资料/轻量备份/Figure1整图调整前_20261115'
C=P/'图件编辑资料/可视化库重绘_当前'
S=P/'独立小图_最新版'
URL='https://app.biorender.com/illustrations/6af93b20061047a1155dc09d'
changed=set()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def backup(p):
    p=Path(p);q=B/p.relative_to(P)
    if p.exists() and not q.exists():q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
def write(p,b):
    p=Path(p);backup(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_bytes(b if isinstance(b,bytes) else b.encode());changed.add(p)
def copy(src,dst):write(dst,Path(src).read_bytes())
def save(d,p):
    d.set_metadata({});d.del_xml_metadata();write(p,d.tobytes(garbage=4,deflate=True))
def merge(paths,p):
    d=fitz.open()
    for q in paths:
        with fitz.open(q) as s:d.insert_pdf(s)
    save(d,p)
def tsvread(p):
    with Path(p).open() as h:
        r=csv.DictReader(h,delimiter='\t');return r.fieldnames,list(r)
def tsvwrite(p,fields,rows):
    o=io.StringIO();w=csv.DictWriter(o,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows);write(p,o.getvalue())

protected=json.loads((B/'unrelated_hashes.json').read_text())
for q in [P/'Manuscript.docx',P/'Supplementary_Tables.xlsx']:
    if q.exists():protected[str(q.relative_to(P))]=sha(q)
assert all(sha(P/q)==v for q,v in protected.items())
geom=json.loads((W/'panel_placement.json').read_text())
pg=fitz.open(W/'Figure1.pdf')[0]
assert abs(pg.rect.width*25.4/72-183)<.001 and abs(pg.rect.height*25.4/72-170)<.001

# Preserve the original native BioRender illustration as the reproducible asset input.
assets=W/'source_assets';assets.mkdir(exist_ok=True)
native=assets/'Figure1a_BioRender_native.svg'
if not native.exists():shutil.copy2(P/'Figures/独立示意图/Figure1a.svg',native)
code=W/'build_figure1.py'
text=code.read_text().replace("native=E.parse(str(P/'Figures/独立示意图/Figure1a.svg'))","native=E.parse(str(O/'source_assets/Figure1a_BioRender_native.svg'))")
code.write_text(text)

copy(W/'Figure1.pdf',P/'Figures/Figure1.pdf')
copy(W/'Figure1.svg',P/'图件编辑资料/editable_svg/Figure1.svg')
for root in (C,S):
    for name in 'abcdefg':
        for ext in ('pdf','svg','png'):copy(W/f'panels/Figure1{name}.{ext}',root/f'Figure1/Figure1{name}.{ext}')
for ext in ('pdf','svg','png'):copy(W/f'panels/Figure1a.{ext}',P/f'Figures/独立示意图/Figure1a.{ext}')
copy(W/'Figure1.pdf',C/'review/Figure1_整图审阅.pdf')
copy(W/'review/Figure1.png',C/'review/Figure1_整图审阅.png')
copy(W/'review/Figure1.png',P/'图件编辑资料/BioRender_当前素材/Figure1_合入预览.png')
merge([C/f'review/Figure{k}_整图审阅.pdf' for k in range(1,6)],C/'主图1–5_整图审阅预览.pdf')
merge([C/'主图1–5_整图审阅预览.pdf',C/'ExtendedData2–5_整图审阅预览.pdf'],C/'主图及ED_整图审阅预览.pdf')
formal=[P/f'Figures/Figure{k}.pdf'for k in range(1,6)]+[P/f'Extended_Data/ExtendedData_Figure{k}.pdf'for k in range(1,11)]
merge(formal,P/'全部图件_15页审阅.pdf')
merge([P/'Figures/独立示意图'/f'{n}.pdf'for n in ['Central_illustration','Figure1a','Figure2a']],P/'BioRender三图_审核.pdf')
combined=P/'Nature_Submission/Manuscript_with_15_Figures.pdf';backup(combined)
with fitz.open(B/combined.relative_to(P)) as old:
    body=len(old)-15;d=fitz.open();d.insert_pdf(old,from_page=0,to_page=body-1)
    for f in formal:
        with fitz.open(f)as s:d.insert_pdf(s)
    save(d,combined)

f=P/'图件采用版本.tsv';fields,rr=tsvread(f)
for r in rr:
    if r['figure']=='Figure1':
        r.update(version='Figure1 native-size label refinement; independent BioRender panel letters 20261115',status='AUTHOR_REVIEW_PENDING',width_mm=183,height_mm=170,source_pdf=str(P/'Figures/Figure1.pdf'),source_svg=str(P/'图件编辑资料/editable_svg/Figure1.svg'),provenance_source_pdf=str(W/'Figure1.pdf'),provenance_source_svg=str(W/'Figure1.svg'))
        for k in ('sha256','package_pdf_sha256'):r[k]=sha(P/'Figures/Figure1.pdf')
        for k in ('svg_sha256','package_svg_sha256'):r[k]=sha(P/'图件编辑资料/editable_svg/Figure1.svg')
tsvwrite(f,fields,rr)

for root,filename,key in [(C,'逐小图_来源与修改.tsv','panel'),(S,'逐panel来源与修改清单.tsv','panel_name')]:
    f=root/filename;fields,rr=tsvread(f)
    for r in rr:
        n=r.get(key,'')
        if n in ['Figure1'+x for x in 'abcdefg']:
            g=geom[n[-1]]
            r.update(operation='native_size_layout_from_preserved_values',method='native_size_layout_from_preserved_values',width_mm=g['width_mm'],height_mm=g['height_mm'],font_p10_pt=6,change='Current layout and palettes retained; b/e/g labels enlarged; lowercase panel letters removed from standalone files and added as independent BioRender text objects; fixed values and coordinates unchanged',pdf_sha256=sha(root/f'Figure1/{n}.pdf'),svg_sha256=sha(root/f'Figure1/{n}.svg'),data_integrity='Stored numerical inputs retained; no scientific recomputation',layout_legibility='Native-size text and render reviewed',visual_review='Agent render review completed; author review pending')
            if n=='Figure1a':r['source_files']=str(native)
    tsvwrite(f,fields,rr)
f=C/'review/preview_placements.tsv';fields,rr=tsvread(f);rr=[r for r in rr if r['figure']!='Figure1']
for n,g in geom.items():rr.append(dict(figure='Figure1',panel='Figure1'+n,source_pdf=str(C/f'Figure1/Figure1{n}.pdf'),scale_x=1,scale_y=1,x_mm=g['origin_mm'][0],y_mm=g['origin_mm'][1],width_mm=g['width_mm'],height_mm=g['height_mm']))
tsvwrite(f,fields,rr)

# Figure1-only display map; never overwrite shared identity or palette files.
f=P/'图件编辑资料/Source_Data_INDEX_current.tsv';fields,rr=tsvread(f)
for r in rr:
    if r['panel'].replace(' ','') in ['Figure1'+x for x in 'abcdefg']:
        refs=[s.strip() for s in r.get('supporting_source','').split(';') if s.strip()]
        for ref in [str(W/'build_figure1.py'),str(W/'L1_display_palette.tsv')]:
            if ref not in refs:refs.append(ref)
        r['supporting_source']='; '.join(refs)
tsvwrite(f,fields,rr)

buf=io.BytesIO()
with zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED)as z:
    for f in formal[:5]+[P/'Figures/Central_illustration.pdf']:
        info=zipfile.ZipInfo('Figures/'+f.name,(2000,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;z.writestr(info,f.read_bytes())
write(P/'Figures.zip',buf.getvalue())

note=f'''# Figure1 当前排版审阅

BioRender整图：{URL}

- 当前Figure1为183×170 mm。保留现有布局和配色，本轮仅放大b细胞名称、e视窗标题和g的marker标签，修正f的边缘标签留白。
- BioRender中有7个可独立移动和等比例缩放的900 dpi图块，以及7个独立的Arial粗体8 pt面板字母。小图PDF/SVG/PNG不含a–g；e内的A/B/C是视窗对应标记，保留。用户可直接单独移动面板字母。图内标签和数据元素需在本地SVG/绘图代码修改后替换相应图块。
- SVG导入试验出现色条图像缺失，未采用。正文图的本地PDF/SVG保留矢量文字和图形；UMAP点为900 dpi图层。
- 原细胞坐标、标签、窗口、计数、样本组成、14×28 marker统计和当前色表不变。面板位置、宽高与字母锚点见panel_placement.json。
- 本轮无科学分析或注释验证。其余14幅图、正文、表格保持。是否最终采用仍待作者审核。
- BioRender实际导出预览：review/Figure1_BioRender_export.pdf和PNG。配套原始图标来源、显示色表、标签锚点和物理尺寸保留在本目录；本轮轻量备份位于图件编辑资料/轻量备份/Figure1标签独立化前_20261115。
'''
write(W/'README_审阅说明.md',note)
f=P/'BioRender图件更新说明.md';old=f.read_text()
if not old.startswith('> 当前Figure1已进一步更新为183×170 mm'):
    write(f,'> 当前Figure1已进一步更新为183×170 mm紧凑审阅稿，整图可在BioRender中按小图编辑。见图件编辑资料/BioRender_当前素材/Figure1整图排版/README_审阅说明.md。以下为此前合入记录；其中Figure1旧高度已被本轮替代。\n\n'+old)
else:changed.add(f)
f=P/'README_当前状态.txt';old=f.read_text()
if '2026-11-15 Figure1紧凑排版：' not in old:
    write(f,old+'\n2026-11-15 Figure1紧凑排版：Figures/Figure1.pdf与15页图集已同步，183×170 mm；BioRender整图链接见图件编辑资料/BioRender_当前素材/Figure1整图排版/README_审阅说明.md。作者审核中。\n')
else:changed.add(f)

checks={'scientific_analysis_run':False,'biorender_url':URL,'width_mm':183,'height_mm':170,'status':'AUTHOR_REVIEW_PENDING','protected_files_unchanged':all(sha(P/q)==v for q,v in protected.items()),'panels':geom}
assert checks['protected_files_unchanged']
with fitz.open(P/'全部图件_15页审阅.pdf')as m:
    assert len(m)==15
    checks['merged_render_matches_masters']=[]
    for i,src in enumerate(formal):
        with fitz.open(src)as d:
            a=m[i].get_pixmap(matrix=fitz.Matrix(1,1));b=d[0].get_pixmap(matrix=fitz.Matrix(1,1))
            assert (a.width,a.height,a.n)==(b.width,b.height,b.n)
            av=np.frombuffer(a.samples,dtype=np.uint8).reshape(a.height,a.width,a.n)
            bv=np.frombuffer(b.samples,dtype=np.uint8).reshape(b.height,b.width,b.n)
            diff=np.abs(av.astype(int)-bv.astype(int))
            exact=bool(np.max(diff)==0)
            interior_exact=bool(np.max(diff[2:-2,2:-2])==0)
            # PDF rewriting can change page-edge antialiasing by one 8-bit level.
            # This is recorded, not described as byte-identical rendering.
            ok=exact or (interior_exact and np.max(diff)<=1)
            assert ok and m[i].get_text()==d[0].get_text(),(i,src)
            checks['merged_render_matches_masters'].append({'figure':src.name,'exact_pixel_match':exact,'interior_pixel_match':interior_exact,'max_channel_difference':int(diff.max()),'page_geometry_equal':m[i].rect==d[0].rect,'text_equal':True,'note':'' if exact else 'Only page-edge antialiasing differs by 1/255; content pixels identical'})
write(W/'installation_checks.json',json.dumps(checks,ensure_ascii=False,indent=2))
for filename in ['SUBMISSION_FILE_LIST.tsv','文件来源与校验.tsv']:
    f=P/filename;fields,rr=tsvread(f)
    for r in rr:
        target=P/r['path']
        if target in changed:
            r.update(bytes=target.stat().st_size,sha256=sha(target))
            if 'status'in r:r['status']='CURRENT_AUTHOR_REVIEW_PENDING'
    tsvwrite(f,fields,rr)
f=P/'PACKAGE_SHA256SUMS.txt';lines=[]
for line in f.read_text().splitlines():
    a=line.split('  ',1)
    if len(a)==2 and P/a[1]in changed:line=sha(P/a[1])+'  '+a[1]
    lines.append(line)
write(f,'\n'.join(lines)+'\n')
print(json.dumps({'changed_files':len(changed),'protected_files':len(protected),'merged_pages':15,'native_size_mm':[183,170],'biorender_url':URL},ensure_ascii=False))
