from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
import fitz, shutil, hashlib, json, csv, io, datetime, copy

P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');R=Path(__file__).resolve().parent;B=R/'backup_before_adoption'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit=json.loads((R/'review/adoption_audit.json').read_text());changed=[Path(p) for p in audit['changed_files']]
def backup(p):
    q=B/p.relative_to(P)
    if p.exists() and not q.exists():q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
def install(s,p):
    backup(p);t=p.with_name(p.name+'.figure3_tmp');shutil.copy2(s,t);t.replace(p);changed.append(p)
def write(p,t):
    backup(p);p.write_text(t);changed.append(p)

for name in ['Figure_Legends','Manuscript']:
    d=fitz.open(R/'review'/f'{name}.pdf');text='\n'.join(p.get_text() for p in d)
    assert 'Paired stacked bars show' in text
    assert 'Circles/squares identify' not in text
    install(R/'review'/f'{name}.pdf',P/'Nature_Submission'/f'{name}.pdf')

# Replace only page 3 in the 15-page figure review, retaining all other page content.
review=P/'全部图件_15页审阅.pdf';backup(review)
old=fitz.open(review);new=fitz.open();f3=fitz.open(P/'Figures/Figure3.pdf')
assert len(old)==15
for i in range(15):new.insert_pdf(f3 if i==2 else old,from_page=0 if i==2 else i,to_page=0 if i==2 else i)
for i in range(15):
    if i!=2:assert hashlib.sha256(old[i].get_pixmap().samples).digest()==hashlib.sha256(new[i].get_pixmap().samples).digest()
tmp=R/'review/全部图件_15页审阅.pdf';new.save(tmp);install(tmp,review)
combined=fitz.open();combined.insert_pdf(fitz.open(P/'Nature_Submission/Manuscript.pdf'));combined.insert_pdf(new)
tmp=R/'review/Manuscript_with_15_Figures.pdf';combined.save(tmp);install(tmp,P/'Nature_Submission/Manuscript_with_15_Figures.pdf')
audit['other_14_review_pages_render_identical']=True
audit['combined_review_pages']=len(combined)

# Synchronize the existing shared-row d/e review with the adopted pair.
pair=fitz.open();page=pair.new_page(width=183*72/25.4,height=105*72/25.4)
for k,x,w in [('d',2,114),('e',117,65)]:page.show_pdf_page(fitz.Rect(x*72/25.4,0,(x+w)*72/25.4,105*72/25.4),fitz.open(R/'panels'/f'Figure3{k}.pdf'),0)
pair.save(R/'review/Figure3de_shared_rows.pdf');page.get_pixmap(dpi=300).save(R/'review/Figure3de_shared_rows.png')
ns='http://www.w3.org/2000/svg';root=E.Element('{'+ns+'}svg',nsmap={None:ns},width=f'{183*72/25.4}pt',height=f'{105*72/25.4}pt',viewBox=f'0 0 {183*72/25.4} {105*72/25.4}')
whole=E.parse(str(R/'review/Figure3.svg'))
for k in 'de':
    el=copy.deepcopy(whole.xpath('//*[@id=$v]',v='adopted_panel_Figure3'+k)[0]);el.set('y','0');root.append(el)
E.ElementTree(root).write(str(R/'review/Figure3de_shared_rows.svg'),encoding='utf-8',xml_declaration=True)
for ext in ['pdf','png','svg']:install(R/'review'/f'Figure3de_shared_rows.{ext}',P/'图件编辑资料/可视化库重绘_当前/review'/f'Figure3de_shared_rows.{ext}')

# Replace only the Figure3 entry inside the existing submission figure archive.
zp=P/'Figures.zip';backup(zp)
with ZipFile(zp) as z:infos=z.infolist();contents={i.filename:z.read(i.filename) for i in infos}
old_contents=contents.copy();contents['Figures/Figure3.pdf']=(P/'Figures/Figure3.pdf').read_bytes()
tmp=zp.with_name('Figures.figure3_tmp.zip')
with ZipFile(tmp,'w') as z:
    for info in infos:z.writestr(info,contents[info.filename])
with ZipFile(tmp) as z:
    assert z.testzip() is None
    for n,b in old_contents.items():
        if n!='Figures/Figure3.pdf':assert z.read(n)==b
tmp.replace(zp);changed.append(zp);audit['other_figure_zip_entries_byte_identical']=True

url='https://app.biorender.com/illustrations/6af9a1433d7c6083ec2f747f'
status={'url':url,'title':'Panvascular Figure 3 | Approved d-f','status':'SAVED_RELOADED_AND_EXPORTED','page_mm':[183,247],
        'independently_arrangeable_panels':6,'biorender_import_formats':{'a':'PNG 600 dpi from original formal PDF','b':'PNG 600 dpi from original formal PDF','c':'PNG 600 dpi from original formal PDF','d':'SVG','e':'SVG','f':'SVG'},
        'native_export_png':str(R/'review/Figure3_BioRender.png'),'native_export_pdf':str(R/'review/Figure3_BioRender.pdf'),
        'local_vector_pdf':str(P/'Figures/Figure3.pdf'),'render_checked':True,'biorender_panel_internals_are_not_native_editable_text':True}
write(R/'review/biorender_status.json',json.dumps(status,ensure_ascii=False,indent=2)+'\n')
coord=P/'图件编辑资料/BioRender_当前素材/browser_coordination.json';obj=json.loads(coord.read_text());obj.update(figure3_using_chrome=False,figure3_browser_done=True,figure3_url=url,active_browser_task=None,chrome_available_for_figure1=True,chrome_available_for_figure2=True,chrome_available_for_figure4=True,updated_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat());write(coord,json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
readme=f'''# Figure 3 已采用版 · 2026-11-15

BioRender：{url}

用户已选择 d 森林图、e 共享行证据符号矩阵、f 双情境后验构成图。六个面板分别导入，均可单独移动、缩放、排版；内部文字并非 BioRender 原生文本。a–c 来自原正式 PDF 的 600 dpi 渲染，d–f 为 SVG 导入。已保存、刷新复核并完成 PNG/PDF 原生导出。

正式矢量 PDF：`Figures/Figure3.pdf`；正式 SVG：`图件编辑资料/editable_svg/Figure3.svg`。正式复合矢量文件保留原 a–c 对象、坐标和比例；未用不同版本的独立 a–c 重绘图替换。独立小图目录内的 a–c 文件没有改动。

d 保留原 29 个效应及 CI；e 保留原 145 个证据类别；f 保留原 10 条情境记录、H0–H4、共享 SNP 数和分类。没有重跑科学分析。采用数据由“泛血管预处理和L1”任务确认，来源、SHA256 和交接说明见 `source_data/`。这是当前本地采用链的确认，不是服务器全盘最新数据审计。

d–f 多余标题、脚注已去除，仅保留面板字母、基因/细胞标签、坐标轴和必要图例。图注相应同步至 Manuscript.docx 和 Figure_Legends.docx/PDF；a–c 图注及其他段落不变。15 页审阅图、稿件合并审阅 PDF、Figures.zip 和版本/校验索引已同步。

`review/Figure3_BioRender.png` 是 BioRender 实际导出预览；`review/Figure3_BioRender.pdf` 是实际导出存档。BioRender PDF 导出可带极细面板边缘线；正式矢量 PDF 无该边缘，正式投稿优先使用本地矢量文件。BioRender 导出尚未申请 publication license；当前为编辑/审阅画布。

原文件见 `backup_before_adoption/`。重建采用版请用 `build_figure3.py`；不要直接运行旧 `redraw_main.py` 覆盖此版。`adopt_figure3.py` 和 `finalize_delivery.py` 是本轮采用记录，不应作为批量覆盖其他任务改动的常规脚本。尚未向期刊提交，也未发布公共模板。
'''
write(R/'README.md',readme)

# Update only changed records in submission and checksum manifests; retain unrelated entries.
changed_map={str(p.relative_to(P)):p for p in changed}
def update_tsv(p,pathkey,add=False):
    rows=list(csv.DictReader(io.StringIO(p.read_text()),delimiter='\t'));keys=list(rows[0]);seen=set()
    for row in rows:
        key=row[pathkey];seen.add(key)
        if key in changed_map:
            f=changed_map[key];row['bytes']=str(f.stat().st_size);row['sha256']=sha(f)
    if add:
        for key,f in changed_map.items():
            if key not in seen:rows.append({pathkey:key,'bytes':str(f.stat().st_size),'sha256':sha(f)})
    out=io.StringIO();w=csv.DictWriter(out,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows);write(p,out.getvalue());changed_map[str(p.relative_to(P))]=p
update_tsv(P/'SUBMISSION_FILE_LIST.tsv','path')
sumfile=P/'PACKAGE_SHA256SUMS.txt';lines=[]
for line in sumfile.read_text().splitlines():
    parts=line.split('  ',1)
    if len(parts)==2 and parts[1] in changed_map:line=sha(changed_map[parts[1]])+'  '+parts[1]
    lines.append(line)
write(sumfile,'\n'.join(lines)+'\n');changed_map[str(sumfile.relative_to(P))]=sumfile

audit['changed_files']=[str(p) for p in dict.fromkeys(changed)]
audit['biorender']=status
audit['local_replacements_completed']=True
for p,h in audit['protected_files_unchanged'].items():assert sha(p)==h
for item in audit['source_data_hashes_unchanged'].values():assert sha(item['path'])==item['sha256']
(R/'review/adoption_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
changed_map[str((R/'review/adoption_audit.json').relative_to(P))]=R/'review/adoption_audit.json'
update_tsv(P/'文件来源与校验.tsv','path',add=True)
print(json.dumps({'status':'COMPLETE','modified_paths':len(changed_map),'review_pages':len(combined),'url':url},ensure_ascii=False))
