from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
import hashlib,json,shutil,csv,io
import fitz
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
A=P/'图件编辑资料/BioRender_当前素材/Figure3整图排版'
R=Path(__file__).resolve().parent;S=R/'stage';B=R/'backup_before_C'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit=json.loads((R/'adoption_audit.json').read_text());changed=[]
assert audit['status']=='C_APPROVED_STAGED'
for p,h in audit['baseline_formal'].items():assert sha(p)==h,p
for v in audit['source_hashes'].values():assert sha(v['path'])==v['sha256']
protected={str(p):sha(p) for p in (P/'Figures').glob('*.pdf') if p.stem!='Figure3'}
def backup(p):
 q=B/p.relative_to(P)
 if p.exists() and not q.exists():q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
def install(src,dst):
 backup(dst);dst.parent.mkdir(parents=True,exist_ok=True);tmp=dst.with_name(dst.name+'.f3c_tmp');shutil.copy2(src,tmp);tmp.replace(dst);changed.append(dst)
def write(p,text):
 backup(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.f3c_tmp');tmp.write_text(text);tmp.replace(p);changed.append(p)
def table(p,fn):
 rows=list(csv.DictReader(io.StringIO(p.read_text()),delimiter='\t'));keys=list(rows[0])
 for row in rows:fn(row)
 o=io.StringIO();w=csv.DictWriter(o,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows);write(p,o.getvalue())
for src,dst in [('review/Figure3.pdf',P/'Figures/Figure3.pdf'),('review/Figure3.svg',P/'图件编辑资料/editable_svg/Figure3.svg'),('review/Figure3.svg',P/'Public_Release/Figure_Code/frozen_svg/Figure3.svg'),('review/Figure3.pdf',P/'图件编辑资料/可视化库重绘_当前/review/Figure3_整图审阅.pdf'),('review/Figure3.png',P/'图件编辑资料/可视化库重绘_当前/review/Figure3_整图审阅.png')]:install(S/src,dst)
for ext in ['pdf','svg','png']:
 install(S/'review'/f'Figure3.{ext}',A/'review'/f'Figure3.{ext}')
 for dstbase in [A/'review',P/'图件编辑资料/可视化库重绘_当前/review']:install(S/'review'/f'Figure3de_shared_rows.{ext}',dstbase/f'Figure3de_shared_rows.{ext}')
for k in 'ef':
 for base in [A/'panels',P/'独立小图_最新版/Figure3',P/'图件编辑资料/可视化库重绘_当前/Figure3']:
  for ext in ['pdf','svg','png']:install(S/'panels'/f'Figure3{k}.{ext}',base/f'Figure3{k}.{ext}')
 install(S/'panels'/f'Figure3{k}.svg',A/'upload_ready'/f'Figure3{k}.svg')
for name,rel in [('Manuscript','Manuscript.docx'),('Figure_Legends','Nature_Submission/Figure_Legends.docx')]:
 install(S/'documents'/f'{name}.docx',P/rel)
 doc=fitz.open(S/'review'/f'{name}.pdf');txt='\n'.join(p.get_text() for p in doc)
 assert 'Blue indicates' in txt and 'horizontal stems' in txt
 install(S/'review'/f'{name}.pdf',P/'Nature_Submission'/f'{name}.pdf')
 install(S/'review'/f'{name}.pdf',A/'review'/f'{name}.pdf')
 install(S/'review'/f'{name}_Figure3_caption.md',A/'review'/('Figure3_caption_manuscript.md' if name=='Manuscript' else 'Figure3_caption_legends.md'))

old=fitz.open(P/'全部图件_15页审阅.pdf');f3=fitz.open(P/'Figures/Figure3.pdf');new=fitz.open();assert len(old)==15
for i in range(15):new.insert_pdf(f3 if i==2 else old,from_page=0 if i==2 else i,to_page=0 if i==2 else i)
for i in range(15):
 if i!=2:assert old[i].get_pixmap().samples==new[i].get_pixmap().samples
new.save(S/'review/全部图件_15页审阅.pdf');install(S/'review/全部图件_15页审阅.pdf',P/'全部图件_15页审阅.pdf')
combined=fitz.open();combined.insert_pdf(fitz.open(P/'Nature_Submission/Manuscript.pdf'));combined.insert_pdf(new);combined.save(S/'review/Manuscript_with_15_Figures.pdf');install(S/'review/Manuscript_with_15_Figures.pdf',P/'Nature_Submission/Manuscript_with_15_Figures.pdf')
zp=P/'Figures.zip';backup(zp)
with ZipFile(zp) as z:infos=z.infolist();entries={i.filename:z.read(i.filename) for i in infos}
assert 'Figures/Figure3.pdf' in entries
entries['Figures/Figure3.pdf']=(P/'Figures/Figure3.pdf').read_bytes();tmp=R/'Figures.tmp.zip'
with ZipFile(tmp,'w') as z:
 for info in infos:z.writestr(info,entries[info.filename])
with ZipFile(tmp) as z:
 assert z.testzip() is None
 for k,v in entries.items():assert z.read(k)==v
install(tmp,zp);tmp.unlink()

pdf=P/'Figures/Figure3.pdf';svg=P/'图件编辑资料/editable_svg/Figure3.svg'
def figure(row):
 if row['figure']=='Figure3':row.update(source_pdf=str(pdf),source_svg=str(svg),version='Approved C navy/gold 3e tiles and 3f regional GWAS display 20261115; a-d retained',status='LOCAL_ADOPTED_BIORENDER_UPDATE_PENDING',sha256=sha(pdf),svg_sha256=sha(svg),package_pdf=str(pdf),package_svg=str(svg),package_pdf_sha256=sha(pdf),package_svg_sha256=sha(svg),provenance_source_pdf=str(S/'review/Figure3.pdf'),provenance_source_svg=str(S/'review/Figure3.svg'))
table(P/'图件采用版本.tsv',figure)
changes={'e':'Approved C navy/pale-gold categorical tiles: 87 variable cells shown; 58 common primary cells explained in legend; categories and row order unchanged','f':'Paired posterior compositions plus -log10 stored minimum regional GWAS P for same 10 records; classification and shared SNP text removed; original classifications retained in source'}
def panel(row):
 if row['figure']=='Figure3' and row['panel'] in 'ef':
  k=row['panel'];base=P/'独立小图_最新版/Figure3'
  row.update(method='approved_display_revision_existing_data',source_files=audit['source_hashes'][k]['path'],source_objects='145 frozen categories; 87 displayed and 58 constant cells documented' if k=='e' else '10 original gene-context-scenario records',change=changes[k],pdf_sha256=sha(base/f'Figure3{k}.pdf'),svg_sha256=sha(base/f'Figure3{k}.svg'),data_integrity='Source SHA256 verified; original categories and values preserved; no analysis rerun',layout_legibility='Shared d/e row alignment; essential axes and legends only',visual_review='C selected by user; local render checked; BioRender update pending')
table(P/'独立小图_最新版/逐panel来源与修改清单.tsv',panel)
def current(row):
 if row['panel'] in ['Figure3e','Figure3f']:row.update(operation='approved_display_revision_existing_data',source_files=audit['source_hashes'][row['panel'][-1]]['path'],change=changes[row['panel'][-1]])
table(P/'图件编辑资料/可视化库重绘_当前/逐小图_来源与修改.tsv',current)
meta=json.loads((A/'panel_placement.json').read_text());meta.update(status='LOCAL_ADOPTED_BIORENDER_UPDATE_PENDING',latest_update_audit=str(R/'adoption_audit.json'))
for row in meta['panels']:
 if row['panel'] in ['Figure3e','Figure3f']:row['source']=changes[row['panel'][-1]]
write(A/'panel_placement.json',json.dumps(meta,ensure_ascii=False,indent=2)+'\n')
status=json.loads((A/'review/biorender_status.json').read_text());status['status']='PREVIOUS_CANVAS_SAVED_CURRENT_EF_UPDATE_PENDING';status['current_local_version']='20261115_C_navy_gold';status['latest_update_audit']=str(R/'adoption_audit.json');write(A/'review/biorender_status.json',json.dumps(status,ensure_ascii=False,indent=2)+'\n')
readme='''# Figure 3 当前采用版 · C 深蓝＋浅金 · 2026-11-15

用户确认采用 3e 配色 C。3e 显示 29 行×3 列证据色块；原 2 列共同条件（各29项均支持）在图注说明。深蓝表示原标准下的支持，浅金加向下三角表示反向空间关联，浅灰表示原标准下无支持。原分类与行序不变。

3f 保留原五个基因—组织背景、两情境的后验构成；右侧改为储存的最小区域 GWAS P 值的 −log10 点图，删除 Classification 与 Shared SNPs 文字列。线段从零出发，不是置信区间。原分类及 SNP 数仍保留在源数据中。

a–d 与源数据均原样保留。已同步本地正式 PDF/SVG、独立小图、图注、审阅 PDF、Figures.zip 及相关索引。没有重跑科学分析。

本轮脚本、备份与核对记录：`updates/20261115_C_navy_gold/`。原 adopted d–f 脚本是历史记录，不可直接重跑覆盖此版。当前源数据与保护检查见该目录 `adoption_audit.json`。

BioRender 画布：https://app.biorender.com/illustrations/6af9a1433d7c6083ec2f747f 。当前本地已采用，等待画布更新；旧 BioRender 导出图已备份，不能代表本轮 C 版。
'''
write(A/'README.md',readme)
for p,h in audit['protected_panel_hashes'].items():assert sha(p)==h,p
for p,h in protected.items():assert sha(p)==h,p
for v in audit['source_hashes'].values():assert sha(v['path'])==v['sha256']
audit.update(status='LOCAL_ADOPTED_BIORENDER_UPDATE_PENDING',other_14_review_pages_render_identical=True,other_figure_pdf_hashes_unchanged=protected,combined_review_pages=len(combined),source_hashes_unchanged=True,modified_files=[str(p) for p in dict.fromkeys(changed)],new_figure3_sha256=sha(pdf),new_figure3_svg_sha256=sha(svg))
(R/'adoption_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'adopted':True,'changed_files':len(set(changed)),'source_records':[29,145,10],'combined_pages':len(combined),'biorender_pending':True},ensure_ascii=False))
