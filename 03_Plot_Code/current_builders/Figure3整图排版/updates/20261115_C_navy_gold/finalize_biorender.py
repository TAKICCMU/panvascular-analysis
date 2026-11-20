"""Run only after actual BioRender reload and downloaded exports are inspected."""
from pathlib import Path
import argparse,json,hashlib,shutil,datetime,csv,io
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');A=P/'图件编辑资料/BioRender_当前素材/Figure3整图排版';R=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--png',required=True);parser.add_argument('--pdf');parser.add_argument('--verified-reload',action='store_true');args=parser.parse_args();assert args.verified_reload
assert Path(args.png).is_file()
audit=json.loads((R/'adoption_audit.json').read_text());changed=[Path(p) for p in audit['modified_files']]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def backup(p):
 q=R/'backup_before_C'/p.relative_to(P)
 if p.exists() and not q.exists():q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
def write(p,v):backup(p);p.write_text(v);changed.append(p)
for ext,src in [('png',args.png),('pdf',args.pdf)]:
 if src:
  dest=A/'review'/f'Figure3_BioRender.{ext}';backup(dest);shutil.copy2(src,dest);changed.append(dest)
  shutil.copy2(src,R/'stage/review'/f'Figure3_BioRender_C.{ext}')
status=json.loads((A/'review/biorender_status.json').read_text());status.update(status='C_PALETTE_SAVED_RELOADED_AND_EXPORTED',current_local_version='20261115_C_navy_gold',render_checked=True,independently_arrangeable_panels=6,updated_at=datetime.datetime.now().astimezone().isoformat(),latest_update_audit=str(R/'adoption_audit.json'))
status['biorender_import_formats']['e']='SVG, C navy/pale-gold tiles'
status['biorender_import_formats']['f']='PNG 600 dpi rendered from exact vector PDF; SVG formula rendering was incompatible'
if not args.pdf:status['native_export_pdf']='Previous export only; not updated in this revision'
write(A/'review/biorender_status.json',json.dumps(status,ensure_ascii=False,indent=2)+'\n')
meta=json.loads((A/'panel_placement.json').read_text());meta['status']='C_PALETTE_SAVED_RELOADED_AND_EXPORTED';write(A/'panel_placement.json',json.dumps(meta,ensure_ascii=False,indent=2)+'\n')
readme=(A/'README.md').read_text().replace('当前本地已采用，等待画布更新；旧 BioRender 导出图已备份，不能代表本轮 C 版。','C 版已同步原画布，刷新复核并导出。六个小图可分别移动、缩放；小图内部文字不是 BioRender 原生文本。旧版导出已备份。3e 为 SVG；3f 使用原矢量 PDF 的 600 dpi 渲染，避免 BioRender SVG 公式排版错误。本地 PDF/SVG 仍为矢量源文件。')
write(A/'README.md',readme)
def table(path,fn):
 rows=list(csv.DictReader(io.StringIO(path.read_text()),delimiter='\t'));keys=list(rows[0])
 for row in rows:fn(row)
 o=io.StringIO();w=csv.DictWriter(o,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows);write(path,o.getvalue())
def fig(row):
 if row['figure']=='Figure3':row['status']='C_PALETTE_ADOPTED_BIORENDER_SAVED'
table(P/'图件采用版本.tsv',fig)
def pan(row):
 if row['figure']=='Figure3' and row['panel'] in 'ef':row['visual_review']='C selected by user; local data and render verified; BioRender saved, reloaded and exported'
table(P/'独立小图_最新版/逐panel来源与修改清单.tsv',pan)
audit.update(status='COMPLETE_C_ADOPTED_BIORENDER_SAVED',biorender_updated=True,biorender=status,modified_files=[str(p) for p in dict.fromkeys(changed)])
(R/'adoption_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
coord=P/'图件编辑资料/BioRender_当前素材/browser_coordination.json';obj=json.loads(coord.read_text());obj.update(figure3_using_chrome=False,figure3_browser_done=True,chrome_available_for_figure3=True,chrome_available_for_figure4=True,active_browser_task=None,updated_at=datetime.datetime.now().astimezone().isoformat());coord.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'complete':True,'url':status['url'],'png_sha256':sha(A/'review/Figure3_BioRender.png')}))
