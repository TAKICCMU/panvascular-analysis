from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
import hashlib,json,csv,io
import pandas as pd
import fitz
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');A=P/'图件编辑资料/BioRender_当前素材/Figure3整图排版';R=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit=json.loads((R/'adoption_audit.json').read_text());changed={str(Path(p).relative_to(P)):Path(p) for p in audit['modified_files']}
for v in audit['source_hashes'].values():assert sha(v['path'])==v['sha256']
for p,h in audit['protected_panel_hashes'].items():assert sha(p)==h,p
T=pd.read_csv(audit['source_hashes']['e']['path'],sep='\t');M=T.pivot(index='row',columns='col',values='source_value').to_numpy()
rt=E.parse(str(A/'panels/Figure3e.svg'));ids=rt.xpath('//*[@id]/@id')
for i in range(29):
 for j in [2,3,4]:
  v=int(M[i,j]);gid=f'evidence_row{i}_col{j}_code{v}';assert gid in ids
  group=rt.xpath('//*[@id=$gid]',gid=gid)[0];styles=' '.join(x.get('style','') for x in group.iter()).lower();assert {1:'#f1f3f5',2:'#d5b47b',3:'#45698f'}[v] in styles
assert sum(x.startswith('evidence_row') for x in ids)==87
assert (M[:,:2]==3).all()
whole=E.parse(str(P/'图件编辑资料/editable_svg/Figure3.svg'));old=E.parse(str(R/'backup_before_C/图件编辑资料/editable_svg/Figure3.svg'))
for gid in audit['a_to_d_svg_subtrees_preserved']:assert E.tostring(whole.xpath('//*[@id=$gid]',gid=gid)[0])==E.tostring(old.xpath('//*[@id=$gid]',gid=gid)[0])
for p in [P/'Manuscript.docx',P/'Nature_Submission/Figure_Legends.docx']:
 with ZipFile(p) as z:root=E.fromstring(z.read('word/document.xml'))
 txt=' '.join(root.xpath('//*[local-name()="t"]/text()'));assert 'Blue indicates' in txt and 'pale gold' in txt and 'horizontal stems' in txt
for k in 'ef':
 for ext in ['svg','pdf','png']:
  assert sha(A/'panels'/f'Figure3{k}.{ext}')==sha(P/'独立小图_最新版/Figure3'/f'Figure3{k}.{ext}')
assert sha(P/'Figures/Figure3.pdf')==sha(A/'review/Figure3.pdf')
assert sha(P/'图件编辑资料/editable_svg/Figure3.svg')==sha(P/'Public_Release/Figure_Code/frozen_svg/Figure3.svg')
with ZipFile(P/'Figures.zip') as z:assert z.read('Figures/Figure3.pdf')==(P/'Figures/Figure3.pdf').read_bytes()
txt=fitz.open(P/'Figures/Figure3.pdf')[0].get_text();assert 'Classification' not in txt and 'Shared SNP' not in txt
for k in 'ef':
 d=fitz.open(A/'panels'/f'Figure3{k}.pdf');d[0].get_pixmap(dpi=600).save(R/'stage/upload_ready'/f'Figure3{k}_C_600dpi.png')
audit['C_colour_and_all_87_categories_verified']=True
audit['formal_figure_and_archive_consistent']=True
audit['a_to_d_svg_byte_identical']=True
(R/'adoption_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
changed[str((R/'adoption_audit.json').relative_to(P))]=R/'adoption_audit.json'
def update_tsv(p,key,add=False):
 rows=list(csv.DictReader(io.StringIO(p.read_text()),delimiter='\t'));keys=list(rows[0]);seen=set()
 for row in rows:
  k=row[key];seen.add(k)
  if k in changed:row['bytes']=str(changed[k].stat().st_size);row['sha256']=sha(changed[k])
 if add:
  for k,v in changed.items():
   if k not in seen:rows.append({key:k,'bytes':str(v.stat().st_size),'sha256':sha(v)})
 o=io.StringIO();w=csv.DictWriter(o,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 backup=R/'backup_before_C'/p.relative_to(P)
 if not backup.exists():backup.parent.mkdir(parents=True,exist_ok=True);backup.write_bytes(p.read_bytes())
 p.write_text(o.getvalue());changed[str(p.relative_to(P))]=p
update_tsv(P/'SUBMISSION_FILE_LIST.tsv','path')
p=P/'PACKAGE_SHA256SUMS.txt';backup=R/'backup_before_C'/p.relative_to(P)
if not backup.exists():backup.write_bytes(p.read_bytes())
lines=[]
for line in p.read_text().splitlines():
 parts=line.split('  ',1)
 if len(parts)==2 and parts[1] in changed:line=sha(changed[parts[1]])+'  '+parts[1]
 lines.append(line)
p.write_text('\n'.join(lines)+'\n');changed[str(p.relative_to(P))]=p
update_tsv(P/'文件来源与校验.tsv','path',True)
print(json.dumps({'checks':'passed','C_category_cells':87,'protected_panels':'a-d','indices_updated':3,'biorender':audit['status']},ensure_ascii=False))
