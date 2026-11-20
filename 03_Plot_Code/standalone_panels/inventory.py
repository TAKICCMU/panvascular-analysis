from pathlib import Path
from lxml import etree as E
import json,hashlib,fitz,pandas as pd
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'独立小图_最新版';ns={'s':'http://www.w3.org/2000/svg'}
rows=[]; texts=[]
for row in pd.read_csv(P/'图件采用版本.tsv',sep='\t').to_dict('records'):
 n=row['figure'];p=Path(row['package_pdf']);s=Path(row['package_svg']);r=E.parse(str(s));page=fitz.open(p)[0];page.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(O/'review'/f'baseline_{n}.png')
 item={'figure':n,'pdf':str(p),'svg':str(s),'hash_match':hashlib.sha256(p.read_bytes()).hexdigest()==row['package_pdf_sha256'],'svg_hash_match':hashlib.sha256(s.read_bytes()).hexdigest()==row['package_svg_sha256'],'width_pt':page.rect.width,'height_pt':page.rect.height};rows.append(item)
 texts.append('\n## '+n)
 for el in r.xpath('//*[@id]'):
  if el.getparent().tag.endswith('svg') or el.getparent().get('id')=='figure_1':
   text=' | '.join(''.join(t.itertext()).strip() for t in el.xpath('.//s:text',namespaces=ns))
   if el.get('id') not in ['figure_1','patch_1']:
    texts.append(el.get('id')+' :: '+text[:850])
(O/'review/BASELINE_MANIFEST.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2));(O/'review/SVG_COMPONENTS.txt').write_text('\n'.join(texts)); print('figures',len(rows),'hash matches',sum(r['hash_match'] and r['svg_hash_match'] for r in rows))
