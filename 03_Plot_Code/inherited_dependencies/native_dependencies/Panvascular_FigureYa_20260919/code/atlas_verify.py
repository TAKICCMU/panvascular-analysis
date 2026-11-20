from pathlib import Path
import xml.etree.ElementTree as E,json,csv,hashlib
O=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');B=Path('/Users/taki/Desktop/投稿包_最新版_20260919/图件编辑资料/editable_svg');N={'s':'http://www.w3.org/2000/svg'}
rows=[]
for name in ['Figure1']+[f'ExtendedData_Figure{i}' for i in range(2,6)]:
 a=E.parse(B/(name+'.svg'));b=E.parse(O/'editable_svg'/(name+'.svg'))
 ga={g.get('id'):g for g in a.findall('.//s:g',N)};gb={g.get('id'):g for g in b.findall('.//s:g',N)}
 for key,g in ga.items():
  if not (key or '').startswith('PathCollection_'):continue
  pa=g.findall('s:path',N);pb=gb[key].findall('s:path',N)
  assert [p.get('d') for p in pa]==[p.get('d') for p in pb]
  tiles=gb[key][0].findall('s:rect',N)
  assert len(pa)==len(tiles)
  for p,t in zip(pa,tiles):assert t.get('fill') in p.get('style','')
  rows.append(dict(figure=name,collection=key,glyph_count=len(pa),exact_paths_and_colours='PASS'))
 # Every original text remains unchanged except declared size-legend description.
 ta=[t.text for t in a.findall('.//s:text',N)];tb=[t.text.replace('Circle: detected cells (%)','Detected cells (%)') if t.text else t.text for t in b.findall('.//s:text',N)];assert ta==tb
json.dump(rows,open(O/'review/atlas_exact_geometry_check.json','w'),indent=2)
A=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Full_Layout_Correction_20260918/panels/atlas')
mapping=[('ExtendedData_Figure3','b','PathCollection_1'),('ExtendedData_Figure4','b','PathCollection_1'),('ExtendedData_Figure4','d','PathCollection_2'),('ExtendedData_Figure5','b','PathCollection_1')]
source=[]
for fig,p,col in mapping:
 f=A/f'{fig}_{p}_full_matrix.tsv';n=sum(1 for _ in csv.DictReader(open(f),delimiter='\t'));assert n==next(r['glyph_count'] for r in rows if r['figure']==fig and r['collection']==col)
 source.append(dict(figure=fig,panel=p,source=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),matrix_entries=n,rendered_entries=n))
json.dump(source,open(O/'notes/atlas_source_mapping.json','w'),indent=2)
print('PASS: all original text, bubble geometry, circle area, exact rendered colour and full matrix entry counts retained')
