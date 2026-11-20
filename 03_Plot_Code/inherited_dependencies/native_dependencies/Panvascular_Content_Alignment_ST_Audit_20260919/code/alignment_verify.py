from pathlib import Path
import xml.etree.ElementTree as E,json,hashlib,collections,fitz
B=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');O=Path('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919');N={'s':'http://www.w3.org/2000/svg'};out=[]
for n in ['Figure1','Figure2','Figure4']+[f'ExtendedData_Figure{i}' for i in range(1,11)]:
 a=E.parse(B/'editable_svg'/f'{n}.svg');b=E.parse(O/'editable_svg'/f'{n}.svg')
 for tag,attr in [('path','d'),('image','{http://www.w3.org/1999/xlink}href')]:
  assert collections.Counter(e.get(attr) for e in a.findall('.//s:'+tag,N))==collections.Counter(e.get(attr) for e in b.findall('.//s:'+tag,N)),(n,tag)
 assert collections.Counter(e.text for e in a.findall('.//s:text',N))==collections.Counter(e.text for e in b.findall('.//s:text',N)),n
 assert a.getroot().get('viewBox')==b.getroot().get('viewBox')
 for e in b.findall('.//s:g',N):
  assert 'scale(' not in e.get('transform','') or e.get('id','').startswith('image'),(n,e.get('id'),e.get('transform'))
 out.append(dict(figure=n,paths='EXACT_UNCHANGED',embedded_images='EXACT_UNCHANGED',text_content='EXACT_UNCHANGED',page='UNCHANGED',operation='translation only; legend anchor normalized to start'))
(O/'audit/alignment_integrity.json').write_text(json.dumps(out,indent=2));print('PASS',len(out),'figures')
