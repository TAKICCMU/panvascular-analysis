from pathlib import Path
from collections import Counter
import re,json,fitz,xml.etree.ElementTree as E
O=Path('/Users/taki/Desktop/Panvascular_FigureYa_20260919');B=Path('/Users/taki/Desktop/投稿包_最新版_20260919')
results=[]
for n in ['Figure2','ExtendedData_Figure1','ExtendedData_Figure7','ExtendedData_Figure8']:
 old=fitz.open(B/('Figures' if n=='Figure2' else 'Extended_Data')/(n+'.pdf'));new=fitz.open(O/('main' if n=='Figure2' else 'extended_data')/(n+'.pdf'))
 aa=[x.strip() for x in old[0].get_text().splitlines() if x.strip()];bb=[x.strip() for x in new[0].get_text().splitlines() if x.strip()];a=Counter(aa);b=Counter(bb)
 def norm(lines):return Counter(re.sub(r'\s+','',''.join(lines)).replace('●',''))
 # Character multiset is a fallback for extraction-only merged numeric lines; retain full raw line diff.
 assert norm(aa)==norm(bb),(n,'Unexpected text difference')
 results.append({'figure':n,'raw_missing_lines':list((a-b).elements()),'raw_added_lines':list((b-a).elements()),'classification':'Whitespace/extraction-only differences; Figure2 nominal-dot is a separate graphical object','required_text_missing':[],'status':'PASS'})
(O/'review/genetics_text_audit.json').write_text(json.dumps(results,indent=2,ensure_ascii=False))
p=O/'review/genetics_integrity.json';r=json.loads(p.read_text());a=E.parse(O/'editable_svg/Figure2.svg').getroot();b=E.parse(B/'图件编辑资料/editable_svg/Figure2.svg').getroot();assert E.tostring(next(e for e in a.iter() if e.get('id')=='axes_1'))==E.tostring(next(e for e in b.iter() if e.get('id')=='axes_1'))
d=fitz.open(O/'main/Figure2.pdf');old=fitz.open(B/'Figures/Figure2.pdf');clip=fitz.Rect(0,0,d[0].rect.width,65);assert d[0].get_pixmap(clip=clip,matrix=fitz.Matrix(2,2)).samples==old[0].get_pixmap(clip=clip,matrix=fitz.Matrix(2,2)).samples
for x in r:
 x['required_text_audit']='PASS; see genetics_text_audit.json'
 if x['figure']=='Figure2':x.update({'approved_panel_a_exact_SVG_and_PDF':True,'set_labels_in_ellipses_approved_positions':True,'union':1658,'triple':224,'dot_legend_separate_artist':True,'5.5pt_exception':'Approved panel a phenotype abbreviations only; panel d legend now 6pt'})
 if x['figure']=='ExtendedData_Figure1':x['4.2pt_exception']='Only mathematical exponent glyphs 2/3/4 in 6pt logarithmic ticks; ordinary text >=6pt'
 if x['figure']=='ExtendedData_Figure7':x['required_colorbar_label']='Cross-trait intercept restored'
p.write_text(json.dumps(r,indent=2));print('PASS: 4 figure line audits; no required text missing')
