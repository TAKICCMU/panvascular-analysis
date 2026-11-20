from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from lxml import etree as ET
import json, re, hashlib, csv

ROOT = Path('/Users/taki/Desktop/投稿包_最新版_20260919')
OUT = Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
WORK = Path(__file__).parent
OUT.mkdir(exist_ok=True)
(OUT/'修订说明').mkdir(exist_ok=True)
REMOVED = {38,39,49,62}
TABLE_MAP = {i:i-sum(j<i for j in REMOVED) for i in range(1,65) if i not in REMOVED}
S='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
R='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
NS={'s':S}
q=lambda x:'{'+S+'}'+x
def val(cell):
    if cell is None:return ''
    return ''.join(cell.itertext())
def setstr(cell,text):
    for x in list(cell):cell.remove(x)
    cell.set('t','inlineStr');ET.SubElement(ET.SubElement(cell,q('is')),q('t')).text=str(text)
def setnum(cell,num):
    for x in list(cell):cell.remove(x)
    cell.set('t','n');ET.SubElement(cell,q('v')).text=str(num)
def compact_rows(root,deleted):
    sd=root.find(q('sheetData'));oldrows=list(sd)
    removed_numbers=[int(r.get('r')) for r in oldrows if int(r.get('r')) in deleted]
    for row in oldrows:
        old=int(row.get('r'))
        if old in deleted:sd.remove(row);continue
        new=old-sum(i<old for i in removed_numbers);row.set('r',str(new))
        for c in row:c.set('r',re.sub(r'\d+$',str(new),c.get('r')))
    last=max(int(x.get('r')) for x in sd)
    for el in root.xpath('./s:dimension|./s:autoFilter',namespaces=NS):
        ref=el.get('ref');el.set('ref',re.sub(r'\d+$',str(last),ref))
    # Only removed rows below the header; its merges and frozen panes stay intact.
    return removed_numbers
def serialize(root):return ET.tostring(root,encoding='utf-8')

with ZipFile(ROOT/'Supplementary_Tables.xlsx') as z:
    parts={n:z.read(n) for n in z.namelist()}
wb=ET.fromstring(parts['xl/workbook.xml']);rels=ET.fromstring(parts['xl/_rels/workbook.xml.rels'])
targets={e.get('Id'):e.get('Target').lstrip('/') for e in rels}
sheet_nodes=list(wb.find(q('sheets')))
sheet_by_name={s.get('name'):targets[s.get('{'+R+'}id')] for s in sheet_nodes}
deleted_parts=set();deleted_relids=set();old_to_new_sheetidx={};idx=0
for oldidx,s in enumerate(sheet_nodes):
    m=re.fullmatch(r'Table (\d+)',s.get('name',''))
    if m and int(m[1]) in REMOVED:
        deleted_parts.add(targets[s.get('{'+R+'}id')]);deleted_relids.add(s.get('{'+R+'}id'))
        wb.find(q('sheets')).remove(s);continue
    old_to_new_sheetidx[oldidx]=idx;idx+=1
    if m:s.set('name',f'Table {TABLE_MAP[int(m[1])]:02d}')
for e in list(rels):
    if e.get('Id') in deleted_relids:rels.remove(e)
for dn in wb.xpath('./s:definedNames/s:definedName',namespaces=NS):
    if dn.get('localSheetId'):
        oldidx=int(dn.get('localSheetId'))
        if oldidx not in old_to_new_sheetidx:dn.getparent().remove(dn);continue
        dn.set('localSheetId',str(old_to_new_sheetidx[oldidx]))
    if dn.text:
        dn.text=re.sub(r'Table (\d{2})',lambda m:f'Table {TABLE_MAP[int(m[1])]:02d}' if int(m[1]) in TABLE_MAP else m[0],dn.text)
parts['xl/workbook.xml']=serialize(wb);parts['xl/_rels/workbook.xml.rels']=serialize(rels)
ct=ET.fromstring(parts['[Content_Types].xml'])
for e in list(ct):
    if e.get('PartName','').lstrip('/') in deleted_parts:ct.remove(e)
parts['[Content_Types].xml']=serialize(ct)
for f in deleted_parts:parts.pop(f)

changes=[]
for old,new in TABLE_MAP.items():
    path=sheet_by_name[f'Table {old:02d}'];root=ET.fromstring(parts[path])
    cell=root.find('.//'+q('c')+'[@r="A1"]')
    if old!=new:
        setstr(cell,re.sub(r'Supplementary Table \d+',f'Supplementary Table {new}',val(cell),count=1))
    if old in {24,47,50,51}:
        deleted=[]
        for row in root.find(q('sheetData')):
            if int(row.get('r'))<4:continue
            text=' '.join(val(c) for c in row)
            if re.search(r'communication|LIANA|CellChat',text,re.I):
                deleted.append(int(row.get('r')));changes.append({'old_table':old,'old_row':int(row.get('r')),'removed_record':text})
        if deleted:compact_rows(root,set(deleted))
    if old!=new or old in {24,47,50,51}:parts[path]=serialize(root)

for name in ['Information','Data dictionary']:
    path=sheet_by_name[name];root=ET.fromstring(parts[path]);deleted=[]
    for row in root.find(q('sheetData')):
        if int(row.get('r'))==1:continue
        cells=list(row);first=cells[0]
        old=int(float(val(first)))
        if old in REMOVED:deleted.append(int(row.get('r')));continue
        setnum(first,TABLE_MAP[old])
        if name=='Information':setstr(cells[4],f'Table {TABLE_MAP[old]:02d}')
    compact_rows(root,set(deleted));parts[path]=serialize(root)

# Update package metadata only. Do not resave scientific cells through a numeric parser.
for name in ['docProps/core.xml','docProps/app.xml']:
    if name not in parts:continue
    root=ET.fromstring(parts[name])
    if name.endswith('core.xml'):
        for e in list(root):
            if ET.QName(e).localname in {'creator','lastModifiedBy','created','modified','lastPrinted','revision'}:root.remove(e)
    else:
        for e in root.iter():
            if ET.QName(e).localname in {'Company','Manager'}:e.text=None
        # Worksheet list is optional; remove stale extended-property caches.
        for e in list(root):
            if ET.QName(e).localname in {'HeadingPairs','TitlesOfParts'}:root.remove(e)
    parts[name]=serialize(root)
target=OUT/'Supplementary_Tables.xlsx'
with ZipFile(target,'w',ZIP_DEFLATED,compresslevel=6) as z:
    for n,b in parts.items():z.writestr(n,b)

with (OUT/'修订说明/补充表旧新编号.tsv').open('w') as h:
    writer=csv.writer(h,delimiter='\t');writer.writerow(['原表号','新表号','处理'])
    for old in range(1,65):writer.writerow([old,TABLE_MAP.get(old,''),'移出本次论文；原结果保留在原包' if old in REMOVED else '保留，数据不变'])
(WORK/'table_map.json').write_text(json.dumps(TABLE_MAP,indent=2))
(WORK/'supplement_removals.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2))

# Exact retained numeric XML values and cell styles are unchanged outside index metadata.
validation=[]
with ZipFile(ROOT/'Supplementary_Tables.xlsx') as a, ZipFile(target) as b:
    for old,new in TABLE_MAP.items():
        path=sheet_by_name[f'Table {old:02d}']
        ar=ET.fromstring(a.read(path));br=ET.fromstring(b.read(path))
        def rows(root):
            return [[(c.get('t'),c.get('s'),val(c)) for c in row] for row in root.find(q('sheetData')) if int(row.get('r'))>=4]
        aa=rows(ar);bb=rows(br)
        if old in {24,47,50,51}:aa=[r for r in aa if not re.search(r'communication|LIANA|CellChat',' '.join(c[2] for c in r),re.I)]
        ok=aa==bb
        validation.append({'old_table':old,'new_table':new,'data_rows':len(bb),'exact_retained_cells_and_styles':ok})
        assert ok,old
    assert not b.testzip()
(WORK/'spreadsheet_validation.json').write_text(json.dumps(validation,indent=2))
print('Saved',target,'60 tables; exact retained-cell comparison passed. Removed metadata rows:',len(changes))
