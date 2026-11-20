from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from lxml import etree as E
from copy import deepcopy
import json,re,hashlib,csv

W=Path(__file__).resolve().parents[1]
SRC=W/'tables/Supplementary_Tables.xlsx'
OUT=W/'staged/Supplementary_Tables.xlsx'
MAP={int(k):v for k,v in json.loads((W/'table_map_current_to_final.json').read_text()).items()}
assert set(MAP)==set(range(1,61)) and set(MAP.values())==set(range(1,61))
S='http://schemas.openxmlformats.org/spreadsheetml/2006/main';R='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
N={'s':S};q=lambda x:'{'+S+'}'+x
def val(c):return ''.join(c.itertext())
def serial(r):return E.tostring(r,encoding='UTF-8')
def setstr(c,t):
    for e in list(c):c.remove(e)
    c.set('t','inlineStr');E.SubElement(E.SubElement(c,q('is')),q('t')).text=t
def setnum(c,n):
    for e in list(c):c.remove(e)
    c.set('t','n');E.SubElement(c,q('v')).text=str(n)
def table_refs(s):return re.sub(r'Table (\d{2})',lambda m:f'Table {MAP[int(m[1])]:02d}',s)
def csignature(c):return c.get('t'),c.get('s'),val(c)
def data_signature(root):
    return [[csignature(c) for c in row] for row in root.find(q('sheetData')) if int(row.get('r'))>=3]

with ZipFile(SRC) as z:parts={n:z.read(n) for n in z.namelist()}
original_parts=parts.copy()
wb=E.fromstring(parts['xl/workbook.xml']);sheets=wb.find(q('sheets'))
rels=E.fromstring(parts['xl/_rels/workbook.xml.rels'])
targets={e.get('Id'):e.get('Target').lstrip('/') for e in rels}
old_nodes=list(sheets);paths={s.get('name'):targets[s.get('{'+R+'}id')] for s in sheets}
old_index={s.get('{'+R+'}id'):i for i,s in enumerate(old_nodes)}
table_validation=[];layout=[]
for s in sheets:
    m=re.fullmatch(r'Table (\d{2})',s.get('name'))
    if not m:continue
    old=int(m[1]);new=MAP[old];path=targets[s.get('{'+R+'}id')]
    a=E.fromstring(parts[path]);before=data_signature(a)
    cell=a.find('s:sheetData/s:row[@r="1"]/s:c[@r="A1"]',N)
    title=val(cell);assert title.startswith(f'Supplementary Table {old} |')
    setstr(cell,re.sub(r'^Supplementary Table \d+',f'Supplementary Table {new}',title))
    for h in a.findall('.//s:hyperlink',N):
        if h.get('location'):h.set('location',table_refs(h.get('location')))
    parts[path]=serial(a);s.set('name',f'Table {new:02d}')
    after=data_signature(a);assert before==after
    table_validation.append({'old_table':old,'new_table':new,'retained_cells':sum(len(r) for r in before),'values_types_styles_identical':True})
    layout.append({'old_table':old,'new_table':new,'sheet':f'Table {new:02d}','title':val(cell),'header_row':2,'first_data_row':3,'xml_path':path})

# Sheet order follows citation order; dictionary and abbreviations remain accessible.
def sortkey(s):
    name=s.get('name')
    if name=='Information':return -2
    if name=='Abbreviations':return -1
    if name=='Data dictionary':return 61
    return int(name.split()[1])
ordered=sorted(list(sheets),key=sortkey)
for s in list(sheets):sheets.remove(s)
sheets.extend(ordered)
new_index={s.get('{'+R+'}id'):i for i,s in enumerate(ordered)}
indexmap={old_index[k]:new_index[k] for k in old_index}
for dn in wb.findall('s:definedNames/s:definedName',N):
    if dn.get('localSheetId') is not None:dn.set('localSheetId',str(indexmap[int(dn.get('localSheetId'))]))
    if dn.text:dn.text=table_refs(dn.text)
for v in wb.findall('s:bookViews/s:workbookView',N):
    if v.get('activeTab') is not None:v.set('activeTab',str(indexmap[int(v.get('activeTab'))]))
parts['xl/workbook.xml']=serial(wb)

def reorder_index(name):
    path=paths[name];root=E.fromstring(parts[path]);sd=root.find(q('sheetData'));header=sd[0];rows=list(sd)[1:]
    for r in rows:
        cells=list(r);old=int(float(val(cells[0])));setnum(cells[0],MAP[old])
        if name=='Information':setstr(cells[4],f'Table {MAP[old]:02d}')
    rows.sort(key=lambda r:int(float(val(r[0]))))
    for r in list(sd):sd.remove(r)
    sd.append(header)
    for i,r in enumerate(rows,2):
        r.set('r',str(i))
        for c in r:c.set('r',re.sub(r'\d+$',str(i),c.get('r')))
        sd.append(r)
    if name=='Information':
        # Rebuild true internal links. The previous package retained 64 stale external-mode link receipts.
        old=root.find(q('hyperlinks'))
        if old is not None:root.remove(old)
        hs=E.Element(q('hyperlinks'))
        for i in range(1,61):E.SubElement(hs,q('hyperlink'),ref=f'E{i+1}',location=f"'Table {i:02d}'!A1",display=f'Table {i:02d}')
        # OOXML order: hyperlinks precede page margins/setup.
        page=root.find(q('pageMargins'))
        root.insert(list(root).index(page) if page is not None else len(root),hs)
        relpath=str(Path(path).parent/'_rels'/(Path(path).name+'.rels'))
        if relpath in parts:
            rel=E.fromstring(parts[relpath])
            assert all(e.get('Type','').endswith('/hyperlink') for e in rel)
            parts.pop(relpath)
    parts[path]=serial(root)
reorder_index('Information');reorder_index('Data dictionary')

OUT.parent.mkdir(exist_ok=True)
with ZipFile(OUT,'w',ZIP_DEFLATED,compresslevel=6) as z:
    for n,b in parts.items():z.writestr(n,b)
with ZipFile(OUT) as z:
    assert z.testzip() is None
    outwb=E.fromstring(z.read('xl/workbook.xml'))
    names=[s.get('name') for s in outwb.find(q('sheets'))]
    assert names==['Information','Abbreviations']+[f'Table {i:02d}' for i in range(1,61)]+['Data dictionary']
    for old,new in MAP.items():
        path=paths[f'Table {old:02d}']
        assert data_signature(E.fromstring(original_parts[path]))==data_signature(E.fromstring(z.read(path)))
    info=E.fromstring(z.read(paths['Information']))
    for i,row in enumerate(info.findall('s:sheetData/s:row',N)[1:],1):
        assert int(val(row[0]))==i and val(row[4])==f'Table {i:02d}'
    assert len(info.findall('s:hyperlinks/s:hyperlink',N))==60

validation={'source':str(SRC),'output':str(OUT),'source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(OUT.read_bytes()).hexdigest(),'tables':60,'sheets':63,'all_retained_scientific_cells_exact':True,'table_headers_row':2,'data_start_row':3,'numbering_and_navigation_only':True,'internal_links_rebuilt':60,'old_stale_link_count':64,'sheet_order_matches_final_numbers':True,'definedNames_localSheetId_rebased':True,'per_table':sorted(table_validation,key=lambda r:r['new_table'])}
(OUT.parent/'renumber_validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
(W/'tables/final_worksheet_layout.json').write_text(json.dumps(sorted(layout,key=lambda r:r['new_table']),ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in validation.items() if k!='per_table'},ensure_ascii=False,indent=2))
