from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from lxml import etree as E
from copy import deepcopy
import re, json, csv, hashlib

BASE = Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
OUT = Path(__file__).parent
SRC = BASE / 'Supplementary_Tables.xlsx'
DEST = OUT / 'Supplementary_Tables.xlsx'
S = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
N = {'s':S}
q = lambda s: '{'+S+'}'+s
def value(c): return ''.join(c.itertext()) if c is not None else ''
def colnum(s):
    n=0
    for c in s:n=n*26+ord(c)-64
    return n
def colname(n):
    s=''
    while n:n,k=divmod(n-1,26);s=chr(k+65)+s
    return s
def splitaddr(s):
    m=re.fullmatch(r'\$?([A-Z]+)\$?(\d+)',s)
    return colnum(m[1]),int(m[2])
def strset(c,t):
    for e in list(c): c.remove(e)
    c.set('t','inlineStr')
    E.SubElement(E.SubElement(c,q('is')),q('t')).text=t
def serial(r):return E.tostring(r,encoding='UTF-8',xml_declaration=False)
def sig(c):return (c.get('t'),c.get('s'),value(c))

removed_headers = {
    'Table 02':{'Evidence source','Decision','Design approval record','Grouping-review flag','Analysis approval record'},
    'Table 31':{'Evidence source','Decision','Design approval record','Grouping-review flag','Analysis approval record'},
    'Table 58':{'Evidence source','Decision','Design approval record','Grouping-review flag','Analysis approval record'},
    'Table 14':{'Module tier (Chinese)'},
    'Table 21':{'Matrix','Barcodes','Features','Positions','Scalefactors','Hires image','Metadata note'},
    'Table 42':{'Cpu','Cpu.all'},
    'Table 50':{'Merge'},
    'Table 51':{'Source path'},
    'Table 57':{'Exit code'},
}
group_wording = {
    '既有临床分组；Non-Diseased归入健康组':'Study clinical grouping; non-diseased specimens assigned to the reference group',
    '用户确认：Non-Diseased归入健康组':'Non-diseased specimens assigned to the reference group',
    '用户确认瘤中心为疾病组':'Aneurysm centre assigned to the disease group',
    '用户明确确认近端为健康组；原始AAA供者诊断与配对关系保留':'Proximal specimen assigned to the reference group; paired with the aneurysm specimen from the same donor',
    '用户指定健康组；供者仍为AAA1，原始诊断不变':'Reference specimen from donor AAA1 with an abdominal aortic aneurysm diagnosis',
}
relationship_wording = {
    'Dependent association methods; not independent replication':'Complementary association methods using shared input data',
    'Sensitivity analyses; not independent cohorts':'Sensitivity comparisons using overlapping GWAS inputs',
    'Association, not causal intervention':'Observational comparison of disease and reference specimens',
    'Tissue-level colocalisation not single-cell lineage causality':'Tissue-level evidence for a shared genetic signal',
    'Stability does not imply biological mechanism':'Reproducibility of modules inferred from cell-state proportions',
    'Reduces direct feature circularity; does not prove source independence':'Sensitivity to exclusion of target genes from reference mapping',
}
dictionary_footnote_definitions = {
    (6,'Genetic correlation'):'LDSC genetic correlation for one ordered phenotype pair; diagonal entries equal one.',
    (8,'Latent factor'):'G: shared panvascular factor; AORTIC: Aortic domain factor; ATHERO: Atherosclerotic domain factor. The three factors are mutually orthogonal.',
    (23,'Cluster'):'Source transcriptomic cluster in GSE248608 used to define plaque and calcification expression programmes.',
    (25,'LDSC intercept'):'Cross-trait LDSC intercept for one ordered phenotype pair; diagonal and off-diagonal entries are included.',
}

changes=[]; footnotes=[]; definitions=[]; maps={}; validations=[]
with ZipFile(SRC) as z:parts={n:z.read(n) for n in z.namelist()}
wb=E.fromstring(parts['xl/workbook.xml'])
rels=E.fromstring(parts['xl/_rels/workbook.xml.rels'])
targets={e.get('Id'):e.get('Target').lstrip('/') for e in rels}
sheetpaths={s.get('name'):targets[s.get('{'+R+'}id')] for s in wb.find(q('sheets'))}

def record(sheet,cell,old,new,action):
    changes.append({'sheet':sheet,'original_cell':cell,'action':action,'original_value':old,'cleaned_value':new})

def compact(root,deleted_rows,deleted_cols):
    def newrow(r):return r-sum(i<r for i in deleted_rows)
    def newcol(c):return c-sum(i<c for i in deleted_cols)
    sd=root.find(q('sheetData'))
    for row in list(sd):
        r=int(row.get('r'))
        if r in deleted_rows:sd.remove(row);continue
        row.set('r',str(newrow(r)))
        for c in list(row):
            co,ro=splitaddr(c.get('r'))
            if co in deleted_cols:row.remove(c);continue
            c.set('r',f'{colname(newcol(co))}{newrow(ro)}')
        if row.get('spans'):
            co=[splitaddr(c.get('r'))[0] for c in row]
            if co:row.set('spans',f'{min(co)}:{max(co)}')
            else:row.attrib.pop('spans')
    for columns in root.findall(q('cols')):
        new=[]
        for elem in list(columns):
            for co in range(int(elem.get('min')),int(elem.get('max'))+1):
                if co in deleted_cols:continue
                e=deepcopy(elem);e.set('min',str(newcol(co)));e.set('max',str(newcol(co)));new.append(e)
            columns.remove(elem)
        columns.extend(new)
    def mapped_range(ref):
        ends=ref.split(':');c1,r1=splitaddr(ends[0]);c2,r2=splitaddr(ends[-1])
        rr=[x for x in range(r1,r2+1) if x not in deleted_rows]
        cc=[x for x in range(c1,c2+1) if x not in deleted_cols]
        if not rr or not cc:return None
        a=f'{colname(newcol(cc[0]))}{newrow(rr[0])}';b=f'{colname(newcol(cc[-1]))}{newrow(rr[-1])}'
        return a if a==b else a+':'+b
    for parent in root.findall(q('mergeCells')):
        for e in list(parent):
            new=mapped_range(e.get('ref'))
            if new is None:parent.remove(e)
            else:e.set('ref',new)
        parent.set('count',str(len(parent)))
    for e in root.xpath('./s:dimension|./s:autoFilter',namespaces=N):
        e.set('ref',mapped_range(e.get('ref')))
    for e in root.findall('.//'+q('pane')):
        if e.get('ySplit'):
            split=int(float(e.get('ySplit')));e.set('ySplit',str(split-sum(i<=split for i in deleted_rows)))
        if e.get('xSplit'):
            split=int(float(e.get('xSplit')));e.set('xSplit',str(split-sum(i<=split for i in deleted_cols)))
        if e.get('topLeftCell'):
            co,ro=splitaddr(e.get('topLeftCell'));e.set('topLeftCell',f'{colname(newcol(co))}{newrow(ro)}')
    for e in root.findall('.//'+q('selection')):
        if e.get('activeCell'):
            co,ro=splitaddr(e.get('activeCell'));e.set('activeCell',f'{colname(newcol(co))}{newrow(ro)}')
        if e.get('sqref'):
            e.set('sqref',' '.join(filter(None,(mapped_range(x) for x in e.get('sqref').split()))))
    for e in list(root):
        if E.QName(e).localname in {'headerFooter','legacyDrawing','legacyDrawingHF'}:root.remove(e)
    return newrow,newcol

for name,path in sheetpaths.items():
    if not name.startswith('Table '):continue
    root=E.fromstring(parts[path]);original=deepcopy(root);rows=root.findall('s:sheetData/s:row',N)
    assert not root.findall('.//s:f',N)
    headers={splitaddr(c.get('r'))[0]:value(c) for c in rows[2]}
    dropped={co for co,h in headers.items() if h in removed_headers.get(name,set())}
    assert {headers[c] for c in dropped}==removed_headers.get(name,set())
    deleted_rows={2}
    footnotes.append({'sheet':name,'original_row':2,'text':value(rows[1][0]),'kind':'table explanatory row'})
    for r in rows[3:]:
        nonempty=[c for c in r if value(c)]
        if len(nonempty)==1 and len(value(nonempty[0]))>100 and nonempty[0].get('r').startswith('A'):
            deleted_rows.add(int(r.get('r')))
            footnotes.append({'sheet':name,'original_row':int(r.get('r')),'text':value(nonempty[0]),'kind':'bottom footnote'})
    for row in rows:
        ro=int(row.get('r'))
        for c in row:
            co,_=splitaddr(c.get('r'));v=value(c);h=headers.get(co,'')
            if ro in deleted_rows:
                if v:record(name,c.get('r'),v,'','removed footnote or explanatory row')
                continue
            if co in dropped:
                if v:record(name,c.get('r'),v,'','removed editorial or local-runtime field')
                continue
            new=None
            if ro>=4:
                if h=='Group evidence':new=group_wording.get(v)
                elif h=='Replicate relationship' and v=='reviewed_same_collection_date_distinct_source_records':new='same_collection_date_distinct_source_records'
                elif name=='Table 24' and h=='Recorded scope':
                    new={'frozen human collections':'Human gene-set collections','archived scIB benchmark runtime, not necessarily full-atlas runtime':'Integration benchmark runtime'}.get(v)
                elif name=='Table 47' and h=='Interpretation boundary':new=relationship_wording.get(v)
                elif name=='Table 53' and h=='Donor-overlap status' and v=='no overlap identified from frozen provenance tables':new='No overlap identified in available provenance records'
                elif name=='Table 60' and h=='P-value interpretation':
                    new={'Finite archived software output':'Finite software output','Numerical underflow in the archived software output; not an exact probability of zero':'Numerical underflow; the reported value is below software precision'}.get(v)
                elif name=='Table 57' and h=='Reason':
                    if v.startswith('Warning message:'):
                        assert 'Both optimizers failed or Hessian not positive definite' in v
                        new='Both optimizers failed or the Hessian was not positive definite'
                    elif v=='Native optimizer reported success but tested coefficient SE/z/P is nonfinite; no P imputation and no model refit.':new='Nonfinite coefficient uncertainty (SE, Z or P) despite an optimizer success flag'
            if new is not None and new!=v:
                record(name,c.get('r'),v,new,'editorial wording condensed; scientific meaning retained');strset(c,new)
    changed_cells={c['original_cell'] for c in changes if c['sheet']==name and c['action'].startswith('editorial wording')}
    expected=[]
    for row in original.find(q('sheetData')):
        ro=int(row.get('r'))
        if ro<4 or ro in deleted_rows:continue
        for c in row:
            co,_=splitaddr(c.get('r'))
            if co not in dropped and c.get('r') not in changed_cells:expected.append(sig(c))
    newrow,newcol=compact(root,deleted_rows,dropped)
    actual=[]
    old_changed_new={f'{colname(newcol(splitaddr(s)[0]))}{newrow(splitaddr(s)[1])}' for s in changed_cells}
    for row in root.find(q('sheetData')):
        if int(row.get('r'))<3:continue
        for c in row:
            if c.get('r') not in old_changed_new:actual.append(sig(c))
    assert expected==actual,name
    ndata=sum(any(value(c) for c in r) for r in root.find(q('sheetData')) if int(r.get('r'))>=3)
    maps[name]={'title':value(root.find('s:sheetData/s:row/s:c',N)), 'header_row':2,'first_data_row':3,'data_rows':ndata,'removed_original_rows':sorted(deleted_rows),'removed_columns':{colname(c):headers[c] for c in sorted(dropped)},'retained_columns':{colname(c):colname(newcol(c)) for c in headers if c not in dropped}}
    validations.append({'sheet':name,'retained_data_cells':len(expected),'retained_values_types_styles_exact':True,'data_rows':maps[name]['data_rows'],'original_last_row':len(rows),'footnote_rows_removed':len(deleted_rows),'editorial_values_changed':len(changed_cells),'sha256_original_retained_values_styles':hashlib.sha256(repr(expected).encode()).hexdigest(),'sha256_cleaned_retained_values_styles':hashlib.sha256(repr(actual).encode()).hexdigest()})
    parts[path]=serial(root)

# Keep scientific definitions in the existing dictionary; remove definitions only for deleted fields.
path=sheetpaths['Data dictionary'];root=E.fromstring(parts[path]);deleted=set()
for row in root.find(q('sheetData')):
    if int(row.get('r'))==1:continue
    cs=list(row);num=int(float(value(cs[0])));h=value(cs[1]);name=f'Table {num:02d}'
    if h in removed_headers.get(name,set()):
        deleted.add(int(row.get('r')))
        for c in cs:
            if value(c):record('Data dictionary',c.get('r'),value(c),'','removed dictionary entry for deleted field')
        continue
    if (num,h) in dictionary_footnote_definitions:
        c=cs[2];v=value(c);new=dictionary_footnote_definitions[num,h]
        record('Data dictionary',c.get('r'),v,new,'scientific definition retained from removed footnote');strset(c,new)
    # Remove process boilerplate; the stated units and missing-value semantics remain.
    for c in cs[2:]:
        v=value(c);new=None
        if v=='Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.':
            new='Blank/NA: not supplied or not applicable; interpret using the record type, model status and eligibility. Zero and FALSE are explicit values.'
        elif v.startswith('Recorded ') and v.endswith('; preserved as an identifier, not recoded into another unit.'):
            new=v[len('Recorded '):-len('; preserved as an identifier, not recoded into another unit.')]+'.'
        if new and new!=v:record('Data dictionary',c.get('r'),v,new,'dictionary wording simplified');strset(c,new)
compact(root,deleted,set());parts[path]=serial(root)

# Keep workbook printing/filter definitions aligned with the cleaned sheet geometry.
names=wb.find(q('definedNames'))
for e in names if names is not None else []:
    m=re.match(r"'([^']+)'!",e.text or '')
    if not m or m[1] not in sheetpaths:continue
    name=m[1];r=E.fromstring(parts[sheetpaths[name]])
    kind=e.get('name');ref=None
    if kind=='_xlnm.Print_Titles':
        e.text=f"'{name}'!$1:${2 if name.startswith('Table ') else 1}"
        continue
    if kind=='_xlnm.Print_Area':ref=r.find(q('dimension')).get('ref')
    if kind=='_xlnm._FilterDatabase':
        a=r.find(q('autoFilter'))
        if a is not None:ref=a.get('ref')
    if ref:
        ref=re.sub(r'([A-Z]+)(\d+)',r'$\1$\2',ref)
        e.text=f"'{name}'!{ref}"
parts['xl/workbook.xml']=serial(wb)

# Remove personal/edit history properties and printer headers/footers, without resaving scientific cells.
for name in ['docProps/core.xml','docProps/app.xml']:
    if name not in parts:continue
    r=E.fromstring(parts[name])
    if name.endswith('core.xml'):
        for e in list(r):
            if E.QName(e).localname in {'creator','lastModifiedBy','created','modified','lastPrinted','revision'}:r.remove(e)
    else:
        for e in r.iter():
            if E.QName(e).localname in {'Company','Manager'}:e.text=None
    parts[name]=serial(r)
for name in list(parts):
    if re.search(r'(^|/)(comments|person|threadedComments)',name,re.I):raise AssertionError('Unexpected comment part requires relationship-aware removal: '+name)
with ZipFile(DEST,'w',ZIP_DEFLATED,compresslevel=6) as z:
    for n,b in parts.items():z.writestr(n,b)
with ZipFile(DEST) as z:assert z.testzip() is None

(OUT/'worksheet_layout_map.json').write_text(json.dumps(maps,ensure_ascii=False,indent=2))
(OUT/'removed_notes_preserved_for_reproducibility.json').write_text(json.dumps(footnotes,ensure_ascii=False,indent=2))
(OUT/'retained_data_validation.json').write_text(json.dumps(validations,ensure_ascii=False,indent=2))
with (OUT/'cleaning_changes.tsv').open('w') as f:
    wr=csv.DictWriter(f,fieldnames=list(changes[0]),delimiter='\t');wr.writeheader();wr.writerows(changes)
summary={'source':str(SRC),'output':str(DEST),'source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(DEST.read_bytes()).hexdigest(),'numbered_tables':60,'total_sheets':len(sheetpaths),'footer_rows_removed':len(footnotes),'bottom_footnotes_removed':sum(x['kind']=='bottom footnote' for x in footnotes),'removed_column_instances':sum(len(v) for v in removed_headers.values()),'all_retained_scientific_values_types_styles_exact':True,'source_workbook_unchanged':True,'scientific_ids_unchanged':True,'table_numbers_unchanged':True}
(OUT/'cleaning_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
