from pathlib import Path
from zipfile import ZipFile
from collections import Counter
import hashlib, json, csv, re, math, datetime
import openpyxl
from lxml import etree as E

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
O=P/'图件编辑资料/补充表核查_20261118'
X=P/'Supplementary_Tables.xlsx'
O.mkdir(exist_ok=True)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(name,obj): (O/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str))
def tsv(name, rows):
    if not rows: return
    with (O/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
baseline={'workbook':str(X),'sha256':sha(X),'size_bytes':X.stat().st_size,
          'protected':{str(p.relative_to(P)):sha(p) for p in P.rglob('*') if p.is_file() and ('图件编辑资料' not in p.parts) and (p.suffix.lower() in {'.xlsx','.docx','.pdf','.tsv','.csv','.gz','.svg'})}}
dump('input_identity.json',baseline)
W=openpyxl.load_workbook(X,read_only=True,data_only=False)
tables={}; sheetmeta=[]; colmeta=[]; anomalies=[]; previews={}
numeric_string=re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$')
for s in W:
    rows=list(s.iter_rows())
    numbered=bool(re.fullmatch(r'Table \d{2}',s.title))
    start=2 if numbered or s.title=='Abbreviations' else 1
    hs=[c.value for c in rows[start-1]]
    vals=[tuple(c.value for c in rr) for rr in rows[start:] if any(c.value is not None for c in rr)]
    tables[s.title]={'headers':hs,'rows':vals}
    previews[s.title]={'title':rows[0][0].value,'headers':hs,'first':vals[:3],'last':vals[-2:]}
    dup=Counter(vals); duprows=sum(n-1 for n in dup.values() if n>1)
    sheetmeta.append({'sheet':s.title,'title':rows[0][0].value,'data_rows':len(vals),'columns':len(hs),'exact_duplicate_rows':duprows,
                      'duplicate_headers':';'.join(str(k) for k,v in Counter(hs).items() if v>1), 'blank_headers':sum(h is None for h in hs)})
    for j,h in enumerate(hs):
        cells=[r[j] for r in rows[start:] if j<len(r)]; present=[c for c in cells if c.value is not None]
        nums=[float(c.value) for c in present if isinstance(c.value,(int,float)) and not isinstance(c.value,bool)]
        txtnums=[c for c in present if isinstance(c.value,str) and numeric_string.fullmatch(c.value.strip())]
        meta={'sheet':s.title,'column':openpyxl.utils.get_column_letter(j+1),'header':h,'values':len(present),
              'types':json.dumps(Counter(type(c.value).__name__ for c in present)),
              'formats':json.dumps(Counter(c.number_format for c in present)),
              'min':min(nums) if nums else '', 'max':max(nums) if nums else '',
              'numeric_strings':len(txtnums),'numeric_string_examples':'; '.join(f'{c.coordinate}={c.value}' for c in txtnums[:3]),
              'na_markers':json.dumps(Counter(str(c.value) for c in present if str(c.value).strip().lower() in {'na','nan','n/a','nr','not reported','none','null'}))}
        colmeta.append(meta)
        for c in present:
            if c.data_type in {'e','f'}:anomalies.append({'sheet':s.title,'cell':c.coordinate,'field':h,'kind':c.data_type,'value':c.value})
            if isinstance(c.value,(datetime.datetime,datetime.date)):anomalies.append({'sheet':s.title,'cell':c.coordinate,'field':h,'kind':'date_cell','value':str(c.value)})
        if numbered:
            ht=str(h).lower()
            prob=bool(re.fullmatch(r'(?:meta |magma |scdrs |empirical |minimum gwas |minimum eqtl |raw software |original-reference |target-excluded |fallback within-spot |fallback lag |parent-matched within-spot |parent-matched lag )?p(?: value)?',ht) or ht.startswith('pp.h') or ht in {'adjusted p','target-set bh fdr','bh fdr','magma bh fdr','scdrs fdr','fdr','primary posterior inclusion probability (pip)'} or re.fullmatch(r'bh fdr \(\d+ tests\)',ht))
            unitfrac=('detected fraction' in ht or ht in {'direction agreement','top study share','top tissue share','positive pair fraction','positive fdr pair fraction'})
            for c in present:
                v=c.value
                if not isinstance(v,(int,float)) or isinstance(v,bool):continue
                if not math.isfinite(v):anomalies.append({'sheet':s.title,'cell':c.coordinate,'field':h,'kind':'nonfinite','value':v})
                if (prob or unitfrac) and not (-1e-12<=v<=1+1e-12):anomalies.append({'sheet':s.title,'cell':c.coordinate,'field':h,'kind':'outside_0_1','value':v})
                if ht in {'se','meta se','magma se','τ²','original-reference se','target-excluded se'} and v<0:anomalies.append({'sheet':s.title,'cell':c.coordinate,'field':h,'kind':'negative_se_or_tau2','value':v})
                if ht.endswith('(n)') and (v<0 or abs(v-round(v))>1e-8):anomalies.append({'sheet':s.title,'cell':c.coordinate,'field':h,'kind':'noninteger_count','value':v})
    print(s.title,len(vals),flush=True)
dump('sheet_previews.json',previews);tsv('worksheet_inventory.tsv',sheetmeta);tsv('column_profiles.tsv',colmeta);dump('initial_anomalies.json',anomalies)

ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with ZipFile(X) as z:
    wb=E.fromstring(z.read('xl/workbook.xml'))
    rel=E.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    targets={r.get('Id'):r.get('Target').lstrip('/') for r in rel}
    structure=[]
    for sn in wb.find('s:sheets',ns):
        target=targets[sn.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')]
        if not target.startswith('xl/'):target='xl/'+target
        root=E.fromstring(z.read(target))
        cols=root.findall('s:cols/s:col',ns)
        structure.append({'sheet':sn.get('name'),'state':sn.get('state','visible'),
                          'pane':[dict(c.attrib) for c in root.findall('s:sheetViews/s:sheetView/s:pane',ns)],
                          'filter':[dict(c.attrib) for c in root.findall('s:autoFilter',ns)],
                          'merged':[c.get('ref') for c in root.findall('s:mergeCells/s:mergeCell',ns)],
                          'column_dimensions':[dict(c.attrib) for c in cols],
                          'hidden_rows':[c.get('r') for c in root.findall('s:sheetData/s:row',ns) if c.get('hidden')=='1'],
                          'hidden_columns':[dict(c.attrib) for c in cols if c.get('hidden')=='1'],
                          'row_heights':Counter(c.get('ht','default') for c in root.findall('s:sheetData/s:row',ns)),
                          'hyperlinks':[dict(c.attrib) for c in root.findall('s:hyperlinks/s:hyperlink',ns)]})
    dump('workbook_structure.json',{'zip_integrity':z.testzip(),'sheets':structure,'external_parts':[n for n in z.namelist() if any(q in n for q in ['externalLink','vba','comments','person'])],
                                    'defined_names':[dict(c.attrib, value=c.text) for c in wb.findall('s:definedNames/s:definedName',ns)]})

checks=[]
def check(name,ok,detail):checks.append({'check':name,'pass':bool(ok),'detail':detail})
check('numbered_sheets',W.sheetnames==['Information','Abbreviations']+[f'Table {i:02d}' for i in range(1,61)]+['Data dictionary'],W.sheetnames)
info=tables['Information']['rows']; dd=tables['Data dictionary']['rows']; ddkeys=Counter((int(r[0]),r[1]) for r in dd)
missing=[];orphan=[]
for n in range(1,61):
    name=f'Table {n:02d}'; headers=tables[name]['headers']; title=previews[name]['title']
    ir=[r for r in info if r[0]==n]
    check('table_index_'+str(n),len(ir)==1 and title==f'Supplementary Table {n} | {ir[0][1]}' and ir[0][4]==name,ir)
    for h in headers:
        if (n,h) not in ddkeys: missing.append((n,h))
for n,h in ddkeys:
    if n not in range(1,61) or h not in tables[f'Table {n:02d}']['headers']:orphan.append((n,h))
check('dictionary_coverage',not missing and not orphan and all(v==1 for v in ddkeys.values()),{'missing':missing,'orphan':orphan,'duplicates':[k for k,v in ddkeys.items() if v>1]})
for n in [1,2,3]:
    tab=tables[f'Table {n:02d}']; j=tab['headers'].index('Cells (n)'); total=sum(r[j] for r in tab['rows'])
    check(f'cell_total_table{n}',total==1175115,total)
for n,expected in [(1,26),(2,247),(3,247),(6,116),(7,116),(11,8),(24,29),(26,15639),(30,14),(31,89),(35,49),(36,87),(37,348)]:
    check(f'rowcount_table{n}',len(tables[f'Table {n:02d}']['rows'])==expected,len(tables[f'Table {n:02d}']['rows']))
check('gene_state_population_table26',len({r[0] for r in tables['Table 26']['rows']})==1730,len({r[0] for r in tables['Table 26']['rows']}))
check('GWAS_total_N',all(r[1]+r[2]==r[3] for r in tables['Table 11']['rows']),tables['Table 11']['rows'])
for n in [8]:
    failures=[(i+3,r[0],sum(r[1:])) for i,r in enumerate(tables[f'Table {n:02d}']['rows']) if not math.isclose(sum(r[1:]),1,abs_tol=1e-8)]
    check(f'composition_sum_table{n}',not failures,failures)
dump('structural_and_total_checks.json',checks)
dump('table_dictionary.json',[dict(zip(tables['Data dictionary']['headers'],r)) for r in dd])

from docx import Document
docs={}
for p in [P/'Manuscript.docx',P/'Figure_Legends.docx',P/'Supplementary_Information.docx',P/'Figure/Figure_Legends.docx']:
    doc=Document(p);pars=[x.text for x in doc.paragraphs]+[c.text for t in doc.tables for r in t.rows for c in r.cells]
    docs[str(p.relative_to(P))]=[{'paragraph':i+1,'text':t} for i,t in enumerate(pars) if re.search(r'Supplementary|Table |Figure |Fig\.',t)]
dump('document_table_references.json',docs)
check('workbook_unchanged',sha(X)==baseline['sha256'],sha(X))
dump('structural_and_total_checks.json',checks)
print('INITIAL AUDIT COMPLETE',json.dumps({'checks':len(checks),'failed':[q for q in checks if not q['pass']],'anomalies':len(anomalies),'numeric_text_columns':sum(bool(c['numeric_strings']) for c in colmeta)},ensure_ascii=False),flush=True)
