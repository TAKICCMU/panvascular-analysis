from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
from copy import deepcopy
import json, hashlib, csv
import openpyxl

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
O=P/'图件编辑资料/补充表核查_20261118'; S=O/'staging'
spec=json.loads((O/'correction_spec.json').read_text())
manifest=json.loads((O/'correction_input_hashes.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for f,h in manifest['affected'].items():assert sha(P/f)==h,('Concurrent original change',f)
for f,h in manifest['protected'].items():assert sha(P/f)==h,('Protected file change',f)

old=openpyxl.load_workbook(P/'Supplementary_Tables.xlsx',read_only=True,data_only=False)
new=openpyxl.load_workbook(S/'Supplementary_Tables.xlsx',read_only=True,data_only=False)
assert old.sheetnames==new.sheetnames
edits={(c['sheet'],c['cell']):c for c in spec['changes']}
changes=[]; checked=0; cells8=0; table24=0
for name in old.sheetnames:
    a=old[name];b=new[name]
    if name in ('Table 59','Table 60'):
        expected=spec['families' if name=='Table 59' else 'versions']
        actual=[list(r) for r in b.iter_rows(min_row=3,values_only=True)]
        assert actual==expected,(name,'data mismatch')
        assert b.max_row==len(expected)+2
        if name=='Table 60':
            assert expected[:13]==[list(r) for r in a.iter_rows(min_row=3,values_only=True)]
        continue
    assert (a.max_row,a.max_column)==(b.max_row,b.max_column),name
    for ra,rb in zip(a.iter_rows(),b.iter_rows()):
        for ca,cb in zip(ra,rb):
            if ca.value is None and cb.value is None:continue
            checked+=1
            key=(name,ca.coordinate)
            if key in edits:
                e=edits[key]
                assert ca.value==e['before'] and cb.value==e['after'],key
                changes.append(key)
            else:
                assert ca.value==cb.value and ca.data_type==cb.data_type,(key,ca.value,cb.value,ca.data_type,cb.data_type)
                assert ca.number_format==cb.number_format,('Unexpected format',key)
            if name=='Table 08' and ca.row>=3 and ca.column>=2:cells8+=1
            if name=='Table 24' and ca.row>=3:table24+=1
assert cells8==192
for e in spec['numeric_cells']:
    c=new[e['sheet']][e['cell']]
    assert c.data_type=='n' and c.value==int(e['stored_value'])
    assert c.number_format=='#,##0'
assert 'Figure 1h' in new['Data dictionary']['C83'].value
assert 'Figure 1g' not in new['Data dictionary']['C83'].value
assert any('Supplementary Table 22' in str(v) for r in spec['families'] for v in r)

N='http://schemas.openxmlformats.org/spreadsheetml/2006/main';T=lambda x:'{'+N+'}'+x
def sheetmap(z):
    w=E.fromstring(z.read('xl/workbook.xml'));r=E.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    rel={x.get('Id'):x.get('Target') for x in r};o={}
    for x in w.find(T('sheets')):
        t=rel[x.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')]
        o[x.get('name')]=t.lstrip('/') if t.startswith('/') else 'xl/'+t
    return o
layout=[]
with ZipFile(P/'Supplementary_Tables.xlsx') as za,ZipFile(S/'Supplementary_Tables.xlsx') as zb:
    assert za.testzip() is None and zb.testzip() is None
    ma=sheetmap(za);mb=sheetmap(zb);assert ma==mb
    assert za.read(ma['Table 08'])==zb.read(mb['Table 08'])
    sa=E.fromstring(za.read('xl/styles.xml'));sb=E.fromstring(zb.read('xl/styles.xml'))
    for typ in ['fonts','fills','borders','cellXfs','cellStyleXfs','numFmts']:
        for x,y in zip(sa.find(T(typ)),sb.find(T(typ))):assert E.tostring(x)==E.tostring(y),typ
    for name in ma:
        ar=E.fromstring(za.read(ma[name]));br=E.fromstring(zb.read(mb[name]))
        for tag in ['sheetViews','mergeCells','sheetProtection','pageMargins','pageSetup','printOptions']:
            x=ar.find(T(tag));y=br.find(T(tag))
            assert (None if x is None else E.tostring(x))==(None if y is None else E.tostring(y)),(name,tag)
        af=ar.find(T('autoFilter'));bf=br.find(T('autoFilter'))
        expected={'Table 59':'A2:D20','Table 60':'A2:C65'}.get(name,af.get('ref') if af is not None else None)
        assert (bf.get('ref') if bf is not None else None)==expected,(name,'filter')
        if name=='Table 24':
            rows={int(x.get('r')):x for x in br.find(T('sheetData'))}
            assert all(float(rows[r].get('ht'))==60 for r in range(3,32))
            assert max(len(str(c.value)) for row in new[name].iter_rows(min_row=3) for c in row)>=1706
        layout.append({'sheet':name,'filter':expected,'frozen_views_unchanged':True})

ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
docreport=json.loads((O/'document_corrections.json').read_text())
for d in docreport:
    f=d['file'];assert sha(P/f)==d['before_sha256'];assert sha(S/f)==d['after_sha256']
    with ZipFile(P/f) as za,ZipFile(S/f) as zb:
        assert za.namelist()==zb.namelist()
        for n in za.namelist():
            if n!='word/document.xml':assert za.read(n)==zb.read(n),(f,n)
        ar=E.fromstring(za.read('word/document.xml'));br=E.fromstring(zb.read('word/document.xml'))
        for x in ar.xpath('//w:t',namespaces=ns):x.text='TEXT'
        for x in br.xpath('//w:t',namespaces=ns):x.text='TEXT'
        assert E.tostring(ar)==E.tostring(br),(f,'non-text XML')

rankpath='Source_Data_Figures/ExtendedData1/benchmark_publication_CPU_only.tsv'
with (P/rankpath).open() as f:ar=list(csv.DictReader(f,delimiter='\t'))
with (S/rankpath).open() as f:br=list(csv.DictReader(f,delimiter='\t'))
assert len(ar)==len(br)==8
for a,b in zip(ar,br):
    a.pop('rank');b.pop('rank');assert a==b
with (S/rankpath).open() as f:br=list(csv.DictReader(f,delimiter='\t'))
assert sorted(int(r['rank']) for r in br)==list(range(1,9))

files={f:sha(S/f) for f in manifest['affected']}
result={'status':'VERIFIED_FOR_IN_PLACE_REPLACEMENT','sheets':len(new.sheetnames),'original_cells_checked_outside_expanded_tables':checked,
        'table8_numeric_values_unchanged':cells8,'table8_worksheet_byte_identical':True,'table24_data_cells_unchanged':table24,
        'table24_data_row_height_pt':60,'numeric_text_conversions':len(spec['numeric_cells']),
        'table59_definitions':len(spec['families']),'table60_version_resource_seed_records':len(spec['versions']),
        'table60_original_13_records_retained':True,'docx_non_text_xml_and_citation_fields_preserved':True,
        'benchmark_score_strings_unchanged':True,'scientific_models_and_p_values_not_recomputed':True,
        'protected_files_unchanged':len(manifest['protected']),'staged_sha256':files,'layout_checks':layout}
(O/'correction_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in result.items() if k not in ('layout_checks','staged_sha256')},ensure_ascii=False))
