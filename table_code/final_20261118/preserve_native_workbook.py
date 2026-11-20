"""Retain untouched OOXML; transplant only the edits authored by Artifact Tool.

Artifact export omits worksheet AutoFilters. This preservation pass retains all
original package parts, filters, metadata and non-target cells, importing only
the authorized edited cells/styles and layout attributes from its export.
"""
from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
from copy import deepcopy,copy
import json,re,hashlib

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116');O=P/'图件编辑资料/补充表核查_20261118';S=O/'staging'
spec=json.loads((O/'correction_spec.json').read_text())
N='http://schemas.openxmlformats.org/spreadsheetml/2006/main';ns={'x':N};T=lambda t:'{'+N+'}'+t
with ZipFile(P/'Supplementary_Tables.xlsx') as z:
    infos=z.infolist();parts={i.filename:z.read(i.filename) for i in infos};comment=z.comment
with ZipFile(S/'artifact_edited.xlsx') as z:ap={n:z.read(n) for n in z.namelist()}
def sheetmap(data):
    root=E.fromstring(data['xl/workbook.xml']);rels=E.fromstring(data['xl/_rels/workbook.xml.rels'])
    rr={x.get('Id'):x.get('Target') for x in rels};out={}
    for x in root.find(T('sheets')):
        target=rr[x.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')]
        out[x.get('name')]=target.lstrip('/') if target.startswith('/') else 'xl/'+target
    return out
om=sheetmap(parts);am=sheetmap(ap)
sst=[]
if 'xl/sharedStrings.xml' in ap:
    sst=E.fromstring(ap['xl/sharedStrings.xml']).findall(T('si'))
os=E.fromstring(parts['xl/styles.xml']);ats=E.fromstring(ap['xl/styles.xml'])
cache={}
def import_component(kind,i):
    key=(kind,i)
    if key in cache:return cache[key]
    src=ats.find(T(kind));dst=os.find(T(kind));node=deepcopy(src[i])
    if kind in ['cellXfs','cellStyleXfs']:
        for key2,kind2 in [('fontId','fonts'),('fillId','fills'),('borderId','borders')]:
            if node.get(key2) is not None:node.set(key2,str(import_component(kind2,int(node.get(key2)))))
        if node.get('numFmtId') is not None:node.set('numFmtId',str(import_numfmt(int(node.get('numFmtId')))))
        if node.get('xfId') is not None:node.set('xfId',str(import_component('cellStyleXfs',int(node.get('xfId')))))
    raw=E.tostring(node)
    for j,el in enumerate(dst):
        if E.tostring(el)==raw:cache[key]=j;return j
    j=len(dst);dst.append(node);dst.set('count',str(len(dst)));cache[key]=j;return j
def import_numfmt(i):
    if i<164:return i
    source=ats.find(T('numFmts'));dst=os.find(T('numFmts'))
    node=next(x for x in source if int(x.get('numFmtId'))==i)
    for el in dst:
        if el.get('formatCode')==node.get('formatCode'):return int(el.get('numFmtId'))
    new=max([163]+[int(x.get('numFmtId')) for x in dst])+1
    el=deepcopy(node);el.set('numFmtId',str(new));dst.append(el);dst.set('count',str(len(dst)));return new
def artifact_cell(c):
    c=deepcopy(c)
    if c.get('s') is not None:c.set('s',str(import_component('cellXfs',int(c.get('s')))))
    if c.get('t')=='s':
        si=sst[int(c.find(T('v')).text)];c.remove(c.find(T('v')));c.set('t','inlineStr');inline=E.SubElement(c,T('is'))
        for x in si:inline.append(deepcopy(x))
    return c
def colindex(addr):
    n=0
    for c in re.match('[A-Z]+',addr)[0]:n=n*26+ord(c)-64
    return n
modified={};log=[]
names=set(x['sheet'] for x in spec['changes'])|{'Table 24','Table 59','Table 60'}
for name in sorted(names):
    root=E.fromstring(parts[om[name]]);ar=E.fromstring(ap[am[name]])
    sd=root.find(T('sheetData'));ad=ar.find(T('sheetData'))
    rows={int(x.get('r')):x for x in sd};acells={x.get('r'):x for x in ad.iter(T('c'))}
    addresses={x['cell'] for x in spec['changes'] if x['sheet']==name}
    if name in ['Table 59','Table 60']:
        n=len(spec['families' if name=='Table 59' else 'versions'])+2;last='D' if name=='Table 59' else 'C'
        addresses|={f'{c}{r}' for r in range(1,n+1) for c in ('ABCD' if last=='D' else 'ABC') if f'{c}{r}' in acells}
        root.find(T('dimension')).set('ref',f'A1:{last}{n}')
        root.find(T('autoFilter')).set('ref',f'A2:{last}{n}')
        original_cols=root.find(T('cols'));newcols=deepcopy(ar.find(T('cols')))
        idx=list(root).index(original_cols);root.remove(original_cols);root.insert(idx,newcols)
    for addr in sorted(addresses,key=lambda a:(int(re.search(r'\d+',a)[0]),colindex(a))):
        r=int(re.search(r'\d+',addr)[0]);row=rows.get(r)
        if row is None:
            row=E.Element(T('row'),r=str(r));sd.append(row);rows[r]=row
        src=acells[addr];new=artifact_cell(src);old=next((x for x in row if x.get('r')==addr),None)
        if old is not None:row.replace(old,new)
        else:
            idx=next((i for i,x in enumerate(row) if colindex(x.get('r'))>colindex(addr)),len(row));row.insert(idx,new)
    # Only formatting, never the scientific cell values, is taken for Table 24.
    if name=='Table 24':
        for r in range(3,32):
            for c in spec['table24_singleline_columns']:
                addr=f'{c}{r}';target=next(x for x in rows[r] if x.get('r')==addr);target.set('s',str(import_component('cellXfs',int(acells[addr].get('s','0')))))
    height_rows=[]
    if name=='Table 24':height_rows=range(3,32)
    elif name in ['Table 59','Table 60']:height_rows=range(2,n+1)
    elif name=='Information':height_rows=[9,60,61]
    elif name=='Data dictionary':height_rows=list(range(114,130))+list(range(904,911))
    arows={int(x.get('r')):x for x in ad}
    for r in height_rows:
        rows[r].set('ht',arows[r].get('ht'));rows[r].set('customHeight','1')
    modified[om[name]]=E.tostring(root,encoding='utf-8')
    log.append({'sheet':name,'cells_from_artifact':len(addresses),'layout_rows':list(height_rows),'native_filter':root.find(T('autoFilter')).get('ref')})
modified['xl/styles.xml']=E.tostring(os,encoding='utf-8')
out=S/'Supplementary_Tables.xlsx'
with ZipFile(out,'w') as z:
    z.comment=comment
    for i in infos:z.writestr(copy(i),modified.get(i.filename,parts[i.filename]))
with ZipFile(out) as z:
    assert z.testzip() is None
    for k,v in parts.items():
        if k not in modified:assert z.read(k)==v,k
report={'artifact_authored_edits':True,'native_filters_retained':True,'untouched_package_parts_byte_identical':len(parts)-len(modified),'changed_parts':list(modified),'sheets':log}
(O/'native_preservation_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
print(json.dumps(report,ensure_ascii=False)[:2500])
