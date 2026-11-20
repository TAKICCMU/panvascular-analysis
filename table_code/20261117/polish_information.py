from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from lxml import etree as E
import json,csv,hashlib,os

W=Path(__file__).resolve().parents[1]
P=W/'staged/Supplementary_Tables.xlsx'
S='http://schemas.openxmlformats.org/spreadsheetml/2006/main';R='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
N={'s':S};q=lambda t:'{'+S+'}'+t
def val(c):return ''.join(c.itertext())
def setstr(c,t):
    for e in list(c):c.remove(e)
    c.set('t','inlineStr');E.SubElement(E.SubElement(c,q('is')),q('t')).text=t

replacements={
    (2,'C'):'Technical-library, specimen and donor crosswalk with clinical fields and sample grouping.',
    (3,'C'):'Cell counts after quality control and the corresponding sample, specimen and donor identifiers.',
    (17,'C'):'Lead-variant annotations and evidence for novelty relative to the component diseases.',
    (21,'C'):'MAGMA GSEA and ORA results. ORA includes pathways with overlap ≥3 before testing, with BH correction within each database and selected-gene definition.',
    (21,'D'):'One pathway and enrichment analysis, with fields specific to GSEA or ORA.',
    (22,'C'):'Analysis index and significant GSEA/ORA records. Complete results are in Supplementary_Table22_full_results/.',
    (26,'C'):'All 15,639 evaluable gene–state pairs with REML effects and target-set FDR.',
    (27,'C'):'FUMA eQTL mapping and colocalisation summaries across primary and leave-CAD scenarios.',
    (32,'C'):'Structural criteria and contextual evidence used to classify module confidence.',
    (35,'C'):'The 49 spatial sections, their 35 analysis units, source accessions and inclusion metadata.',
    (39,'C'):'Dataset medians after section-to-unit pooling; P-value columns contain descriptive medians of the unit-level P values.',
    (42,'C'):'Colocalisation models allowing multiple signals for the gene–QTL contexts listed in the table, with fit status.',
    (42,'D'):'One tested gene–QTL context and sensitivity fit.',
    (43,'C'):'Disease-associated cell-proportion estimates, uncertainty and inference status.',
    (45,'C'):'NMF quality metrics and automatic rank-selection criteria for ranks 2–25.',
    (49,'C'):'Numerical diagnostics and model output availability by spatial dataset.',
    (50,'C'):'Correspondence between spatial data sources and single-cell reference sources.',
    (52,'C'):'Within-unit pooled coefficients and standard errors used in the candidate meta-analysis.',
    (55,'C'):'Differential-expression results from the source study used to define plaque and calcification expression programmes.',
    (60,'C'):'Software and resource versions and the workflows to which they apply.',
}
with ZipFile(P) as z:parts={n:z.read(n) for n in z.namelist()}
before=parts.copy();before_sha=hashlib.sha256(P.read_bytes()).hexdigest()
tar={e.get('Id'):e.get('Target').lstrip('/') for e in E.fromstring(parts['xl/_rels/workbook.xml.rels'])}
w=E.fromstring(parts['xl/workbook.xml']);node=next(s for s in w.find(q('sheets')) if s.get('name')=='Information');path=tar[node.get('{'+R+'}id')]
r=E.fromstring(parts[path]);changes=[]
for row in r.findall('s:sheetData/s:row',N)[1:]:
    num=int(float(val(row[0])))
    for c in row:
        col=''.join(s for s in c.get('r') if s.isalpha())
        if (num,col) in replacements:
            new=replacements[num,col];old=val(c)
            if new==old:continue
            changes.append({'final_table':num,'cell':c.get('r'),'field':'Contents' if col=='C' else 'Row unit','before':old,'after':new})
            setstr(c,new)
assert len(changes)==len(replacements)
parts[path]=E.tostring(r,encoding='UTF-8')
tmp=P.with_suffix('.tmp.xlsx')
with ZipFile(tmp,'w',ZIP_DEFLATED,compresslevel=6) as z:
    for n,b in parts.items():z.writestr(n,b)
os.replace(tmp,P)
with ZipFile(P) as z:
    assert z.testzip() is None
    diff=[n for n in z.namelist() if z.read(n)!=before[n]]
assert diff==[path]
with (W/'tables/Information_wording_changes.tsv').open('w') as f:
    wr=csv.DictWriter(f,fieldnames=list(changes[0]),delimiter='\t');wr.writeheader();wr.writerows(changes)
(W/'tables/Information_wording_changes.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2))
vpath=W/'staged/renumber_validation.json';v=json.loads(vpath.read_text());v['output_sha256']=hashlib.sha256(P.read_bytes()).hexdigest();v['information_prose_changes']=len(changes);v['information_only_xml_changed']=True;v['information_prose_before_sha256']=before_sha
vpath.write_text(json.dumps(v,ensure_ascii=False,indent=2))
print(json.dumps({'changed_cells':len(changes),'changed_zip_parts':diff,'sha256':v['output_sha256']},ensure_ascii=False))
