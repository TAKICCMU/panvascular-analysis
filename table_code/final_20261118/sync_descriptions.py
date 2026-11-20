from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
from xml.parsers import expat
from xml.sax.saxutils import escape
import copy,hashlib,json

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
O=P/'图件编辑资料/补充表核查_20261118'; S=O/'staging'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
sha=lambda b:hashlib.sha256(b).hexdigest()
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def spans(raw):
    parser=expat.ParserCreate();out=[];stack=[]
    def start(n,a):
        if n=='w:t':
            i=parser.CurrentByteIndex;j=raw.index(b'>',i)+1;stack.append((i,j,raw[j-2:j]==b'/>'))
    def end(n):
        if n=='w:t':
            a,b,empty=stack.pop();c=parser.CurrentByteIndex;d=raw.index(b'>',c)+1
            out.append((a,b,b if empty else c,b if empty else d,empty))
    parser.StartElementHandler=start;parser.EndElementHandler=end;parser.Parse(raw,True)
    return sorted(out)
def replace_phrase(p,old,new):
    original=text(p);a=original.index(old);b=a+len(old);pos=0;nodes=[]
    for node in p.xpath('.//w:t',namespaces=NS):
        t=node.text or '';nodes.append((node,pos,pos+len(t)));pos+=len(t)
    hits=[x for x in nodes if x[1]<b and x[2]>a]
    for i,(node,start,end) in enumerate(hits):
        t=node.text or '';prefix=t[:max(0,a-start)] if i==0 else '';suffix=t[max(0,b-start):] if i==len(hits)-1 else ''
        node.text=prefix+(new if i==0 else '')+suffix
    assert text(p)==original[:a]+new+original[b:]

legend_old='e, Mean Level 1 cell fractions across samples within each tissue.'
legend_new='e, Mean Level 1 cell fractions across samples within each tissue. Supplementary Table 8 reports averages across donor identifiers.'
results_old='The atlas also captures variation in cell composition across vascular beds and the available age and sex metadata (Fig. 1e,f and Supplementary Table 8).'
results_new='The atlas also captures variation in cell composition across vascular beds, summarised across samples in Fig. 1e and across donor identifiers in Supplementary Table 8. Available age and sex metadata are shown in Fig. 1f.'
method_old='Study counts per Level 1 category are shown in Fig. 1d; recurrence of the 116 named states is reported in Supplementary Table 7, and coverage of all 128 combinations is provided in Source Data.'
method_new=method_old+' Cell fractions use all atlas cells in each sample as the denominator. Figure 1e averages these fractions across samples within each tissue. Supplementary Table 8 first averages them within each tissue and recorded donor identifier, then across donor identifiers.'
spec=json.loads((O/'correction_spec.json').read_text())
t8=next(x['after'] for x in spec['changes'] if x['sheet']=='Information' and x['cell']=='C9')
common=[(legend_old,legend_new)]
rules={
 'Manuscript.docx':common+[(results_old,results_new),(method_old,method_new),('Workflow-specific versions and random seeds are listed in Supplementary Table 60.','Workflow-specific software versions and recorded random seeds are listed in Supplementary Table 60.')],
 'Figure_Legends.docx':common,
 'Figure/Figure_Legends.docx':common,
 'Supplementary_Information.docx':[
  ('Sample-averaged L1 composition within tissue categories. One tissue category; columns are L1 fractions.',t8+' One tissue category; columns are mean L1 fractions (0–1).'),
  ('Prespecified test families and correction boundaries across analyses. One multiple-testing family definition.','Multiple-testing families and correction boundaries for the retained genetic, expression, composition, CoVarNet and spatial analyses. One family definition per row.'),
  ('Software versions','Software versions and random seeds'),
  ('Software and resource versions and the workflows to which they apply. One software or resource record.','Workflow-specific software versions, resources and random seeds from retained environment, parameter and execution records. One software, resource or seed record.')]
}
report=[]
for name,edits in rules.items():
    p=P/name;data=p.read_bytes()
    with ZipFile(p) as z:
        entries=z.infolist();parts={i.filename:z.read(i.filename) for i in entries};comment=z.comment
    raw=parts['word/document.xml'];root=E.fromstring(raw);oldroot=copy.deepcopy(root);changes=[]
    for old,new in edits:
        hits=[q for q in root.xpath('//w:p',namespaces=NS) if old in text(q)]
        assert len(hits)==1,(name,old,len(hits))
        q=hits[0];before=text(q);replace_phrase(q,old,new);changes.append({'old':old,'new':new})
    aa=oldroot.xpath('//w:t',namespaces=NS);bb=root.xpath('//w:t',namespaces=NS);ss=spans(raw)
    assert len(aa)==len(bb)==len(ss)
    replacements=[]
    for a,b,s in zip(aa,bb,ss):
        assert a.attrib==b.attrib
        if a.text!=b.text:
            assert not s[4];replacements.append((s[1],s[2],escape(b.text or '').encode()))
    newraw=raw
    for a,b,v in reversed(replacements):newraw=newraw[:a]+v+newraw[b:]
    parsed=E.fromstring(newraw)
    for tree in (parsed,root):
        for q in tree.xpath('//w:t',namespaces=NS):
            if not q.text:q.text=None
    assert E.tostring(parsed)==E.tostring(root)
    a=copy.deepcopy(oldroot);b=copy.deepcopy(parsed)
    for q in a.xpath('//w:t',namespaces=NS):q.text='__TEXT__'
    for q in b.xpath('//w:t',namespaces=NS):q.text='__TEXT__'
    assert E.tostring(a)==E.tostring(b),'formatting/field structure altered'
    out=S/name;out.parent.mkdir(parents=True,exist_ok=True)
    with ZipFile(out,'w') as z:
        z.comment=comment
        for i in entries:z.writestr(copy.copy(i),newraw if i.filename=='word/document.xml' else parts[i.filename])
    with ZipFile(out) as z:
        assert z.testzip() is None
        for k,v in parts.items():
            if k!='word/document.xml':assert z.read(k)==v
    report.append({'file':name,'before_sha256':sha(data),'after_sha256':sha(out.read_bytes()),'changes':changes,'all_non_text_xml_and_citation_fields_preserved':True,'unchanged_other_parts':len(parts)-1})

p=P/'Figure/Figure_Legends.md';s=p.read_text();assert s.count(legend_old)==1
(S/'Figure/Figure_Legends.md').write_text(s.replace(legend_old,legend_new))

p=P/'Source_Data_Figures/ExtendedData1/benchmark_publication_CPU_only.tsv';lines=p.read_text().splitlines(keepends=True)
headers=lines[0].rstrip('\r\n').split('\t');ri=headers.index('rank');mi=headers.index('method')
import openpyxl
w=openpyxl.load_workbook(P/'Supplementary_Tables.xlsx',read_only=True,data_only=True)
tab=list(w['Table 04'].values);h=tab[1];rank_by={r[h.index('Method')]:r[h.index('Rank')] for r in tab[2:]}
updated=[lines[0]];ranking=[]
for line in lines[1:]:
    newline='\r\n' if line.endswith('\r\n') else '\n';fields=line.rstrip('\r\n').split('\t');old=fields[ri]
    fields[ri]=str(int(rank_by[fields[mi]]));updated.append('\t'.join(fields)+newline);ranking.append({'method':fields[mi],'before':old,'after':fields[ri]})
out=S/'Source_Data_Figures/ExtendedData1/benchmark_publication_CPU_only.tsv';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(''.join(updated))
for a,b in zip(lines[1:],updated[1:]):
    aa=a.rstrip('\r\n').split('\t');bb=b.rstrip('\r\n').split('\t');aa.pop(ri);bb.pop(ri);assert aa==bb
(O/'document_corrections.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
(O/'benchmark_rank_correction.json').write_text(json.dumps(ranking,ensure_ascii=False,indent=2))
print(json.dumps({'documents':len(report),'legend_markdown_synced':True,'rank_rows':len(ranking),'statistical_values_changed':False},ensure_ascii=False))
