from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
from collections import Counter
import json,csv,hashlib,shutil,re,gzip
from copy import copy
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/五项字典纠正_20260922';B=O/'轻量备份';B.mkdir(exist_ok=True);X=P/'Supplementary_Tables.xlsx';tmp=O/'Supplementary_Tables_corrected.xlsx';ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def backup(p):
 q=B/p.relative_to(P);q.parent.mkdir(parents=True,exist_ok=True)
 if p.exists() and not q.exists():shutil.copy2(p,q)
backup(X);w=load_workbook(X,read_only=True);dd=list(w['Data dictionary'].values);ddidx={(int(r[0]),r[1]):i for i,r in enumerate(dd[1:],2)}
fixes=[]
def fix(n,old,new,unit,definition,source):fixes.append(dict(table=n,old=old,new=new,unit=unit,definition=definition,source=source))
base='Public_Release/Analysis_Code/native/downstream_final/overnight_manuscript_integration_20260917/'
fix(13,'Primary FDR','Primary MAGMA FDR <0.05','boolean','Boolean indicator that the primary MAGMA BH-adjusted P value is below 0.05. Numeric FDR is retained in column Y.',base+'candidate_prioritization/build_candidates.py; source original_GWAS_master; AE versus Y checked')
fix(13,'Retained leave one CAD FDR','Leave-CAD MAGMA FDR <0.05','boolean','Boolean indicator that the leave-CAD MAGMA BH-adjusted P value is below 0.05. Numeric FDR is retained in column X. This is distinct from the nominal-P eligibility flag in AF.',base+'candidate_prioritization/build_candidates.py; source original_GWAS_master; AG versus X checked')
fix(13,'Strict three context shared','Meets three-context criteria','boolean','Boolean result of the original rule: n_disease_contexts >= 3 AND direction_consistency >= 0.75. This field does not itself impose additional GWAS or meta-FDR criteria.',base+'candidate_prioritization/build_candidates.py:22; original_05b.R:507')
fix(40,'Top shared variant','Top variant under H4','variant identifier','Variant identifier with the highest variant-level posterior conditional on H4 in the tested region. It is not itself a probability and does not imply that the region meets strong-colocalisation criteria.','Public_Release/Analysis_Code/retained_genetics/final_updates/02_run_full_coloc.R:126-138')
fix(50,'Shared evidence','Shared inputs','text','Inputs, participant sources or datasets shared between the named analyses.','Existing Table50 text records; values unchanged')
for old,new in [('GSEA FDR05','GSEA pathways, FDR<0.05 (n)'),('ORA FDR05 overlap3','ORA pathways, FDR<0.05 and overlap≥3 (n)'),('N GSEA','GSEA pathways (n)'),('N ORA','ORA pathways (n)'),('Native NA','Unavailable enrichment results (n)')]:
 d=dd[ddidx[(20,old)]-1];fix(20,old,new,'count',d[2],'Existing Table20 index-row definition and integer-valued counts; structural blanks retained')
fix(39,'Both measured','Interactions measured by both methods (n)','count','Number of interaction records matched between LIANA and CellChat on disease_context, threshold_scheme, interaction_key, source, target, ligand and receptor within this tissue/background.',base+'review_communication.py:86-92')
fix(39,'Same direction','Interactions with same effect direction (n)','count','Number of matched LIANA/CellChat interactions whose disease-minus-control pooled effects have equal signs. The source uses equality of numpy.sign, including equal zero signs if present.',base+'review_communication.py:86-92')
fix(39,'Both supported same direction','Interactions supported by both, same direction (n)','count','Number of matched interactions with equal pooled-effect signs and context_meta_supported TRUE for both LIANA and CellChat. Original method-specific support rules are retained.',base+'review_communication.py:86-92')
headers={};values={}
for n in sorted({f['table'] for f in fixes}):
 rows=list(w[f'Table {n:02}'].values);headers[n]=list(rows[2]);values[n]=rows[3:]
for r in values[13]:
 if not any(v is not None for v in r):continue
 assert isinstance(r[30],bool) and r[30]==(r[24]<.05)
 assert isinstance(r[32],bool) and r[32]==(r[23]<.05)
 assert r[34]==(r[3]>=3 and r[13]>=.75)
assert Counter(r[34] for r in values[13])=={False:27,True:2}
# Check hard type/scale conflicts across all current column-type exports, without treating vague categories as errors.
fmt=list(csv.DictReader((P/'Spreadsheet_Format_Check.tsv').open(),delimiter='\t'));types={(int(r['Table']),r['Actual header']):json.loads(r['Actual cell type distribution']) for r in fmt};v=json.loads((P/'图件编辑资料/三处文稿修正与Excel文本核查_20260922/excel_export_validation.json').read_text());assert sha(X)==v['workbook_sha256_before_after']
conflicts=[]
for r in dd[1:]:
 n=int(r[0]);typ=set(types[(n,r[1])]) - {'blank'};unit=str(r[3]);bad=(typ=={'boolean'} and ('probability' in unit or 'fraction' in unit or unit=='count')) or (typ=={'text'} and any(z in unit for z in ['probability','fraction']) )
 if bad:conflicts.append((n,r[1]))
assert set(conflicts)=={(f['table'],f['old']) for f in fixes[:5]},conflicts
# Edit explicit OOXML text cells only. All scientific data and styles remain byte-identical.
z=ZipFile(X);rels={e.get('Id'):e.get('Target').lstrip('/') for e in E.fromstring(z.read('xl/_rels/workbook.xml.rels'))};paths={e.get('name'):rels[e.get('{'+ns['r']+'}id')] for e in E.fromstring(z.read('xl/workbook.xml')).findall('s:sheets/s:sheet',ns)};paths={k:(v if v.startswith('xl/') else 'xl/'+v) for k,v in paths.items()};changes={};report=[]
def textpatch(sheet,cell,text):
 path=paths[sheet];b=changes.get(path,z.read(path));pat=rb'(<c\b[^>]*\br="'+cell.encode()+rb'"[^>]*>)(.*?)(</c>)';m=re.search(pat,b,re.S);assert m,(sheet,cell)
 c=E.fromstring(m.group(),E.XMLParser());assert c.get('t')=='inlineStr';t=c.find('is/t');assert t is not None and len(c.findall('is/t'))==1
 # Only replace the existing t contents; leave style IDs and all surrounding XML bytes unchanged.
 start=m.start()+m.group().index(b'<t>')+3 if b'<t>' in m.group() else None
 if start is None:
  tm=re.search(rb'<t\b[^>]*>',m.group());assert tm;start=m.start()+tm.end()
 end=b.index(b'</t>',start);safe=text.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;').encode();changes[path]=b[:start]+safe+b[end:]
for f in fixes:
 n=f['table'];sheet=f'Table {n:02}';j=headers[n].index(f['old'])+1;cell=get_column_letter(j)+'3';di=ddidx[(n,f['old'])];original=dd[di-1]
 nonblank=[r[j-1] for r in values[n] if r[j-1] is not None]
 if f['unit']=='boolean':assert all(isinstance(v,bool) for v in nonblank)
 elif f['unit']=='count':assert all(isinstance(v,(int,float)) and not isinstance(v,bool) and v>=0 and float(v).is_integer() for v in nonblank)
 else:assert all(isinstance(v,str) for v in nonblank)
 textpatch(sheet,cell,f['new']);textpatch('Data dictionary',f'B{di}',f['new']);textpatch('Data dictionary',f'C{di}',f['definition']);textpatch('Data dictionary',f'D{di}',f['unit'])
 report.append({'Table':n,'Worksheet':sheet,'Header cell':cell,'Dictionary row':di,'Old header':f['old'],'New header':f['new'],'Old definition':original[2],'New definition':f['definition'],'Old unit':original[3],'New unit':f['unit'],'Source':f['source']})
with ZipFile(tmp,'w') as out:
 for item in z.infolist():out.writestr(copy(item),changes.get(item.filename,z.read(item.filename)))
# Assert every numbered-table data row, including scientific text/NA/bool/formula, is byte-identical; all formatting unchanged.
with ZipFile(tmp) as new:
 for n in range(1,65):
  path=paths[f'Table {n:02}'];a=z.read(path);b=new.read(path)
  # Mask only authorised header cells: everything else in each sheet must match.
  for f in [r for r in report if r['Table']==n]:
   pat=rb'<c\b[^>]*\br="'+f['Header cell'].encode()+rb'"[^>]*>.*?</c>';a=re.sub(pat,b'HEADER',a,flags=re.S);b=re.sub(pat,b'HEADER',b,flags=re.S)
  assert a==b,n
 for item in z.infolist():
  if item.filename not in changes:assert z.read(item)==new.read(item.filename)
 # Per-cell style preserved, including changed headers and dictionary cells.
 for path in changes:
  a=E.fromstring(z.read(path));b=E.fromstring(new.read(path));assert [(c.get('r'),c.get('s'),c.get('t')) for c in a.findall('.//s:c',ns)]==[(c.get('r'),c.get('s'),c.get('t')) for c in b.findall('.//s:c',ns)]
w.close();z.close()
with (O/'Changed_dictionary_and_headers.tsv').open('w') as f:wr=csv.DictWriter(f,list(report[0]),delimiter='\t');wr.writeheader();wr.writerows(report)
# Mapping updates only: original machine fields remain unchanged.
mp=P/'column_name_map.tsv';backup(mp);maps=list(csv.DictReader(mp.open(),delimiter='\t'));fields=list(maps[0]);lookup={(f'Table {f["table"]:02}',f['old']):f for f in fixes}
for r in maps:
 f=lookup.get((r['Worksheet'],r['Display column']))
 if f:r['Display column']=f['new'];r['Definition']=f['definition'];r['Change reason']='Source-verified type/meaning correction; scientific values unchanged'
with mp.open('w') as f:wr=csv.DictWriter(f,fields,delimiter='\t');wr.writeheader();wr.writerows(maps)
shutil.copy2(tmp,X)
R=P/'Public_Release/Supplementary_Tables_TSV';backup(R/'column_name_map.tsv');shutil.copy2(mp,R/'column_name_map.tsv')
ww=load_workbook(X,read_only=True);backup(R/'Data dictionary.tsv.gz')
with gzip.open(R/'Data dictionary.tsv.gz','wt',newline='') as f:csv.writer(f,delimiter='\t').writerows(ww['Data dictionary'].values)
ww.close()
(O/'fixes.json').write_text(json.dumps(fixes,ensure_ascii=False,indent=2));(O/'validation.json').write_text(json.dumps({'numbered_tables_data_rows_byte_identical':64,'styles_freezes_filters_merges_unchanged':True,'header_changes':len(fixes),'dictionary_rows_changed':len(fixes),'AE_equals_Y_lt_0_05':True,'AG_equals_X_lt_0_05':True,'AI_original_rule_and_counts':{'FALSE':27,'TRUE':2},'additional_clear_type_conflicts':0,'xlsx_after_sha256':sha(X),'scientific_reanalysis':False},indent=2));print('Corrected 13 headers/dictionary rows; all 64 data regions and formatting byte-identical.')
