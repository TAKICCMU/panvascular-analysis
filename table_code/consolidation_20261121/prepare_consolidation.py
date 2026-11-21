from pathlib import Path
import json,csv,collections,hashlib,datetime,gzip
from openpyxl import load_workbook
O=Path(__file__).parent; P=O.parent.parent
D=json.loads((O/'SOURCE_TABLES_DATA.json').read_text() if (O/'SOURCE_TABLES_DATA.json').exists() else gzip.decompress((O/'SOURCE_TABLES_DATA.json.gz').read_bytes()).decode()); checks=json.loads((O/'REDUNDANCY_CHECKS.json').read_text())
# Verified wide-to-long candidate evidence; no recalculation.
h=D['24']['h']; cand={tuple(r[h.index(x)] for x in ['Gene','Annotation level','Cell state']):dict(zip(h,r)) for r in D['24']['r']}
for r in D['25']['r']:
 c=cand[tuple(r[:3])]; scenario={'primary_r0.1':'primary r0.1','primary_r0.6':'primary r0.6','leave_one_CAD_r0.1':'leave one CAD r0.1'}[r[3]]
 assert [r[4],r[5],r[6]]==[c['P '+scenario],c['FDR BH '+scenario],c['Bonferroni 0.05 '+scenario]],r[:4]
checks['25_in_24']={'identical_values':True,'rows':87}
# Join one Unknown count onto the same sample metadata, proving the common fields.
base=D['2']; unknown={r[D['10']['h'].index('Sample ID')]:dict(zip(D['10']['h'],r)) for r in D['10']['r']}
assert len(unknown)==len(base['r'])==247
for r in base['r']:
 a=dict(zip(base['h'],r)); b=unknown[a['Sample ID']]
 for k,v in a.items():assert v==b[{'Cells (n)':'Total sample cells (n)'}.get(k,k)],(a['Sample ID'],k)
 assert b['L1 category']==b['State']=='Unknown'
merged2={'h':base['h']+['Unknown cells (n)'],'r':[r+[unknown[r[0]]['Unknown cells (n)']] for r in base['r']]}
# Join the same 64 phenotype pairs without changing correlation or intercept.
inter={tuple(r[:2]):r[2] for r in D['41']['r']}; assert len(inter)==64
merged8={'h':D['12']['h']+['LDSC intercept'],'r':[r+[inter[tuple(r[:2])]] for r in D['12']['r']]}
checks['2_10_sample_join']={'all_common_fields_exact':True,'rows':247};checks['12_41_pair_join']={'all_pair_keys_match':True,'rows':64}
(O/'REDUNDANCY_CHECKS.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
# One final table per coherent analysis theme; detailed records remain separately deliverable.
groups=[[1],[2],[4],[5],[6],[8,9],[11],[12],[13,14,15,16],[17,18],[19,20],[21],[22,23],[24],[26],[28,29,42],[30,31],[45,33,34],[35,48,49,50],[36,37,51],[40],[43,44],[55,56],[57,58,59,60]]
titles='''Single cell transcriptomic studies and source publications
Sample and donor clinical metadata, quality control cell counts and Unknown cell coverage
Single cell integration benchmark results
Complete marker gene statistics
Cell state annotations and study coverage
Donor mean cell proportions by tissue and within study sample mean cell proportions
GWAS data sources
Pairwise LDSC genetic correlations and intercepts
Exploratory factor analysis loadings, candidate Genomic SEM fits, bifactor loadings and leave one trait correlations
Novel shared loci and literature or GWAS Catalog correspondence
MAGMA and scDRS cell population enrichment results and summaries
MAGMA pathway enrichment from GSEA and original ORA
Pseudobulk pathway GSEA and ORA results index, significant results and cross contrast recurrence
Complete evidence for the 29 candidate gene–cell population pairs
Complete pseudobulk meta-analysis results across all 15,639 pairs
eQTL coverage, colocalisation and multi-signal sensitivity analyses
CoVarNet module evidence and node weights
CoVarNet rank selection, matched null comparisons, resampling and leave one study refits
Spatial section metadata, regional annotation availability, cell2location diagnostics and reference overlap
Primary spatial association results for candidate genes, four sensitivity settings and target exclusion controls
Spatial coherence of CoVarNet modules
Disease cell proportion results and model eligibility
Plaque and calcification programme marker data and module–programme directional associations
Analysis units, evidence source relationships, multiple testing families, software versions and random seeds'''.splitlines()
# Native formats/definitions retained as source metadata for the artifact-tool writer.
w=load_workbook(O/'60表基线备份/Supplementary_Tables.xlsx',data_only=False)
fmt={}
for k,x in D.items():
 s=w[f'Table {int(k):02d}'];valid=[row for row in s.iter_rows(min_row=3,max_row=s.max_row,max_col=len(x['h'])) if any(c.value is not None for c in row)]
 assert [[c.value for c in row] for row in valid]==x['r'],k
 fmt[k]={'widths':[s.column_dimensions[c.column_letter].width for c in s[2]],'formats':[s.cell(3,j+1).number_format for j in range(len(x['h']))], 'row_height':s.row_dimensions[3].height or 24,'cell_formats':[[c.number_format for c in row] for row in valid]}
info=[list(r) for r in w['Information'].values]; infoMap={r[0]:r for r in info[1:]}
dictionary=[list(r) for r in w['Data dictionary'].values];abbr=[list(r) for r in w['Abbreviations'].values]
out=[];mapping={};positions=[]
for new,(olds,title) in enumerate(zip(groups,titles),1):
 blocks=[]
 # Do not duplicate 87 primary rows: four-setting table includes them exactly.
 actual=[x for x in olds if x!=36]
 for old in actual:
  x=D[str(old)]; content=merged2 if old==2 else merged8 if old==12 else x
  block=chr(65+len(blocks)) if len(actual)>1 else ''
  b={'old':old,'block':block,'title':x['title'].split('|',1)[-1].strip(),'h':content['h'],'r':content['r'],'formats':fmt[str(old)]['formats'],'widths':fmt[str(old)]['widths'],'row_height':fmt[str(old)]['row_height'],'row_unit':infoMap[old][3],'cell_formats':fmt[str(old)]['cell_formats']}
  if old==2:b['formats']+=['#,##0'];b['widths']+=[17];b['cell_formats']=[r+['#,##0'] for r in b['cell_formats']]
  if old==12:b['formats']+=[fmt['41']['formats'][2]];b['widths']+=[16];b['cell_formats']=[r+[fmt['41']['formats'][2]] for r in b['cell_formats']]
  if old==8:b['title']='Cell-type fractions by tissue (donor means)'
  if old==9:b['title']='Cell-type composition by study (sample means)'
  if old==37:b['title']='Candidate spatial associations across the four posterior and QC settings'
  blocks.append(b);mapping[str(old)]={'file':'Supplementary_Tables.xlsx','sheet':f'Table {new:02d}','new_table':new,'block':block,'action':'retained' if len(actual)==1 else 'merged as distinct block'}
 out.append({'number':new,'sheet':f'Table {new:02d}','title':title,'blocks':blocks})
details={38:'Module_section_statistics',39:'Module_dataset_statistics',47:'CoVarNet_sample_metadata',52:'Candidate_unit_coefficients',53:'Candidate_section_coefficients',54:'Spatial_model_eligibility'}
sd=[]
for old,name in details.items():
 x=D[str(old)];sd.append({'number':None,'sheet':name,'title':x['title'].split('|',1)[-1].strip(),'blocks':[{'old':old,'block':'','title':x['title'].split('|',1)[-1].strip(),'h':x['h'],'r':x['r'],**fmt[str(old)],'row_unit':infoMap[old][3],'cell_formats':fmt[str(old)]['cell_formats']}]});mapping[str(old)]={'file':'Supplementary_Data.xlsx','sheet':name,'new_table':None,'block':'','action':'moved to detailed data workbook'}
for old,target,reason in [(3,2,'identical subset'),(7,6,'identical subset'),(10,2,'Unknown count joined by exact sample ID; common metadata identical'),(25,24,'exact wide-to-long equivalent'),(27,24,'identical subset'),(32,30,'identical subset'),(36,37,'exact primary subset: q05/all_sections'),(41,12,'LDSC intercept joined by exact phenotype pair'),(46,34,'identical subset')]:
 mapping[str(old)]={**mapping[str(target)],'action':reason,'represented_in_old_table':target}
assert set(mapping)==set(map(str,range(1,61)))
# Build field dictionaries for the new blocks, retaining original definitions and scales.
dictBy=collections.defaultdict(list)
for r in dictionary[1:]:dictBy[int(r[0])].append(r[1:])
for wbSheets in [out,sd]:
 for s in wbSheets:
  s['dictionary']=[]
  for b in s['blocks']:
   fields={r[0]:r for r in dictBy[b['old']]}
   if b['old']==2:fields['Unknown cells (n)']=next(r for r in dictBy[10] if r[0]=='Unknown cells (n)')
   if b['old']==12:fields['LDSC intercept']=next(r for r in dictBy[41] if r[0]=='LDSC intercept')
   for h in b['h']:
    assert h in fields,(b['old'],h)
    s['dictionary'].append([s['sheet'],b['block'],*fields[h]])
plan={'core':out,'data':sd,'abbreviations':abbr,'mapping':mapping,'titles':titles,'source_sha256':hashlib.sha256((O/'60表基线备份/Supplementary_Tables.xlsx').read_bytes()).hexdigest()}
(O/'CONSOLIDATION_PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,separators=(',',':')))
with (P/'Supplementary_Table_Crosswalk.tsv').open('w') as f:
 wr=csv.writer(f,delimiter='\t',lineterminator='\n');wr.writerow(['Original table','Original title','Current file','Current worksheet','Current table','Block','Disposition'])
 for old in range(1,61):m=mapping[str(old)];wr.writerow([old,D[str(old)]['title'].split('|',1)[-1].strip(),m['file'],m['sheet'],m['new_table'],m['block'],m['action']])
mem=Path('/Users/taki/.codex/memories/extensions/ad_hoc/notes');mem.mkdir(parents=True,exist_ok=True)
(mem/'20261121-1900-panvascular-60-table-baseline.md').write_text(f'''User explicitly requested saving and remembering the current 60-table baseline before simplification on 2026-11-21.\nCurrent package: {P}\nImmutable affected-file baseline: {O/'60表基线备份'}\nBaseline workbook SHA256: {plan['source_sha256']}\nExact 60-table inventory with original headers/counts: {O/'CURRENT_60_TABLE_INVENTORY.json'}\nDisposition map (verify latest completion before treating as adopted): {P/'Supplementary_Table_Crosswalk.tsv'}\nSimplification is editorial only; no analysis refit or scientific-value changes. Original Table 8 is donor-mean tissue composition, different from Figure 1e sample means. Preserve zero/NA and failed-model records. Formal final workbooks must be live-verified; this note records the preserved baseline, not a completion claim.\n''')
print('PLAN_READY',len(out),'numbered tables;',len(sd),'detailed sheets;',len(mapping),'original tables accounted for')
