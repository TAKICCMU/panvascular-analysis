from pathlib import Path
import json, re, csv, math
import pandas as pd
import openpyxl
from docx import Document

P = Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
O = P/'图件编辑资料/补充表核查_20261118'
OLD = Path('/Users/taki/Desktop/投稿包_最新版_20260919/图件编辑资料/全流程代码审核_20260922/Panvascular_code_review/analysis_code/native/downstream_final')
W = openpyxl.load_workbook(P/'Supplementary_Tables.xlsx', read_only=True, data_only=False)

def table(n):
    a = list(W[f'Table {n:02d}'].values)
    return pd.DataFrame([r for r in a[2:] if any(v is not None for v in r)], columns=a[1])

def dump(name, data):
    (O/name).write_text(json.dumps(data, ensure_ascii=False, indent=2, default=str))

# Compare saved Table 8 values to both saved weighting outputs, without recomputing either.
t8 = table(8).set_index('Tissue')
donor = pd.read_csv(OLD/'descriptive_updated_v2/L1_tissue_donor_mean_fraction.tsv', sep='\t').set_index('tissue')
sample = pd.read_csv(P/'Source_Data_Figures/Figure1/L1_tissue_sample_mean_fraction.tsv', sep='\t').set_index('tissue')
rows=[]
for ri, (tissue, row) in enumerate(t8.iterrows(),3):
    for ci, field in enumerate(t8.columns,2):
        a, b, c = row[field], donor.loc[tissue,field], sample.loc[tissue,field]
        rows.append(dict(cell=f'{openpyxl.utils.get_column_letter(ci)}{ri}',tissue=tissue,cell_type=field,workbook=a,saved_donor_average=b,saved_sample_average=c,difference_from_sample_percentage_points=100*(a-c),matches_donor=math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-14),matches_sample=math.isclose(a,c,rel_tol=1e-10,abs_tol=1e-14)))
pd.DataFrame(rows).to_csv(O/'Table08_averaging_comparison.tsv',sep='\t',index=False)

t9=table(9)
res9=[]
for kind in ['sample','donor']:
    native=pd.read_csv(OLD/f'descriptive_updated_v2/L1_study_{kind}_mean_fraction.tsv',sep='\t').set_index('study')
    bad=[]
    for i,row in t9.iterrows():
        a=row['Mean fraction across samples']; b=native.loc[row['Study ID'],row['L1 lineage']]
        if not math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-14):
            bad.append({'excel_row':i+3,'study':row['Study ID'],'L1':row['L1 lineage'],'workbook':a,'source':b})
    res9.append({'saved_average':kind,'rows':len(t9),'mismatch_rows':len(bad),'first_examples':bad[:3]})

# Locate textual references in the actual stored cells, excluding first-row titles.
refs=[]
dictionary=[]
numeric=[]
for ws in W:
    for row in ws.iter_rows():
        for cell in row:
            v=cell.value
            if not isinstance(v,str): continue
            if cell.row>1 and re.search(r'(?:Table\s*\d{1,2}\b|Figure\s*[1-9][a-h]?\b|Fig\.\s*[1-9][a-h]?\b)',v):
                refs.append({'sheet':ws.title,'cell':cell.coordinate,'value':v})
            if ws.title=='Data dictionary' and ('min–max' in v or 'min-max' in v or 'Record type' in v or 'Display overlap' in v):
                dictionary.append({'cell':cell.coordinate,'value':v})
            if (ws.title=='Table 17' and cell.column==3 and cell.row>=3) or (ws.title=='Table 18' and cell.column in [5,7,8,15] and cell.row>=3):
                if re.fullmatch(r'[+-]?\d+(?:\.\d+)?',v.strip()):
                    numeric.append({'sheet':ws.title,'cell':cell.coordinate,'stored_value':v,'stored_type':cell.data_type,'number_format':cell.number_format})
dump('workbook_textual_references.json',refs)
dump('numeric_text_cells.json',numeric)

# Confirm Supplementary Information's 60-table index using actual DOCX table rows.
si=Document(P/'Supplementary_Information.docx')
si_tables=[]
for ti,tab in enumerate(si.tables):
    for ri,row in enumerate(tab.rows):
        texts=[c.text for c in row.cells]
        if texts and re.fullmatch(r'\d{1,2}',texts[0].strip()):
            si_tables.append({'doc_table':ti+1,'row':ri+1,'number':int(texts[0]),'cells':texts})
info=list(W['Information'].values)
index={int(r[0]):r for r in info[1:] if isinstance(r[0],int)}
index_errors=[]
for r in si_tables:
    n=r['number']
    if n not in index: index_errors.append({'table':n,'reason':'not in workbook'})
    elif r['cells'][1] != index[n][1]: index_errors.append({'table':n,'reason':'title differs','doc':r['cells'][1],'workbook':index[n][1]})

docrefs=json.loads((O/'document_table_references.json').read_text())
invalid=[]
for file, paras in docrefs.items():
    for para in paras:
        for m in re.finditer(r'Supplementary\s+Tables?\s+(\d+)(?:\s*[–—-]\s*(\d+))?',para['text']):
            nums=[int(x) for x in m.groups() if x]
            if any(n<1 or n>60 for n in nums): invalid.append({'file':file,'paragraph':para['paragraph'],'reference':m.group(0)})

dump('final_focused_checks.json',{
 'Table8':{'cells':len(rows),'matches_donor':all(r['matches_donor'] for r in rows),'cells_differing_from_sample':sum(not r['matches_sample'] for r in rows),'max_absolute_difference_percentage_points':max(abs(r['difference_from_sample_percentage_points']) for r in rows)},
 'Table9':res9,
 'dictionary_reference_cells':dictionary,
 'quantitative_numbers_stored_as_text':len(numeric),
 'Supplementary_Information_catalog_rows':len(si_tables),
 'Supplementary_Information_catalog_title_errors':index_errors,
 'out_of_range_current_document_table_refs':invalid,
 'status':'READ_ONLY_AUDIT; no statistical results or workbook values modified'
})
W.close()
print((O/'final_focused_checks.json').read_text())
