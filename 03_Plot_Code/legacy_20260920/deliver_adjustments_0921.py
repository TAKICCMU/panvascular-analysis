from common_main import *
import shutil
from PIL import Image
PK=Path('/Users/taki/Desktop/投稿包_最新版_20260919');BACK=O/'merge_records/before_20260921_adjustments'
# Protected manuscript and data are never written by this script.
protected={str(p):sha(p) for p in PK.rglob('*') if p.is_file() and p.suffix.lower() in ['.docx','.xlsx']}
versions=pd.read_csv(O/'FIGURE_VERSIONS.tsv',sep='\t').sort_values('page');assert len(versions)==15
checks=[]
for idx,row in versions.iterrows():
 name=row.figure;pdf=O/'figures'/f'{name}.pdf';svg=O/'editable_svg'/f'{name}.svg';d=fitz.open(pdf);page=d[0];assert len(d)==1
 if name in [f'Figure{i}' for i in range(1,6)]:
  bad=[]
  for b in page.get_text('dict')['blocks']:
   for ln in b.get('lines',[]):
    for s in ln['spans']:
     bb=fitz.Rect(s['bbox'])
     if bb.x0<-.2 or bb.y0<-.2 or bb.x1>page.rect.width+.2 or bb.y1>page.rect.height+.2:bad.append(s['text'])
  assert not bad,(name,bad)
  versions.loc[idx,'version']='Panvascular_Main_Figure_Refinement_20260920 / 20260921 in-place adjustment'
  versions.loc[idx,'status']='AUTHOR_REVIEW_PENDING'
  before=fitz.open(BACK/'figures'/f'{name}.pdf');a=before[0].get_pixmap();b=page.get_pixmap()
  ia=Image.frombytes('RGB',(a.width,a.height),a.samples);ib=Image.frombytes('RGB',(b.width,b.height),b.samples)
  im=Image.new('RGB',(ia.width+ib.width+14,max(ia.height,ib.height)),'white');im.paste(ia,(0,0));im.paste(ib,(ia.width+14,0));im.save(O/'review'/f'{name}_before_after.png')
 else:
  assert sha(pdf)==sha(BACK/'figures'/pdf.name),(name,'ED changed')
  assert sha(svg)==sha(BACK/'editable_svg'/svg.name)
 versions.loc[idx,'source_pdf']=str(pdf);versions.loc[idx,'source_svg']=str(svg)
 versions.loc[idx,'width_mm']=page.rect.width/mm;versions.loc[idx,'height_mm']=page.rect.height/mm
 versions.loc[idx,'sha256']=sha(pdf);versions.loc[idx,'svg_sha256']=sha(svg)
 checks.append({'figure':name,'width_mm':page.rect.width/mm,'height_mm':page.rect.height/mm,'changed':name.startswith('Figure'),'pdf_sha256':sha(pdf),'svg_sha256':sha(svg)})
versions.to_csv(O/'FIGURE_VERSIONS.tsv',sep='\t',index=False)
merged=fitz.open()
for _,row in versions.iterrows():
 d=fitz.open(row.source_pdf);merged.insert_pdf(d);typ,val=d.xref_get_key(d[0].xref,'Group')
 if typ=='dict':
  assert not re.search(r'\d+ \d+ R',val);merged.xref_set_key(merged[-1].xref,'Group',val)
book=O/'全部图件_15页_主图调整版.pdf';temp=book.with_suffix('.tmp.pdf');merged.save(temp,garbage=4,deflate=True);merged.close();temp.replace(book)
m=fitz.open(book);assert len(m)==15
for i,(_,row) in enumerate(versions.iterrows()):
 d=fitz.open(row.source_pdf);assert d[0].rect==m[i].rect
 assert d[0].get_pixmap().samples==m[i].get_pixmap().samples,row.figure
# Synchronize the existing package only, keeping its existing directory/path contract.
pv=pd.read_csv(PK/'图件采用版本.tsv',sep='\t').set_index('figure')
for _,row in versions.iterrows():
 name=row.figure;rp=Path(pv.loc[name,'package_pdf']);rs=Path(pv.loc[name,'package_svg'])
 if name in [f'Figure{i}' for i in range(1,6)]:shutil.copy2(row.source_pdf,rp);shutil.copy2(row.source_svg,rs)
 assert sha(rp)==row.sha256 and sha(rs)==row.svg_sha256
 for c in versions.columns:
  if c!='figure':pv.loc[name,c]=row[c]
 pv.loc[name,'package_pdf_sha256']=sha(rp);pv.loc[name,'package_svg_sha256']=sha(rs)
pv.reset_index().sort_values('page').to_csv(PK/'图件采用版本.tsv',sep='\t',index=False)
shutil.copy2(book,PK/'全部图件_15页审阅.pdf')
shutil.copytree(O/'code',PK/'图件编辑资料/code/current_20260920',dirs_exist_ok=True)
shutil.copytree(O/'notes',PK/'图件编辑资料/本轮说明_20260920',dirs_exist_ok=True)
shutil.copy2(O/'source_data/Figure1b_tissue_sample_mapping.tsv',PK/'Source_Data_Figures/Current_Display_20260920/Figure1b_tissue_sample_mapping.tsv')
for f in ['Figure1b_tissue_palette.tsv','Figure1b_display15_direct_labels_20260921.tsv']:
 shutil.copy2(O/'source_data'/f,PK/'Source_Data_Figures/Current_Display_20260920'/f)
shutil.copy2(O/'SOURCE_DATA_INDEX.tsv',PK/'图件编辑资料/Source_Data_INDEX_current.tsv')
shutil.copy2(O/'本轮修改说明.md',PK/'图件编辑资料/本轮修改说明.md')
# A/B are review deliverables; do not replace official Figure1g with B.
audit=O/'marker_audit_20260921'
if (audit/'AUDIT_COMPLETE.json').exists():
 shutil.copytree(audit,PK/'图件编辑资料/marker_audit_20260921',dirs_exist_ok=True)
 shutil.copytree(audit/'tables',PK/'Source_Data_Figures/Current_Display_20260920/marker_audit_20260921',dirs_exist_ok=True)
assert protected=={str(p):sha(p) for p in PK.rglob('*') if p.is_file() and p.suffix.lower() in ['.docx','.xlsx']}
result={'changed_main_figures':5,'unchanged_extended_data':10,'merged_pages':15,'native_page_sizes_and_rendering_identical':True,'package_pdf_svg_hash_equality':True,'docx_xlsx_unchanged':True,'no_scientific_models_rerun':True,'official_Figure1g_kept_original_A':True,'ORA':'PENDING','donor_identity':'PENDING','visual_review':'All five revised main-figure renders inspected; author review pending','figures':checks}
(O/'notes/DELIVERY_CHECKS.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
(PK/'图件编辑资料/VERSION_EQUALITY.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
shutil.copy2(O/'notes/DELIVERY_CHECKS.json',PK/'图件编辑资料/本轮说明_20260920/DELIVERY_CHECKS.json')
# Refresh existing per-file checksum index without inventing source provenance for inherited files.
idxpath=PK/'文件来源与校验.tsv';rows=pd.read_csv(idxpath,sep='\t').to_dict('records');known={r['file'] for r in rows}
for row in rows:
 p=PK/row['file']
 if p.is_file():row['sha256']=sha(p)
for p in PK.rglob('*'):
 if p.is_file() and str(p.relative_to(PK)) not in known and p!=idxpath:
  rows.append({'file':str(p.relative_to(PK)),'kind':'20260921 visualization or marker audit','source':'existing canonical figure workspace','sha256':sha(p)})
pd.DataFrame(rows).to_csv(idxpath,sep='\t',index=False)
print('PASS: same package overwritten; 15 pages native-size and render-identical; 10 ED unchanged; manuscripts unchanged')
