from common import *
from PIL import Image,ImageDraw
import shutil,csv
base=Path('/Users/taki/Desktop/Panvascular_Figure1b_Complete_Labels_20260920/FIGURE_VERSIONS.tsv');versions=pd.read_csv(base,sep='\t')
updated={'Figure1','Figure3','ExtendedData_Figure2','ExtendedData_Figure3','ExtendedData_Figure4','ExtendedData_Figure5'}
rows=[];merged=fitz.open();checks=[]
for _,r in versions.iterrows():
 name=r.figure;old=Path(r.source_pdf);pdf=O/'figures'/f'{name}.pdf' if name in updated else old;svg=O/'editable_svg'/f'{name}.svg' if name in updated else Path(r.source_svg)
 if name not in updated:
  shutil.copy2(pdf,O/'figures'/f'{name}.pdf');shutil.copy2(svg,O/'editable_svg'/f'{name}.svg');pdf=O/'figures'/f'{name}.pdf';svg=O/'editable_svg'/f'{name}.svg'
 d=fitz.open(pdf);assert len(d)==1
 merged.insert_pdf(d);typ,value=d.xref_get_key(d[0].xref,'Group')
 if typ=='dict':assert not re.search(r'\d+ \d+ R',value);merged.xref_set_key(merged[-1].xref,'Group',value)
 rows.append(dict(page=int(r.page),figure=name,source_version=O.name if name in updated else r.adopted_version,status='AUTHOR_REVIEW_PENDING' if name in updated else r.adoption_scope,source_pdf=str(pdf if name in updated else old),source_svg=str(svg if name in updated else Path(r.source_svg)),delivered_pdf=str(pdf),delivered_svg=str(svg),width_mm=d[0].rect.width/mm,height_mm=d[0].rect.height/mm,sha256=sha(pdf),svg_sha256=sha(svg)))
 if name in updated:
  pics=[]
  for f in [old,pdf]:
   pg=fitz.open(f)[0];pm=pg.get_pixmap(matrix=fitz.Matrix(1.5,1.5));pics.append(Image.frombytes('RGB',(pm.width,pm.height),pm.samples))
  cn=Image.new('RGB',(pics[0].width*2+18,max(p.height for p in pics)+25),'white');cn.paste(pics[0],(0,25));cn.paste(pics[1],(pics[0].width+18,25));dr=ImageDraw.Draw(cn);dr.text((5,5),'BEFORE',fill='black');dr.text((pics[0].width+23,5),'AFTER',fill='black');cn.save(O/'review'/f'{name}_before_after.png')
merged.save(O/'全部图件_15页_标注Marker与富集更新.pdf',garbage=4,deflate=True);merged.close()
m=fitz.open(O/'全部图件_15页_标注Marker与富集更新.pdf');assert len(m)==15
for row,p in zip(rows,m):
 src=fitz.open(row['delivered_pdf'])[0];a=src.get_pixmap(matrix=fitz.Matrix(1,1)).samples;b=p.get_pixmap(matrix=fitz.Matrix(1,1)).samples;assert a==b,row['figure'];assert p.rect==src.rect
 checks.append({'figure':row['figure'],'native_page_and_merged_render':'IDENTICAL','size_mm':[round(p.rect.width/mm,3),round(p.rect.height/mm,3)]})
pd.DataFrame(rows).to_csv(O/'FIGURE_VERSIONS.tsv',sep='\t',index=False)
# Marker value verification against the copied complete source matrices.
dis=pd.read_csv(O/'source_data/display_marker_values.tsv',sep='\t');sources=[]
l1=pd.read_csv(O/'source_data/ED2_L1_complete_matrix.tsv',sep='\t').rename(columns={'level1_annotation':'state_id'});sources.append(l1[['state_id','gene_id','mean_log1p','detection_fraction']])
for f in (O/'source_data').glob('*full_matrix.tsv'):
 z=pd.read_csv(f,sep='\t').rename(columns={'level2_annotation':'state_id'});sources.append(z[['state_id','gene_id','mean_log1p','detection_fraction']])
ref=pd.concat(sources);assert not ref.duplicated(['state_id','gene_id']).any();j=dis.merge(ref,on=['state_id','gene_id'],suffixes=['_display','_original'],validate='one_to_one');assert len(j)==len(dis)
for c in ['mean_log1p','detection_fraction']:assert np.allclose(j[c+'_display'],j[c+'_original'],equal_nan=True,rtol=0,atol=1e-12)
assert (dis.loc[dis.detection_fraction.eq(0),'area_pt2']==0).all()
for st,z in dis.groupby('state_id'):
 q=z[z.detection_fraction>0];assert np.allclose(q.area_pt2/q.detection_fraction,(q.area_pt2/q.detection_fraction).iloc[0])
# Verify unchanged frozen atlas maps and untouched Figure3 evidence objects.
for n in [2,3,4,5]:
 old=ET.parse(B/'editable_svg'/f'ExtendedData_Figure{n}.svg').getroot();new=ET.parse(O/'editable_svg'/f'ExtendedData_Figure{n}.svg').getroot()
 aa=json.load(open('/Users/taki/Desktop/Panvascular_Content_Alignment_ST_Audit_20260919/audit/alignment_panels.json'))
 for region in aa:
  if region['figure']==f'ExtendedData_Figure{n}' and region['panel'] in (['a','c'] if n==4 else ['a']):
   for key in region['members']:
    e1=next(e for e in old.iter() if e.get('id')==key);e2=next(e for e in new.iter() if e.get('id')==key);assert ET.tostring(e1)==ET.tostring(e2),(n,key)
old=ET.parse(B/'editable_svg/Figure3.svg').getroot();new=ET.parse(O/'editable_svg/Figure3.svg').getroot()
for key in ['figure3_b_content','axes_5','axes_6','axes_7','figure3_f_original_fill_key']:
 e1=next(e for e in old.iter() if e.get('id')==key);e2=next(e for e in new.iter() if e.get('id')==key);assert ET.tostring(e1)==ET.tostring(e2),key
r5=next(r for r in rows if r['figure']=='Figure5');assert r5['sha256']==versions.loc[versions.figure.eq('Figure5'),'pdf_sha256'].iloc[0]
# Communication review 3 pages, no new formal ED number.
book=fitz.open()
for f in ['Communication_chord_all176','Communication_ligand_receptor_detail','Communication_parent_counts']:
 d=fitz.open(O/'companion'/f'{f}.pdf');book.insert_pdf(d)
book.save(O/'细胞通讯_弦图与配体受体.pdf');book.close()
# Exact marker comparison pages.
book=fitz.open()
for no in [2,3,4,5]:
 name=f'ExtendedData_Figure{no}';before=fitz.open(B/'extended_data'/f'{name}.pdf');after=fitz.open(O/'figures'/f'{name}.pdf');page=book.new_page(width=before[0].rect.width*2,height=before[0].rect.height);page.show_pdf_page(fitz.Rect(0,0,before[0].rect.width,before[0].rect.height),before,0);page.show_pdf_page(fitz.Rect(before[0].rect.width,0,before[0].rect.width*2,before[0].rect.height),after,0)
book.save(O/'review/Marker_before_after_same_size.pdf')
# Gene/source index + dimensions of all companion pages.
pd.DataFrame([dict(file=str(f.relative_to(O)),sha256=sha(f),bytes=f.stat().st_size) for f in (O/'source_data').iterdir()]).to_csv(O/'SOURCE_DATA_INDEX.tsv',sep='\t',index=False)
dims=[]
for f in (O/'companion').glob('*.pdf'):
 for i,p in enumerate(fitz.open(f)):dims.append({'file':f.name,'page':i+1,'width_mm':p.rect.width/mm,'height_mm':p.rect.height/mm,'role':'review companion; not formally numbered main/ED'})
pd.DataFrame(dims).to_csv(O/'notes/companion_page_dimensions.tsv',sep='\t',index=False)
(O/'notes/VALIDATION.json').write_text(json.dumps({'data_integrity':{'displayed_marker_cells':len(dis),'exact_original_values':'PASS','zero_area':'PASS','full_source_matrices_retained':'PASS','original_maps_untouched':'PASS','Figure3_b_d_e_f_unchanged':'PASS','adopted_Figure5_byte_identical':'PASS','no_statistical_reanalysis':True},'layout_legibility':{'dot_max_diameter_rule':'0.73 × min(row pitch, column pitch)','all_15_native_sizes_and_merged_render':checks,'Figure1b':'16 full cell names with leader lines; no numeric key'},'visual_review':'Rendered Figure1, Figure3, ED2-5, chord, gene detail and enrichment companions inspected; author acceptance pending','ORA_statistical_validity':'PENDING'},indent=2))
print('PASS',len(rows),'pages;',len(dis),'marker entries; Figure5 unchanged')
