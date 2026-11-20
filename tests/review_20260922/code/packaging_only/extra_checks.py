from pathlib import Path
import json,hashlib,shutil
import pandas as pd,numpy as np,fitz
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/全流程代码审核_20260922';B=O/'Panvascular_code_review';s=B/'evidence/source_data';rows=[]
def add(n,ok,x):rows.append({'check':n,'status':'PASS' if ok else 'REVIEW_REQUIRED','detail':str(x)})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
# original release indexes vs actual current files
for root,index in [(P/'Public_Release/Analysis_Code','CODE_FILE_INDEX.tsv'),(P/'Source_Data_Figures','SOURCE_FILE_INDEX.tsv')]:
 rr=[]
 for r in pd.read_csv(root/index,sep='\t').itertuples():
  f=root/r.relative_path;ok=f.exists() and sha(f)==r.sha256
  rr.append({'index':index,'file':r.relative_path,'sha_match':ok})
 pd.DataFrame(rr).to_csv(B/'audit'/('INHERITED_'+index),sep='\t',index=False)
 add('inherited_'+index,all(x['sha_match'] for x in rr),{'files':len(rr),'mismatches':[r['file'] for r in rr if not r['sha_match']]})
coords=pd.read_csv(s/'Figure1/full_cells_frozen_umap.tsv.gz',sep='\t');m=json.loads((B/'analysis_code/native/atlas_finalization/PREPARED_RELEASE_MANIFEST.json').read_text())
add('all_L1_counts_match_release',coords.level1_annotation.value_counts().to_dict()==m['level1_counts_verified'],coords.level1_annotation.value_counts().to_dict())
add('cell_ids_unique',coords['_index'].is_unique,len(coords))
add('sample_count247',coords.sample_uid.nunique()==247,coords.sample_uid.nunique())
old=pd.read_csv(B/'recovered_sources/scib/benchmark_CPU_Harmony_additional_v1/benchmark_CPU_Harmony_comparison.tsv',sep='\t').set_index('method');pub=pd.read_csv(s/'ExtendedData1/benchmark_publication_CPU_only.tsv',sep='\t').set_index('method')
# publication naming Harmony denotes CPU branch
pd.DataFrame({'publication_method':pub.index}).to_csv(B/'audit/benchmark_names.tsv',sep='\t',index=False)
mapper={x:('Harmony_CPU' if x.lower() in ['harmony','harmony_cpu'] else x) for x in pub.index}
for x in pub.index:
 y=mapper[x];ok=y in old.index
 if ok:ok=np.allclose(pub.loc[x,['biology_score','batch_score','overall_score']].astype(float),old.loc[y,['biology_score','batch_score','overall_score']].astype(float),atol=1e-12)
 add('benchmark_scores_'+x,ok,{'native_method':y})
sec=pd.read_csv(s/'Figure5/section_input_audit.tsv',sep='\t')
add('spatial_approved_grouping',len(sec)==49 and sec.biological_unit_id.nunique()==35,{'sections':len(sec),'units':sec.biological_unit_id.nunique(),'approval':'Latest user reconfirmed original grouping; not a new identity investigation'})
current=pd.read_csv(s/'Current_Display_20260920/Figure5a_current_sections.tsv',sep='\t');spotfiles=[s/'Figure5_final_companion/CA1_existing_spot_values.tsv',s/'Current_Display_20260920/source_data/Figure5/T1077_existing_spot_values.tsv']
for (_,r),f in zip(current.iterrows(),spotfiles):
 sp=pd.read_csv(f,sep='\t');add('Figure5a_'+r.section,len(sp)==r.n_existing_spots and sha(f)==r.spot_source_sha256,{'n_spots':len(sp),'source_hash_match':sha(f)==r.spot_source_sha256,'study':r.study})
# compare regenerated masters in same native geometry, no overwrite
render=[]
for f in sorted((B/'evidence/current_figures').glob('*.pdf')):
 new=O/'render_smoke'/f.name
 if not new.exists():continue
 with fitz.open(f) as a,fitz.open(new) as b:
  shape=tuple(a[0].rect)==tuple(b[0].rect);x=a[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False);y=b[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False)
  xa=np.frombuffer(x.samples,dtype=np.uint8);ya=np.frombuffer(y.samples,dtype=np.uint8)
  equal=xa.shape==ya.shape and np.array_equal(xa,ya);delta=np.abs(xa.astype(float)-ya.astype(float)) if xa.shape==ya.shape else np.array([255])
  render.append({'figure':f.stem,'page_geometry_equal':shape,'pixel_identical':equal,'mean_abs_channel_delta':float(delta.mean()),'fraction_changed_channels':float((delta>0).mean()),'old_text_length':len(a[0].get_text()),'new_text_length':len(b[0].get_text()),'note':'Raster equality at 108 dpi; not fresh author visual approval'})
pd.DataFrame(render).to_csv(B/'audit/SVG_EXPORT_RENDER_COMPARISON.tsv',sep='\t',index=False)
shutil.copy2(O/'render_smoke/EXPORT_VALIDATION.json',B/'audit/SVG_EXPORT_VALIDATION.json')
add('SVG_export_native_dimensions',len(render)==15 and all(x['page_geometry_equal'] for x in render),{'figures':len(render),'pixel_identical':sum(x['pixel_identical'] for x in render)})
pd.DataFrame(rows).to_csv(B/'audit/ADDITIONAL_RESULT_CHECKS.tsv',sep='\t',index=False)
print(pd.DataFrame(rows)[['check','status']].to_string(index=False));print('render',json.dumps(render,ensure_ascii=False))
