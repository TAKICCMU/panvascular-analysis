"""Read-only input/posterior checks; writes only this branch validation artifacts."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,time,traceback,csv
JOB=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def require_complete(status):
 if status.get('status')!='COMPLETE':raise RuntimeError('Production is not COMPLETE: '+str(status))
def main():
 if '--self-test' in sys.argv:
  for state in ['PARTIAL_FAILED','FAILED','RUNNING',None]:
   try:require_complete({'status':state})
   except RuntimeError:pass
   else:raise AssertionError('Failure state incorrectly accepted')
  require_complete({'status':'COMPLETE'});print('PASS: non-complete states rejected');return
 import numpy as np,pandas as pd,anndata as ad
 from exclusion import token,excluded_mask,assert_excluded
 out=JOB/'outputs';v=JOB/'validation';v.mkdir(exist_ok=True)
 report={'status':'FAIL','started_utc_epoch':time.time(),'validation_code_sha256':sha(__file__),'claim':'technical integrity and exclusion validation; not proof of scientific convergence'}
 hashes=[]
 def bind(path,role,expected=None):
  path=Path(path);actual=sha(path)
  if expected is not None:assert actual==expected,'Hash mismatch '+str(path)
  hashes.append({'path':str(path),'role':role,'size_bytes':path.stat().st_size,'sha256':actual});return actual
 try:
  require_complete(json.loads((JOB/'RUN_STATUS.json').read_text()))
  assert not list(out.glob('*_FAILED.json')) and not list(out.glob('*_BLOCKED.json')),'Unresolved failure/block markers'
  c=json.loads((JOB/'configs/config.json').read_text());r=json.loads((JOB/'configs/PRODUCTION_RECEIPT.json').read_text())
  receipt_sha=bind(JOB/'configs/PRODUCTION_RECEIPT.json','receipt');code_sha=bind(JOB/'code/spatial_target_excluded.py','runner')
  for p,h in r['bound_files'].items():bind(p,'bound_input',h)
  exc=set(c['excluded_feature_tokens']);assert len(c['target_ids'])==29
  rc=json.loads((out/'reference_COMPLETE.json').read_text())
  assert rc['reference_retrained'] is True and rc['target_excluded_before_model_setup'] is True
  assert rc['code_sha']==code_sha and rc['receipt_sha']==receipt_sha
  sig=pd.read_csv(out/'reference_signatures.tsv.gz',sep='\t',index_col=0)
  sigsha=bind(out/'reference_signatures.tsv.gz','reference_signature',rc['signature_sha'])
  assert sig.shape==(20401,128) and sig.index.is_unique and sig.columns.is_unique
  assert np.isfinite(sig.to_numpy()).all() and (sig.to_numpy()>=0).all()
  assert_excluded(pd.DataFrame(index=sig.index),exc)
  ref=ad.read_h5ad(out/'reference_counts.h5ad',backed='r');assert_excluded(ref.var,exc)
  assert ref.shape==(167660,20401) and set(ref.var_names)==set(sig.index)
  original_ref=ad.read_h5ad(c['reference_counts'],backed='r')
  assert ref.obs_names.equals(original_ref.obs_names)
  for col in ['reference_state','reference_study']:assert ref.obs[col].astype(str).equals(original_ref.obs[col].astype(str))
  assert int(excluded_mask(original_ref.var,exc).sum())==27
  ref.file.close();original_ref.file.close()
  bind(out/'reference_counts.h5ad','excluded_reference_counts');bind(out/'reference_model/model.pt','new_reference_model');bind(out/'reference_COMPLETE.json','reference_completion')
  manifest=pd.read_csv(c['spatial_manifest'],sep='\t');manifest=manifest[manifest.include_mapping.astype(str).str.lower().isin(['true','1','yes'])]
  assert manifest.section_id.nunique()==49 and manifest.GSE.nunique()==7
  datasets=[];sections=[];all_spots=[];fallback=[];factor_names=list(sig.columns)
  old=Path(c['reuse_completed_root']);resume=JOB.parents[1]/'spatial_gpu_resume_v2'
  for gse,rows in manifest.groupby('GSE',sort=True):
   cp=out/f'{gse}_COMPLETE.json';z=json.loads(cp.read_text())
   assert z['code_sha']==code_sha and z['receipt_sha']==receipt_sha and z['reference_sha']==sigsha
   assert z['excluded_targets_present']==[]
   posterior=out/f'{gse}_posterior.h5ad';bind(posterior,'target_excluded_posterior',z['posterior_sha']);bind(cp,'dataset_completion')
   olddir=next((d for d in [resume,old] if (d/'outputs'/f'{gse}_provenance.json').exists()),None)
   assert olddir is not None,'Missing original completed reference dataset '+gse
   origpath=olddir/'outputs'/f'{gse}_posterior.h5ad';origprov=json.loads((olddir/'outputs'/f'{gse}_provenance.json').read_text());bind(origpath,'original_posterior',origprov['posterior_sha'])
   a=ad.read_h5ad(posterior,backed='r');orig=ad.read_h5ad(origpath,backed='r')
   assert a.obs_names.equals(orig.obs_names),'Exact spot order differs '+gse
   assert set(a.obs.section_id)==set(rows.section_id)
   for col in ['sample_id','section_id','biological_unit_id','QC_fallback','total_counts_original','n_genes_original']:
    assert a.obs[col].astype(str).equals(orig.obs[col].astype(str)),'Original spot metadata differs '+gse+' '+col
   assert np.array_equal(a.obsm['spatial'],orig.obsm['spatial'])
   assert_excluded(a.var,exc);assert set(a.var_names)<=set(sig.index)
   assert list(map(str,a.uns['mod']['factor_names']))==factor_names
   vals={stat:np.asarray(a.obsm[stat+'_cell_abundance_w_sf']) for stat in ['means','q05','q95']}
   for stat,x in vals.items():
    assert x.shape==(a.n_obs,128) and np.isfinite(x).all() and (x>=0).all(),gse+' invalid '+stat
   tol=1e-7+1e-6*np.abs(vals['q95'])
   assert np.all(vals['q05']<=vals['q95']+tol),'Quantile order violation '+gse
   outside=int(((vals['means']<vals['q05']-tol)|(vals['means']>vals['q95']+tol)).sum())
   f=sorted(a.obs.loc[a.obs.QC_fallback.astype(str).str.lower().eq('true'),'section_id'].unique().tolist())
   assert f==sorted(origprov['QC_fallback_sections']);fallback+=f
   all_spots+=list(a.obs_names)
   for sec,sobs in a.obs.groupby('section_id',observed=True):
    rawpath=out/'full_umi'/(str(sec).replace(':','_')+'.h5ad');raw=ad.read_h5ad(rawpath,backed='r')
    assert raw.obs_names.equals(sobs.index),'Full source UMI spot mismatch '+sec
    row=rows[rows.section_id==sec].iloc[0];gm=pd.read_csv(JOB/'evidence'/f'{gse}__{str(sec).replace(":","_")}__gene_map.tsv.gz',sep='\t',keep_default_na=False)
    assert raw.n_vars==len(gm) and raw.var.index.tolist()==gm.source_row.astype(str).tolist()
    for col in ['source_id','source_symbol','atlas_gene']:assert raw.var[col].fillna('').astype(str).tolist()==gm[col].astype(str).tolist()
    sections.append({'GSE':gse,'section_id':sec,'n_spots':len(sobs),'source_features':raw.n_vars,'QC_fallback':sec in f});raw.file.close();bind(rawpath,'full_source_feature_UMI')
   fullpath=out/f'{gse}_full_atlas_gene_UMI.h5ad';full=ad.read_h5ad(fullpath,backed='r');assert full.obs_names.equals(a.obs_names);full.file.close();bind(fullpath,'full_atlas_gene_UMI')
   datasets.append({'GSE':gse,'n_sections':len(rows),'n_spots':a.n_obs,'n_model_genes':a.n_vars,'n_factors':128,'quantile_order':'PASS','finite_nonnegative':'PASS','mean_outside_q05_q95_entries':outside,'exact_original_spots':'PASS','original_source':str(olddir)})
   a.file.close();orig.file.close()
  assert len(datasets)==7 and len(sections)==49 and len(all_spots)==45348 and len(set(all_spots))==45348
  assert len(fallback)==2 and len(set(fallback))==2
  pd.DataFrame(datasets).to_csv(v/'dataset_QC.tsv',sep='\t',index=False);pd.DataFrame(sections).to_csv(v/'section_QC.tsv',sep='\t',index=False)
  report.update(status='PASS',n_datasets=7,n_sections=49,n_spots=45348,n_targets=29,reference_excluded_features=27,fallback_sections=sorted(fallback),mean_quantile_policy='Require q05<=q95. Mean is not a quantile; its location outside q05/q95 is counted descriptively, not assumed impossible.',reference_retrained=True,source_bytes_total=sum(x['size_bytes'] for x in hashes if x['role']=='bound_input'),dataset_QC=datasets)
  (v/'QC报告.md').write_text('空间双阶段目标排除完整性验证通过：7数据集、49切片、45348 spots；29目标身份均已排除，reference确认重训；两个fallback与原样本/spot顺序完全一致。所有mean/q05/q95有限非负、q05≤q95。完整原始feature UMI保存且spot及feature轴核对通过。mean并非分位数，区间外个数仅作诊断。技术通过不代替收敛及科学验证。\n')
 except Exception as e:
  report.update(status='FAIL',error=repr(e),traceback=traceback.format_exc());raise
 finally:
  report['finished_utc_epoch']=time.time();(v/'VALIDATION.json').write_text(json.dumps(report,indent=2,ensure_ascii=False));pd.DataFrame(hashes).to_csv(v/'file_hashes_and_bytes.tsv',sep='\t',index=False)
if __name__=='__main__':main()
