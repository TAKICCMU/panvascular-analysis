"""Validate cache signatures and freeze copied metric evidence without recomputing."""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
import hashlib,json,shutil
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907')
W=R/'02_scib_benchmark/studywise_sampleonly_v1';O=R/'06_reanalysis_takeover_20260908/scib/audit';O.mkdir(parents=True,exist_ok=True)
def fh(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(16*1024*1024):h.update(b)
 return h.hexdigest()
def nh(v):
 h=hashlib.sha256()
 for x in v:h.update(str(x).encode());h.update(b'\0')
 return h.hexdigest()
records=[];populations=[]
for stage,cap in [('pilot',500),('primary',1500)]:
 inp=W/f'inputs/scib_{stage}_sampleonly_cap{cap}_hvg2000.h5ad';a=ad.read_h5ad(inp,backed='r');obs=a.obs.copy();ch=nh(a.obs_names);shape=list(a.shape);a.file.close()
 baseline=W/f'omicverse_pca/{stage}/omicverse_pca.npy';bm=baseline.with_name('omicverse_pca.metadata.json');bh=fh(baseline);bmh=fh(bm)
 population={'stage':stage,'input_sha256':fh(inp),'shape':shape,'obs_names_sha256':ch,'sample_uid':int(obs.sample_uid.nunique()),
             'study_file':int(obs.study_file.nunique()),'study_canonical':int(obs.study_canonical.nunique()),
             'maximum_cells_per_sample':int(obs.sample_uid.value_counts().max()),'sample_cap':cap,
             'excluded_source_representations_present':bool(obs.study_file.isin(['GSE198750','GSE216860']).any())}
 assert population['sample_uid']==247 and population['study_file']==29 and population['study_canonical']==26
 assert not population['excluded_source_representations_present'] and population['maximum_cells_per_sample']<=cap
 populations.append(population)
 # Pilot cache headers predate the v3 name: validate the actual saved members,
 # rather than treating a version string alone as proof of sampling behavior.
 rng=np.random.default_rng(0);picked=[];frame=obs[['sample_uid']].astype(str).copy();frame['_p']=np.arange(len(obs))
 cap_metric=25000//247
 for _,g in frame.groupby('sample_uid',sort=True,observed=True):
  pos=g['_p'].to_numpy();picked.append(rng.choice(pos,cap_metric,replace=False) if len(pos)>cap_metric else pos)
 expected_subset=np.sort(np.concatenate(picked))
 saved=pd.read_csv(W/stage/'metrics/silhouette_pcr_metric_subset.tsv.gz',sep='\t')
 if not np.array_equal(saved.source_position.to_numpy(),expected_subset):
  raise RuntimeError(f'{stage} metric subset is not the frozen sample-only draw')
 population['sample_only_metric_membership_verified']=True
 table=pd.read_csv(W/stage/'metrics/scib_metrics_by_K.tsv',sep='\t')
 for (method,config),rows in table.groupby(['method','config_id']):
  meta=W/stage/'embeddings'/f'{method}__{config}.json';arr=meta.with_suffix('.npy')
  cache=W/stage/'metrics/per_embedding_cache'/f'{method}__{config}.json';ct=cache.with_suffix('.tsv')
  c=json.loads(cache.read_text());m=json.loads(meta.read_text());sig=c['signature'];sig_expected={
   'evaluator_contract_version':'2026-09-08.sample-only-subset.primary-allowlist.v3','stage':stage,'input_obs_names_sha256':ch,
   'embedding_metadata_sha256':fh(meta),'embedding_file_size':arr.stat().st_size,'baseline_pca_sha256':bh,'baseline_pca_metadata_sha256':bmh}
  content_keys=['input_obs_names_sha256','embedding_metadata_sha256','embedding_file_size','baseline_pca_sha256','baseline_pca_metadata_sha256']
  if not all(sig.get(k)==sig_expected[k] for k in content_keys):
   raise RuntimeError(f'Content signature mismatch: {stage}/{method}/{config}: {sig}')
  assert c['status']=='completed' and m['obs_names_sha256']==ch
  if stage=='primary':assert sig==sig_expected
  cached=pd.read_csv(ct,sep='\t');cols=[x for x in cached.columns if x in rows]
  pd.testing.assert_frame_equal(cached[cols].sort_values('K').reset_index(drop=True),rows[cols].sort_values('K').reset_index(drop=True),check_dtype=False,rtol=1e-10,atol=1e-12)
  records.append({'stage':stage,'method':method,'config_id':config,'cache_signature_valid':True,'cached_rows_equal_source_table':True,
                  'cache_contract_version':sig.get('evaluator_contract_version'),
                  'embedding_file':str(arr),'embedding_sha256_current':fh(arr),'metadata_sha256':fh(meta),'cache_table_sha256':fh(ct),
                  'previous_embedding_content_hash_available':'embedding_sha256' in m,
                  'retrospective_limit':'Old cache bound embedding size and metadata hash, not full embedding content; current content hash is frozen now.'})
pd.DataFrame(records).to_csv(O/'REUSED_EMBEDDING_METRIC_SIGNATURES.tsv',sep='\t',index=False)
(O/'REUSED_METRIC_AUDIT.json').write_text(json.dumps({'status':'passed','populations':populations,'n_validated_embeddings':len(records),
 'current_content_hashes_frozen':True,'historical_cache_limitation':'No content hash existed for all original embeddings; size+metadata+input/PCA signatures are verified and current hashes now frozen.'},indent=2))
print(json.dumps(populations,indent=2),flush=True)
