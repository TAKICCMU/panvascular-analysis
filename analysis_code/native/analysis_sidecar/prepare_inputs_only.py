import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,sys,gzip,shutil,hashlib
import pandas as pd,numpy as np
from scipy import sparse,io
W=Path(__file__).resolve().parent
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
CODE=R/'downstream_CPU_L1L2_v2_1_release_compat_20260914'
sys.path.insert(0,str(CODE))
from src.covarnet import prepare,write_prepared
cfg=json.loads((R/'pseudobulk/descriptive_final_atlas_v1/config.json').read_text())
phase=Path(cfg['cohort_output']).parent.parent
meta=pd.read_csv(W/'sample_profile_decisions.tsv',sep='\t',keep_default_na=False)
meta=meta[meta.proposed_core_include.astype(str).str.lower().eq('true')]
counts=pd.read_csv(phase/'counts/sample_state_counts.tsv',sep='\t')
counts=counts[counts.sample_uid.isin(meta.sample_uid)]
ont=pd.read_csv(phase/'counts/ontology.tsv',sep='\t',keep_default_na=False)
prepared=prepare(counts,meta,ont,{'approve_named_state_frequency_rule':False,'workers':1,'threads_per_worker':4})
write_prepared(prepared,W/'covarnet_prepared_PROPOSED',{'atlas_sha256':cfg['atlas_sha256'],'core_selection':'sample_profile_decisions.tsv; proposed only','frozen_sources_unchanged':True})
PB=R/'pseudobulk/descriptive_final_atlas_v1/outputs'
dest=W/'pseudobulk_prepared_PROPOSED';dest.mkdir()
genes=pd.read_csv(PB/'genes.tsv',sep='\t',keep_default_na=False)
genes.to_csv(dest/'genes.tsv',sep='\t',index=False)
roster=[]
for level in ['L1','L2']:
 u=pd.read_csv(PB/('biological_'+level)/'units.tsv.gz',sep='\t',keep_default_na=False)
 matrix=sparse.load_npz(PB/('biological_'+level)/'counts_units_by_genes.npz')
 assert matrix.shape==(len(u),len(genes))
 z=u[u.study.eq('GSE222318')&(u.n_cells>=20)].copy()
 keys=['level1_annotation']+(['level2_annotation'] if level=='L2' else [])
 for ident,v in z.groupby(keys):
  labels=ident if isinstance(ident,tuple) else (ident,)
  if any(x in ['Unknown','Not_assessed','Pending'] for x in labels):continue
  arms=v.groupby('comparison_group').donor_uid.nunique()
  if min(arms.get('case',0),arms.get('reference',0))<5:continue
  assert v.donor_uid.is_unique
  design=np.column_stack([np.ones(len(v)),v.sex.eq('male').astype(int),v.comparison_group.eq('case').astype(int)])
  assert np.linalg.matrix_rank(design)==3 and len(v)-3>=2
  name=level+'_'+hashlib.sha256(json.dumps(labels).encode()).hexdigest()[:12]
  p=dest/name;p.mkdir()
  sub=matrix[v.index,:].T.tocsc()
  assert (sub.data>=0).all() and (sub.data==np.floor(sub.data)).all()
  with gzip.open(p/'counts.mtx.gz','wb') as h:io.mmwrite(h,sub,field='integer')
  v.to_csv(p/'units.tsv',sep='\t',index=False)
  spec=dict(contrast_id='GSE222318_sex_adjusted_PROPOSED',cell_state=' | '.join(labels),level=level,meta_context='aortic dissection',case='case',reference='reference',formula='~ sex + group',min_tested_genes=100,rows=v.to_dict('records'))
  (p/'selected_design.json').write_text(json.dumps(spec,indent=2))
  roster.append(dict(state=' | '.join(labels),level=level,model_input=str(p),n_units=len(v),n_case=int(arms['case']),n_reference=int(arms['reference']),raw_count_total=int(sub.sum()),model_status='PROPOSED_NOT_RUN'))
 del matrix
pd.DataFrame(roster).to_csv(dest/'model_roster.tsv',sep='\t',index=False)
(W/'PREPARATION_COMPLETE.json').write_text(json.dumps(dict(status='REAL_INPUTS_PREPARED_MODELS_NOT_RUN',n_models=len(roster),covarnet_blocking=prepared['manifest']['blocking_reasons'],n_covarnet_nodes=prepared['manifest']['n_states']),indent=2))
print((W/'PREPARATION_COMPLETE.json').read_text())
