from common import *
import pandas as pd,numpy as np,shutil
old=R/'production_launches/communication_L1_20260914T015040Z';out=W/'communication';(out/'inputs/specimen_rows').mkdir(parents=True,exist_ok=False);(out/'code').mkdir()
m=pd.read_csv(AUDIT/'sample_sidecar_updated.tsv',sep='\t',keep_default_na=False).set_index('sample_uid')
a=pd.read_csv(old/'inputs/cell_membership.tsv.gz',sep='\t',index_col=0);a['row_index']=np.arange(len(a));assert len(a)==1175115 and a.index.is_unique and set(a.sample_uid)==set(m.index)
oldplans=pd.read_csv(old/'inputs/selected_specimens.tsv',sep='\t');oldby={}
for r in oldplans.to_dict('records'):
 z=pd.read_csv(r['membership'],sep='\t',index_col=0);assert z.sample_uid.nunique()==1;oldby[z.sample_uid.iloc[0]]=r
plans=[];elig=[]
for sample,sub in a.groupby('sample_uid',sort=True):
 counts=sub.level1_annotation.value_counts();allowed=sorted(counts[(counts>=50)&~counts.index.isin(['Unknown','Not_assessed','Pending'])].index)
 eligible=len(allowed)>=2;elig.append(dict(sample_uid=sample,n_types=len(allowed),eligible=eligible))
 if not eligible:continue
 if sample in oldby:
  r=oldby[sample];rid=r['run_id'];seed=int(r['seed']);rows=pd.read_csv(r['membership'],sep='\t',index_col=0)
  assert set(rows.index)<=set(sub.index)
 else:
  rid='specimen_'+hashlib.sha256(sample.encode()).hexdigest()[:14];seed=1+((20260820+int(hashlib.sha256(f'20260820|{rid}'.encode()).hexdigest()[:8],16))%(2**31-2));rng=np.random.default_rng(seed);chosen=[]
  for t in allowed:
   ids=sub.index[sub.level1_annotation.eq(t)];chosen.extend(rng.choice(ids,min(len(ids),200),replace=False))
  rows=sub.loc[chosen].sort_values('row_index').copy()
 rows['original_biological_sample_uid']=rows.biological_sample_uid;rows['biological_sample_uid']=sample;rows['analysis_unit_semantics']='independent sample/capture; donor block retained for inference'
 for k in ['donor_uid','group_binary','tissue','cohort_disease_context','source_diagnosis']:rows[k]=m.loc[sample,k]
 rows['aggregate_approved']=True
 p=out/'inputs/specimen_rows'/(rid+'.tsv.gz');rows.to_csv(p,sep='\t')
 plans.append(dict(run_id=rid,seed=seed,n_cells=len(rows),n_types=len(allowed),biological_sample_uid=sample,sample_uid=sample,donor_uid=m.loc[sample,'donor_uid'],membership=str(p),membership_sha256=sha(p),old_selected_cells_reused=sample in oldby,recalculate_native_scores=True))
assert len(plans)==229 and sum(x['old_selected_cells_reused'] for x in plans)==203
pd.DataFrame(plans).to_csv(out/'inputs/selected_specimens.tsv',sep='\t',index=False);pd.DataFrame(elig).to_csv(out/'inputs/eligibility_all247.tsv',sep='\t',index=False);m.to_csv(out/'inputs/sample_metadata.tsv',sep='\t')
shutil.copy2(old/'inputs/feature_mapping.tsv.gz',out/'inputs/feature_mapping.tsv.gz');shutil.copy2(old/'code/cellchat_native.R',out/'code/cellchat_native.R')
s=(old/'code/communication_worker.py').read_text().replace('W=Path(__file__).resolve().parents[1];R=W.parents[1]',f'W=Path(__file__).resolve().parents[1];R=Path({str(R)!r})')
s=s.replace('from scipy import sparse,io','from scipy import sparse,io\nfrom scipy.io._fast_matrix_market import _fmm_core\nfrom threadpoolctl import threadpool_limits\nthreadpool_limits(limits=1)')
(out/'code/communication_worker.py').write_text(s)
(out/'code/communication_runner.py').write_text('import sys\nfrom pathlib import Path\nsys.path.insert(0,str(Path(__file__).resolve().parents[2]/"code"))\nfrom common import verify as check\ndef verify(receipt):\n return check()\n')
archive=out/'archive_203';archive.mkdir();p=D/'communication_data_v1/parallel4_r2/ALL_VALID_RESULTS_INDEX.json';shutil.copy2(p,archive/p.name)
dump(archive/'ARCHIVE.json',dict(original_root=str(D/'communication_data_v1/parallel4_r2'),index_sha256=sha(p),policy='Original 203 results retained immutable; replacement current index published only after all 229 validate'))
dump(out/'PREPARATION.json',dict(eligible=229,low_representation=18,prior=203,new=26,recomputed=229,unchanged_sampling_for_old203=True,min_cells=50,max_cells_per_type=200,liana_expr_prop=.1,liana_permutations=1000,CellChat_nboot=100,CellChat_population_size=False,communication_level='L1',L2='not in original manuscript design',donor_handling='Samples separate; donor-aware inference downstream'))
print('Prepared 229 independent sample communication runs',flush=True)
