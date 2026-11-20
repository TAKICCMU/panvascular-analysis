from common import *
import pandas as pd,numpy as np,gzip,shutil
from scipy import sparse,io
from scipy.io._fast_matrix_market import _fmm_core
from threadpoolctl import threadpool_limits
out=W/'pseudobulk';out.mkdir(exist_ok=True)
meta=pd.read_csv(AUDIT/'sample_sidecar_updated.tsv',sep='\t',keep_default_na=False).set_index('sample_uid')
full=pd.read_csv(AUDIT.parent/'full_review_v1/all247_sample_aggregation_audit_source_updated.tsv',sep='\t',keep_default_na=False).set_index('sample_uid')
specs=pd.read_csv(C/'configs/contrasts.tsv',sep='\t',keep_default_na=False)
aud=[]
for s in specs.to_dict('records'):
 z=full.copy()
 for k in ['study','vascular_bed','cohort_disease_context','assay_type']:z=z[z[k].eq(s[k])]
 z=z.join(meta[['group_binary','donor_uid']].rename(columns={'donor_uid':'updated_donor'}))
 case=set(z.loc[z.group_binary.eq('disease'),'updated_donor']);control=set(z.loc[z.group_binary.eq('healthy'),'updated_donor']);paired=str(s['paired']).lower()=='true'
 n=int(s['min_independent_per_arm']);enough=len(case&control)>=n if paired else min(len(case),len(control))>=n
 aud.append({**s,'n_case_donors':len(case),'n_reference_donors':len(control),'n_pairs':len(case&control),'passes_donor_gate':enough})
audit=pd.DataFrame(aud);audit.to_csv(out/'all_contrast_readiness.tsv',sep='\t',index=False)
passed=audit[audit.passes_donor_gate]
assert set(passed.study)=={'GSE222318'},'Changed estimable study set requires model-specific design review'
pb=R/'pseudobulk/descriptive_final_atlas_v1/outputs';genes=pd.read_csv(pb/'genes.tsv',sep='\t',keep_default_na=False);genes.to_csv(out/'genes.tsv',sep='\t',index=False)
roster=[]
for level in ['L1','L2']:
 u=pd.read_csv(pb/('sample_'+level)/'units.tsv.gz',sep='\t',keep_default_na=False);x=sparse.load_npz(pb/('sample_'+level)/'counts_units_by_genes.npz')
 assert x.shape==(len(u),len(genes))
 z=u[u.sample_uid.isin(meta.index[meta.study.eq('GSE222318')])&(u.n_cells>=20)].copy()
 keys=['level1_annotation']+(['level2_annotation'] if level=='L2' else [])
 for label,v in z.groupby(keys):
  labels=label if isinstance(label,tuple) else (label,)
  if any(a in ['Unknown','Not_assessed','Pending'] for a in labels):continue
  v=v.copy();v['donor_uid']=v.sample_uid.map(meta.donor_uid);v['comparison_group']=v.sample_uid.map(meta.group_binary).map({'disease':'case','healthy':'reference'})
  v['sex']=v.sample_uid.map(meta.sex)
  arms=v.groupby('comparison_group').donor_uid.nunique()
  if min(arms.get('case',0),arms.get('reference',0))<5:continue
  assert v.donor_uid.is_unique
  design=np.column_stack([np.ones(len(v)),v.sex.eq('male').astype(int),v.comparison_group.eq('case').astype(int)])
  assert np.linalg.matrix_rank(design)==3 and len(v)-3>=2
  name=level+'_'+hashlib.sha256(json.dumps(labels).encode()).hexdigest()[:12];p=out/name;p.mkdir()
  sub=x[v.index,:].T.tocsc();assert (sub.data>=0).all() and (sub.data==np.floor(sub.data)).all()
  with gzip.open(p/'counts.mtx.gz','wb') as h:
   with threadpool_limits(limits=1):io.mmwrite(h,sub,field='integer')
  v.to_csv(p/'units.tsv',sep='\t',index=False)
  spec=dict(contrast_id='GSE222318_sex_adjusted_v2',cell_state=' | '.join(labels),level=level,meta_context='aortic dissection',case='case',reference='reference',formula='~ sex + group',min_tested_genes=100,rows=v.to_dict('records'))
  dump(p/'selected_design.json',spec);roster.append(dict(model_id=name,cell_state=' | '.join(labels),level=level,n_case=int(arms['case']),n_reference=int(arms['reference']),input=str(p)))
 del x
pd.DataFrame(roster).to_csv(out/'model_roster.tsv',sep='\t',index=False)
shutil.copy2(C/'scripts/pseudobulk_edgeR.R',W/'code/pseudobulk_edgeR.R')
dump(out/'PREPARATION.json',dict(status='PREPARED_FOR_AUTHORIZED_EXECUTION',models=len(roster),level_counts=pd.DataFrame(roster).level.value_counts().to_dict(),all_contrasts_reviewed=len(audit),cross_disease_meta='NOT_ESTIMABLE: only one study passes original independent-donor threshold',age='entirely missing in GSE222318; not included',sex='case9 male; reference3 male2female; identifiable but imbalanced',reused='sample-level raw pseudobulk only; no old biological aggregation reused'))
print('Pseudobulk prepared',len(roster),'models',flush=True)
