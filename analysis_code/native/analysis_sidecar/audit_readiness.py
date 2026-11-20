import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,collections,sys
import pandas as pd,numpy as np
W=Path(__file__).resolve().parent
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
PBC=json.loads((R/'pseudobulk/descriptive_final_atlas_v1/config.json').read_text())
cohort=Path(PBC['cohort_output']);phase=cohort.parent.parent
m=pd.read_csv(cohort,sep='\t',keep_default_na=False)
c=pd.read_csv(phase/'counts/sample_state_counts.tsv',sep='\t',keep_default_na=False)
ont=pd.read_csv(phase/'counts/ontology.tsv',sep='\t',keep_default_na=False)
assert m.sample_uid.is_unique and len(m)==247
keys=['level1_annotation','level2_annotation']
assert not c.duplicated(['sample_uid']+keys).any()
assert c['count'].sum()==1175115
def truth(s):return s.astype(str).str.lower().isin(['true','1'])
m['condition']=m.comparison_group
m['person']=np.where(truth(m.donor_resolved),'DONOR::'+m.donor_uid,'UNRESOLVED_BIO::'+m.biological_sample_uid)
m['profile']=[json.dumps(v,separators=(',',':')) for v in zip(m.study,m.person,m.tissue,m.condition)]
bio_n=m.groupby('biological_sample_uid').sample_uid.transform('size')
bad_technical=(bio_n>1)&~truth(m.aggregate_approved)
profile_multi=m.groupby('profile').biological_sample_uid.transform('nunique')>1
bad_profile=profile_multi&~truth(m.profile_aggregate_approved)
m['unresolved_technical_aggregation']=bad_technical
m['unreviewed_distinct_specimen_profile_merge']=bad_profile
m['proposed_core_include']=~bad_technical&~bad_profile
m.to_csv(W/'sample_profile_decisions.tsv',sep='\t',index=False)
core=m[m.proposed_core_include]
x=c.merge(core,on='sample_uid',validate='many_to_one')
totals=x.groupby('profile')['count'].sum();profiles=totals[totals>=100].index
x=x[x.profile.isin(profiles)]
a=x[(x['count']>0)&~x.level1_annotation.isin(['Unknown','Not_assessed','Pending'])&~x.level2_annotation.isin(['Unknown','Not_assessed','Pending'])].groupby(keys).agg(cells=('count','sum'),profiles=('profile','nunique'),studies=('study','nunique')).reset_index()
a['frequency_pass']=(a.cells>=500)&(a.profiles>=10)&(a.studies>=3)
a=ont.merge(a,on=keys,how='left',validate='one_to_one')
a.to_csv(W/'CoVarNet_node_frequency_review.tsv',sep='\t',index=False)
a[a.frequency_pass.eq(True)].to_csv(W/'CoVarNet_candidate_nodes.tsv',sep='\t',index=False)
# Current counts and metadata only; does not change inference_eligible.
units=pd.read_csv(R/'pseudobulk/descriptive_final_atlas_v1/outputs/biological_L2/units.tsv.gz',sep='\t',keep_default_na=False)
(W/'pseudobulk_units_schema.json').write_text(json.dumps(dict(columns=list(units),rows=len(units)),indent=2))
g=m[m.study.eq('GSE222318')].copy()
assert len(g)==14 and g.donor_uid.nunique()==14 and truth(g.aggregate_approved).all()
sex=pd.get_dummies(g.sex,drop_first=True).to_numpy(dtype=float)
design=np.column_stack([np.ones(len(g)),sex,g.comparison_group.eq('case').astype(int)])
assert np.linalg.matrix_rank(design)==design.shape[1]
g.to_csv(W/'GSE222318_contrast_membership.tsv',sep='\t',index=False)
# Coverage after current min20 cells for each biological sample/state.
cx=c.merge(g[['sample_uid','biological_sample_uid','donor_uid','comparison_group','sex']],on='sample_uid',validate='many_to_one')
state_rows=[]
for level,col in [('L1','level1_annotation'),('L2','level2_annotation')]:
 groupcols=['biological_sample_uid','donor_uid','comparison_group','sex','level1_annotation']+(['level2_annotation'] if level=='L2' else [])
 z=cx.groupby(groupcols,as_index=False)['count'].sum()
 for ident,v in z.groupby(['level1_annotation']+(['level2_annotation'] if level=='L2' else [])):
  labels=ident if isinstance(ident,tuple) else (ident,)
  if any(t in ['Unknown','Not_assessed','Pending'] for t in labels):continue
  v=v[v['count']>=20]
  arms=v.groupby('comparison_group').donor_uid.nunique()
  state_rows.append(dict(level=level,state=str(ident),case_donors=int(arms.get('case',0)),reference_donors=int(arms.get('reference',0)),passes_min5=(arms.get('case',0)>=5 and arms.get('reference',0)>=5)))
pd.DataFrame(state_rows).to_csv(W/'GSE222318_state_coverage_review.tsv',sep='\t',index=False)
summary=dict(status='REAL_INPUT_REVIEW_COMPLETE_NO_MODELS_RUN',atlas_cells=1175115,total_samples=len(m),unapproved_aggregate_records=int((~truth(m.aggregate_approved)).sum()),technical_merge_blocked_records=int(bad_technical.sum()),distinct_specimen_merge_blocked_records=int(bad_profile.sum()),proposed_core_samples=len(core),proposed_core_profiles=len(profiles),proposed_core_cells=int(x['count'].sum()),frequency_candidate_nodes=int(a.frequency_pass.eq(True).sum()),named_ontology_unreviewed=116,GSE222318_case_donors=9,GSE222318_reference_donors=5,GSE222318_age_known=int(g.age.ne('').sum()),GSE222318_sex_by_arm=g.groupby(['comparison_group','sex']).size().to_dict(),sex_adjusted_design_rank=int(np.linalg.matrix_rank(design)),sex_adjusted_design_columns=design.shape[1],eligible_L1_states=sum(v['passes_min5'] for v in state_rows if v['level']=='L1'),eligible_L2_states=sum(v['passes_min5'] for v in state_rows if v['level']=='L2'),formal_models_started=False)
summary['GSE222318_sex_by_arm']={str(k):int(v) for k,v in summary['GSE222318_sex_by_arm'].items()}
(W/'READINESS.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,default=lambda v:int(v)))
sources=[cohort,phase/'counts/sample_state_counts.tsv',phase/'counts/ontology.tsv',R/'pseudobulk/descriptive_final_atlas_v1/outputs/biological_L2/units.tsv.gz']
(W/'SOURCE_SHA256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},indent=2))
print(json.dumps(summary,indent=2,ensure_ascii=False,default=lambda v:int(v)))
