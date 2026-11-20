from common import *
import pandas as pd,numpy as np,shutil
from scipy import sparse
P=Path(json.loads((R/'pseudobulk/descriptive_final_atlas_v1/config.json').read_text())['cohort_output']).parent.parent
out=W/'covarnet';out.mkdir(exist_ok=False)
m=pd.read_csv(AUDIT/'sample_sidecar_updated.tsv',sep='\t',keep_default_na=False)
m=m[m.covarnet_include.astype(str).str.lower().eq('true')&(m.n_cells>=100)].copy()
review=pd.read_csv(AUDIT/'L2_node_eligibility_review.tsv',sep='\t')
core=review[review.recommendation.eq('RECOMMEND_CORE')]
assert len(m)==223 and len(core)==89
c=pd.read_csv(P/'counts/sample_state_counts.tsv',sep='\t',keep_default_na=False);c=c[c.sample_uid.isin(m.sample_uid)]
k=['level1_annotation','level2_annotation'];c=c.merge(review[k+['level2_state_id']],on=k,how='left',validate='many_to_one')
ids=sorted(m.sample_uid);states=sorted(core.level2_state_id)
counts=c[c.level2_state_id.isin(states)].pivot(index='level2_state_id',columns='sample_uid',values='count').fillna(0).reindex(index=states,columns=ids,fill_value=0).astype('int64')
counts.index.name='cell_state'
total=c.groupby('sample_uid')['count'].sum().reindex(ids)
parent=c.groupby(['level1_annotation','sample_uid'])['count'].sum().unstack(fill_value=0).reindex(columns=ids)
mapping=core[['level2_state_id','level1_annotation','level2_annotation']].rename(columns={'level2_state_id':'cell_state','level1_annotation':'parent_lineage'}).set_index('cell_state').loc[states].reset_index()
assert not counts.sum(axis=0).eq(0).any() and int(total.sum())==1100160
within=counts.astype(float).copy()
for x in mapping.itertuples():within.loc[x.cell_state]=counts.loc[x.cell_state]/parent.loc[x.parent_lineage].replace(0,np.nan)
meta=m.set_index('sample_uid').loc[ids].reset_index().rename(columns={'sample_uid':'covarnet_sample_uid','study':'study_canonical','tissue':'tissue_v2','group_binary':'condition_v2','donor_block':'covarnet_donor_key'})
meta['n_total_atlas_cells']=total.to_numpy();meta['n_technical_samples']=1
for name,df,index in [('CoVarNet_node_counts.tsv',counts,True),('CoVarNet_node_proportion_of_all_atlas_cells.tsv',counts.divide(total,axis=1),True),('CoVarNet_node_proportion_among_stable_nodes_sensitivity.tsv',counts.divide(counts.sum(axis=0),axis=1),True),('CoVarNet_proportion_within_full_parent.tsv',within,True),('CoVarNet_full_parent_counts.tsv',parent,True),('CoVarNet_sample_metadata.tsv',meta,False),('CoVarNet_node_parent_map.tsv',mapping,False)]:df.to_csv(out/name,sep='\t',index=index,na_rep='NA')
core.to_csv(out/'node_eligibility_applied.tsv',sep='\t',index=False)
for name in ['covarnet_nmf.R','covarnet_parent_null.R','covarnet_review_helpers.R']:shutil.copy2(C/'src'/name,W/'code'/name)
shutil.copytree(C/'src/covarnet_legacy',W/'code/covarnet_legacy')
dump(out/'INPUT_MANIFEST.json',dict(status='PREPARED_NOT_EXECUTED',n_profiles=223,n_states=89,n_cells=1100160,profile='separate sample_uid; donor-block resampling',source_counts_sha256=sha(P/'counts/sample_state_counts.tsv'),source_node_review_sha256=sha(AUDIT/'L2_node_eligibility_review.tsv'),known_program_validation='descriptive eligibility; NMF stability not yet established',files={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
print('CoVarNet inputs validated: 223 profiles x 89 nodes',flush=True)
