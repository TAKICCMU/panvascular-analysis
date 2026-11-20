import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json
import pandas as pd,numpy as np
from scipy import sparse
W=Path(__file__).resolve().parent
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
PB=R/'pseudobulk/descriptive_final_atlas_v1/outputs'
genes=pd.read_csv(PB/'genes.tsv',sep='\t',keep_default_na=False)
assert genes.gene_id.is_unique
gi={g:i for i,g in enumerate(genes.gene_id)}
nodes=pd.read_csv(W/'all116_CoVarNet_node_evidence_review.tsv',sep='\t',keep_default_na=False)
m=pd.read_csv(W/'all247_sample_aggregation_audit_source_updated.tsv',sep='\t',keep_default_na=False)
meta=m.set_index('sample_uid')
u1=pd.read_csv(PB/'sample_L1/units.tsv.gz',sep='\t',keep_default_na=False)
u2=pd.read_csv(PB/'sample_L2/units.tsv.gz',sep='\t',keep_default_na=False)
x1=sparse.load_npz(PB/'sample_L1/counts_units_by_genes.npz')
x2=sparse.load_npz(PB/'sample_L2/counts_units_by_genes.npz')
assert x1.shape==(len(u1),len(genes)) and x2.shape==(len(u2),len(genes))
lib1=np.asarray(x1.sum(axis=1)).ravel();lib2=np.asarray(x2.sum(axis=1)).ravel()
parent={(x.sample_uid,x.level1_annotation):i for i,x in enumerate(u1.itertuples())}
details=[];summaries=[]
for n in nodes.itertuples():
 anchor=str(n.marker_anchor); gene=gi.get(anchor)
 if gene is None:continue
 ix=np.flatnonzero(u2.level1_annotation.eq(n.level1_annotation)&u2.level2_annotation.eq(n.level2_annotation)&(u2.n_cells>=20))
 local=[]
 for i in ix:
  unit=u2.iloc[i];j=parent[(unit.sample_uid,unit.level1_annotation)]
  rest_cells=int(u1.iloc[j].n_cells-unit.n_cells)
  rest_lib=float(lib1[j]-lib2[i])
  if rest_cells<20 or rest_lib<=0 or lib2[i]<=0:continue
  val=float(x2[i,gene]);rest=float(x1[j,gene])-val
  assert rest>=0
  case_rate=val/float(lib2[i])*1e6;rest_rate=rest/rest_lib*1e6
  mm=meta.loc[unit.sample_uid]
  row=dict(level1_annotation=n.level1_annotation,level2_annotation=n.level2_annotation,marker_anchor=anchor,sample_uid=unit.sample_uid,study=mm.study,donor_uid=mm.donor_uid,biological_sample_uid=mm.biological_sample_uid,safe_aggregation_status=mm.audit_status,analysis_unit='technical_sample_descriptive_only',state_cells=int(unit.n_cells),rest_parent_cells=rest_cells,state_anchor_CPM=case_rate,rest_parent_anchor_CPM=rest_rate,anchor_difference_CPM=case_rate-rest_rate,anchor_positive=case_rate>rest_rate)
  local.append(row);details.append(row)
 if local:
  z=pd.DataFrame(local);safe=z[z.safe_aggregation_status.eq('NO_CURRENT_AGGREGATION_BLOCK')]
  bs=z.groupby('study').anchor_difference_CPM.median()
  summaries.append(dict(level1_annotation=n.level1_annotation,level2_annotation=n.level2_annotation,anchor=anchor,n_sample_comparisons=len(z),n_study_comparisons=len(bs),fraction_samples_anchor_positive=float(z.anchor_positive.mean()),fraction_studies_median_positive=float((bs>0).mean()),n_unambiguous_sample_comparisons=len(safe),fraction_unambiguous_samples_positive=float(safe.anchor_positive.mean()) if len(safe) else None,interpretation='Descriptive consistency only; no independent test, no significance, no automatic node eligibility'))
pd.DataFrame(details).to_csv(W/'anchor_consistency_all_samples_descriptive.tsv',sep='\t',index=False)
pd.DataFrame(summaries).to_csv(W/'anchor_consistency_by_state_descriptive.tsv',sep='\t',index=False)
(W/'ANCHOR_AUDIT_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE_DESCRIPTIVE_ONLY',states=len(summaries),sample_state_comparisons=len(details),matrices_reused=True,network_models_run=False,DE_models_run=False),indent=2))
print((W/'ANCHOR_AUDIT_COMPLETE.json').read_text())
