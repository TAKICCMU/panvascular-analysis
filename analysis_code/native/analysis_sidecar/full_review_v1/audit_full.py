import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
W=Path(__file__).resolve().parent
cfg=json.loads((R/'pseudobulk/descriptive_final_atlas_v1/config.json').read_text())
phase=Path(cfg['cohort_output']).parent.parent
prov=R/'L2_targeted_followup_metadata_v4_20260913_154850_913989/audit/metadata_provenance_v4.csv'
m=pd.read_csv(W.parent/'sample_profile_decisions.tsv',sep='\t',keep_default_na=False)
e=pd.read_csv(prov,keep_default_na=False)
assert e.sample_uid.is_unique
efields=['sample_uid','sample_name','source_sample_accession','sra_experiment_accession','biosample_accession','tissue_detail','replicate_type','replicate_evidence_status','primary_inference_constraint','donor_evidence','source_diagnosis','user_decision_text','metadata_sources']
m=m.merge(e[efields],on='sample_uid',validate='one_to_one')
def tr(s):return s.astype(str).str.lower().isin(['true','1'])
audit=[]
for x in m.to_dict('records'):
 tech=str(x['unresolved_technical_aggregation']).lower()=='true';distinct=str(x['unreviewed_distinct_specimen_profile_merge']).lower()=='true'
 if not (tech or distinct):status='NO_CURRENT_AGGREGATION_BLOCK';rec='Retain existing identity; normal profile rules'
 elif x['replicate_type']=='independent_capture_same_donor_user_reported':
  status='USER_CONFIRMED_INDEPENDENT_CAPTURE';rec='Propose technical count aggregation within recorded biological sample; one donor, not two independent people; user-reported evidence retained'
 elif 'HOLD_PAIR' in x['primary_inference_constraint']:
  status='HOLD_CELL_IDENTITY_RECONCILIATION';rec='No sum or automatic deletion; trace capture/library and possible same-cell sequencing overlap'
 elif any(t in x['replicate_type'] for t in ['tissue_layers','anatomical_regions','valve_layers']):
  status='DISTINCT_REGION_OR_LAYER';rec='Not technical duplicates; retain specimen/layer; donor-block dependence; pooled profile changes weighting and requires a stated estimand'
 elif distinct:
  status='DISTINCT_SPECIMENS_SAME_DONOR_UNRESOLVED_POOLING';rec='Trace specimen relation before donor-profile pooling; no invented independent donor'
 else:
  status='INHERITED_TECHNICAL_OR_SEQUENCING_RELATION_REQUIRES_SOURCE';rec='Current flags alone do not prove disjoint capture; review source library/cell identity before count summation'
 audit.append({**x,'audit_status':status,'recommended_handling':rec})
a=pd.DataFrame(audit);a.to_csv(W/'all247_sample_aggregation_audit.tsv',sep='\t',index=False)
a[a.audit_status.ne('NO_CURRENT_AGGREGATION_BLOCK')].to_csv(W/'aggregation_70_records_detailed.tsv',sep='\t',index=False)
keys=['level1_annotation','level2_annotation']
c=pd.read_csv(phase/'counts/sample_state_counts.tsv',sep='\t',keep_default_na=False)
ont=pd.read_csv(phase/'counts/ontology.tsv',sep='\t',keep_default_na=False)
named=ont[~ont.level1_annotation.isin(['Unknown','Not_assessed','Pending'])&~ont.level2_annotation.isin(['Unknown','Not_assessed','Pending'])].copy()
x=c[c['count']>0].merge(m[['sample_uid','study','donor_uid','profile','n_cells']],on='sample_uid',validate='many_to_one')
nodes=[]
for ident,v in x.groupby(keys):
 if ident not in set(zip(named.level1_annotation,named.level2_annotation)):continue
 bysample=v.groupby('sample_uid')['count'].sum();bystudy=v.groupby('study')['count'].sum()
 nodes.append(dict(zip(keys,ident),n_cells=int(v['count'].sum()),n_samples=v.sample_uid.nunique(),n_studies=v.study.nunique(),n_recorded_donors=v.donor_uid.nunique(),n_provisional_profiles=v.profile.nunique(),max_sample_fraction=float(bysample.max()/bysample.sum()),max_study_fraction=float(bystudy.max()/bystudy.sum()),frequency_pass_full_provisional=bool(v['count'].sum()>=500 and v.profile.nunique()>=10 and v.study.nunique()>=3)))
node=named.merge(pd.DataFrame(nodes),on=keys,how='left',validate='one_to_one',suffixes=('_ontology',''))
anchors={(r.level1_annotation,r.level2_annotation):str(r.marker_anchor) for r in named.itertuples()}
extract=[];top=[]
path=phase/'markers/L2_within_parent_marker_evidence.tsv.gz'
for z in pd.read_csv(path,sep='\t',chunksize=70000):
 mask=[str(g)==anchors.get((l1,l2),'') for g,l1,l2 in zip(z.gene_id,z.level1_annotation,z.level2_annotation)]
 extract.append(z.loc[mask])
 top.append(z[z.descriptive_rank<=20])
anchor=pd.concat(extract,ignore_index=True);tops=pd.concat(top,ignore_index=True)
anchor.to_csv(W/'L2_anchor_within_parent_evidence.tsv',sep='\t',index=False)
tops.to_csv(W/'L2_top20_within_parent_descriptive_markers.tsv',sep='\t',index=False)
node=node.merge(anchor[keys+['gene_id','mean_log1p','detection_fraction','mean_log1p_rest','detection_fraction_rest','mean_difference','descriptive_rank']],on=keys,how='left',validate='one_to_one')
node['anchor_positive_within_parent']=node.mean_difference>0
node['review_status']=np.where(~node.frequency_pass_full_provisional,'FREQUENCY_INSUFFICIENT',np.where(node.label_scope.eq('marker_parent_group'),'BROAD_MARKER_PARENT_GROUP_REVIEW',np.where(node.mean_difference.isna()|node.mean_difference.le(0),'ANCHOR_SPECIFICITY_REVIEW','MULTIGENE_AND_REPRODUCIBILITY_REVIEW')))
node['eligible_approved']=False
node['evidence_limit']='Descriptive all-cell marker means do not establish donor-level reproducibility or coherent per-cell programs; anchor alone does not determine eligibility'
node.to_csv(W/'all116_CoVarNet_node_evidence_review.tsv',sep='\t',index=False)
specs=pd.read_csv(R/'downstream_CPU_L1L2_v2_1_release_compat_20260914/configs/contrasts.tsv',sep='\t',keep_default_na=False)
design=[]
for spec in specs.to_dict('records'):
 g=m.copy()
 for k in ['study','vascular_bed','cohort_disease_context','assay_type']:g=g[g[k].eq(spec[k])]
 g=g[g.comparison_group.isin(['case','reference'])]
 ca=g[g.comparison_group.eq('case')];re=g[g.comparison_group.eq('reference')]
 cd=set(ca.donor_uid);rd=set(re.donor_uid);pairs=cd&rd
 paired=str(spec['paired']).lower()=='true'
 nmin=int(spec['min_independent_per_arm'])
 both=bool(cd and rd);enough=len(pairs)>=nmin if paired else len(cd)>=nmin and len(rd)>=nmin
 repeated=g.groupby(['donor_uid','comparison_group']).biological_sample_uid.nunique().gt(1).sum()
 status='NO_WITHIN_STRATUM_TWO_ARMS' if not both else 'BELOW_CURRENT_DONOR_THRESHOLD' if not enough else 'MODEL_COVARIATES_AND_STATE_COVERAGE_REVIEW'
 notes=[]
 if repeated:notes.append('Multiple specimens per donor/arm; not independent replicates')
 if g.age.eq('').all():notes.append('Age entirely missing; age-adjusted model not estimable')
 if ca.sex.nunique()==1 and re.sex.nunique()>1:notes.append('Case/control sex imbalance; adjust only under prespecified identifiable model')
 if paired:notes.append('Paired tissue contrast, not automatically disease-person vs healthy-person')
 if g.audit_status.isin(['HOLD_CELL_IDENTITY_RECONCILIATION']).any() if 'audit_status' in g else False:notes.append('Cell-identity hold')
 design.append(dict(**spec,n_case_donors_actual=len(cd),n_reference_donors_actual=len(rd),n_complete_pairs_actual=len(pairs),n_case_records=len(ca),n_reference_records=len(re),repeated_donor_arm_groups=int(repeated),case_diagnoses=json.dumps(sorted(set(ca.source_diagnosis)),ensure_ascii=False),reference_diagnoses=json.dumps(sorted(set(re.source_diagnosis)),ensure_ascii=False),case_age_available=int(ca.age.ne('').sum()),reference_age_available=int(re.age.ne('').sum()),case_sex_distribution=json.dumps(ca.sex.value_counts().to_dict()),reference_sex_distribution=json.dumps(re.sex.value_counts().to_dict()),current_audit_status=status,limitations='; '.join(notes),models_run=False))
pd.DataFrame(design).to_csv(W/'all38_clinical_contrast_design_audit.tsv',sep='\t',index=False)
summary=dict(status='CURRENT_EVIDENCE_AUDIT_COMPLETE_OPEN_ITEMS_REMAIN',models_started=False,samples=247,aggregation_review_counts=a.audit_status.value_counts().to_dict(),nodes=116,node_review_counts=node.review_status.value_counts().to_dict(),clinical_contrasts=len(design),clinical_status_counts=pd.DataFrame(design).current_audit_status.value_counts().to_dict(),no_full_atlas_expression_loaded=True,small_exploratory_model_authorization=False)
(W/'AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
sources=[prov,phase/'counts/sample_state_counts.tsv',phase/'counts/ontology.tsv',path,Path(cfg['cohort_output'])]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
(W/'SOURCE_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in sources},indent=2))
print(json.dumps(summary,indent=2,ensure_ascii=False))
