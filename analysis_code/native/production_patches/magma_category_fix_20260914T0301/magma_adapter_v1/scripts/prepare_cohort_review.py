"""Only current observed metadata; does not use unfinished L2 for scientific outputs."""
from pathlib import Path
import sys,json
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.cohort import annotate_units,propose_contrasts
root=Path(__file__).resolve().parents[1]
m=pd.read_csv(root/'evidence/working_sample_metadata_v4.tsv',sep='\t')
c=annotate_units(m);c['identity_eligibility_suggestion']=c['inference_eligible'];c['inference_eligible']=False;c['condition']=c['comparison_group'];c['profile_aggregate_approved']=False
c['review_evidence']='Current working h5ad imported metadata v4; final release reconciliation and downstream unit review still required'
c.to_csv(root/'configs/cohort_design.tsv',sep='\t',index=False)
propose_contrasts(c).to_csv(root/'configs/contrasts.tsv',sep='\t',index=False)
models=[]
# Actual known L1 names only; L2 models are generated after release from approved ontology.
for state in json.load(open(root/'evidence/READ_ONLY_INPUT_CHECK.json'))['l1_counts']:
 for pred in ['age','sex']:
  models.append({'model_id':f'L1_{len(models):03d}_{pred}','level':'L1','state':state,'approved':False,'formula':f'cbind(n_state,n_parent-n_state) ~ {"age_z" if pred=="age" else "sex"} + comparison_group + (1|study) + (1|donor_uid)','tested_terms':['age_z'] if pred=='age' else ['sexmale'],'reference_levels':{'sex':'female','comparison_group':'reference'},'filters':{},'min_denominator':100,'min_units':8,'min_donors':8,'FDR_family':'L1_'+pred+'_primary','review_note':'Review study/assay/bed confounding and independent replication; sex study>=5male/5female; no demographic imputation'})
(root/'configs/composition_models.json').write_text(json.dumps(models,indent=2,ensure_ascii=False))
# ID field references known original gene symbols but matching must validate var mapping at release.
genes=['PECAM1','VWF','KDR','DCN','LUM','COL1A1','COL1A2','MYH11','ACTA2','TAGLN','RGS5','ABCC9','KCNJ8','C1QA','C1QC','TYROBP','FCN1','S100A8','S100A9','FCGR3B','CSF3R','CD3D','CD3E','TRAC','NKG7','GNLY','MS4A1','CD79A','JCHAIN','MZB1','TPSAB1','KIT','FCER1A','CD1C','PLP1','S100B']
pd.DataFrame({'gene_id':genes,'panel':'canonical_L1_PROPOSED','gene_id_policy':'exact var_names; validate symbol mapping; L2 use release ontology marker_anchor'}).to_csv(root/'configs/marker_panel.tsv',sep='\t',index=False)
print(json.dumps({'samples':len(c),'donor_keys':c.donor_uid.nunique(),'repeat_partition_rows':int((~c.aggregate_approved).sum()),'contrasts_proposed':len(propose_contrasts(c)),'annotation_source':'working metadata only','no_inference_executed':True}))
