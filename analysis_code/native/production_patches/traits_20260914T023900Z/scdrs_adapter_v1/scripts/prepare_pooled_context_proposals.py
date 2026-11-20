"""Explicit alternative preserving original study-adjusted pooled-context edgeR design."""
import json,sys
from pathlib import Path
import pandas as pd
p=Path(__file__).resolve().parents[1]
m=pd.read_csv(p/'configs/cohort_design.tsv',sep='\t');out=[]
for (context,bed,assay),g in m.groupby(['cohort_disease_context','vascular_bed','assay_type'],dropna=False):
 studies=[]
 for study,s in g.groupby('study'):
  counts=s.groupby('comparison_group').donor_uid.nunique()
  if counts.get('case',0)>=2 and counts.get('reference',0)>=2:studies.append(study)
 if len(studies)<2:continue
 out.append({'contrast_id':f'pooled_context_{len(out):03d}','study':'','study_allowlist':json.dumps(studies),'vascular_bed':bed,'cohort_disease_context':context,'assay_type':assay,'case':'case','reference':'reference','formula':'~ study + group','covariates':'study','min_independent_per_arm':5,'min_per_arm_per_study':2,'paired':False,'approved':False,'meta_context':context,'review_reason':'Original within-context study-adjusted model; repeated donors require explicitly reviewed block formula, no condition-only fallback'})
pd.DataFrame(out,columns=['contrast_id','study','study_allowlist','vascular_bed','cohort_disease_context','assay_type','case','reference','formula','covariates','min_independent_per_arm','min_per_arm_per_study','paired','approved','meta_context','review_reason']).to_csv(p/'configs/pooled_context_contrasts_PROPOSED.tsv',sep='\t',index=False)
print('Candidate pooled contexts:',len(out))
