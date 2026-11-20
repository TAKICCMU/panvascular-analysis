"""Metadata-only proposal expansion after approved ontology release; never fits models."""
import argparse,json
from pathlib import Path
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--ontology',required=True);p.add_argument('--out',required=True);a=p.parse_args()
out=Path(a.out)
if out.exists():raise SystemExit('Output already exists')
o=pd.read_csv(a.ontology,sep='\t',keep_default_na=False)
if not {'level1_annotation','level2_annotation','approved'}.issubset(o):raise SystemExit('Wrong ontology schema')
o=o[~o.level2_annotation.isin(['Unknown','Not_assessed'])]
if o.level2_annotation.duplicated().any():raise SystemExit('Named state must have one parent')
if not o.approved.astype(str).str.lower().isin(['true','1']).all():raise SystemExit('Ontology not approved')
models=[]
for row in o.itertuples():
 for pred in ['age','sex']:
  models.append({'model_id':f'L2_{len(models):04d}_{pred}','level':'L2','state':row.level2_annotation,'parent':row.level1_annotation,'approved':False,'formula':f'cbind(n_state,n_parent-n_state) ~ {"age_z" if pred=="age" else "sex"} + comparison_group + (1|study) + (1|donor_uid)','tested_terms':['age_z'] if pred=='age' else ['sexmale'],'reference_levels':{'sex':'female','comparison_group':'reference'},'filters':{},'min_denominator':50,'min_units':8,'min_donors':8,'FDR_family':'L2_'+pred+'_primary','note':'Review fixed-effect confounding and minimum sex representation before approval'})
out.write_text(json.dumps(models,indent=2,ensure_ascii=False))
print('Proposals only; no model run:',len(models))
