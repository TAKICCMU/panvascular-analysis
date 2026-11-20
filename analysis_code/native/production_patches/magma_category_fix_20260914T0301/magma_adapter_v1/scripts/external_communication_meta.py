#!/usr/bin/env python3
"""Future reviewed sample-level communication contrasts and donor CR2 sensitivity."""
import argparse,json,sys,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_gate import bootstrap_cli,approved_preflight,reviewed_cohort
if __name__=='__main__':bootstrap_cli()
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_statistics import communication_effects_allow_repeated,communication_estimator_namespace

def calculate(long,design,out,rscript):
    # Designer must provide verified person IDs; no sample-name donor imputation.
    required=['sample_unit','donor_uid','donor_verified','study_canonical','tissue_v2','disease_context','condition_class']
    if design[required].isna().any().any() or design.sample_unit.duplicated().any():raise ValueError('Invalid specimen design')
    truth=design.donor_verified.astype(str).str.lower().isin(['true','1'])
    if not truth.all():raise ValueError('Donor inference BLOCKED: unresolved or provisional people')
    if not set(design.condition_class)<= {'disease','control'}:raise ValueError('Explicit condition mapping required')
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    effects,_=communication_effects_allow_repeated(long,design)
    if effects.empty:
        effects.to_csv(out/'cohort_effects.tsv',sep='\t',index=False);return
    x=long.merge(design,on='sample_unit',validate='many_to_one')
    patterns=[];lookup={}
    for i,row in effects.iterrows():
        z=x[(x.method==row.method)&(x.threshold_scheme==row.threshold_scheme)&(x.interaction_key==row.interaction_key)&
            (x.study_canonical==row.study_canonical)&(x.tissue_v2==row.tissue_v2)&(x.disease_context==row.disease_context)]
        if not z.donor_uid.duplicated().any():continue
        pid=f'p{i}';lookup[pid]=i
        for j,r in enumerate(z.itertuples()):patterns.append(dict(pattern_id=pid,sample_order=j,condition_class=r.condition_class,resolved_person_key=r.donor_uid,validation_y=r.score))
    effects['donor_model']='independent_biological_specimens'
    if patterns:
        p=out/'CR2_actual_scores.tsv';pd.DataFrame(patterns).to_csv(p,sep='\t',index=False)
        q=out/'CR2_estimates.tsv';script=Path(__file__).resolve().parents[1]/'external_legacy/communication/g08_cr2_pattern_reference.R'
        subprocess.run([rscript,str(script),str(p),str(q)],check=True)
        res=pd.read_csv(q,sep='\t')
        if set(res.pattern_id)!=set(lookup):raise ValueError('Incomplete CR2 results')
        for r in res.itertuples():
            i=lookup[r.pattern_id]
            if r.cr2_reference_status!='passed':
                effects.loc[i,'se']=np.nan;effects.loc[i,'donor_model']='CR2_not_estimable';continue
            effects.loc[i,'effect_disease_minus_control']=r.reference_estimate;effects.loc[i,'se']=r.reference_se
            effects.loc[i,'donor_model']='CR2_Satterthwaite';effects.loc[i,'df_Satt']=r.reference_df_satterthwaite;effects.loc[i,'cohort_p']=r.reference_p
    effects.to_csv(out/'cohort_effects_with_CR2.tsv',sep='\t',index=False)
    valid=effects[np.isfinite(effects.se)&(effects.se>0)]
    # A donor occurring in multiple cohort effects invalidates the independent DL meta.
    reuse=design[['donor_uid','study_canonical','tissue_v2']].drop_duplicates().groupby('donor_uid').size()
    if (reuse>1).any():raise ValueError('Cross-cohort repeated donors: cohort results saved; meta BLOCKED pending reviewed dependence model')
    communication_estimator_namespace()['context_meta'](valid).to_csv(out/'context_DL_meta.tsv',sep='\t',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--long',required=True);p.add_argument('--design',required=True);p.add_argument('--out',required=True);p.add_argument('--rscript',required=True);p.add_argument('--config',required=True);a=p.parse_args()
    cfg,obs,ontology,info=approved_preflight(a.config,Path(a.out).with_name(Path(a.out).name+'.authorization'))
    units=reviewed_cohort(cfg,obs)
    units=units[units.inference_eligible & units.aggregate_approved]
    checked=units.drop_duplicates('biological_sample_uid').set_index('biological_sample_uid')
    design=pd.read_csv(a.design,sep='\t',keep_default_na=False)
    if not set(design.sample_unit)<=set(checked.index):raise RuntimeError('Design contains unapproved inference/aggregation units')
    for col,source in [('donor_uid','donor_uid'),('study_canonical','study'),('tissue_v2','tissue'),('disease_context','cohort_disease_context')]:
        if not np.array_equal(design[col].astype(str),checked.loc[design.sample_unit,source].astype(str)):
            raise RuntimeError('Supplied communication design conflicts with reviewed final metadata: '+col)
    expected=checked.loc[design.sample_unit,'comparison_group'].map({'case':'disease','reference':'control'})
    if expected.isna().any() or not np.array_equal(expected.to_numpy(),design.condition_class):raise RuntimeError('Unapproved communication contrast class')
    design['donor_verified']=True
    calculate(pd.read_csv(a.long,sep='\t'),design,a.out,a.rscript)
