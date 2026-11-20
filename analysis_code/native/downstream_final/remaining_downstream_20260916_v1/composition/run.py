import os,json,hashlib,subprocess,sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control
N=Path(__file__).resolve().parent; W=N.parent.parent
SIDE=W/'disease_background_3groups_review_v2/sample_background_sidecar_v2.tsv'
SHA='852c08353c53a1b4d32b0b5928c9d1ac6a0278e365748912adfaf2ff30cf2ce4'
RS='/home/anzhen_vhd/miniconda3/envs/proportion_r/bin/Rscript'
def dump(p,x): p.write_text(json.dumps(x,indent=2,ensure_ascii=False))
def gate(d):
 d=d.copy(); d['comparison_group']=d.group_binary.map({'healthy':'reference','disease':'disease','reference':'reference'})
 keep=d.comparison_group.notna() & d.donor_uid.notna() & d.biological_sample_uid.notna() & d.donor_resolved.astype(str).str.lower().eq('true') & (d.n_parent>=100)
 for c in ['study','tissue','assay_type']: keep &= d[c].notna() & d[c].astype(str).ne('')
 d=d[keep].copy()
 d['stratum']=d.study.astype(str)+' | '+d.tissue.astype(str)+' | '+d.assay_type.astype(str)
 if d.biological_sample_uid.duplicated().any():
  duplicated=d[d.biological_sample_uid.duplicated(False)]
  allowed={'GSE210248__GSM6424655_rep1','GSE210248__GSM6424655_rep2'}
  if 'sample_uid' not in d or not set(duplicated.sample_uid)<=allowed or not duplicated.replicate_relationship.eq('independent_capture_same_donor').all():return d,'UNKNOWN_DUPLICATE_BIOLOGICAL_UNIT'
 if 'analysis_profile_uid' not in d:d['analysis_profile_uid']=d.biological_sample_uid
 if d.analysis_profile_uid.isna().any() or d.analysis_profile_uid.duplicated().any():return d,'INVALID_OBSERVATION_UID'
 support=d.groupby(['stratum','comparison_group']).donor_uid.nunique().unstack(fill_value=0).reindex(columns=['reference','disease'],fill_value=0)
 good=support.index[(support>=3).all(axis=1)];d=d[d.stratum.isin(good)].copy()
 if d.stratum.nunique()<1:return d,'NO_COMPARABLE_STUDY_TISSUE_ASSAY_STRATUM'
 if d.donor_uid.nunique()<8:return d,'FEWER_THAN_EIGHT_DONORS_TOTAL'
 return d,None

def spec(d,mid,level,state,family):
 f='cbind(n_state,n_parent-n_state) ~ comparison_group'
 if d.stratum.nunique()>1:f+=' + stratum'
 if d.donor_uid.duplicated().any():f+=' + (1|donor_uid)'
 return dict(approved=True,formula=f,tested_terms=['comparison_groupdisease'],reference_levels={'comparison_group':'reference'},min_denominator=100,min_units=8,min_donors=8,model_id=mid,FDR_family=family,state=state,level=level)

def fit(d,s,out):
 inp=out.with_suffix('.tsv'); js=out.with_suffix('.json');d.to_csv(inp,sep='\t',index=False);dump(js,s)
 with out.with_suffix('.log').open('w') as log:r=subprocess.run([RS,str(N/'composition_glmm_capture_safe.R'),str(inp),str(js),str(out)],stdout=log,stderr=subprocess.STDOUT)
 return r.returncode

def main():
 assert hashlib.sha256(SIDE.read_bytes()).hexdigest()==SHA
 mode=sys.argv[1] if len(sys.argv)>1 else 'prepare'
 out=N/('production' if mode=='run' else 'preparation');out.mkdir(exist_ok=(mode!='run'))
 m=pd.read_csv(SIDE,sep='\t');assert not m.sample_uid.duplicated().any()
 assert len(m)==247
 mainmask=m.main_analysis_scope_candidate.astype(str).str.lower().eq('true');assert mainmask.sum()==241
 m['composition_scope']=mainmask;m.to_csv(out/'scope_all247.tsv',sep='\t',index=False)
 audit=[]; effects=[];i=0;hashes={str(SIDE):SHA}
 for level in ['L1','L2']:
  p=W/'descriptive_updated_v2'/f'{level}_sample_composition.tsv.gz';hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
  original=pd.read_csv(p,sep='\t');counts=['sample_uid','state','n_cells']+(['n_all_qc','n_parent','level1_annotation','level2_annotation'] if level=='L2' else ['level1_annotation'])
  d=original[counts].merge(m,on='sample_uid',validate='many_to_one',suffixes=('','_clinical'))
  assert d.n_cells.notna().all() and (d.n_cells>=0).all()
  if level=='L1':
   totals=d.groupby('sample_uid').n_cells.sum();assert (totals==m.set_index('sample_uid').n_cells.reindex(totals.index)).all()
   d['n_all_qc']=d.sample_uid.map(totals);d['n_parent']=d.n_all_qc
  else:
   assert (d.groupby('sample_uid').n_cells.sum()==d.groupby('sample_uid').n_all_qc.first()).all()
  d['fraction_all_cells']=d.n_cells/d.n_all_qc;d['fraction_within_parent']=d.n_cells/d.n_parent.replace(0,np.nan)
  d.to_csv(out/f'{level}_sample_counts_latest.tsv.gz',sep='\t',index=False)
  for fields,name in [(['proposed_background_id','clinical_group_CN'],'three_backgrounds'),(['diagnosis_subtype_preserved','clinical_group_CN'],'fine_context_updated')]:
   z=d[d.composition_scope].groupby(fields+['donor_uid','state'],dropna=False).fraction_all_cells.mean().reset_index()
   z.groupby(fields+['state'],dropna=False).agg(mean_donor_fraction=('fraction_all_cells','mean'),n_donors=('donor_uid','nunique')).to_csv(out/f'{level}_{name}_donor_mean.tsv',sep='\t')
  for denominator in (['all_cells'] if level=='L1' else ['within_parent','all_cells_sensitivity']):
   for (bg,state),s in d[d.composition_scope].groupby(['proposed_background_id','state']):
    s=s.copy();s['n_state']=s.n_cells
    if denominator!='within_parent':s['n_parent']=s.n_all_qc
    eligible,reason=gate(s); mid=f'model_{i:05d}';i+=1
    rec=dict(model_id=mid,level=level,background=bg,state=state,denominator=denominator,n_samples=len(eligible),n_donors=eligible.donor_uid.nunique(),status='NOT_ESTIMABLE' if reason else 'PREPARED',reason=reason)
    specx=spec(eligible,mid,level,state,f'{level}_{bg}_{denominator}');rec['formula']=specx['formula'];rec['scope_kind']='single_homogeneous_stratum' if eligible.stratum.nunique()==1 else 'adjusted_multiple_strata'
    excluded=s[~s.sample_uid.isin(eligible.sample_uid)].copy();excluded.to_csv(out/f'{mid}_excluded.tsv',sep='\t',index=False)
    if not reason:
     if mode=='run':
      code=fit(eligible,specx,out/mid);rec['exit_code']=code
      if code==0:effects.append(pd.read_csv(out/mid/'effects.tsv',sep='\t'));rec['status']='COMPLETE'
      else:rec['status']='MODEL_FAILED';rec['reason']=(out/f'{mid}.log').read_text()[-4000:]
     else:eligible.to_csv(out/f'{mid}.tsv',sep='\t',index=False);dump(out/f'{mid}.json',specx)
    audit.append(rec)
 pd.DataFrame(audit).to_csv(out/'model_audit.tsv',sep='\t',index=False)
 if effects:
  e=pd.concat(effects);e['FDR_BH']=e.groupby('FDR_family').P.transform(lambda x:false_discovery_control(x));e.to_csv(out/'effects_BH.tsv',sep='\t',index=False)
 dump(out/'receipt.json',dict(status='PREPARED_NOT_RUN' if mode!='run' else ('PARTIAL_WITH_FAILED_MODELS' if any(a['status']=='MODEL_FAILED' for a in audit) else 'COMPLETE_WITH_MODEL_AUDIT'),source_sha256=hashes,model_status=pd.Series([a['status'] for a in audit]).value_counts().to_dict(),patched_R_sha256=hashlib.sha256((N/'composition_glmm_capture_safe.R').read_bytes()).hexdigest(),scope_n=241,excluded_unmatched=6))
if __name__=='__main__':main()
