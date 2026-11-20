"""Entity separation, model-specific inclusion and explicit contrast proposals."""
import json
import numpy as np
import pandas as pd
from .core import Blocked,missing

FIELDS=['sample_uid','biological_sample_uid','donor_uid','study','atlas_source_id','tissue','vascular_bed','vascular_bed_group','comparison_group','cohort_disease_context','assay_type','age','sex']
def cohort_from_obs(obs,fields=FIELDS):
    absent=set(fields)-set(obs)
    if absent:raise Blocked(f'Cohort metadata absent: {absent}')
    data=obs[fields].copy()
    for c in fields:
        if c=='sample_uid':continue
        n=data.assign(_value=data[c].map(lambda v:'<MISSING>' if missing(v) else str(v))).groupby('sample_uid')['_value'].nunique()
        if (n>1).any():raise Blocked(f'{c} conflicts within samples {list(n[n>1].index)}')
    result=data.drop_duplicates('sample_uid').copy()
    result['n_cells']=result.sample_uid.map(obs.sample_uid.value_counts())
    return annotate_units(result)

def annotate_units(meta):
    m=meta.reset_index(drop=True).copy()
    if m.sample_uid.duplicated().any() or m.sample_uid.map(missing).any():raise Blocked('Invalid sample keys')
    m['donor_resolved']=~m.donor_uid.map(missing) & ~m.donor_uid.astype(str).str.contains('PROVISIONAL|UNRESOLVED|PLACEHOLDER',case=False,regex=True)
    m['biological_sample_resolved']=~m.biological_sample_uid.map(missing)
    m['age']=pd.to_numeric(m.age,errors='coerce');m['sex']=m.sex.map(lambda v:np.nan if missing(v) else v)
    n=m.groupby('biological_sample_uid',dropna=False).sample_uid.transform('size')
    m['n_sample_partitions_per_biosample']=n
    m['aggregate_approved']=n.eq(1) & m.biological_sample_resolved
    # Multiple capture/specimen records require independent evidence of non-overlap.
    m['inference_eligible']=m.donor_resolved & m.biological_sample_resolved & m.aggregate_approved
    m['inference_reason']=np.select([~m.donor_resolved,~m.biological_sample_resolved,~m.aggregate_approved],['MISSING_DONOR','MISSING_BIOLOGICAL_SAMPLE','REVIEW_TECHNICAL_PARTITIONS'],default='IDENTITY_AVAILABLE_MODEL_CHECKS_REQUIRED')
    m['design_status']='PROPOSED_REQUIRES_REVIEW'
    return m

def read_unit_review(path):
    # Preserve leading-zero IDs while retaining parsed missing masks.
    return pd.read_csv(path,sep='\t',dtype={k:'string' for k in ['sample_uid','donor_uid','biological_sample_uid']})

def identical_identity(left,right):
    """Strict value and NA-mask comparison on disposable Series only."""
    a=left.reset_index(drop=True);b=right.reset_index(drop=True)
    if len(a)!=len(b) or not np.array_equal(a.isna().to_numpy(),b.isna().to_numpy()):return False
    keep=~a.isna()
    return np.array_equal(a.loc[keep].astype('string').to_numpy(),b.loc[keep].astype('string').to_numpy())

def apply_unit_review(meta,review):
    """Review is external and keyed exactly; it cannot change supplied donor identity."""
    if meta.sample_uid.isna().any() or meta.sample_uid.duplicated().any() or review.sample_uid.isna().any() or review.sample_uid.duplicated().any() or set(review.sample_uid)!=set(meta.sample_uid):raise Blocked('Unit review must cover each sample exactly once')
    m=meta.set_index('sample_uid');r=review.set_index('sample_uid').reindex(m.index)
    for c in ['donor_uid','biological_sample_uid']:
        if c not in r or not identical_identity(m[c],r[c]):raise Blocked(f'Review cannot silently rewrite {c}')
    for c in ['aggregate_approved','inference_eligible'] + (['profile_aggregate_approved'] if 'profile_aggregate_approved' in r else []):
        m[c]=r[c].astype(str).str.lower().isin(['true','1'])
    for c in ['review_evidence','design_status']:
        if c not in r:raise Blocked(f'Missing review column {c}')
        m[c]=r[c]
    needs_approval=m.inference_eligible | (m.aggregate_approved & m.n_sample_partitions_per_biosample.gt(1)) | m.get('profile_aggregate_approved',False)
    if (needs_approval & (m.review_evidence.map(missing) | m.design_status.ne('APPROVED'))).any():raise Blocked('Positive inference/aggregation requires approved evidence')
    if (m.inference_eligible & (~m.donor_resolved | ~m.biological_sample_resolved | ~m.aggregate_approved)).any():raise Blocked('Unresolved/unsafe units cannot be inferentially included')
    for key,g in m.reset_index().groupby('biological_sample_uid'):
        if len(g)>1 and g.aggregate_approved.any() and not g.aggregate_approved.all():raise Blocked(f'Inconsistent aggregation review for {key}')
        for col in ['donor_uid','study','tissue','vascular_bed','comparison_group','assay_type','age','sex']:
            vals=g[col].map(lambda v:'<MISSING>' if missing(v) else str(v))
            if g.aggregate_approved.all() and vals.nunique()>1:raise Blocked(f'Approved biosample merge has conflicting {col}: {key}')
    return m.reset_index()

def propose_contrasts(m):
    rows=[]
    for (study,bed,context,assay),g in m.groupby(['study','vascular_bed','cohort_disease_context','assay_type'],dropna=False,sort=True):
        case=g[g.comparison_group=='case'];ref=g[g.comparison_group=='reference']
        shared=set(case.donor_uid.dropna()) & set(ref.donor_uid.dropna())
        reason='MODEL_AND_REFERENCE_DEFINITION_REVIEW' if len(case) and len(ref) else 'NO_WITHIN_STUDY_CASE_REFERENCE'
        rows.append({'contrast_id':f'contrast_{len(rows)+1:03d}','study':study,'vascular_bed':bed,'cohort_disease_context':context,'assay_type':assay,'case':'case','reference':'reference','n_case_samples':len(case),'n_reference_samples':len(ref),'n_case_donors':case.loc[~case.donor_uid.map(missing),'donor_uid'].nunique(),'n_reference_donors':ref.loc[~ref.donor_uid.map(missing),'donor_uid'].nunique(),'n_paired_donors':len(shared),'formula':'~ donor_uid + group' if shared else '~ group','covariates':'','min_independent_per_arm':5,'paired':bool(shared),'approved':False,'review_reason':reason,'reference_interpretation':'reference tissue; not assumed healthy donor','meta_context':context})
    return pd.DataFrame(rows)

def select_contrast(meta,spec):
    if str(spec.get('approved')).lower() not in ['true','1']:raise Blocked('Contrast not human approved')
    m=meta.copy()
    for c in ['study','vascular_bed','cohort_disease_context','assay_type']:
        value=spec.get(c)
        if not missing(value):m=m[m[c].astype(str)==str(value)]
    if not missing(spec.get('study_allowlist','')):
        allowed=json.loads(spec['study_allowlist'])
        if not isinstance(allowed,list) or not allowed:raise Blocked('Explicit study allowlist must be a nonempty JSON array')
        m=m[m.study.isin(allowed)]
    m=m[m['comparison_group'].isin([spec['reference'],spec['case']]) & m.inference_eligible.astype(bool)]
    if int(spec.get('min_per_arm_per_study',0) or 0)>0:
        threshold=int(spec['min_per_arm_per_study']);support=m.groupby(['study','comparison_group']).donor_uid.nunique().unstack(fill_value=0)
        if not {spec['reference'],spec['case']}.issubset(support.columns):raise Blocked('No within-study support')
        good=support.index[(support[spec['reference']]>=threshold)&(support[spec['case']]>=threshold)]
        m=m[m.study.isin(good)]
    m['group']=pd.Categorical(m['comparison_group'],categories=[spec['reference'],spec['case']],ordered=True)
    cov=[c for c in str(spec.get('covariates','')).split(';') if c and not missing(c)]
    required=['donor_uid','biological_sample_uid','group']+cov
    before=len(m);m=m.loc[~m[required].apply(lambda s:s.map(missing)).any(axis=1)].copy()
    nmin=int(spec.get('min_independent_per_arm',5))
    paired=str(spec.get('paired')).lower() in ['true','1']
    if paired:
        good=m.groupby('donor_uid',observed=True)['group'].nunique();m=m[m.donor_uid.isin(good[good==2].index)]
        if m.donor_uid.nunique()<nmin:raise Blocked('Insufficient complete independent donor pairs')
        if m.duplicated(['donor_uid','group']).any():raise Blocked('Multiple specimens per donor/arm: explicit aggregation or repeated-measures model required')
    elif m.donor_uid.duplicated().any():raise Blocked('Repeated donor in unpaired design')
    for value in [spec['reference'],spec['case']]:
        if m.loc[m.comparison_group==value,'donor_uid'].nunique()<nmin:raise Blocked(f'Insufficient independent donors: {value}')
    return m,{'n_before_covariates':before,'n_used':len(m),'n_donors':m.donor_uid.nunique(),'missing_strategy':'contrast-specific complete cases; no covariate dropping','covariates':cov,'paired':paired}
