"""Observed zeros and undefined within-parent proportions remain distinct."""
import numpy as np
import pandas as pd
from .core import Blocked,missing

def state_counts(obs,unit='sample_uid',l1='level1_annotation',l2='level2_annotation'):
    pairs=obs[[l1,l2]].drop_duplicates().sort_values([l1,l2])
    units=pd.DataFrame({unit:sorted(obs[unit].astype(str).unique())})
    grid=units.merge(pairs,how='cross')
    n=obs.groupby([unit,l1,l2],observed=True).size().rename('n_cells').reset_index()
    out=grid.merge(n,on=[unit,l1,l2],how='left',validate='one_to_one');out.n_cells=out.n_cells.fillna(0).astype('int64')
    out['count']=out.n_cells
    total=obs.groupby(unit,observed=True).size();parent=obs.groupby([unit,l1],observed=True).size()
    out['n_all_qc']=out[unit].map(total).astype('int64')
    out['n_parent']=[int(parent.get((u,p),0)) for u,p in zip(out[unit],out[l1])]
    out['prop_all_qc']=out.n_cells/out.n_all_qc
    out['prop_within_parent']=out.n_cells/out.n_parent.replace(0,np.nan)
    out['within_parent_estimable']=out.n_parent.gt(0)
    out['state_observed']=out.n_cells.gt(0)
    out['biological_state']=~out[l2].isin(['Unknown','Not_assessed'])
    assigned=obs.loc[~obs[l2].isin(['Unknown','Not_assessed'])].groupby(unit,observed=True).size()
    out['n_assigned_only']=out[unit].map(assigned).fillna(0).astype('int64')
    out['prop_assigned_only_sensitivity']=np.where(out.biological_state,out.n_cells/out.n_assigned_only.replace(0,np.nan),np.nan)
    return out

def biological_counts(obs,cohort):
    eligible=cohort.loc[cohort.aggregate_approved.astype(bool)]
    sub=obs.loc[obs.sample_uid.isin(eligible.sample_uid)].copy()
    sub['biological_sample_uid']=sub.sample_uid.map(eligible.set_index('sample_uid').biological_sample_uid)
    if sub.empty:raise Blocked('No approved biological-sample units')
    out=state_counts(sub,unit='biological_sample_uid')
    out['scope']='APPROVED_NONOVERLAPPING_PARTITIONS_ONLY'
    return out

def coverage_table(obs,cohort,l1='level1_annotation',l2='level2_annotation'):
    m=obs[[l1,l2,'sample_uid']].merge(cohort[['sample_uid','donor_uid','study','tissue']],on='sample_uid',validate='many_to_one')
    m['donor_uid']=m.donor_uid.mask(m.donor_uid.map(missing))
    return m.groupby([l1,l2],observed=True).agg(n_cells=('sample_uid','size'),n_samples=('sample_uid','nunique'),n_donors=('donor_uid','nunique'),n_studies=('study','nunique'),n_tissues=('tissue','nunique')).reset_index()
