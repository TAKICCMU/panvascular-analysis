"""New frozen-state count ledger -> original CoVarNet profile matrices. No NMF."""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
# Reuse validated donor parsing, not a second donor identity inference.
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'composition/code'))
from composition_inputs import unique_meta,boolean,filehash

def build(ann,overlay,samples,min_total=100,min_node_cells=500,min_node_samples=10,min_node_studies=3):
    if 'Level2_CoVarNet_eligible' not in ann.columns:raise ValueError('Independent per-cell Level2_CoVarNet_eligible is required; primary eligibility or node label must not substitute')
    if ann.cell_id.duplicated().any() or overlay.cell_id.duplicated().any() or not ann.cell_id.equals(overlay.cell_id):raise ValueError('Annotation/overlay IDs not unique and exactly ordered')
    if samples.sample_uid.duplicated().any():raise ValueError('Nonunique technical sample metadata')
    meta=unique_meta(samples);lookup=samples.set_index('sample_uid');bio=meta.set_index('analysis_unit_uid')
    cells=ann.copy();cells['sample_uid']=overlay.sample_uid.to_numpy()
    if not cells.sample_uid.isin(lookup.index).all():raise ValueError('Unmapped technical sample')
    actual=cells.sample_uid.value_counts()
    if not samples.sample_uid.map(actual).fillna(0).astype(int).eq(samples.n_cells.astype(int)).all():raise ValueError('Real cell counts differ from sample contract')
    cells['biological_sample_uid']=cells.sample_uid.map(lookup.biological_sample_uid)
    for col in ['study_canonical','condition_v2','tissue_v2','sex_v2','age_v2','donor_uid_analysis','donor_identity_resolved','resolved_person_key','donor_resolution']:
        cells[col]=cells.biological_sample_uid.map(bio[col])
    for col in ['study_canonical','condition_v2','tissue_v2']:
        if cells[col].isna().any():raise ValueError('Missing identity field '+col)
    cells['covarnet_donor_key']=cells.resolved_person_key
    cells['covarnet_sample_uid']=cells.study_canonical+'::DONOR::'+cells.covarnet_donor_key+'::TISSUE::'+cells.tissue_v2+'::CONDITION::'+cells.condition_v2
    rows=[]
    for uid,g in cells.groupby('covarnet_sample_uid',sort=True,observed=True):
        row={'covarnet_sample_uid':uid,'n_total_atlas_cells':len(g),'biological_sample_uids':';'.join(sorted(g.biological_sample_uid.unique())),'n_biological_sample_uids_collapsed':g.biological_sample_uid.nunique(),'technical_sample_uids':';'.join(sorted(g.sample_uid.unique()))}
        for col in ['study_canonical','covarnet_donor_key','condition_v2','tissue_v2','sex_v2','age_v2','donor_uid_analysis','donor_resolution']:
            if g[col].nunique(dropna=False)!=1:raise ValueError('Profile identity/covariate conflict '+col)
            row[col]=g[col].iloc[0]
        row['donor_uid_v2']=row['donor_uid_analysis'] # G08 adapter receives cleaned key; original values remain sample contract.
        row['covarnet_donor_resolution']='reported_donor_uid' if row['donor_resolution']=='resolved' else 'biological_sample_uid_fallback'
        row['condition_metadata_source_v2']='current_sample_metadata_contract'
        rows.append(row)
    all_profiles=pd.DataFrame(rows);profiles=all_profiles[all_profiles.n_total_atlas_cells>=min_total].copy()
    primary=boolean(cells.Level2_primary_eligible)
    covarnet=boolean(cells.Level2_CoVarNet_eligible)
    if (covarnet & ~primary).any():raise ValueError('CoVarNet eligible cells must also be primary eligible; inconsistent per-cell eligibility')
    if (primary & (cells.Level2_rebuilt.isna()|cells.Level2_rebuilt.eq('Unknown'))).any():raise ValueError('Eligible Unknown L2')
    if (covarnet & (cells.Level1_rebuilt.isna()|cells.Level1_rebuilt.eq('Unknown'))).any():raise ValueError('CoVarNet eligible cells require a resolved L1 parent')
    flags=primary & covarnet
    cells['Level2_primary_eligible']=primary;cells['Level2_CoVarNet_eligible']=covarnet
    eligible=cells[flags & cells.covarnet_sample_uid.isin(profiles.covarnet_sample_uid)].copy()
    audit=eligible.groupby('Level2_rebuilt',observed=True).agg(n_cells=('cell_id','size'),n_samples=('covarnet_sample_uid','nunique'),n_studies=('study_canonical','nunique'),n_parents=('Level1_rebuilt','nunique'),Level1=('Level1_rebuilt','first')).reset_index().rename(columns={'Level2_rebuilt':'cell_state'})
    if audit.n_parents.gt(1).any():raise ValueError('Node has multiple parents')
    audit['passes_frequency_thresholds']=(audit.n_cells>=min_node_cells)&(audit.n_samples>=min_node_samples)&(audit.n_studies>=min_node_studies)
    states=sorted(audit.loc[audit.passes_frequency_thresholds,'cell_state']);valid=eligible[eligible.Level2_rebuilt.isin(states)]
    counts=pd.crosstab(valid.Level2_rebuilt,valid.covarnet_sample_uid).reindex(index=states,columns=profiles.covarnet_sample_uid,fill_value=0).astype(int)
    counts.index.name='cell_state';den=profiles.set_index('covarnet_sample_uid').n_total_atlas_cells.reindex(counts.columns)
    prop=counts.divide(den,axis=1);closed=counts.divide(counts.sum(axis=0).replace(0,np.nan),axis=1).fillna(0)
    if not (counts.sum(axis=0)<=den).all():raise ValueError('Node counts exceed all-atlas denominator')
    assert int(all_profiles.n_total_atlas_cells.sum())==len(cells)
    ledger=cells.groupby(['Level1_rebuilt','Level2_rebuilt'],observed=True).size().rename('n_cells').reset_index()
    eligibility_ledger=cells.groupby(['Level1_rebuilt','Level2_rebuilt','Level2_primary_eligible','Level2_CoVarNet_eligible'],observed=True).size().rename('n_cells').reset_index()
    mapping=audit.loc[audit.passes_frequency_thresholds,['cell_state','Level1']].rename(columns={'Level1':'parent_lineage'}).sort_values('cell_state')
    return dict(counts=counts,prop=prop,closed=closed,profiles=profiles,all_profiles=all_profiles,node_audit=audit,mapping=mapping,ledger=ledger,eligibility_ledger=eligibility_ledger)

def write_inputs(result,out,source_hashes):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    names={'counts':'CoVarNet_node_counts.tsv','prop':'CoVarNet_node_proportion_of_all_atlas_cells.tsv','closed':'CoVarNet_node_proportion_among_stable_nodes_sensitivity.tsv','profiles':'CoVarNet_sample_metadata.tsv','all_profiles':'all_profile_count_ledger.tsv','node_audit':'CoVarNet_node_replication_audit.tsv','mapping':'CoVarNet_node_parent_map.tsv','ledger':'all_cell_state_count_ledger.tsv','eligibility_ledger':'per_state_dual_eligibility_count_ledger.tsv'}
    for key,name in names.items():result[key].to_csv(out/name,sep='\t',index=key in ['counts','prop','closed'])
    prop=result['prop'];variable=prop.max(axis=1).gt(prop.min(axis=1)) if len(prop) else pd.Series(dtype=bool)
    manifest=dict(stage='inputs_prepared',n_all_cells=int(result['all_profiles'].n_total_atlas_cells.sum()),n_profiles=len(result['profiles']),n_candidate_nodes=len(result['node_audit']),n_retained_nodes=len(prop),n_variable_nodes=int(variable.sum()),rank_max=min(25,int(variable.sum())-1,len(result['profiles'])-1),nmf_dimension_ready=len(result['profiles'])>=20 and int(variable.sum())>=5,parameters=dict(min_profile_cells=100,min_node_cells=500,min_node_profiles=10,min_node_studies=3),all_atlas_denominator=True,unknown_donors='sample-specific biological placeholder, never same NR person',source_sha256=source_hashes,files={name:filehash(out/name) for name in names.values()})
    manifest.update(per_cell_eligibility='strict Level2_primary_eligible AND independent Level2_CoVarNet_eligible before any state/profile aggregation',missing_independent_covarnet_flag='hard error; never infer from primary flag or state name')
    (out/'INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2));return manifest
