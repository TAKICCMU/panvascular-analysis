"""Final-atlas CoVarNet preparation and inherited CPU nsNMF adapter.

No h5ad is opened here. Inputs are audited sample-state counts and cohort tables.
`execute` must be called only by the release/approval-gated orchestrator.
"""
from pathlib import Path
import hashlib, json, os, subprocess, sys
import numpy as np
import pandas as pd
from .state_identity import build_state_identity_map

MISSING = {'', 'nan', 'none', 'na', 'nr', 'not_available', 'unknown'}
EXCLUDED = {'Unknown', 'Not_assessed', 'Pending'}
DEFAULT = dict(min_profile_cells=100, min_node_cells=500, min_node_profiles=10,
               min_node_studies=3, rank_min=2, rank_max=25, seed=20260806,
               rank_nrun=50, final_nrun=200, stability_nrun=30,
               n_bootstrap=30, parent_null_permutations=2000, parent_null_seed=20260826,
               threads_per_worker=20, workers=2)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024), b''): h.update(b)
    return h.hexdigest()

def boolean(s):
    x=s.astype(str).str.lower()
    if not x.isin(['true','false','1','0']).all(): raise ValueError('Explicit boolean required')
    return x.isin(['true','1'])

def require(frame, cols, label):
    missing=set(cols)-set(frame.columns)
    if missing: raise ValueError(f'{label}: missing columns {sorted(missing)}')

def prepare(counts, cohort, ontology, options=None):
    """Return matrices and an explicit eligibility review; preserves zero vs NA.

    ontology: level1_annotation, level2_annotation; covarnet_eligible optional.
    Missing eligibility produces suggestions and BLOCKED status unless the user
    explicitly approves the documented named-state+replication rule in options.
    Donor fallback is a biological-sample placeholder, never an inferred person.
    """
    cfg={**DEFAULT, **(options or {})}; c=counts.copy(); m=cohort.copy(); ont=ontology.copy()
    keys=['level1_annotation','level2_annotation']
    require(c,['sample_uid',*keys,'count'],'counts')
    require(m,['sample_uid','biological_sample_uid','donor_uid','donor_resolved','study','tissue','condition','aggregate_approved'],'cohort')
    require(ont, keys,'ontology')
    if m.sample_uid.duplicated().any(): raise ValueError('Nonunique sample metadata')
    if c.duplicated(['sample_uid',*keys]).any(): raise ValueError('Duplicate sample-state rows')
    if ont.duplicated(keys).any(): raise ValueError('Duplicate ontology identity')
    if c[keys+['sample_uid']].isna().any().any(): raise ValueError('Missing state/sample identity')
    numeric=pd.to_numeric(c['count'],errors='raise').to_numpy(float)
    if not np.isfinite(numeric).all() or (numeric<0).any() or (numeric!=np.floor(numeric)).any(): raise ValueError('Counts must be finite nonnegative integers')
    c['count']=numeric.astype(np.int64)
    if not c.sample_uid.isin(m.sample_uid).all(): raise ValueError('Unmapped sample in counts')
    if not m.sample_uid.isin(c.sample_uid).all(): raise ValueError('Cohort sample absent from cell ledger')
    if 'Pending' in set(c.level2_annotation): raise ValueError('Pending final L2 forbidden')
    for k in ['sample_uid','biological_sample_uid','study','tissue','condition']:
        if m[k].isna().any() or m[k].astype(str).str.lower().isin(MISSING).any(): raise ValueError('Missing profile identity '+k)
    joined=c.merge(ont[keys].assign(_ontology_present=True),on=keys,how='left',validate='many_to_one')
    if joined._ontology_present.isna().any(): raise ValueError('Count state missing from ontology')
    # Reviewed state_id takes priority; without explicit IDs each L2 must name one parent.
    state_identity_map=build_state_identity_map(ont)
    resolved=boolean(m.donor_resolved)
    if (resolved & (m.donor_uid.isna()|m.donor_uid.astype(str).str.lower().isin(MISSING))).any(): raise ValueError('Resolved donor marked missing')
    m['resolved_person_key']=np.where(resolved,'DONOR::'+m.donor_uid.astype(str),'UNRESOLVED_BIO::'+m.biological_sample_uid.astype(str))
    # donor_uid must already be canonical/global; do not split a known person
    # into separate bootstrap blocks merely because two studies report them.
    m['covarnet_donor_key']=m.resolved_person_key
    # JSON tuple avoids accidental delimiter collisions, age/sex do not define units.
    m['covarnet_sample_uid']=[json.dumps(v,ensure_ascii=False,separators=(',',':')) for v in zip(m.study,m.resolved_person_key,m.tissue,m.condition)]
    m['aggregate_approved']=boolean(m.aggregate_approved)
    bio_identity=m.groupby('biological_sample_uid')[['study','resolved_person_key','tissue','condition']].nunique(dropna=False)
    if bio_identity.gt(1).any().any(): raise ValueError('Biological sample identity conflicts across technical partitions')
    repeated=m.groupby('biological_sample_uid').sample_uid.transform('size').gt(1)
    if (repeated & ~m.aggregate_approved).any(): raise ValueError('Multiple technical partitions require explicit disjoint-cell aggregation approval')
    multi_bio=m.groupby('covarnet_sample_uid').biological_sample_uid.transform('nunique').gt(1)
    profile_approved=boolean(m.profile_aggregate_approved) if 'profile_aggregate_approved' in m else pd.Series(False,index=m.index)
    if (multi_bio & ~profile_approved).any(): raise ValueError('Distinct biological samples require separate profile_aggregate_approved; technical aggregation approval is insufficient')
    x=c.merge(m,on='sample_uid',validate='many_to_one')
    total_sample=c.groupby('sample_uid')['count'].sum()
    if 'n_all_qc' in c:
        if c.groupby('sample_uid').n_all_qc.nunique(dropna=False).gt(1).any(): raise ValueError('Conflicting all-QC denominators')
        if not c.groupby('sample_uid').n_all_qc.first().eq(total_sample).all(): raise ValueError('All-QC counts omit states or denominator mismatch')
    parent_counts=x.groupby(['level1_annotation','covarnet_sample_uid'])['count'].sum().unstack(fill_value=0)
    profiles=m.groupby('covarnet_sample_uid',sort=True).agg(study_canonical=('study','first'),tissue_v2=('tissue','first'),condition_v2=('condition','first'),covarnet_donor_key=('covarnet_donor_key','first'),donor_resolved=('donor_resolved','first'),n_technical_samples=('sample_uid','size')).reset_index()
    profiles['n_total_atlas_cells']=profiles.covarnet_sample_uid.map(x.groupby('covarnet_sample_uid')['count'].sum())
    all_profiles=profiles.copy(); profiles=profiles[profiles.n_total_atlas_cells>=cfg['min_profile_cells']].copy()
    profile_ids=profiles.covarnet_sample_uid.tolist()
    parent_counts=parent_counts.reindex(columns=profile_ids,fill_value=0)
    valid=x[x.covarnet_sample_uid.isin(profile_ids)&~x.level2_annotation.isin(EXCLUDED)&~x.level1_annotation.isin(EXCLUDED)&x['count'].gt(0)]
    audit=valid.groupby(keys).agg(n_cells=('count','sum'),n_profiles=('covarnet_sample_uid','nunique'),n_studies=('study','nunique')).reset_index()
    audit['passes_frequency_thresholds']=(audit.n_cells>=cfg['min_node_cells'])&(audit.n_profiles>=cfg['min_node_profiles'])&(audit.n_studies>=cfg['min_node_studies'])
    # counts adds explicit False for reserved states. If the incoming ontology
    # had no eligibility column, this creates blank values for EVERY named state;
    # preserve the existing opt-in named-state rule instead of treating that
    # transport artifact as a partially reviewed eligibility table.
    named_mask=~ont.level2_annotation.isin(EXCLUDED)&~ont.level1_annotation.isin(EXCLUDED)
    explicit_named_eligibility=False
    if 'covarnet_eligible' in ont:
        values=ont.loc[named_mask,'covarnet_eligible']
        absent=values.isna()|values.astype(str).str.strip().eq('')
        if absent.any() and not absent.all():
            raise ValueError('Partial named covarnet_eligible is missing; explicit boolean required for every named state')
        if len(values) and not absent.all():
            ont.loc[named_mask,'covarnet_eligible']=boolean(values)
            explicit_named_eligibility=True
        ont.loc[~named_mask,'covarnet_eligible']=False
    if explicit_named_eligibility:
        audit=audit.merge(ont[keys+['covarnet_eligible']],on=keys,validate='one_to_one')
        audit['covarnet_eligible']=boolean(audit.covarnet_eligible)
        eligibility_approved=True; eligibility_basis='explicit_ontology_review'
    else:
        audit['covarnet_eligible']=True
        eligibility_approved=cfg.get('approve_named_state_frequency_rule',False) is True
        eligibility_basis='named_state_frequency_rule_PROPOSED' if not eligibility_approved else 'named_state_frequency_rule_APPROVED'
    audit=audit.merge(state_identity_map[keys+['state_id']],on=keys,validate='one_to_one')
    valid=valid.merge(state_identity_map[keys+['state_id']],on=keys,validate='many_to_one')
    audit['include_node']=audit.passes_frequency_thresholds & audit.covarnet_eligible
    states=sorted(audit.loc[audit.include_node,'state_id'])
    counts_mat=valid.groupby(['state_id','covarnet_sample_uid'])['count'].sum().unstack(fill_value=0).reindex(index=states,columns=profile_ids,fill_value=0)
    counts_mat.index.name='cell_state'
    mapping=audit.loc[audit.include_node,keys+['state_id']].rename(columns=dict(state_id='cell_state',level1_annotation='parent_lineage')).set_index('cell_state').reindex(states).reset_index()
    den=profiles.set_index('covarnet_sample_uid').n_total_atlas_cells
    prop=counts_mat.divide(den,axis=1)
    closed=counts_mat.divide(counts_mat.sum(axis=0).replace(0,np.nan),axis=1)
    # Keep closure undefined when no retained state is observed; NMF adapter blocks those profiles.
    within=counts_mat.astype(float).copy()
    for row in mapping.itertuples(index=False): within.loc[row.cell_state]=counts_mat.loc[row.cell_state]/parent_counts.loc[row.parent_lineage].replace(0,np.nan)
    variable=prop.max(axis=1)>prop.min(axis=1)
    blocks=[]
    if not eligibility_approved: blocks.append('REVIEW_COVARNET_NODE_ELIGIBILITY')
    if len(profiles)<20: blocks.append('FEWER_THAN_20_PROFILES')
    if int(variable.sum())<5: blocks.append('FEWER_THAN_5_VARIABLE_STATES')
    if counts_mat.sum(axis=0).eq(0).any(): blocks.append('RETAINED_NODE_DENOMINATOR_ZERO_REVIEW_PROFILE_EXCLUSION')
    return dict(counts=counts_mat,prop=prop,closed=closed,within_parent=within,parents=parent_counts,
                profiles=profiles,all_profiles=all_profiles,node_audit=audit,mapping=mapping,state_identity_map=state_identity_map,
                sample_profile_map=m[['sample_uid','biological_sample_uid','covarnet_sample_uid','covarnet_donor_key','donor_resolved']],
                manifest=dict(status='BLOCKED' if blocks else 'PREPARED_NOT_EXECUTED',phase='INPUTS_BLOCKED' if blocks else 'INPUTS_READY',blocking_reasons=blocks,parameters=cfg,eligibility_basis=eligibility_basis,n_cells=int(c['count'].sum()),n_profiles=len(profiles),n_states=len(states),n_variable_states=int(variable.sum()),unknown_donor_policy='biological_sample_placeholder_not_inferred_person'))

FILES=dict(state_identity_map='CoVarNet_state_identity_map.tsv',counts='CoVarNet_node_counts.tsv',prop='CoVarNet_node_proportion_of_all_atlas_cells.tsv',closed='CoVarNet_node_proportion_among_stable_nodes_sensitivity.tsv',within_parent='CoVarNet_proportion_within_full_parent.tsv',parents='CoVarNet_full_parent_counts.tsv',profiles='CoVarNet_sample_metadata.tsv',all_profiles='all_profile_count_ledger.tsv',node_audit='CoVarNet_node_replication_audit.tsv',mapping='CoVarNet_node_parent_map.tsv',sample_profile_map='sample_to_covarnet_profile.tsv')

def write_prepared(result, out, provenance):
    out=Path(out); out.mkdir(parents=True,exist_ok=False)
    for key,name in FILES.items(): result[key].to_csv(out/name,sep='\t',index=key in {'counts','prop','closed','within_parent','parents'},na_rep='NA')
    manifest={**result['manifest'],'source_identity':provenance,'code_sha256':sha(__file__),'files':{n:sha(out/n) for n in FILES.values()}}
    (out/'INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    return manifest

def run_nmf(prepared, options=None):
    cfg={**DEFAULT,**(options or {})}; p=Path(prepared)
    manifest=json.loads((p/'INPUT_MANIFEST.json').read_text())
    if manifest['status']!='PREPARED_NOT_EXECUTED' or manifest.get('phase')!='INPUTS_READY': raise ValueError('CoVarNet blocked: '+str(manifest['blocking_reasons']))
    for name,h in manifest['files'].items():
        if sha(p/name)!=h: raise ValueError('Prepared input hash mismatch: '+name)
    workers=int(cfg['workers']); threads=int(cfg['threads_per_worker'])
    if workers<1 or workers>2 or threads<1 or workers*threads>40: raise ValueError('CPU budget exceeds2workers/40threads')
    env=os.environ.copy()
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']: env[name]=str(threads)
    env['CUDA_VISIBLE_DEVICES']=''; env['COVARNET_WORKERS']=str(workers)
    for k in ['rank_nrun','final_nrun','stability_nrun','n_bootstrap','parent_null_permutations','parent_null_seed','rank_min','rank_max','seed']: env['COVARNET_'+k.upper()]=str(cfg[k])
    rscript=cfg.get('rscript','Rscript'); script=Path(__file__).with_name('covarnet_nmf.R')
    with (p/'NMF.log').open('x') as log: subprocess.run([rscript,str(script),str(p)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    with (p/'PARENT_NULL.log').open('x') as log: subprocess.run([rscript,str(script.with_name('covarnet_parent_null.R')),str(p)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    result={'status':'COMPLETE','output':str(p),'input_manifest_sha256':sha(p/'INPUT_MANIFEST.json'),'code_sha256':{str(x.name):sha(x) for x in [Path(__file__),script,script.with_name('covarnet_parent_null.R'),script.with_name('covarnet_review_helpers.R')]},'parameters':cfg,'files':{str(x.relative_to(p)):sha(x) for x in p.rglob('*') if x.is_file()}}
    (p/'COVARNET_COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def execute(config, run_dir):
    r=Path(run_dir); opts=config.get('covarnet',{}); source=r/'counts/sample_state_counts.tsv'; cohort=r/'cohort/cohort_design.tsv'
    ontology=Path(opts.get('ontology_path') or r/'counts/ontology.tsv')
    if not ontology.exists(): raise FileNotFoundError('Reviewed final ontology required: '+str(ontology))
    m=pd.read_csv(cohort,sep='\t',keep_default_na=False)
    aliases=dict(opts.get('cohort_columns') or {})
    if 'condition' not in m and 'comparison_group' in m and 'condition' not in aliases:
        aliases['condition']='comparison_group'  # root cohort contract alias, recorded below
    for canonical,actual in aliases.items():
        if actual not in m: raise ValueError('Configured metadata alias absent: '+actual)
        m[canonical]=m[actual]
    result=prepare(pd.read_csv(source,sep='\t'),m,pd.read_csv(ontology,sep='\t',keep_default_na=False),opts)
    target=r/'covarnet'; manifest=write_prepared(result,target,{'files':{str(p):sha(p) for p in [source,cohort,ontology]},'cohort_column_aliases':aliases})
    if manifest['status']=='BLOCKED':
        raise RuntimeError('BLOCKED CoVarNet: '+','.join(manifest['blocking_reasons'])+'; full audit saved in '+str(target/'INPUT_MANIFEST.json'))
    if not opts.get('run_nmf',False): return {**manifest,'status':'PREPARED_NOT_EXECUTED'}
    if not config.get('execution',{}).get('run_full_analysis',False):
        raise PermissionError('Full downstream execution not approved')
    return run_nmf(target,opts)
