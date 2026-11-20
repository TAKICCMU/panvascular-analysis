"""Inherited communication/spatial estimators loaded by AST, without top-level effects."""
import ast,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

def inherited_functions(relative,names,globals_extra=None):
    p=Path(__file__).resolve().parents[1]/'external_legacy'/relative
    tree=ast.parse(p.read_text());nodes=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in names]
    if {n.name for n in nodes}!=set(names):raise ValueError('Inherited source API mismatch')
    ns={'np':np,'pd':pd,'math':math,'Path':Path,'sparse':sparse,'cKDTree':cKDTree}
    ns.update(globals_extra or {})
    # Postponed annotations avoid importing full old analysis dependencies.
    module=ast.Module(body=[ast.ImportFrom(module='__future__',names=[ast.alias(name='annotations')],level=0)]+nodes,type_ignores=[])
    exec(compile(ast.fix_missing_locations(module),str(p),'exec'),ns)
    return ns

def communication_effects(long,design):
    """Unpaired biological-sample contrast, inherited disease-context DL meta.
    design: sample_unit, donor_uid, donor_verified, study_canonical, tissue_v2,
            disease_context, condition_class ('disease'/'control').
    long: sample_unit, interaction_key, source,target,ligand,receptor,score,
          threshold_scheme,method. Exactly one row per specimen/edge/method.
    Repeated-donor cohorts must use the CR2 runner; not independent observations.
    """
    required=['sample_unit','donor_uid','donor_verified','study_canonical','tissue_v2','disease_context','condition_class']
    if design[required].isna().any().any():raise ValueError('Missing reviewed communication design')
    if design.sample_unit.duplicated().any():raise ValueError('Duplicate biological specimen in design')
    if not set(design.condition_class)<= {'disease','control'}:raise ValueError('Explicit disease/control contrast required')
    dup=design.groupby(['study_canonical','tissue_v2','disease_context','donor_uid'],observed=True).size()>1
    if dup.any():raise ValueError('Repeated donors: use external_communication_meta.py with CR2')
    return communication_effects_allow_repeated(long,design)

def communication_effects_allow_repeated(long,design):
    if long.duplicated(['sample_unit','interaction_key','method','threshold_scheme']).any():
        raise ValueError('Repeated edge rows would inflate sample N')
    x=long.merge(design,on='sample_unit',how='inner',validate='many_to_one')
    if len(x)!=len(long):raise ValueError('Unmapped communication specimen')
    if 'signaling_category' not in x:x['signaling_category']='Not_applicable'
    keys=['method','disease_context','study_canonical','tissue_v2','condition_class','threshold_scheme',
          'interaction_key','source','target','ligand','receptor','signaling_category']
    x=x[np.isfinite(x.score)].copy();x['score_sq']=x.score**2
    stats=x.groupby(keys,observed=True,dropna=False).agg(n=('score','size'),score_sum=('score','sum'),score_sq_sum=('score_sq','sum')).reset_index()
    ns=communication_estimator_namespace()
    effects=ns['cohort_effects'](stats)
    return effects,ns['context_meta'](effects)

def communication_estimator_namespace():
    return inherited_functions('communication/disease_context_meta.py',['bh','normal_p','dl_meta','cohort_effects','context_meta'],
      {'MIN_SAMPLES_PER_CONDITION':2,'MIN_COHORTS_FOR_CONTEXT_META':2,'MIN_DIRECTION_CONSISTENCY':0.70,
       'MAX_SINGLE_COHORT_ABS_EFFECT_SHARE':0.75,'FDR_THRESHOLD':0.05})

def spatial_estimators(n_permutations=999,n_matched=100,min_spots=40,min_positive_spots=10):
    import statsmodels.api as sm
    ns=inherited_functions('spatial/spatial_autocorrelation_matched_gene.py',
      ['block_groups','zscore','neighbor_edges','moran_statistic','moran_permutation','matched_gene_roster','linear_association_z'],
      {'sm':sm,'MIN_SPOTS':min_spots,'MIN_POSITIVE_SPOTS':min_positive_spots,'K_NEIGHBORS':6,'N_MORAN_PERMUTATIONS':n_permutations,'N_MATCHED_GENES':n_matched})
    if __package__:
        from .spatial_contract import strict_fit_count_models
    else:
        from spatial_contract import strict_fit_count_models
    ns['fit_count_models']=strict_fit_count_models
    return ns

def spatial_module_projection(abundance,weights,state_map=None,abundance_key='state_id',weight_key='state_id',allow_missing=False):
    """Inherited weighted CLR projection with dynamic final state/module universe."""
    abundance=abundance.copy();weights=weights.copy()
    if state_map is not None:
        if __package__:
            from .state_identity import map_state_axis
        else:
            from state_identity import map_state_axis
        abundance.columns=map_state_axis(abundance.columns,state_map,abundance_key)
        weights.index=map_state_axis(weights.index,state_map,weight_key)
    if not np.isfinite(abundance.to_numpy()).all() or (abundance.to_numpy()<0).any():raise ValueError('Invalid spatial abundance')
    if weights.index.has_duplicates or abundance.columns.has_duplicates:raise ValueError('Ambiguous state identity')
    if not np.isfinite(weights.to_numpy()).all() or (weights.to_numpy()<0).any():raise ValueError('Invalid NMF weights')
    ns=inherited_functions('spatial/03_spatial_v4_targets_CoVarNet_validation.py',['module_scores'])
    audit=[]
    for module in weights:
        positive=weights.index[weights[module]>0];matched=positive.intersection(abundance.columns)
        total=weights[module].sum()
        audit.append({'module':module,'n_weighted_states':len(positive),'n_mapped_states':len(matched),
                      'mapped_weight_fraction':weights.loc[matched,module].sum()/total if total>0 else np.nan,
                      'missing_states':';'.join(positive.difference(abundance.columns))})
    audit=pd.DataFrame(audit)
    audit['missing_weight_fraction']=1-audit.mapped_weight_fraction
    audit['full_weight_coverage']=audit.n_weighted_states==audit.n_mapped_states
    audit['projection_mode']='descriptive_partial_explicit' if allow_missing else 'full_weight_coverage_required'
    if not audit.full_weight_coverage.all() and not allow_missing:
        exc=ValueError('Missing positive CoVarNet weight states; no silent partial projection')
        exc.audit=audit
        raise exc
    scores=ns['module_scores'](abundance,weights)
    return scores,audit
