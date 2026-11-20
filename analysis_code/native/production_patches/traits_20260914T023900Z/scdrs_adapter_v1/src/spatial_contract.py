"""Spatial provenance and statistical decisions are independent hard gates."""
import hashlib,json
from pathlib import Path
import numpy as np

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()).hexdigest()
def file_digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()

def validate_spatial_provenance(record,expected,reference_genes,targets,spot_ids,sample_ids,mode='target_validation',posterior_path=None,spatial_units=None):
    """Verify content-bound sidecar against independently supplied current identities.
    Descriptive mode relaxes ONLY target exclusion, never atlas/ontology identity.
    No attestation-only shortcut: actual posterior hash and reference gene list required.
    """
    if mode not in {'target_validation','descriptive'}:raise ValueError('Explicit spatial mode required')
    if record.get('schema')!='spatial_provenance_v2':raise ValueError('Unverified/legacy posterior provenance')
    if 'locked71' in json.dumps(record).lower():raise ValueError('Old locked71 spatial input prohibited')
    for key in ['atlas_sha256','ontology_sha256','parameters_sha256','reference_artifact_sha256']:
        if not expected.get(key) or record.get(key)!=expected[key]:raise ValueError('Spatial identity mismatch: '+key)
    if posterior_path is None or record.get('posterior_sha256')!=file_digest(posterior_path):raise ValueError('Posterior content digest mismatch')
    genes=list(map(str,reference_genes));spots=list(map(str,spot_ids));samples=list(map(str,sample_ids));tg=sorted(set(map(str,targets)))
    if len(set(genes))!=len(genes) or len(set(spots))!=len(spots) or len(spots)!=len(samples):raise ValueError('Ambiguous spatial identities')
    for key,value in [('reference_genes_sha256',genes),('spot_ids_sha256',spots),('sample_ids_sha256',samples),('targets_sha256',tg)]:
        if record.get(key)!=digest(value):raise ValueError('Spatial identity mismatch: '+key)
    if spatial_units is None or record.get('spatial_units_sha256')!=digest(spatial_units):raise ValueError('Spatial section/person identity not bound')
    if mode=='target_validation':
        if not tg or record.get('reference_mode')!='target_excluded' or set(genes)&set(tg):raise ValueError('Target exclusion not proven; use explicit descriptive mode')
        if sorted(record.get('excluded_targets',[]))!=tg:raise ValueError('Target exclusion roster mismatch')
    return {'status':'COMPLETE','mode':mode,'inferential_target_validation_allowed':mode=='target_validation'}

def spatial_design_gate(design):
    required={'min_spots','min_positive_spots','interaction_model','section_to_person','bh_family','variance_method'}
    if design.get('approved') is not True or not required<=set(design):raise ValueError('BLOCKED: spatial statistical decisions need approval')
    # This implementation exposes the inherited additive section model. An interaction
    # proposal is not silently converted to the additive model.
    if design['interaction_model']!='none_additive_section':raise ValueError('BLOCKED: interaction model requires separate reviewed implementation')
    if design['section_to_person']!='mean_section_z_same_person_complete_support':raise ValueError('BLOCKED: section-to-person strategy not implemented/approved')
    if design['bh_family']!='all_target_state_pairs_person_empirical':raise ValueError('BLOCKED: BH family not implemented/approved')
    if design['variance_method']!='spatial_block_cluster':raise ValueError('BLOCKED: primary covariance must be spatial_block_cluster')
    if int(design['min_spots'])<1 or int(design['min_positive_spots'])<1:raise ValueError('Invalid spatial support thresholds')
    return dict(design)

def strict_fit_count_models(y,context,library,coords,min_positive_spots=10,min_spots=40):
    """Inherited Poisson/NB point model; cluster covariance failures are NOT HC1."""
    import statsmodels.api as sm
    if __package__:
        from .external_statistics import inherited_functions
    else:
        from external_statistics import inherited_functions
    ns=inherited_functions('spatial/spatial_autocorrelation_matched_gene.py',['block_groups','zscore'])
    y=np.asarray(y);groups=ns['block_groups'](coords);n_blocks=len(np.unique(groups))
    x=ns['zscore'](np.log1p(np.maximum(context,0)))
    if len(y)<min_spots or (y>0).sum()<min_positive_spots or np.std(x)<=0:raise ValueError('Insufficient target expression/context support')
    if n_blocks<2:raise ValueError('spatial_block_cluster failed: n_blocks='+str(n_blocks))
    design=sm.add_constant(x,has_constant='add');offset=np.log(np.maximum(library,1))
    alpha=min(max((float(np.var(y,ddof=1))-float(np.mean(y)))/max(float(np.mean(y))**2,1e-8),1e-4),10.)
    out={'variance_method':'spatial_block_cluster','n_blocks':n_blocks,'covariance_failure_reason':'','nb_alpha':alpha}
    for name,family in [('poisson',sm.families.Poisson()),('nb',sm.families.NegativeBinomial(alpha=alpha))]:
        try:
            fit=sm.GLM(y,design,family=family,offset=offset).fit(maxiter=100,disp=0,cov_type='cluster',cov_kwds={'groups':groups,'use_correction':True})
            if not np.isfinite(fit.bse).all() or (np.asarray(fit.bse)<=0).any():raise ValueError('nonfinite/nonpositive cluster SE')
        except Exception as exc:raise ValueError(f'{name} spatial_block_cluster failed; n_blocks={n_blocks}; {type(exc).__name__}: {exc}') from exc
        out.update({name+'_beta':float(fit.params[1]),name+'_se':float(fit.bse[1]),name+'_p':float(fit.pvalues[1])})
        if name=='poisson':out['poisson_dispersion']=float(np.sum(fit.resid_pearson**2)/max(fit.df_resid,1))
        else:out['pearson_residuals']=np.asarray(fit.resid_pearson)
    return out

def matched_person_empirical(controls,design):
    """Same complete sections/people for target and all null match ranks.
    Gene identity may vary by dataset at the same match rank; never pool different
    gene sets as if every person provided all of their union. Decision remains gated.
    """
    import pandas as pd
    from statsmodels.stats.multitest import multipletests
    spatial_design_gate(design)
    if controls.duplicated(['pair_id','dataset','section_id','match_rank']).any():raise ValueError('Duplicate section/rank observation')
    rows=[]
    for pair,g in controls.groupby('pair_id',observed=True):
        if not (g.is_target==(g.match_rank==0)).all():raise ValueError('Target/match rank mismatch')
        sections=g.pivot(index=['dataset','section_id','biological_unit_id'],columns='match_rank',values='z')
        sections=sections.replace([np.inf,-np.inf],np.nan).dropna(axis=0,how='any')
        if 0 not in sections or sections.shape[1]<2:continue
        person=sections.groupby(level='biological_unit_id').mean()
        if len(person)<2:continue
        observed=abs(person[0].mean());null=abs(person.drop(columns=0).mean()).to_numpy()
        rows.append({'pair_id':pair,'n_biological_units':len(person),'n_sections_complete':len(sections),'n_control_ranks':len(null),'empirical_p':(1+int((null>=observed).sum()))/(1+len(null)),'method':'dataset-specific match ranks; identical complete sections and persons; section mean within person'})
    out=pd.DataFrame(rows)
    if len(out):
        out['BH_q']=multipletests(out.empirical_p,method='fdr_bh')[1];out['BH_family']=design['bh_family']
    return out
