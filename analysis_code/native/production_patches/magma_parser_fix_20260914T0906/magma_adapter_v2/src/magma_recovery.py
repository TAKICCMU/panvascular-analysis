"""Strict MAGMA conditional-model parsing and content-certified limited restart."""
from pathlib import Path
import json,shutil,time
import numpy as np
import pandas as pd
from external_modules import sha256,specificity_from_counts,magma_command
from external_statistics import communication_estimator_namespace

def parse_conditional_results(result,mapping,family):
    required={'VARIABLE','TYPE','MODEL','NGENES','BETA','BETA_STD','SE','P'}
    if not required<=set(result):raise ValueError('MAGMA output schema incomplete')
    if not {'variable','group'}<=set(mapping) or mapping.variable.isna().any() or mapping.variable.duplicated().any() or mapping.group.duplicated().any():raise ValueError('MAGMA target mapping must be one-to-one')
    if 'Average' in set(mapping.variable):raise ValueError('Average cannot be a target')
    r=result.copy();targets=set(mapping.variable)
    if r.VARIABLE.isna().any() or not set(r.VARIABLE)<=targets|{'Average'}:raise ValueError('Unexpected MAGMA output variable')
    if not r.TYPE.eq('COVAR').all():raise ValueError('Expected only covariate results')
    for col in ['MODEL','NGENES','BETA','BETA_STD','SE','P']:
        r[col]=pd.to_numeric(r[col],errors='raise')
        if not np.isfinite(r[col].to_numpy()).all():raise ValueError('Nonfinite MAGMA output')
    if not ((r.MODEL>=1)&(r.MODEL==np.floor(r.MODEL))).all():raise ValueError('Invalid model identifier')
    r['MODEL']=r.MODEL.astype(int)
    if not ((r.P>=0)&(r.P<=1)&(r.SE>0)&(r.NGENES>0)).all():raise ValueError('Invalid MAGMA numeric result')
    if r.duplicated(['VARIABLE','MODEL']).any():raise ValueError('Duplicate variable-model rows')
    t=r.loc[r.VARIABLE.isin(targets)].copy();n=r.loc[r.VARIABLE.eq('Average')].copy()
    if set(t.VARIABLE)!=targets:raise ValueError('Missing mapped target variable')
    if t.VARIABLE.duplicated().any():raise ValueError('Target appears in repeated models')
    if t.MODEL.duplicated().any():raise ValueError('Multiple target variables in one model')
    if n.MODEL.duplicated().any() or set(n.MODEL)!=set(t.MODEL):raise ValueError('Every target model requires exactly one Average condition')
    # Now, and only now, target-level one-to-one joining is well-defined.
    t=t.merge(mapping,left_on='VARIABLE',right_on='variable',how='left',validate='one_to_one')
    t['FDR_BH']=communication_estimator_namespace()['bh'](t.P.to_numpy())
    t['BH_family']=family;t['n_tests_in_family']=len(t)
    model=t[['MODEL','VARIABLE','group']].merge(n[['MODEL','VARIABLE']].rename(columns={'VARIABLE':'condition_variable'}),on='MODEL',validate='one_to_one')
    n['role']='nuisance_condition_not_an_enrichment_target'
    return t,n,model

def _read_and_validate(path,mapping,family):
    headers=[s for s in Path(path).read_text().splitlines() if s.startswith('#')]
    if not any(s.startswith('# CONDITIONED_VARIABLES') and 'Average (covar)' in s for s in headers):raise ValueError('Native output does not declare Average conditioning')
    return parse_conditional_results(pd.read_csv(path,sep=r'\s+',comment='#'),mapping,family)

def run_recovery(a,ac,cfg,ext,out,call):
    out=Path(out);mp=Path(cfg['recovery_manifest'])
    if sha256(mp)!=cfg['recovery_manifest_sha256']:raise ValueError('Recovery manifest drift')
    reuse=json.loads(mp.read_text());assert reuse['atlas_sha256']==ac['input_sha256']
    for name,item in reuse['files'].items():
        if sha256(item['path'])!=item['sha256']:raise ValueError('Prior completed artifact drift: '+name)
    expected={'level1_key':ac.get('level1_key','level1_annotation'),'level2_key':ac.get('level2_key','level2_annotation'),'study_key':ac.get('study_key','study')}
    if reuse['annotation_and_study_keys']!=expected:raise ValueError('Prior annotation/study semantics differ')
    if sha256(Path(__file__).with_name('external_modules.py'))!=reuse['specificity_source_sha256']:raise ValueError('Specificity computation differs from certified run')
    gm=pd.read_csv(ext['gene_map'],sep='\t');assert gm.feature_id.tolist()==a.var_names.tolist();genes=gm.ensembl_id.str.replace(r'\.\d+$','',regex=True)
    def state(**x):(out/'RECOVERY_PROGRESS.json').write_text(json.dumps({'updated_unix':time.time(),**x},indent=2))
    state(stage='CERTIFIED_REUSE',source=mp.as_posix())
    for key in ['L1_cell_weighted_specificity.tsv','L1_variables.tsv','L1_eligibility.tsv','L1_gene_mapping_audit.tsv','primary_r0.1_L1_cell_weighted.gsa.out','primary_r0.1_L1_cell_weighted.log']:
        item=reuse['files'][key];target=out/key
        if target.exists():raise FileExistsError(target)
        shutil.copyfile(item['path'],target)
    reused=[]
    for level,key in [('L1',expected['level1_key']),('L2','_state')]:
        state(stage='RAW_SPECIFICITY',level=level,note='L1 raw pass is required for previously unsaved study-balanced output; copied cell-weighted table is verified, not republished from the recomputation')
        cell,study,mapping,audit=specificity_from_counts(a.X,a.obs,genes,a.obs[key],expected['study_key'])
        if level=='L1':
            prior=pd.read_csv(out/'L1_cell_weighted_specificity.tsv',sep='\t')
            copy=cell.copy();copy.attrs={};pd.testing.assert_frame_equal(prior,copy,check_exact=False,rtol=1e-12,atol=1e-12)
            pd.testing.assert_frame_equal(pd.read_csv(out/'L1_variables.tsv',sep='\t'),mapping)
            pd.testing.assert_frame_equal(pd.read_csv(out/'L1_eligibility.tsv',sep='\t'),audit)
        else:
            cell.attrs['gene_mapping_audit'].to_csv(out/f'{level}_gene_mapping_audit.tsv',sep='\t',index=False)
            mapping.to_csv(out/f'{level}_variables.tsv',sep='\t',index=False);audit.to_csv(out/f'{level}_eligibility.tsv',sep='\t',index=False)
        for variant,data in [('cell_weighted',cell),('study_balanced',study)]:
            path=out/f'{level}_{variant}_specificity.tsv'
            if not (level=='L1' and variant=='cell_weighted'):data.to_csv(path,sep='\t',index=False)
            for scenario in cfg['scenarios']:
                if sha256(scenario['gene_results'])!=scenario['sha256']:raise ValueError('MAGMA GWAS input drift')
                name=f'{scenario["name"]}_{level}_{variant}';prefix=out/name;gsa=Path(str(prefix)+'.gsa.out')
                certified=(name=='primary_r0.1_L1_cell_weighted')
                state(stage='PARSE_REUSED_NATIVE_OUTPUT' if certified else 'NATIVE_MAGMA',analysis=name)
                if certified:reused.append({'analysis':name,'raw_output_sha256':sha256(gsa),'source':reuse['files'][gsa.name]['path']})
                else:call(magma_command(cfg['executable'],scenario['gene_results'],path,prefix))
                target,nuisance,models=_read_and_validate(gsa,mapping,name)
                target.to_csv(str(prefix)+'_with_BH.tsv',sep='\t',index=False);nuisance.to_csv(str(prefix)+'_condition_rows.tsv',sep='\t',index=False);models.to_csv(str(prefix)+'_model_mapping.tsv',sep='\t',index=False)
                state(stage='ANALYSIS_COMPLETE',analysis=name,n_targets=len(target),reused_native_output=certified)
    (out/'CERTIFIED_REUSE_APPLIED.json').write_text(json.dumps(reused,indent=2));state(stage='ALL_12_ANALYSES_COMPLETE',native_calls=11,reused_native_outputs=1)
