#!/usr/bin/env python3
"""Future approved external modules. Current design_and_test_only config always blocks."""
import argparse,json,os,sys,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_gate import bootstrap_cli,approved_preflight,reviewed_cohort
if __name__=='__main__':bootstrap_cli()
import numpy as np,pandas as pd,yaml
from scipy import sparse,io
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_modules import (sha256,qualified_states,specificity_from_counts,communication_eligibility,
    validate_sample_metadata,prepare_scdrs_counts,scdrs_commands,magma_command,target_evidence_join)

def run(c,module,out,config_path=None):
    if config_path is None:raise RuntimeError('Exact reviewed config path required')
    c,release_obs,ontology,preflight_info=approved_preflight(config_path,out)
    e=c['execution'];ap=Path(c['atlas']['path']);ac=c['atlas'];ext=c.get('external',{})
    out=Path(out)
    import anndata as ad
    a=ad.read_h5ad(ap,backed='r')
    scdrs_lock=None
    l1=ac.get('level1_key','level1_annotation'); l2=ac.get('level2_key','level2_annotation')
    a.obs['_state']=qualified_states(a.obs,l1,l2)
    cfg=ext.get(module,{})
    def call(cmd):
        with open(out/'commands.jsonl','a') as f:f.write(json.dumps(list(map(str,cmd)))+'\n')
        with open(out/'runtime.log','a') as f:subprocess.run(list(map(str,cmd)),check=True,stdout=f,stderr=subprocess.STDOUT)
    try:
        if module=='magma':
            gm=pd.read_csv(ext['gene_map'],sep='\t'); assert gm.feature_id.tolist()==a.var_names.tolist()
            genes=gm.ensembl_id.str.replace(r'\.\d+$','',regex=True)
            for level,key in [('L1',l1),('L2','_state')]:
                cell,study,mapping,audit=specificity_from_counts(a.X,a.obs,genes,a.obs[key],ac.get('study_key','study'))
                cell.attrs['gene_mapping_audit'].to_csv(out/f'{level}_gene_mapping_audit.tsv',sep='\t',index=False)
                mapping.to_csv(out/f'{level}_variables.tsv',sep='\t',index=False);audit.to_csv(out/f'{level}_eligibility.tsv',sep='\t',index=False)
                for name,data in [('cell_weighted',cell),('study_balanced',study)]:
                    path=out/f'{level}_{name}_specificity.tsv';data.to_csv(path,sep='\t',index=False)
                    for scenario in cfg['scenarios']:
                        if sha256(scenario['gene_results'])!=scenario['sha256']: raise ValueError('MAGMA GWAS hash mismatch')
                        prefix=out/f'{scenario["name"]}_{level}_{name}'
                        call(magma_command(cfg['executable'],scenario['gene_results'],path,prefix))
                        result=pd.read_csv(str(prefix)+'.gsa.out',sep=r'\s+',comment='#')
                        if not {'VARIABLE','P'}<=set(result):raise ValueError('Unexpected MAGMA output schema')
                        result=result.merge(mapping,left_on='VARIABLE',right_on='variable',how='inner',validate='one_to_one')
                        from external_statistics import communication_estimator_namespace
                        result['FDR_BH']=communication_estimator_namespace()['bh'](result.P.to_numpy())
                        result.to_csv(str(prefix)+'_with_BH.tsv',sep='\t',index=False)
        elif module=='scdrs':
            # Native scDRS needs one in-memory count matrix; never run concurrent full scores.
            import fcntl
            scdrs_lock=open(Path(c['outputs']['root'])/'.external_scdrs_single_full.lock','a')
            fcntl.flock(scdrs_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            mem={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if ':' in line}
            if mem['MemAvailable'] < max(220,int(cfg.get('minimum_available_gib',220)))*1024**3:
                raise RuntimeError('Insufficient reviewed native scDRS memory budget; no matrix materialized')
            gm=pd.read_csv(ext['gene_map'],sep='\t');assert gm.feature_id.tolist()==a.var_names.tolist()
            # One raw-count symbol matrix, shared across serial or resource-bounded scenario scoring.
            X,names=prepare_scdrs_counts(a.X[:],a.obs,gm)
            obs=a.obs.copy();obs[l2]=obs['_state'].fillna('Excluded_Not_assessed')
            b=ad.AnnData(X=X,obs=obs,var=pd.DataFrame(index=names))
            group_keys=[]
            for label_key in [l1,l2]:
                support=b.obs.groupby([label_key,ac.get('study_key','study')],observed=True).size().unstack(fill_value=0)
                eligible=support.index[(support.sum(axis=1)>=500)&((support>=20).sum(axis=1)>=3)]
                eligible=[g for g in eligible if g not in ['Unknown','Not_assessed','Excluded_Not_assessed']]
                test_key=label_key+'_scdrs_test';group_keys.append(test_key)
                b.obs[test_key]=b.obs[label_key].astype(str).where(b.obs[label_key].isin(eligible),'Excluded_Not_eligible')
                support.to_csv(out/f'{label_key}_study_support.tsv',sep='\t')
            # Preserve the exact frozen global graph for native scDRS heterogeneity.
            # No reclustering or neighbor recomputation; b retains exact original cell order.
            ng=a.uns['neighbors_cpu_n30']
            ck=ng['connectivities_key']; dk=ng['distances_key']
            if not b.obs_names.equals(a.obs_names): raise ValueError('Frozen graph cell order differs')
            if a.obsp[ck].shape != (b.n_obs,b.n_obs): raise ValueError('Frozen graph shape differs')
            b.obsp['connectivities']=a.obsp[ck].copy()
            b.obsp['distances']=a.obsp[dk].copy()
            b.uns['neighbors']={'connectivities_key':'connectivities','distances_key':'distances','params':dict(ng['params'])}
            hp=out/'scdrs_raw_symbols.h5ad';b.write_h5ad(hp,compression='lzf')
            # Covariates explicitly inherited from old scoring; file rows must match exact cell IDs.
            cov=pd.read_csv(cfg['covariates'],sep='\t',index_col=0)
            if not cov.index.equals(b.obs_names): raise ValueError('scDRS covariate cell IDs/order mismatch')
            # Release parent count copy before native score subprocess allocates its input/controls.
            del b,X
            __import__('gc').collect()
            for scenario in cfg['scenarios']:
                if sha256(scenario['geneset'])!=scenario['sha256']: raise ValueError('scDRS geneset hash mismatch')
                so=out/scenario['name'];so.mkdir()
                for cmd in scdrs_commands(cfg['executable'],hp,cfg['covariates'],scenario['geneset'],so):call(cmd)
                score=list((so/'scores').glob('*.full_score.gz'))
                if len(score)!=1: raise ValueError('Expected one native score file')
                for key in group_keys:
                    (so/key).mkdir(parents=True,exist_ok=False)
                    call([cfg['executable'],'perform-downstream','--h5ad-file',hp,'--score-file',score[0],
                          '--out-folder',so/key,'--group-analysis',key, '--flag-filter-data','True','--flag-raw-count','True'])
                    files=list((so/key).glob('*.scdrs_group.*'))
                    if len(files)!=1:raise ValueError('Native group result schema missing')
                    tab=pd.read_csv(files[0],sep='\t')
                    if 'group' not in tab:raise ValueError('Native group column absent')
                    tab=tab[~tab.group.astype(str).str.startswith('Excluded_')].copy()
                    from external_statistics import communication_estimator_namespace
                    bh=communication_estimator_namespace()['bh']
                    for metric in ['assoc_mcp','hetero_mcp']:
                        if metric in tab:tab[metric+'_BH']=bh(tab[metric].to_numpy())
                    tab.to_csv(so/key/'eligible_native_group_BH.tsv',sep='\t',index=False)
        elif module=='communication':
            bio=ac.get('biological_sample_key','biological_sample_uid'); donor=ac.get('donor_key','donor_uid')
            units=reviewed_cohort(c,release_obs)
            units.to_csv(out/'reviewed_cohort.tsv',sep='\t',index=False)
            safe=units.aggregate_approved & units.biological_sample_resolved
            units.loc[~safe].to_csv(out/'blocked_communication_units.tsv',sep='\t',index=False)
            allowed_samples=set(units.loc[safe,'sample_uid'])
            allowed_obs=a.obs[a.obs[ac.get('sample_key','sample_uid')].isin(allowed_samples)].copy()
            if allowed_obs.empty:raise RuntimeError('No aggregation-approved biological units')
            # Exact biological specimen membership is reviewed before this stage. Unknown donor stays missing.
            validate_sample_metadata(allowed_obs,bio,[donor,ac.get('study_key','study'),ac.get('comparison_key','comparison_group')])
            elig=communication_eligibility(allowed_obs,bio,l1,cfg.get('min_cells',50));elig.to_csv(out/'pair_eligibility.tsv',sep='\t',index=False)
            gm=pd.read_csv(ext['gene_map'],sep='\t');assert gm.feature_id.tolist()==a.var_names.tolist()
            manifest=[]
            for number,(sample,local_idx) in enumerate(allowed_obs.groupby(bio,observed=True).indices.items()):
                idx=a.obs_names.get_indexer(allowed_obs.index[local_idx])
                if (idx<0).any():raise ValueError('Approved specimen cell mapping mismatch')
                counts=a.obs.iloc[idx][l1].value_counts(); types=counts[counts>=cfg.get('min_cells',50)].index
                types=[t for t in types if t not in ['Unknown','Not_assessed']]
                if len(types)<2:continue
                idx=idx[a.obs.iloc[idx][l1].isin(types).to_numpy()]
                # Inherited bounded per-specimen/type sampling, no donor pooling.
                seed=int.from_bytes(__import__('hashlib').sha256(f"20260820|{sample}".encode()).digest()[:4],'little') % 2147483647
                rng=np.random.default_rng(seed); chosen=[]
                for t in types:
                    ti=idx[a.obs.iloc[idx][l1].astype(str).to_numpy()==t]
                    chosen.extend(rng.choice(ti,min(len(ti),cfg.get('max_cells_per_sample_type',200)),replace=False).tolist())
                idx=np.sort(chosen)
                # Collapse duplicate symbols at raw-count level, then CP10k/log1p per cell.
                X,names=prepare_scdrs_counts(a.X[idx,:],a.obs.iloc[idx],gm)
                totals=np.asarray(X.sum(axis=1)).ravel(); assert (totals>0).all()
                X=sparse.diags(1e4/totals)@X;X.data=np.log1p(X.data)
                b=ad.AnnData(X=X,obs=a.obs.iloc[idx].copy(),var=pd.DataFrame(index=names))
                sd=out/f'specimen_{number:04d}';sd.mkdir()
                b.obs['_communication_group']=b.obs[l1].astype('category');b.obs['samples']=str(sample)
                b.layers['log1p']=b.X
                if cfg.get('run_liana',True):
                    import liana as li
                    if sha256(cfg['liana_resource'])!=cfg['liana_resource_sha256']:raise ValueError('LIANA resource hash mismatch')
                    resource=pd.read_csv(cfg['liana_resource'],sep='\t')
                    li.mt.rank_aggregate.by_sample(b,sample_key=bio,groupby='_communication_group',resource=resource,
                        expr_prop=0.10,min_cells=cfg.get('min_cells',50),layer='log1p',use_raw=False,inplace=True,seed=seed)
                    b.uns['liana_res'].to_csv(sd/'liana_interactions.tsv',sep='\t',index=False)
                io.mmwrite(sd/'expression.mtx',b.X.T)
                pd.Series(b.var_names).to_csv(sd/'genes.tsv',sep='\t',header=False,index=False)
                pd.Series(b.obs_names).to_csv(sd/'cells.tsv',sep='\t',header=False,index=False)
                b.obs.to_csv(sd/'obs.tsv',sep='\t')
                manifest.append({'sample':sample,'directory':str(sd),'min_cells':cfg.get('min_cells',50),'seed':seed})
                if cfg.get('run_cellchat',True):
                    if sha256(cfg['cellchat_db'])!=cfg['cellchat_db_sha256']: raise ValueError('CellChat database hash mismatch')
                    call([cfg['rscript'],Path(__file__).with_name('external_cellchat.R'),sd,cfg['cellchat_db'],cfg.get('min_cells',50)])
            pd.DataFrame(manifest).to_csv(out/'communication_sample_runs.tsv',sep='\t',index=False)
        elif module=='targets':
            tables={name:pd.read_csv(path,sep='\t') for name,path in cfg['evidence_tables'].items()}
            target_evidence_join(tables).to_csv(out/'target_evidence_all_genes.tsv',sep='\t',index=False)
        else: raise ValueError('Unsupported module; genetics/coloc/spatial use audited standalone adapters')
        (out/'COMPLETE.json').write_text(json.dumps({'module':module,'atlas_sha256':ac['input_sha256'],'labels_written':False}))
    finally:
        a.file.close()
        if scdrs_lock is not None:scdrs_lock.close()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--module',choices=['magma','scdrs','communication','targets'],required=True);p.add_argument('--out',required=True);v=p.parse_args()
    run(yaml.safe_load(Path(v.config).read_text()),v.module,v.out,v.config)
