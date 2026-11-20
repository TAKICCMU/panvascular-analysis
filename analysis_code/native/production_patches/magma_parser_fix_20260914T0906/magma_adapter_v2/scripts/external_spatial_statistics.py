#!/usr/bin/env python3
"""Inherited count models, Moran test and expression/detection-matched gene controls.
Each dataset manifest row supplies mapped_h5ad AND full_raw_h5ad with exact spot IDs.
No pooling spots across donors; output section effects and biological-unit identifiers.
"""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_gate import bootstrap_cli,approved_preflight
if __name__=='__main__':bootstrap_cli()
import numpy as np,pandas as pd
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_statistics import spatial_estimators,inherited_functions
from spatial_contract import spatial_design_gate,validate_spatial_provenance,digest,file_digest,strict_fit_count_models
from spatial_renderer import effective_parameters
from state_identity import build_state_identity_map,map_state_axis

def calculate(manifest,roster,label_map,out,config,ontology,mode="target_validation"):
    design=spatial_design_gate(config["external"]["spatial"].get("statistical_design",{}))
    state_map=build_state_identity_map(ontology)
    for col in ["state_id","parent_qualified","safe_label"]:
        if col not in label_map or not label_map.set_index("state_id")[col if col!="state_id" else "parent_qualified"].sort_index().equals(state_map.set_index("state_id")[col if col!="state_id" else "parent_qualified"].sort_index()):raise ValueError("Label map does not match final ontology")
    import anndata as ad
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    ns=spatial_estimators(min_spots=design["min_spots"],min_positive_spots=design["min_positive_spots"]);ns.update(inherited_functions('spatial/spatial_autocorrelation_matched_gene.py',
       ['factor_names','abundance_frame']))
    original_to_safe=dict(zip(state_map.state_id,state_map.safe_label))
    spcfg=config['external']['spatial'];params=effective_parameters(spcfg)
    params.update(reference_mode=spcfg.get('reference_mode','descriptive'),excluded_targets=sorted(set(spcfg.get('excluded_targets',[]))),forced_genes=sorted(set(spcfg.get('forced_genes',[]))),only_gses=sorted(set(spcfg.get('only_gses',[]))))
    ontology_sha=digest(ontology.fillna('').astype(str).sort_values(['level1_annotation','level2_annotation']).to_dict('records'))
    (out/'EFFECTIVE_STATISTICAL_PARAMETERS.json').write_text(json.dumps(dict(design,mode=mode),indent=2))
    rng=np.random.default_rng(20260826);rows=[];controls=[];failed=[]
    for resource in manifest.itertuples():
        a=ad.read_h5ad(resource.mapped_h5ad);raw=ad.read_h5ad(resource.full_raw_h5ad,backed='r')
        try:
            ref=pd.read_csv(resource.reference_signatures,sep='\t',index_col=0)
            if file_digest(resource.full_raw_h5ad)!=resource.full_raw_sha256:raise ValueError('Full raw spatial input hash mismatch')
            record=json.loads(Path(resource.provenance_json).read_text())
            expected={'atlas_sha256':config['atlas']['input_sha256'],'ontology_sha256':ontology_sha,'parameters_sha256':digest(params),'reference_artifact_sha256':file_digest(resource.reference_signatures)}
            validate_spatial_provenance(record,expected,ref.index,roster.gene,a.obs_names,a.obs['sample_id'],mode,resource.mapped_h5ad,a.obs[['section_id','biological_unit_id']].astype(str).values.tolist())
            ix=raw.obs_names.get_indexer(a.obs_names)
            if (ix<0).any() or raw.obs_names.has_duplicates:raise ValueError('Spatial spot IDs do not align')
            # Native-backed row indexing must use increasing indices before restoring mapped order.
            order=np.argsort(ix);X=sparse.csr_matrix(raw.X[ix[order],:])[np.argsort(order)]
            genes=raw.var_names
            if genes.has_duplicates:raise ValueError('Audited unique full spatial gene symbols required')
            if (X.data<0).any() or not np.equal(X.data,np.floor(X.data)).all():raise ValueError('Full raw spatial counts required')
            lib=np.asarray(X.sum(axis=1)).ravel();abundance=ns['abundance_frame'](a,'q05_cell_abundance_w_sf')
            matches=ns['matched_gene_roster'](X,genes,set(roster.gene));matches.to_csv(out/f'{resource.dataset}_matched_genes.tsv',sep='\t',index=False)
            for section,idx in a.obs.groupby('section_id',observed=True).indices.items():
                donors=a.obs.iloc[idx]['biological_unit_id'].drop_duplicates()
                if len(donors)!=1 or donors.astype(str).isin(['NR','nan','Unknown','not_available']).any():raise ValueError('Unresolved/multiple biological units within spatial section')
                if len(idx)<design['min_spots']:
                    failed.append(dict(dataset=resource.dataset,section=section,error='Insufficient spots',n_spots=len(idx),required=design['min_spots']))
                    continue
                coords=np.asarray(a.obsm['spatial'])[idx]
                for _,target in roster.iterrows():
                    pair=f'{target.gene} | {target.cell_state}'
                    try:
                        if target.annotation_level=='level1_annotation':
                            safes=state_map.loc[state_map.level1_annotation==target.cell_state,'safe_label'].tolist()
                            if not safes or not set(safes)<=set(abundance):raise ValueError('Incomplete parent state coverage')
                            context=abundance[safes].sum(axis=1).to_numpy()[idx]
                        elif target.annotation_level=='level2_annotation':
                            context=abundance[original_to_safe[target.cell_state]].to_numpy()[idx]
                        else:raise ValueError('Explicit final annotation level required')
                        gi=genes.get_loc(target.gene);y=X[idx,gi].toarray().ravel()
                        identity={'dataset':resource.dataset,'section_id':section,'biological_unit_id':donors.iloc[0],'pair_id':pair,'gene':target.gene,'cell_state':target.cell_state}
                        if mode=='descriptive':
                            rows.append(dict(identity,n_spots=len(idx),mean_target=float(y.mean()),mean_context=float(context.mean()),mode='descriptive_no_target_validation_no_pvalues'))
                            continue
                        result=strict_fit_count_models(y,context,lib[idx],coords,design['min_positive_spots'],design['min_spots'])
                        mi,mp=ns['moran_permutation'](result.pop('pearson_residuals'),coords,rng)
                        identity={'dataset':resource.dataset,'section_id':section,'biological_unit_id':donors.iloc[0],'pair_id':pair,'gene':target.gene,'cell_state':target.cell_state}
                        rows.append(dict(identity,**result,n_spots=len(idx),moran_I=mi,moran_p=mp))
                        selected=[(0,target.gene)]+list(matches.loc[matches.target_gene==target.gene,['match_rank','matched_gene']].itertuples(index=False,name=None))
                        for rank,gene in selected:
                            yy=X[idx,genes.get_loc(gene)].toarray().ravel()
                            z=ns['linear_association_z'](np.log1p(yy/np.maximum(lib[idx],1)*1e4),ns['zscore'](np.log1p(np.maximum(context,0))))
                            controls.append(dict(identity,tested_gene=gene,match_rank=int(rank),is_target=gene==target.gene,z=z))
                    except (ValueError,KeyError,np.linalg.LinAlgError) as exc:failed.append(dict(dataset=resource.dataset,section=section,pair_id=pair,error=str(exc),status='NOT_ESTIMABLE',variance_method='spatial_block_cluster',n_blocks=len(np.unique(ns['block_groups'](coords)))))
        finally:raw.file.close()
    pd.DataFrame(rows).to_csv(out/'section_count_models.tsv',sep='\t',index=False)
    ctrl=pd.DataFrame(controls);ctrl.to_csv(out/'matched_gene_section_statistics.tsv',sep='\t',index=False)
    from spatial_contract import matched_person_empirical
    empirical=matched_person_empirical(ctrl,design) if len(ctrl) else pd.DataFrame()
    empirical.to_csv(out/'matched_gene_person_level_empirical.tsv',sep='\t',index=False)
    pd.DataFrame(failed).to_csv(out/'not_estimable.tsv',sep='\t',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--roster',required=True);p.add_argument('--label-map',required=True);p.add_argument('--out',required=True);p.add_argument('--config',required=True);p.add_argument('--mode',choices=['target_validation','descriptive'],required=True);a=p.parse_args()
    config,obs,ontology,info=approved_preflight(a.config,Path(a.out).with_name(Path(a.out).name+'.authorization'))
    calculate(pd.read_csv(a.manifest,sep='\t'),pd.read_csv(a.roster,sep='\t'),pd.read_csv(a.label_map,sep='\t'),a.out,config,ontology,a.mode)
