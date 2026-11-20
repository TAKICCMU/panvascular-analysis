"""Bind actual CLI tables to approved CoVarNet and spatial artifacts; descriptive only."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
if __package__:
    from .spatial_contract import digest,file_digest,validate_spatial_provenance
    from .state_identity import build_state_identity_map,map_state_axis
    from .spatial_renderer import effective_parameters
    from .external_statistics import inherited_functions
else:
    from spatial_contract import digest,file_digest,validate_spatial_provenance
    from state_identity import build_state_identity_map,map_state_axis
    from spatial_renderer import effective_parameters
    from external_statistics import inherited_functions

WEIGHT_SEMANTICS='column_normalized_nsNMF_W_state_by_module_weights_not_cell_abundance'

def ontology_digest(ont):
    return digest(ont.fillna('').astype(str).sort_values(['level1_annotation','level2_annotation']).to_dict('records'))

def read_record(path):
    if not path or not Path(path).is_file():raise ValueError('Missing source manifest')
    return json.loads(Path(path).read_text())

def same_path(a,b,label):
    if not a or not b or Path(a).resolve()!=Path(b).resolve():raise ValueError('Source path mismatch: '+label)

def verify_hash(path,expected,label):
    if not expected or not Path(path).is_file() or file_digest(path)!=expected:raise ValueError('Content digest mismatch: '+label)

def table(path):
    d=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False);d=d.set_index(d.columns[0]);d.index=d.index.astype(str);d.columns=d.columns.astype(str);d=d.astype(float)
    if d.index.has_duplicates or d.columns.has_duplicates or not np.isfinite(d.to_numpy(float)).all() or (d.to_numpy(float)<0).any():raise ValueError('Ambiguous or invalid artifact table')
    return d

def weights_record(weights,run_dir,ontology,atlas_sha):
    """Prepare a review sidecar from already complete artifacts; grants no approval."""
    run=Path(run_dir).resolve();cp=run/'covarnet/COVARNET_COMPLETE.json';identity=run/'RUN_IDENTITY.json';preflight=run/'preflight.json';status=run/'RUN_STATUS.json';inp=run/'covarnet/INPUT_MANIFEST.json';co=run/'counts/ontology.tsv'
    completion=read_record(cp);ri=read_record(identity);pf=read_record(preflight);rs=read_record(status)
    if completion.get('status')!='COMPLETE' or rs.get('covarnet',{}).get('status')!='COMPLETE':raise ValueError('CoVarNet run not complete')
    if ri.get('atlas_sha256')!=atlas_sha or pf.get('sha256')!=atlas_sha or pf.get('status')!='PASS_RELEASE_PREFLIGHT':raise ValueError('CoVarNet run atlas mismatch')
    if digest(build_state_identity_map(pd.read_csv(co,sep='\t')).sort_values('state_id').to_dict('records'))!=digest(build_state_identity_map(ontology).sort_values('state_id').to_dict('records')):raise ValueError('CoVarNet run ontology mismatch')
    verify_hash(inp,completion.get('input_manifest_sha256'),'CoVarNet inputs')
    inputs=read_record(inp);verify_hash(co,inputs.get('source_identity',{}).get('files',{}).get(str(co)),'CoVarNet input ontology')
    selected=run/'covarnet/stable_NMF/selected_K.txt'
    verify_hash(selected,completion.get('files',{}).get('stable_NMF/selected_K.txt'),'selected primary K')
    k_lines=[v.split('=',1)[1] for v in selected.read_text().splitlines() if v.startswith('K_FINAL=')]
    if len(k_lines)!=1 or not k_lines[0].isdigit() or int(k_lines[0])<1:raise ValueError('Invalid selected primary K')
    k=int(k_lines[0]);expected_w=run/'covarnet/stable_NMF'/f'module_node_weights_K{k:02d}.tsv'
    same_path(weights,expected_w,'selected primary W artifact, not counts/proportions')
    nodes=run/'covarnet/CoVarNet_node_counts.tsv'
    verify_hash(nodes,inputs.get('files',{}).get(nodes.name),'input state roster')
    weights=Path(weights).resolve()
    try:relative=str(weights.relative_to(run/'covarnet'))
    except ValueError:raise ValueError('Weights not inside selected CoVarNet run')
    verify_hash(weights,completion.get('files',{}).get(relative),'complete run W table')
    w=table(weights);state_map=build_state_identity_map(ontology);map_state_axis(w.index,state_map,'state_id')
    if list(w.index)!=list(table(nodes).index):raise ValueError('W state order differs from fitted input nodes')
    if list(w.columns)!=[f'CM{i:02d}' for i in range(1,k+1)]:raise ValueError('W module IDs do not match selected K')
    if not np.allclose(w.sum(axis=0).to_numpy(float),1.,rtol=1e-8,atol=1e-10):raise ValueError('Expected column-normalized native W weights')
    return {'schema':'spatial_module_weights_v2_1','atlas_sha256':atlas_sha,'ontology_sha256':ontology_digest(ontology),'run_dir':str(run),'run_identity_sha256':file_digest(identity),'completion_sha256':file_digest(cp),'preflight_sha256':file_digest(preflight),'run_status_sha256':file_digest(status),'input_manifest_sha256':file_digest(inp),'weights_path':str(weights),'weights_sha256':file_digest(weights),'state_ids':list(w.index),'module_ids':list(w.columns),'weight_semantics':WEIGHT_SEMANTICS,'selected_K':k,'selected_K_sha256':file_digest(selected)}

def abundance_record(abundance,posterior,provenance,reference,ontology,config,abundance_key,source_obsm_key='q05_cell_abundance_w_sf'):
    """Check every exported value/spot against the actual posterior, never just names."""
    import anndata as ad
    if source_obsm_key not in {'q05_cell_abundance_w_sf','means_cell_abundance_w_sf'}:raise ValueError('Unsupported posterior abundance key; embeddings are not abundance')
    ab=table(abundance);m=build_state_identity_map(ontology)
    canonical=map_state_axis(ab.columns,m,abundance_key);safes=map_state_axis(canonical,m,'state_id','safe_label')
    a=ad.read_h5ad(posterior,backed='r')
    try:
        ref=pd.read_csv(reference,sep='\t',index_col=0);record=read_record(provenance)
        sc=config['external']['spatial'];params=effective_parameters(sc)
        params.update(reference_mode=sc.get('reference_mode','descriptive'),excluded_targets=sorted(set(sc.get('excluded_targets',[]))),forced_genes=sorted(set(sc.get('forced_genes',[]))),only_gses=sorted(set(sc.get('only_gses',[]))))
        expected={'atlas_sha256':config['atlas']['input_sha256'],'ontology_sha256':ontology_digest(ontology),'parameters_sha256':digest(params),'reference_artifact_sha256':file_digest(reference)}
        validate_spatial_provenance(record,expected,ref.index,params['excluded_targets'],a.obs_names,a.obs['sample_id'],'descriptive',posterior,a.obs[['section_id','biological_unit_id']].astype(str).values.tolist())
        if list(ab.index)!=list(map(str,a.obs_names)):raise ValueError('Abundance/posterior spot identity or order mismatch')
        ns=inherited_functions('spatial/spatial_autocorrelation_matched_gene.py',['factor_names','abundance_frame'])
        value=a.obsm[source_obsm_key]
        if isinstance(value,pd.DataFrame):
            prefix=source_obsm_key.replace('_cell','cell',1)+'_'
            columns=[str(v)[len(prefix):] if str(v).startswith(prefix) else str(v) for v in value.columns]
            names=a.uns.get('mod',{}).get('factor_names')
            if names is not None and list(map(str,names))!=columns:raise ValueError('Posterior column labels conflict with factor_names order')
        source=ns['abundance_frame'](a,source_obsm_key)
        if not set(safes)<=set(source.columns):raise ValueError('Abundance states absent from posterior')
        if not np.allclose(ab.to_numpy(float),source.loc[:,safes].to_numpy(float),rtol=1e-10,atol=1e-12):raise ValueError('Abundance values differ from actual posterior export')
        return {'schema':'spatial_module_abundance_v2_1','atlas_sha256':config['atlas']['input_sha256'],'ontology_sha256':ontology_digest(ontology),'abundance_path':str(Path(abundance).resolve()),'abundance_sha256':file_digest(abundance),'posterior_path':str(Path(posterior).resolve()),'posterior_sha256':file_digest(posterior),'provenance_path':str(Path(provenance).resolve()),'provenance_sha256':file_digest(provenance),'reference_path':str(Path(reference).resolve()),'reference_sha256':file_digest(reference),'source_obsm_key':source_obsm_key,'abundance_key':abundance_key,'state_ids':list(canonical),'spot_ids_sha256':digest(list(ab.index)),'sample_ids_sha256':digest(list(map(str,a.obs.sample_id))),'semantics':'posterior_abundance_export_for_descriptive_projection'}
    finally:a.file.close()

def validate_projection_inputs(config,ontology,abundance,weights,abundance_key):
    cfg=config.get('external',{}).get('spatial',{}).get('module_projection',{})
    if cfg.get('approved') is not True:raise ValueError('Module projection artifact selection not approved')
    wm=read_record(cfg.get('weights_manifest'));am=read_record(cfg.get('abundance_manifest'))
    verify_hash(cfg['weights_manifest'],cfg.get('weights_manifest_sha256'),'approved W manifest')
    verify_hash(cfg['abundance_manifest'],cfg.get('abundance_manifest_sha256'),'approved abundance manifest')
    same_path(weights,wm.get('weights_path'),'actual CLI weights');same_path(abundance,am.get('abundance_path'),'actual CLI abundance')
    same_path(wm.get('run_dir'),cfg.get('covarnet_run_dir'),'approved CoVarNet run')
    if wm.get('run_identity_sha256')!=cfg.get('covarnet_run_identity_sha256'):raise ValueError('Approved CoVarNet run identity differs')
    same_path(am.get('posterior_path'),cfg.get('posterior_path'),'approved posterior')
    if am.get('posterior_sha256')!=cfg.get('posterior_sha256') or am.get('abundance_key')!=abundance_key:raise ValueError('Approved posterior or abundance namespace differs')
    actual_w=weights_record(weights,cfg['covarnet_run_dir'],ontology,config['atlas']['input_sha256'])
    if actual_w!=wm:raise ValueError('Actual W content/source does not match reviewed manifest')
    actual_ab=abundance_record(abundance,am['posterior_path'],am['provenance_path'],am['reference_path'],ontology,config,abundance_key,am['source_obsm_key'])
    if actual_ab!=am:raise ValueError('Actual abundance content/source does not match reviewed manifest')
    return {'status':'COMPLETE','mode':'DESCRIPTIVE_MODULE_PROJECTION','independent_spatial_validation':False,'weights':actual_w,'abundance':actual_ab}
