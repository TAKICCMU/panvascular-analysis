"""Auditable AST adapter retaining inherited cell2location training algorithms."""
import ast,json
from pathlib import Path
if __package__:
    from .spatial_contract import digest
else:
    from spatial_contract import digest
PARAMETER_DEFAULTS={'reference_max_epochs':250,'mapping_max_epochs':30000,'posterior_samples':1000,'N_cells_per_location':12,'detection_alpha':20.,'max_cells_per_state':1500,'min_cells_per_state':40,'max_studies_per_state':20,'seed':20260814}
GLOBALS={'reference_max_epochs':'REF_MAX_EPOCHS','mapping_max_epochs':'SPATIAL_MAX_EPOCHS','posterior_samples':'POSTERIOR_SAMPLES','N_cells_per_location':'N_CELLS_PER_LOCATION','detection_alpha':'DETECTION_ALPHA','max_cells_per_state':'MAX_CELLS_PER_STATE','min_cells_per_state':'MIN_CELLS_PER_STATE','max_studies_per_state':'MAX_STUDIES_PER_STATE','seed':'RANDOM_SEED'}

def effective_parameters(training_config=None):
    c=dict(training_config or {});p=dict(PARAMETER_DEFAULTS)
    for key in p:
        if key in c:p[key]=c[key]
    for key in p:
        if isinstance(p[key],bool) or not isinstance(p[key],(int,float)) or p[key]<0 or (key!='seed' and p[key]==0):raise ValueError('Invalid spatial parameter: '+key)
        if key!='detection_alpha' and int(p[key])!=p[key]:raise ValueError('Integer spatial parameter required: '+key)
    p.update(CPU_only=True,reference_posterior_batch_size=2500,mapping_batch_size=None,mapping_train_size=1,reference_filter_cell_count_cutoff=5,reference_filter_cell_percentage_cutoff2=.03,reference_filter_nonz_mean_cutoff=1.12)
    return p

def render(source,atlas_path,output_root,manifest_path,l1_key,l2_key,study_key,approved_forced_genes,atlas_sha256,config_path=None,training_config=None):
    p=effective_parameters(training_config);cfg=dict(training_config or {})
    mode=cfg.get('reference_mode','descriptive');targets=sorted(set(cfg.get('excluded_targets',[])))
    if mode not in {'descriptive','target_excluded'}:raise ValueError('Explicit reference_mode required')
    if mode=='target_excluded' and not targets:raise ValueError('Target-excluded reference needs explicit target roster')
    if set(targets)&set(approved_forced_genes):raise ValueError('Targets cannot be forced into reference')
    p.update(reference_mode=mode,excluded_targets=targets,forced_genes=sorted(set(approved_forced_genes)),only_gses=sorted(set(cfg.get('only_gses',[]))))
    tree=ast.parse(source);lines=source.splitlines(keepends=True)
    replacements={'ROOT':f'ROOT = Path({str(output_root)!r})\n','SOURCE_ROOT':f'SOURCE_ROOT = Path({str(Path(manifest_path).parent)!r})\n','ATLAS_PATH':f'ATLAS_PATH = Path({str(atlas_path)!r})\n','MANIFEST_PATH':f'MANIFEST_PATH = Path({str(manifest_path)!r})\n','USE_GPU':'USE_GPU = False\n','FORCED_GENES':f'FORCED_GENES = set({p["forced_genes"]!r})\n'}
    replacements.update({glob:f'{glob} = {p[key]!r}\n' for key,glob in GLOBALS.items()})
    funcs={
      'select_state_key':f'def select_state_key(obs):\n    return {l2_key!r}\n',
      'select_level1_key':f'def select_level1_key(obs):\n    return {l1_key!r}\n',
      'select_study_key':f'def select_study_key(obs):\n    return {study_key!r}\n',
      'gpu_audit':'''def gpu_audit():
    from importlib.metadata import version,PackageNotFoundError
    versions={}
    for package in ['cell2location','scvi-tools','scanpy','omicverse','torch','anndata']:
        try:versions[package]=version(package)
        except PackageNotFoundError:versions[package]='not_installed'
    STATE_IDENTITY_MAP.to_csv(REF_OUT/'state_identity_map.tsv',sep='\\t',index=False)
    audit=dict(EFFECTIVE_PARAMETERS, package_versions=versions, parameter_sha256=PARAMETERS_SHA256, atlas_sha256=EXPECTED_ATLAS_SHA256, backend='CPU', python=sys.version, torch=torch.__version__)
    (AUDIT_OUT/'EFFECTIVE_MODEL_PARAMETERS.json').write_text(json.dumps(audit,indent=2))
    return audit
''',
      'build_spatial_reference_labels':f'''def build_spatial_reference_labels(obs):
    from state_identity import map_state_axis
    state=obs[{l2_key!r}].astype(str);parent=obs[{l1_key!r}].astype(str)
    pair=parent+' | '+state
    known=~state.isin({{'Unknown','Not_assessed','nan','NR','not_available',''}})
    if not pair[known].isin(STATE_IDENTITY_MAP.parent_qualified).all():raise ValueError('Atlas/ontology parent mismatch')
    # Unknown cells remain explicit nuisance reference groups, never CoVarNet nodes.
    labels=pair.where(known,'Other | '+parent)
    return labels, {{'state_key':{l2_key!r},'level1_key':{l1_key!r}}}
''',
      'clean_label':'''def clean_label(value):
    matches=STATE_IDENTITY_MAP.loc[STATE_IDENTITY_MAP.parent_qualified==str(value),'safe_label']
    if len(matches)==1:return matches.iloc[0]
    if str(value).startswith('Other | '):return 'NUISANCE_'+_digest(str(value))[:24]
    raise ValueError('Unmapped spatial label: '+str(value))
'''}
    edits=[];found=set()
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in replacements:
                    edits.append((node.lineno,node.end_lineno,replacements[target.id]));found.add(target.id)
        if isinstance(node,ast.FunctionDef) and node.name in funcs:
            edits.append((node.lineno,node.end_lineno,funcs[node.name]));found.add(node.name)
    if found!=set(replacements)|set(funcs):raise ValueError('Inherited spatial source changed: patch aborted')
    for start,end,new in sorted(edits,reverse=True):lines[start-1:end]=[new]
    result=''.join(lines)
    for old,new in [('atlas.raw.var_names','atlas.var_names'),('atlas.raw.var','atlas.var'),('atlas.raw.X','atlas.X')]:
        if old not in result:raise ValueError('Inherited counts source not found')
        result=result.replace(old,new)
    result=result.replace('for value in os.environ.get("SPATIAL_ONLY_GSES", "").split(",")','for value in EFFECTIVE_PARAMETERS["only_gses"]')
    result=result.replace('    ref = ref[:, sorted(selected_genes)].copy()', '    selected_genes.difference_update(EFFECTIVE_PARAMETERS["excluded_targets"])\n    ref = ref[:, sorted(selected_genes)].copy()')
    result=result.replace('"cell2location.filter_genes plus forced target/pathology genes; no fixed maximum-gene cap"','"cell2location.filter_genes plus reviewed forced genes minus explicit excluded targets; no fixed gene cap"')
    # Never reuse outputs based only on size. Resume needs a new authenticated artifact workflow.
    result=result.replace('        print("Reusing reference signatures:", signature_path)\n        return pd.read_csv(signature_path, sep="\\t", index_col=0)','        raise RuntimeError("Existing reference forbidden without reviewed content-provenance resume")')
    result=result.replace('        print("Reusing mapped output:", output)\n        return output','        raise RuntimeError("Existing posterior forbidden without reviewed content-provenance resume")')
    result=result.replace('    spatial.write_h5ad(output, compression="gzip")','''    spatial.write_h5ad(output, compression="gzip")
    provenance={"schema":"spatial_provenance_v2","atlas_sha256":EXPECTED_ATLAS_SHA256,"ontology_sha256":ONTOLOGY_SHA256,"parameters_sha256":PARAMETERS_SHA256,"reference_artifact_sha256":_file_digest(REF_OUT/'final_atlas_reference_signatures.tsv.gz'),"posterior_sha256":_file_digest(output),"reference_genes_sha256":_digest(list(map(str,signatures.index))),"spot_ids_sha256":_digest(list(map(str,spatial.obs_names))),"sample_ids_sha256":_digest(list(map(str,spatial.obs['sample_id']))),"spatial_units_sha256":_digest(spatial.obs[['section_id','biological_unit_id']].astype(str).values.tolist()),"targets_sha256":_digest(EFFECTIVE_PARAMETERS['excluded_targets']),"reference_mode":EFFECTIVE_PARAMETERS['reference_mode'],"excluded_targets":EFFECTIVE_PARAMETERS['excluded_targets']}
    output.with_suffix('.provenance.json').write_text(json.dumps(provenance,indent=2))''')
    gate=Path(__file__).resolve().parent
    guard=f'''import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0,{str(gate)!r})
from external_gate import approved_preflight as _approved_preflight
from spatial_contract import digest as _digest, file_digest as _file_digest
from state_identity import build_state_identity_map as _build_state_identity_map
EFFECTIVE_PARAMETERS={p!r}
PARAMETERS_SHA256={digest(p)!r}
EXPECTED_ATLAS_SHA256={atlas_sha256!r}
_config_path={str(config_path) if config_path is not None else None!r}
if _config_path is None:raise RuntimeError('WAITING_FOR_FINAL_L2_ATLAS: reviewed config required')
_c,_obs,_ontology,_info=_approved_preflight(_config_path,_Path({str(output_root)!r}).with_name(_Path({str(output_root)!r}).name+'.authorization'))
if _c['atlas']['input_sha256']!=EXPECTED_ATLAS_SHA256:raise RuntimeError('Atlas digest changed since rendering')
from spatial_renderer import effective_parameters as _effective_parameters
_current=_effective_parameters(_c['external']['spatial'])
_current.update(reference_mode=_c['external']['spatial'].get('reference_mode','descriptive'),excluded_targets=sorted(set(_c['external']['spatial'].get('excluded_targets',[]))),forced_genes=sorted(set(_c['external']['spatial'].get('forced_genes',[]))),only_gses=sorted(set(_c['external']['spatial'].get('only_gses',[]))))
if _digest(_current)!=PARAMETERS_SHA256:raise RuntimeError('Spatial parameters changed since rendering')
ONTOLOGY_SHA256=_digest(_ontology.fillna('').astype(str).sort_values(['level1_annotation','level2_annotation']).to_dict('records'))
STATE_IDENTITY_MAP=_build_state_identity_map(_ontology)
del _obs
'''
    pos=result.find('\n',result.find('from __future__ import'))+1;result=result[:pos]+guard+result[pos:]
    # No hidden environment reads allowed in generated spatial program.
    if 'os.environ.get' in result:raise ValueError('Unbound inherited environment setting remains')
    compile(result,'spatial_adapted.py','exec')
    return result
