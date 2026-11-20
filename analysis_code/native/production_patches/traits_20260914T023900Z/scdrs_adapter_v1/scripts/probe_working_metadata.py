"""Read-only metadata probe; never substitutes working file for the final atlas."""
import json, sys, pathlib, importlib.metadata as im
import h5py, numpy as np, pandas as pd
from anndata.io import read_elem
root=pathlib.Path(sys.argv[1]).resolve(); out=pathlib.Path(sys.argv[2]).resolve(); out.mkdir(parents=True,exist_ok=True)
p=root/'panvascular_CPU_level1_L2_partial_metadata_v4_20260913_155051_440821.h5ad'
contract=json.load(open(root/'handoff_L1L2_final_v1/atlas_input_contract_v1.json'))
cols=list(dict.fromkeys(contract['data_layout']['metadata_columns_expected_from_v4']+['donor_identity_resolved','donor_resolution_status','donor_evidence_class','replicate_type','replicate_evidence_status','replicate_group_uid','overlap_class','recommended_analysis_action','primary_inference_constraint','merge_evidence','downstream_repeat_policy','condition_binary','cohort_disease_context']))
stat=p.stat(); audit={'source':str(p),'source_role':'WORKING_METADATA_ONLY_NOT_FINAL_INPUT','final_exists':(root/'panvascular_CPU_theta2_L1L2_final_v1.h5ad').exists(),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns,'missing_fields':[],'conflicts':{}}
with h5py.File(p,'r') as f:
    s=np.asarray(read_elem(f['obs/sample_uid'])).astype(str); keys,first=np.unique(s,return_index=True)
    df=pd.DataFrame({'sample_uid':keys,'n_cells':pd.Series(s).value_counts().reindex(keys).values})
    for c in cols:
        if c=='sample_uid':continue
        if c not in f['obs']:audit['missing_fields'].append(c);continue
        x=np.asarray(read_elem(f['obs'][c])); df[c]=x[first]
        u=pd.DataFrame({'s':s,'v':pd.Series(x).astype('string').fillna('<NA>')}).groupby('s')['v'].nunique(dropna=False)
        if (u>1).any(): audit['conflicts'][c]=u[u>1].to_dict()
    audit['shape']=list(f['X'].attrs['shape']); audit['obsm_keys']=list(f['obsm'].keys()); audit['layers']=list(f['layers'].keys())
    audit['l1_counts']=pd.Series(np.asarray(read_elem(f['obs/level1_annotation'])).astype(str)).value_counts().to_dict()
    for key in ['level2_annotation','L2_annotation']:
        if key in f['obs']:audit[key+'_counts']=pd.Series(np.asarray(read_elem(f['obs'][key])).astype(str)).value_counts().to_dict()
audit['stat_unchanged']=(p.stat().st_mtime_ns==stat.st_mtime_ns and p.stat().st_size==stat.st_size)
audit['versions']={}
for pkg in ['omicverse','scanpy','anndata','h5py','pandas','numpy','scipy','statsmodels','scdrs','liana','cell2location','PyYAML']:
    try:audit['versions'][pkg]=im.version(pkg)
    except im.PackageNotFoundError:audit['versions'][pkg]='NOT_INSTALLED_IN_THIS_ENV'
df.to_csv(out/'working_sample_metadata_v4.tsv',sep='\t',index=False)
(out/'READ_ONLY_INPUT_CHECK.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False,default=int))
print(json.dumps({k:audit[k] for k in ['final_exists','shape','stat_unchanged','missing_fields','conflicts','versions']},ensure_ascii=False,default=int))
