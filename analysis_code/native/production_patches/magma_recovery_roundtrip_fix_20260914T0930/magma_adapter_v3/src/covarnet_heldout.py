"""Explicit historical CD-NMF algorithm sensitivity, never primary nsNMF replacement."""
from pathlib import Path
import importlib.util, json
import numpy as np
import pandas as pd

def run(prepared, *, n_perm=200, min_heldout=5, seed=20260824):
    p=Path(prepared)
    source=Path(__file__).with_name('covarnet_legacy')/'heldout_cd_sensitivity.py'
    spec=importlib.util.spec_from_file_location('audited_legacy_heldout',source)
    legacy=importlib.util.module_from_spec(spec); spec.loader.exec_module(legacy)
    done=json.loads((p/'INPUT_MANIFEST.json').read_text())
    if done['status']!='PREPARED_NOT_EXECUTED' or done.get('phase')!='INPUTS_READY': raise ValueError('Prepared node eligibility not approved')
    for file,h in done['files'].items():
        if legacy.sha(p/file)!=h: raise ValueError('Input hash mismatch')
    stable=p/'stable_NMF'; paths=list(stable.glob('module_node_weights_K*.tsv'))
    if len(paths)!=1: raise ValueError('Unique new primary nsNMF basis required')
    w=legacy.matrix_file(paths[0]); k=w.shape[1]
    counts=legacy.matrix_file(p/'CoVarNet_node_counts.tsv').loc[w.index]
    md=pd.read_csv(p/'CoVarNet_sample_metadata.tsv',sep='\t').set_index('covarnet_sample_uid').loc[counts.columns]
    d=dict(counts_df=counts,counts=legacy.validate_counts(counts),k=k,modules=list(w.columns),reference_W=w.to_numpy(),metadata=md,all_atlas_denominator=md.n_total_atlas_cells.to_numpy(float))
    result=legacy.analyze(d,p/'heldout_CD_sensitivity',n_perm=n_perm,min_heldout=min_heldout,seed=seed,denominator='all_atlas',evidence=None)
    result['source_sha256']=legacy.sha(source)
    result['interpretation']='Additional historical algorithm sensitivity. Global fixed-margin null is supplementary; full-parent fixed-margin null is mandatory separately. Full-data state universe and K mean validation is conditional, not nested independent discovery.'
    (p/'heldout_CD_sensitivity/RESULT.json').write_text(json.dumps(result,indent=2))
    return result
