#!/usr/bin/env python3
"""Collect new specimen outputs into inherited long-table schema; no missing-edge zeros."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_gate import bootstrap_cli,approved_preflight
bootstrap_cli()
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--runs',required=True);p.add_argument('--design',required=True);p.add_argument('--out',required=True);v=p.parse_args()
approved_preflight(v.config,Path(v.out).with_name(Path(v.out).name+'.authorization'))
runs=pd.read_csv(v.runs,sep='\t');design=pd.read_csv(v.design,sep='\t').set_index('sample_unit',verify_integrity=True)
o=Path(v.out);o.mkdir(parents=True,exist_ok=False)
for method,file in [('LIANA_plus','liana_interactions.tsv'),('CellChat','cellchat_interactions.tsv')]:
    parts=[]
    for r in runs.itertuples():
        if r.sample not in design.index:raise ValueError('Specimen missing reviewed design')
        source=Path(r.directory)/file
        if not source.is_file():raise FileNotFoundError(source)
        x=pd.read_csv(source,sep='\t');z=pd.DataFrame(index=x.index)
        for col,choices in {'source':['source','sender'],'target':['target','receiver'],'ligand':['ligand_complex','ligand','ligand.complex'],'receptor':['receptor_complex','receptor','receptor.complex']}.items():
            key=next((k for k in choices if k in x),None)
            if key is None:raise ValueError('Unexpected native LR schema')
            z[col]=x[key].astype(str)
        key='magnitude_rank' if method=='LIANA_plus' else 'prob'
        if key not in x:raise ValueError(f'Missing prespecified native score {key}')
        z['score']=1-pd.to_numeric(x[key]) if method=='LIANA_plus' else pd.to_numeric(x[key])
        z['method']=method;z['sample_unit']=r.sample;z['threshold_scheme']=f'min_{r.min_cells}'
        z['interaction_key']=z[['source','target','ligand','receptor']].agg('|'.join,axis=1)
        z['signaling_category']=x['annotation'] if 'annotation' in x else 'Not_applicable'
        parts.append(z)
    result=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame()
    if len(result) and result.duplicated(['sample_unit','interaction_key','threshold_scheme']).any():raise ValueError('Duplicate LR rows require explicit database-specific resolution')
    result.to_csv(o/f'{method}_long.tsv.gz',sep='\t',index=False)
