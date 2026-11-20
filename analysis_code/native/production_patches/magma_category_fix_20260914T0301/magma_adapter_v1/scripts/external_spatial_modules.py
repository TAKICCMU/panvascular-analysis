#!/usr/bin/env python3
"""Exact final ontology state IDs -> inherited weighted CLR spatial projection."""
import argparse,sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_gate import bootstrap_cli,approved_preflight
if __name__=='__main__':bootstrap_cli()
import pandas as pd
from external_statistics import spatial_module_projection
from state_identity import build_state_identity_map
from spatial_projection_contract import validate_projection_inputs,table,verify_hash

def main():
    p=argparse.ArgumentParser();p.add_argument('--abundance',required=True);p.add_argument('--weights',required=True);p.add_argument('--out',required=True);p.add_argument('--config',required=True)
    p.add_argument('--abundance-key',required=True,choices=['state_id','parent_qualified','safe_label']);a=p.parse_args()
    c,obs,ontology,info=approved_preflight(a.config,Path(a.out).with_name(Path(a.out).name+'.authorization'))
    evidence=validate_projection_inputs(c,ontology,a.abundance,a.weights,a.abundance_key)
    del obs
    ab=table(a.abundance);w=table(a.weights)
    verify_hash(a.abundance,evidence['abundance']['abundance_sha256'],'abundance re-read')
    verify_hash(a.weights,evidence['weights']['weights_sha256'],'W re-read')
    state_map=build_state_identity_map(ontology);out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    (out/'PROJECTION_INPUT_VALIDATION.json').write_text(json.dumps(evidence,indent=2))
    state_map.to_csv(out/'state_identity_map.tsv',sep='\t',index=False)
    try:scores,audit=spatial_module_projection(ab,w,state_map,abundance_key=a.abundance_key)
    except ValueError as exc:
        if hasattr(exc,'audit'):exc.audit.to_csv(out/'state_weight_mapping_audit_BLOCKED.tsv',sep='\t',index=False)
        raise
    scores.to_csv(out/'spatial_module_scores.tsv.gz',sep='\t');audit.to_csv(out/'state_weight_mapping_audit.tsv',sep='\t',index=False)
if __name__=='__main__':main()
