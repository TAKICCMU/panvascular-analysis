#!/usr/bin/env python3
"""Read existing artifacts and prepare unsigned review manifests; no model fitting."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.external_gate import load_config
from src.core import preflight
from src.spatial_projection_contract import weights_record,abundance_record,file_digest
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--weights',required=True);p.add_argument('--covarnet-run',required=True);p.add_argument('--abundance',required=True);p.add_argument('--posterior',required=True);p.add_argument('--posterior-provenance',required=True);p.add_argument('--reference',required=True);p.add_argument('--abundance-key',choices=['state_id','parent_qualified','safe_label'],required=True);p.add_argument('--source-obsm-key',default='q05_cell_abundance_w_sf');p.add_argument('--out',required=True);a=p.parse_args()
c,cp=load_config(a.config);out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
obs,ont,info=preflight(c,out/'release_preflight.json');del obs
w=weights_record(a.weights,a.covarnet_run,ont,c['atlas']['input_sha256']);ab=abundance_record(a.abundance,a.posterior,a.posterior_provenance,a.reference,ont,c,a.abundance_key,a.source_obsm_key)
for name,value in [('weights_manifest.json',w),('abundance_manifest.json',ab)]:
 with (out/name).open('x') as f:json.dump(value,f,indent=2)
candidate={'approved':False,'weights_manifest':str((out/'weights_manifest.json').resolve()),'weights_manifest_sha256':file_digest(out/'weights_manifest.json'),'abundance_manifest':str((out/'abundance_manifest.json').resolve()),'abundance_manifest_sha256':file_digest(out/'abundance_manifest.json'),'covarnet_run_dir':w['run_dir'],'covarnet_run_identity_sha256':w['run_identity_sha256'],'posterior_path':ab['posterior_path'],'posterior_sha256':ab['posterior_sha256']}
(out/'MODULE_PROJECTION_SELECTION_NOT_APPROVED.json').write_text(json.dumps(candidate,indent=2))
print(json.dumps({'status':'PREPARED_NOT_EXECUTED','approved':False,'out':str(out),'instruction':'Review selection and obtain a new content-bound approval; no automatic config writeback'}))
