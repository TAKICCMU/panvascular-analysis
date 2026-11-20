#!/usr/bin/env python3
"""Read-only external dependency audit. Does not load expression or execute tools."""
import argparse,json,sys
from pathlib import Path
import pandas as pd
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from external_modules import audit_resources
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--hash-resources',action='store_true');a=p.parse_args()
c=yaml.safe_load(Path(a.config).read_text());o=Path(a.out);o.mkdir(parents=True,exist_ok=False)
r=Path(__file__).resolve().parents[1]/'external_resource_inventory.tsv'
rows=pd.read_csv(r,sep='\t',keep_default_na=False).to_dict('records')
audit=audit_resources(rows,a.hash_resources)
pd.DataFrame(audit).to_csv(o/'external_resource_audit.tsv',sep='\t',index=False)
(o/'STATUS.json').write_text(json.dumps({'status':'WAITING_FOR_FINAL_L2_ATLAS','production_executed':False,'atlas_exists':Path(c['atlas']['path']).is_file(),'resources':audit},indent=2))
print(o/'STATUS.json')
