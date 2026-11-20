#!/usr/bin/env python3
"""Run only after root's full-SHA certificates are complete. Seals a new attempt; no execution here."""
from pathlib import Path
import json,hashlib,datetime
W=Path(__file__).resolve().parents[1];O=W.parents[1];R=O.parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
old_path=O/'configs/scdrs_APPROVED.json';old=json.loads(old_path.read_text());certdir=R/'maintenance_resume_20260914'
certs=[certdir/'ATLAS_VERIFIED.json',certdir/'SCDRS_DERIVED_VERIFIED.json']
for p in certs:
 z=json.loads(p.read_text())
 if not z['full_sha_verified'] or z['observed_sha256']!=z['expected_sha256']:raise ValueError('Root certificate incomplete '+str(p))
 if z['boot_id']!=Path('/proc/sys/kernel/random/boot_id').read_text().strip():raise ValueError('Wrong boot')
report=O/'outputs/derived/DERIVED_OBJECT.json';dr=json.loads(report.read_text());ledger=certdir/'RESOURCE_RESUME_LEDGER.json'
assert sha(ledger)=='d0a80a5e335f7a93c26f29198db50bf24f3776d6ed574ab4dba9fcfa48b8826e'
c=dict(root=str(R),workdir=str(W),original_config=str(old_path),original_receipt=str(O/'approvals/scdrs_APPROVAL_RECEIPT_r3.json'),adapter=str(O/'code/scdrs_adapter.py'),derived=str(O/'outputs/derived/scdrs_raw_symbols.h5ad'),derived_sha256=dr['derived_h5ad_sha256'],derived_report=str(report),covariates=str(O/'outputs/derived/native_covariates.cov.gz'),covariates_sha256=dr['covariates_sha256'],scenarios=old['scenarios'],native_threads=12,n_ctrl=1000,graph_policy=old['graph_policy'],scdrs_executable=old['scdrs_executable'],group_keys=dr['group_keys'],output=str(W/'outputs'),rebuild_derived=False,old_native_outputs=str(O/'outputs/native'),lock=str(R/'production_launches/execution_control_20260914/.scdrs_full_project.lock'),communication_root=str(R/'production_launches/communication_L1_20260914T015040Z'),gpu_root=str(R/'production_launches/spatial_gpu_20260913T183753Z'),resource_ledger=str(ledger),resource_ledger_sha256=sha(ledger),authenticated_files=[],bound_files=[])
for p in certs:
 z=json.loads(p.read_text());c['authenticated_files'].append(dict(certificate=str(p),selector=[],path=z['path'],sha256=z['observed_sha256']))
assert c['authenticated_files'][0]['sha256']==old['atlas_sha256'] and c['authenticated_files'][1]['sha256']==c['derived_sha256']
files=certs+[certdir/'INPUT_REVALIDATION.json',ledger,old_path,Path(c['original_receipt']),Path(c['adapter']),report,O/'STAGE_STATUS.json',R/'maintenance_pause_20260914/MAINTENANCE_VERIFIED.json',R/'maintenance_pause_20260914/SCDRS_NATIVE_INVENTORY.json',W/'TEST_RESULTS.json',W/'HISTORICAL_REUSE_AUDIT.json']+list((W/'code').glob('*.py'))
for p in files:c['bound_files'].append(dict(path=str(p),sha256=sha(p)))
if (W/'config.json').exists() or (W/'RECEIPT.json').exists():raise FileExistsError('Already sealed attempt')
(W/'config.json').write_text(json.dumps(c,indent=2))
r=dict(mode='approved_execution',allowed_modules=['scdrs_scoring_and_native_downstream'],authorization_source='parent_relay_of_user_server_restarted_continue_request_20260914',is_user_digital_signature=False,config=str(W/'config.json'),config_sha256=sha(W/'config.json'),bound_files=c['bound_files']+[dict(path=str(W/'config.json'),sha256=sha(W/'config.json'))],original_receipt_sha256=sha(c['original_receipt']),rebuild_derived=False,old_outputs_overwritten=False,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(W/'RECEIPT.json').write_text(json.dumps(r,indent=2));print(json.dumps(dict(config_sha256=sha(W/'config.json'),receipt_sha256=sha(W/'RECEIPT.json'))))
