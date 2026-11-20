"""Read-only status of the authorized MAGMA/scDRS repaired branches."""
from pathlib import Path
import json,datetime,shutil,subprocess
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches')
out={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'free_disk_GiB':round(shutil.disk_usage(R).free/2**30,2)}
for role,name in [('MAGMA','magma_recovery_roundtrip_fix_20260914T0930'),('scDRS','scdrs_repair_20260914T0910')]:
    p=R/name;d={'directory':str(p),'exists':p.exists()}
    for f in ['STAGE_STATUS.json','RESOURCE_STATUS.json','EXIT_STATUS.json','LAUNCH_MANIFEST.json']:
        q=p/f
        if q.exists():
            x=json.loads(q.read_text())
            if f=='RESOURCE_STATUS.json':x={k:v for k,v in x.items() if k!='other_project_processes'}
            if f=='LAUNCH_MANIFEST.json':
                x={k:v for k,v in x.items() if k!='science_started'}
                x['note']='Static launch definition; current progress is reported by runtime status files below.'
            d[f]=x
    runs=[]
    for q in p.glob('outputs/runs/*/RUN_STATUS.json'):
        v={'run':str(q.parent),'status':json.loads(q.read_text()),'magma_completed_BH_families':len(list((q.parent/'magma').glob('*_with_BH.tsv')))}
        for z in (q.parent/'magma').glob('*PROGRESS*.json'):v[z.name]=json.loads(z.read_text())
        runs.append(v)
    if runs:d['runs']=runs
    d['native_scdrs_completed_scenarios']=[q.parent.name for q in p.glob('outputs/native/*/SCHEMA_VALIDATION.json')]
    out[role]=d
mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines() if ':' in x}
out['available_RAM_GiB']=round(mem['MemAvailable']/1024**2,2)
try:out['GPU']=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,memory.total','--format=csv,noheader'],text=True).strip()
except Exception as e:out['GPU_error']=str(e)
print(json.dumps(out,ensure_ascii=False,indent=2))
