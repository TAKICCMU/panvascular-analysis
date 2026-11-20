from common import *
from datetime import datetime,timezone
z=processes();answer={'time_utc':datetime.now(timezone.utc).isoformat(),'actual_project_threads':sum(p['threads'] for p in z),'processes':z}
answer['statuses']={}
for name in ['covarnet_main_and_refit_sensitivity','scdrs_third_L2','communication_disease']:
 p=W/'status'/(name+'.json')
 if p.exists():answer['statuses'][name]=json.loads(p.read_text())
for name,p in [('pseudobulk',W/'pseudobulk/COMPLETE.json'),('communication',W/'communication/PROGRESS.json'),('descriptive',W/'descriptive_updated_v2/COMPLETE.json'),('GPU',W/'spatial_gpu_resume_v2/RUN_STATUS.json')]:
 if p.exists():answer['statuses'][name]=json.loads(p.read_text())
gp=W/'spatial_gpu_resume_v2/logs/production.log'
if gp.exists():
 with open(gp,'rb') as f:f.seek(max(0,gp.stat().st_size-4000));answer['GPU_tail']=f.read().decode(errors='replace').replace('\r','\n').splitlines()[-3:]
answer['GPU_usage']=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,memory.total','--format=csv,noheader'],text=True).strip()
answer['memory']=subprocess.check_output(['free','-h'],text=True);answer['disk']=subprocess.check_output(['df','-h','/data','/home'],text=True)
print(json.dumps(answer,indent=2,ensure_ascii=False))
