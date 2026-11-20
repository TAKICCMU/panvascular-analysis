import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import sys,json,csv,time,subprocess,hashlib,fcntl
from pathlib import Path
B=Path(__file__).resolve().parents[1];W=B.parent
sys.path.insert(0,str(W/'code'));from common import processes

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for z in iter(lambda:f.read(8<<20),b''):h.update(z)
 return h.hexdigest()
def dump(p,d):
 p=Path(p);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(d,ensure_ascii=False,indent=2,default=str));t.replace(p)
receipt=json.loads((B/'RECEIPT.json').read_text())
assert json.loads((B/'tests/native_v1/PASS.json').read_text())['status']=='PASS'
for p,h in receipt['bound_files'].items():assert sha(p)==h,p
assert json.loads((W/'DISEASE_CONTEXT_REVIEW_HOLD.json').read_text())['pseudobulk_3background_user_approved'] is True
rows=list(csv.DictReader((B/'tables/model_roster.tsv').open(),delimiter='\t'));priority=['Macrophage','Fibroblast','Smooth muscle cell','Endothelial','Mural transitional cell','Monocyte','T cell','NK cell']
rows.sort(key=lambda r:(r['level']!='L1',priority.index(r['cell_state']) if r['cell_state'] in priority else 100,r['cell_state'],r['context']))
todo=[r for r in rows if r['status']=='READY_FOR_DESIGN_CHECK'];active=[];ended=[]
(B/'results').mkdir(exist_ok=True);(B/'logs').mkdir(exist_ok=True)
if (B/'EARLY_VALIDATED_MODELS.json').exists():
 for e in json.loads((B/'EARLY_VALIDATED_MODELS.json').read_text()):
  rp=Path(e['receipt']);assert sha(rp)==e['receipt_sha256'];er=json.loads(rp.read_text())
  for fp,h in er['bound_files'].items():assert sha(fp)==h,fp
  res=B/'results'/e['model_id'];assert sha(res/'MODEL.json')==e['MODEL_sha256']
  if e['DE_sha256']:assert sha(res/'DE.tsv')==e['DE_sha256']
  row=next(r for r in todo if r['model_id']==e['model_id']);todo.remove(row);ended.append({**row,'exit_code':0,'result_status':e['result_status'],'pid':None,'reused_early_native_fit':True})
started=time.time();dump(B/'LAUNCHED.json',{'pid':os.getpid(),'command':sys.argv,'started':started,'receipt_sha256':sha(B/'RECEIPT.json'),'models_in_queue':len(todo)})

def poll():
 global active
 for a in active[:]:
  rc=a['p'].poll()
  if rc is None:continue
  a['log'].close();a['lock'].close();active.remove(a);p=B/'results'/a['row']['model_id']/'MODEL.json'
  status=json.loads(p.read_text()).get('status') if p.exists() else 'FAILED'
  ended.append({**a['row'],'exit_code':rc,'result_status':status,'pid':a['p'].pid})
  print(a['row']['model_id'],status,'exit',rc,flush=True)
  # Refresh readable CKB results without a separate heavy process.
  ck=[]
  for q in ended:
   pp=B/'results'/q['model_id']/'CKB_result.tsv'
   if not pp.exists():continue
   vals=list(csv.DictReader(pp.open(),delimiter='\t'))
   if vals:ck.append({**{k:q[k] for k in ['model_id','level','context','cell_state']},**vals[0]})
   else:ck.append({**{k:q[k] for k in ['model_id','level','context','cell_state']},'gene_id':'ENSG00000166165','gene_symbol':'CKB','status':'FILTERED_LOW_EXPRESSION'})
  if ck:
   keys=list(dict.fromkeys(k for q in ck for k in q));p=B/'tables/CKB_results_current.tsv';t=p.with_suffix('.tmp')
   with t.open('w') as f:w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(ck)
   t.replace(p)
while todo or active:
 poll();ps=processes();threads=sum(p['threads'] for p in ps)
 other=sum([any('covarnet_nmf.R' in p['cmd'] and not p['state'].startswith('T') for p in ps),any('communication_worker' in p['cmd'] and not p['state'].startswith('T') for p in ps),any('scdrs perform-downstream' in p['cmd'] and not p['state'].startswith('T') for p in ps)])
 if todo and len(active)<min(2,max(0,4-other)) and threads+2<=40:
  handle=None
  for i in range(4):
   h=(W/'locks'/('cpu_heavy_'+str(i)+'.lock')).open('a')
   try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB);handle=h;break
   except BlockingIOError:h.close()
  if handle:
   row=todo.pop(0);p=Path(row['input']);cmd=['/usr/bin/Rscript','--vanilla',str(B/'code/fit.R'),str(p),str(B/'genes.tsv'),str(B/'results'/row['model_id'])]
   assert sha(B/'code/fit.R')==receipt['bound_files'][str(B/'code/fit.R')]
   f=(B/'logs'/(row['model_id']+'.log')).open('x');proc=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
   active.append({'p':proc,'row':row,'log':f,'lock':handle});print('START',row['model_id'],proc.pid,'threads_before',threads,flush=True)
 dump(B/'PROGRESS.json',{'status':'RUNNING' if active else 'WAITING_RESOURCE_BUDGET','controller_pid':os.getpid(),'active':[{'pid':a['p'].pid,**a['row']} for a in active],'completed':len(ended),'pending':len(todo),'ended':ended,'actual_project_threads':threads,'other_CPU_heavy_branches':other,'time':time.time()})
 if active or todo:time.sleep(3)
poll()
dump(B/'COMPLETE.json',{'status':'COMPLETE' if all(x['exit_code']==0 for x in ended) else 'COMPLETE_WITH_FAILURES','models':len(ended),'failed':[q for q in ended if q['exit_code']!=0],'not_estimable':[q for q in ended if q['result_status']=='NOT_ESTIMABLE'],'elapsed_seconds':time.time()-started,'results':{q['model_id']:sha(B/'results'/q['model_id']/'DE.tsv') for q in ended if (B/'results'/q['model_id']/'DE.tsv').exists()}})
print('ALL_PRIMARY_MODELS_ENDED',flush=True)

subprocess.run(["/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python",str(B/"code/report.py")],check=True)
