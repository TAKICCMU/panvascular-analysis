import os,sys,time,json,subprocess,hashlib,csv,signal
from pathlib import Path
D=Path(__file__).parent.resolve();os.chdir(D.parent)
REMOTE='codex-shpc-gpu:/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1/overnight_manuscript_integration_20260917/candidate_prioritization/'
status={'controller_pid':os.getpid(),'status':'STARTED','started_epoch':time.time(),'local_dir':str(D),'remote':REMOTE}
def save(**kwargs):
 status.update(kwargs);(D/'DURABLE_STATUS.json').write_text(json.dumps(status,indent=2))
def run(script,log):
 with open(D/log,'a') as f:
  p=subprocess.Popen([sys.executable,str(D/script)],stdout=f,stderr=subprocess.STDOUT,start_new_session=True);save(stage=script,child_pid=p.pid);return p.wait()
def hashes():
 h={}
 for p in D.rglob('*'):
  if p.is_file() and p.name not in ['DELIVERY_SHA256.json','controller.log','DURABLE_STATUS.json'] and not p.name.endswith('.part'):
   a=hashlib.sha256()
   with p.open('rb') as f:
    for b in iter(lambda:f.read(1048576),b''):a.update(b)
   h[str(p.relative_to(D))]=a.hexdigest()
 (D/'DELIVERY_SHA256.json').write_text(json.dumps(h,indent=2))
def upload():
 return subprocess.run(['/usr/bin/rsync','-a','-e','ssh -o BatchMode=yes',str(D)+'/',REMOTE],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
try:
 save(status='WAITING_EXISTING_QTL_DOWNLOAD',adopted_fetch_pid=54658)
 deadline=time.time()+7200
 while True:
  if (D/'NEW_QTL_COMPLETE.json').exists():break
  if (D/'NEW_QTL_BLOCKED.json').exists():raise RuntimeError('QTL service blocked after three failed requests; no automatic retry')
  if time.time()>deadline:raise RuntimeError('Download wait exceeded 2h; incomplete')
  try:os.kill(54658,0)
  except ProcessLookupError:
   save(status='RESUMING_QTL_DOWNLOAD')
   code=run('fetch_new_qtl.py','new_qtl_fetch.log')
   if code:raise RuntimeError('QTL download failed with exit '+str(code))
   break
  time.sleep(15)
 if not (D/'NEW_QTL_COMPLETE.json').exists():raise RuntimeError('Missing download completion marker')
 save(status='RUNNING_NEW_COLOC')
 code=run('finish_new_qtl_coloc.py','finish_new_qtl_coloc.log')
 if code:raise RuntimeError('New coloc or validation failed with exit '+str(code))
 if not (D/'ALL29_COMPLETE.json').exists():raise RuntimeError('Missing ALL29 validation marker')
 complete=json.loads((D/'ALL29_COMPLETE.json').read_text());save(status=complete['status'],finished_epoch=time.time(),completion=complete)
 report='\n\n## 补充分析最终自动验证状态\n\n'+json.dumps(complete,ensure_ascii=False,indent=2)+'\n\n本次新增完整cis QTL已按原ABF方法分析，结果在all29_complete_cis_coloc_results.tsv与Tier1_integrated_evidence_FINAL.tsv。strong列表为原方法区域阈值分类，组织max H4未经跨组织多重检验，不是细胞类型因果认证。缺失/失败context必须按coverage/status解释。没有运行新的SuSiE。\n'
 with open(D/'REPORT_CN.md','a') as f:f.write(report)
except Exception as e:
 save(status='BLOCKED_NOT_COMPLETE',error=str(e),finished_epoch=time.time())
 with open(D/'REPORT_CN.md','a') as f:f.write('\n\n## 补充分析阻断\n\n'+str(e)+'。已完成候选Tier与17基因ABF不受影响，12新基因不得称全部完成。见DURABLE_STATUS.json与日志。\n')
finally:
 hashes();u=upload();save(upload_exit_code=u.returncode,upload_output_tail=u.stdout[-1500:]);
 # Deliver terminal status after prior rsync so remote also receives final status.
 subprocess.run(['/usr/bin/scp','-q',str(D/'DURABLE_STATUS.json'),REMOTE],check=False)
