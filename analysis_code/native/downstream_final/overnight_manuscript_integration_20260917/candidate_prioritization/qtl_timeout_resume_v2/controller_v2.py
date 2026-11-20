import os,sys,time,json,subprocess,hashlib,csv
from pathlib import Path
V=Path(__file__).parent.resolve();D=V.parent;os.chdir(V)
EXPECTED='a872fc169066e37de51bc5292c10d030094d87ec676bdb6bd6cf2319afe3711b'
REMOTE='codex-shpc-gpu:/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1/overnight_manuscript_integration_20260917/candidate_prioritization/qtl_timeout_resume_v2/'
status={'controller_pid':os.getpid(),'status':'STARTED','started_epoch':time.time(),'local_dir':str(V),'target_interface_sha256':EXPECTED}
def frozen():assert hashlib.sha256((D/'Tier1_spatial_interface.tsv').read_bytes()).hexdigest()==EXPECTED,'Frozen target changed externally'
def save(**kw):status.update(kw);(V/'STATUS.json').write_text(json.dumps(status,indent=2))
def run(args,log):
 with open(V/log,'a') as f:
  p=subprocess.Popen(args,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'));save(child_pid=p.pid,stage=args[-1]);return p.wait()
try:
 frozen();pilot=json.loads((V/'PILOT_RESULT.json').read_text());assert pilot['status'] in ['DOWNLOADED','REUSED','NO_TARGET_GENE_QTL'],'Pilot failed: no batch permitted'
 save(status='FETCHING_WITH_LEGACY_360S_TIMEOUT')
 code=run([sys.executable,str(V/'fetch_v2.py'),'--all'],'fetch_v2.log')
 if code:raise RuntimeError('Bounded QTL resume failed exit '+str(code))
 done=json.loads((V/'NEW_QTL_COMPLETE.json').read_text());assert done['rows']==84,'Incomplete manifest (<84)'
 R=V/'coloc_new12';(R/'00_metadata').mkdir(parents=True,exist_ok=True);(R/'00_metadata/qtl_download_manifest.tsv').write_bytes((V/'new_qtl_manifest.tsv').read_bytes())
 s=(D/'original_run_full_coloc.R').read_text().replace('ROOT <- "/Users/taki/Desktop/Codex任务1/Tier1_full_colocalisation_v1"','ROOT <- '+json.dumps(str(R),ensure_ascii=False));s=s.replace('library(data.table)','library(data.table)\n  setDTthreads(1)');s=s.split('# Plot gene-level maxima')[0]+'\nwriteLines(capture.output(sessionInfo()), file.path(ROOT, "00_metadata", "R_sessionInfo.txt"))\n';script=V/'run_new12_coloc.R';script.write_text(s)
 save(status='RUNNING_INHERITED_ABF');code=run(['/usr/local/bin/Rscript',str(script)],'coloc_v2.log')
 if code:raise RuntimeError('ABF failed exit '+str(code))
 code=run([sys.executable,str(V/'validate_all29_v2.py')],'validation_v2.log')
 if code:raise RuntimeError('Validation failed exit '+str(code))
 frozen();complete=json.loads((V/'ALL29_COMPLETE.json').read_text());save(status=complete['status'],completion=complete,finished_epoch=time.time())
 (V/'REPORT_CN.md').write_text('# QTL超时恢复v2\n\n原45秒上限引起真实请求超时，本版本恢复原360秒，仅对失败/未尝试区域发请求；62个既存成功或无目标QTL区域复用。先行SKI Artery_Aorta真实请求验证通过后才启动批次。\n\n'+json.dumps(complete,ensure_ascii=False,indent=2)+'\n\n29个Tier目标名单与接口字节未改变；所有新增QTL只进入evidence，不改变Tier。ABF严格继承原方法/参数，单线程，不运行SuSiE。原v1 blocked状态保留；本版本结果、报告、hash独立。分类中有失败资源则明确标ANALYSIS_COMPLETE_WITH_QTL_RESOURCE_GAPS，不能称84资源全部完成。\n')
except Exception as e:
 save(status='BLOCKED_NOT_COMPLETE',error=str(e),finished_epoch=time.time());(V/'REPORT_CN.md').write_text('# QTL超时恢复v2未完成\n\n'+str(e)+'\n\n冻结29目标不变。原有候选与17基因ABF仍可使用；新12基因完整补充未完成。查看STATUS.json和日志，不以PID判定完成。\n')
finally:
 h={}
 for p in V.rglob('*'):
  if p.is_file() and p.name not in ['SHA256.json','controller_v2.log','STATUS.json'] and not p.name.endswith('.part'):h[str(p.relative_to(V))]=hashlib.sha256(p.read_bytes()).hexdigest()
 (V/'SHA256.json').write_text(json.dumps(h,indent=2))
 u=subprocess.run(['/usr/bin/rsync','-a','-e','ssh -o BatchMode=yes',str(V)+'/',REMOTE],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);save(upload_exit_code=u.returncode,upload_tail=u.stdout[-1000:]);subprocess.run(['/usr/bin/scp','-q',str(V/'STATUS.json'),REMOTE],check=False)
