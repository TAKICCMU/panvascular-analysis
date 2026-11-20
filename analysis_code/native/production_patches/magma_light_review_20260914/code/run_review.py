import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['MPLBACKEND']='Agg'
import json,hashlib,csv,time,sys,fcntl,math
from pathlib import Path
from review_core import read_tsv,native_rows,validate_family
W=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2))
def threads():return len(list(Path('/proc/self/task').iterdir()))
def gate(c):
 sys.path.insert(0,'/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/communication_data_v1/code')
 from mount_resources import cpu_counts
 procs={}
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   d=dict(x.split(':',1) for x in (p/'status').read_text().splitlines() if ':' in x);procs[int(p.name)]={'ppid':int(d['PPid']),'state':d['State'].split()[0],'threads':int(d['Threads']),'cmd':(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')}
  except (OSError,KeyError,ValueError):pass
 live,ps,stopped=cpu_counts(procs,c['project_roots'],os.getpid());return {'other_live_threads':live,'projected_plus_light1_GPU4':live+5,'pass':live+5<=40,'stopped_processes':stopped}
def main():
 receipt=json.loads((W/'RECEIPT.json').read_text());assert receipt['scope']=='magma_completed_output_review_only'
 for f in receipt['bound_files']:assert sha(f['path'])==f['sha256'],f['path']
 c=json.loads((W/'config.json').read_text());out=Path(c['output']);source=Path(c['source']);assert not out.exists()
 while True:
  g=gate(c);write(W/'RESOURCE_STATUS.json',g)
  if g['pass']:break
  write(W/'STATUS.json',{'status':'WAITING_CPU_LIGHT_SLOT','pid':os.getpid()});time.sleep(20)
 assert threads()==1
 identity=json.loads((source.parent/'RUN_IDENTITY.json').read_text());assert identity['atlas_sha256']==c['atlas_sha256']
 status=json.loads((source.parent/'RUN_STATUS.json').read_text());assert status['magma']['status']=='COMPLETE' and status['magma']['exit_code']==0
 assert json.loads((source/'COMPLETE.json').read_text())['atlas_sha256']==c['atlas_sha256']
 assert json.loads((source/'RECOVERY_PROGRESS.json').read_text())['stage']=='ALL_12_ANALYSES_COMPLETE'
 summary=[];allrows=[];tables={};names={};inputs=[];start=time.time()
 for level in c['levels']:
  variables=read_tsv(source/(level+'_variables.tsv'));names[level]=[x['group'] for x in variables]
  for scenario in c['scenarios']:
   for weight in c['weights']:
    family=f'{scenario}_{level}_{weight}';path=source/(family+'_with_BH.tsv');rows=read_tsv(path)
    result=validate_family(rows,read_tsv(source/(family+'_model_mapping.tsv')),variables,native_rows(source/(family+'.gsa.out')),family)
    summary.append(dict(scenario=scenario,annotation_level=level,weighting=weight,BH_family=family,**result));tables[family]={x['group']:float(x['FDR_BH']) for x in rows}
    for row in rows:allrows.append(dict(scenario=scenario,annotation_level=level,weighting=weight,source_file=str(path),source_sha256=sha(path),**row))
 assert len(summary)==12;assert threads()==1
 out.mkdir(parents=True);(out/'figures').mkdir()
 for filename,rows in [('MAGMA_all_12_families_original_values.tsv',allrows),('MAGMA_scenario_level_weighting_summary.tsv',summary)]:
  with open(out/filename,'w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 write(out/'VALIDATION.json',{'status':'PASS','families':summary,'original_BETA_P_FDR_strings_preserved':True,'Average_nuisance_checked_not_included_in_target_BH':True,'native_MAGMA_rerun':False})
 import numpy as np
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 assert threads()==1
 for level in c['levels']:
  groups=names[level];families=[f'{s}_{level}_{w}' for s in c['scenarios'] for w in c['weights']]
  q=np.array([[tables[f][g] for f in families] for g in groups]);colors=-np.log10(np.maximum(q,np.finfo(float).tiny))
  fig,ax=plt.subplots(figsize=(12,max(6,len(groups)*.245+2.2)))
  im=ax.imshow(colors,aspect='auto',cmap='viridis',vmin=0)
  ax.set_yticks(range(len(groups)),groups,fontsize=7 if level=='L2' else 9)
  ax.set_xticks(range(6),[s.replace('primary_','').replace('leave_one_CAD_','LOO CAD ')+'\n'+w.replace('_',' ') for s in c['scenarios'] for w in c['weights']],fontsize=9)
  for i,j in zip(*np.where(q<.05)):ax.text(j,i,'*',ha='center',va='center',color='white',fontsize=10)
  ax.set_title(f'MAGMA {level}: original within-family FDR\n* FDR < 0.05; no new test or combined BH',fontsize=12)
  bar=fig.colorbar(im,ax=ax,pad=.025,shrink=.45);bar.set_label('-log10(original FDR_BH)');fig.tight_layout()
  fig.savefig(out/'figures'/f'MAGMA_{level}_original_FDR_heatmap.png',dpi=180);fig.savefig(out/'figures'/f'MAGMA_{level}_original_FDR_heatmap.pdf');plt.close(fig)
 assert threads()==1
 write(out/'SOURCE_SHA256.json',receipt['bound_files'])
 write(out/'PARAMETERS.json',{'atlas_sha256':c['atlas_sha256'],'source':str(source),'FDR_source':'existing FDR_BH, checked by recomputation but not replaced','BH_family':'scenario x annotation_level x weighting; unchanged','condition':'Average','direction':'one-sided positive','plot_zero_FDR_floor':float(np.finfo(float).tiny),'scenarios':c['scenarios'],'weighting_not_independent_replication':True,'cpu_threads':threads(),'numpy':np.__version__,'matplotlib':matplotlib.__version__})
 text='''# 已完成MAGMA结果的轻量汇总和绘图\n\n本次只读汇总当前最终atlas的MAGMA r3十二个已完成分析家族，没有重新运行MAGMA、改变BETA/P/FDR，或改变BH检验家族。每家族VARIABLE、MODEL、身份映射、原生目标行和Average条件行均通过核验。原FDR可复算一致（容差1e-12），复算仅用于验证。\n\n全量合并表保留原表BETA、P、FDR等字段的原始字符串数值；绘图仅展示原FDR。热图包含全部受检L1/L2，星号表示原FDR<0.05，不把两个weighting规格视作独立重复证据，也不为汇总重新跨家族BH。\n\nMAGMA十二组合已完成，不代表scDRS原生组检验已完成；本报告不能称全部遗传分析已结束。该图不是MAGMA基因通路GSEA/ORA，也不涉及临床DE/meta、疾病通讯、CoVarNet。\n\n结果存于scdrs_storage_v1/magma_completed_review_v1仅为共享已批准存储预算，方法归属仍是MAGMA，不能据目录名称合并两种方法结论。\n\n'''
 text+='|Scenario|Level|Weighting|Tested|Original FDR<0.05|\n|---|---|---|---:|---:|\n'
 for x in summary:text+=f"|{x['scenario']}|{x['annotation_level']}|{x['weighting']}|{x['n_tested']}|{x['n_FDR_lt_0_05']}|\n"
 (out/'README_CN.md').write_text(text)
 manifest={str(p.relative_to(out)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in out.rglob('*') if p.is_file()};write(out/'FILES.json',manifest)
 result={'status':'COMPLETE_MAGMA_EXISTING_RESULT_REVIEW','pid':os.getpid(),'families':12,'rows':len(allrows),'cpu_threads':threads(),'seconds':time.time()-start,'output':str(out),'output_manifest_sha256':sha(out/'FILES.json'),'no_new_scientific_tests':True};write(W/'COMPLETE.json',result);write(W/'STATUS.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':
 with open(W/'.job.lock','a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  try:main()
  except Exception:
   import traceback;traceback.print_exc();write(W/'FAILED.json',{'error':traceback.format_exc()});raise
