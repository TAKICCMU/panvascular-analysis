"""Complete the isolated CellHint/scIB dependency chain, preserving mainline inputs."""
import os
os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',JAX_PLATFORM_NAME='cpu')
import hashlib,json,subprocess,time,shutil,sys
from pathlib import Path
import pandas as pd

R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907')
W=R/'02_scib_benchmark/studywise_sampleonly_v1'
O=R/'06_reanalysis_takeover_20260908/scib';C=O/'code'
OV='/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python'
SCIB='/home/anzhen_vhd/miniconda3/envs/scib/bin/python'

def fh(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def status(phase,**kw):
 p={'phase':phase,'updated_unix':time.time(),**kw};(O/'TAKEOVER_STATUS.json').write_text(json.dumps(p,indent=2));print(json.dumps(p),flush=True)
def run(cmd,log,envname):
 env=os.environ.copy();env['LD_LIBRARY_PATH']=f'/home/anzhen_vhd/miniconda3/envs/{envname}/lib';env['JAX_PLATFORM_NAME']='cpu'
 with open(log,'ab') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
def wait(p,phase):
 watch={'waiting_pilot_meta3':(2403745,'13_run_cellhint.py'),
        'waiting_pilot_meta5':(2403746,'13_run_cellhint.py'),
        'waiting_primary_shared_harmonization':(2508566,'prepare_cellhint_primary_shared.py')}.get(phase)
 while not p.exists():
  if watch:
   pid,token=watch;cmd=Path(f'/proc/{pid}/cmdline')
   if not cmd.exists() or token not in cmd.read_bytes().decode(errors='replace'):
    raise RuntimeError(f'Expected upstream process {pid} is no longer running and output is absent: {p}')
  status(phase,waiting_for=str(p),upstream_pid=watch[0] if watch else None,
         eta_note='Pilot remaining runtime estimated 2-5 hours or more at handover; not guaranteed. Primary shared harmonization runs concurrently.');time.sleep(60)
def snapshot(src,dst):
 dst.parent.mkdir(parents=True,exist_ok=True)
 if not dst.exists():shutil.copy2(src,dst)
 if fh(src)!=fh(dst):raise RuntimeError(f'Frozen snapshot mismatch: {src}')
def evaluate(stage,meta,config):
 out=O/stage/'cellhint_metrics'/config
 inp=W/f'inputs/scib_{stage}_sampleonly_cap{500 if stage=="pilot" else 1500}_hvg2000.h5ad'
 run([SCIB,str(C/'evaluate_cellhint_graph.py'),'--input',str(inp),'--metadata',str(meta),'--output',str(out),'--stage',stage],O/'logs'/f'{stage}_{config}_metrics.log','scib')
 return out/'cellhint_metrics_by_K.tsv'
def combine_score(stage,extra):
 original=W/stage/'metrics/scib_metrics_by_K.tsv';snap=O/stage/'original_scib_metrics_by_K.tsv';snapshot(original,snap)
 table=pd.concat([pd.read_csv(snap,sep='\t')]+[pd.read_csv(p,sep='\t') for p in extra],ignore_index=True)
 dest=O/stage/'combined_raw_metrics_by_K.tsv';table.to_csv(dest,sep='\t',index=False)
 run([OV,str(C/'recompute_unified_ranking.py'),'--metrics-tsv',str(dest),'--output-dir',str(O/stage/'unified')],O/'logs'/f'{stage}_unified_ranking.log','omicverse')
 return pd.read_csv(O/stage/'unified/unified_scib_configuration_ranking.tsv',sep='\t')
def main():
 for d in ['logs','pilot/embeddings','primary']:(O/d).mkdir(parents=True,exist_ok=True)
 cellhint=[]
 for config in ['meta3','meta5']:
  src=W/'pilot/embeddings'/f'cellhint__{config}.json';wait(src,f'waiting_pilot_{config}')
  meta=O/'pilot/embeddings'/src.name
  snapshot(src,meta);snapshot(src.with_name(src.stem+'.graph.npz'),meta.with_name(meta.stem+'.graph.npz'))
  status(f'evaluating_pilot_{config}');cellhint.append(evaluate('pilot',meta,config))
 ranks=combine_score('pilot',cellhint)
 candidates=ranks[(ranks.method=='cellhint')&ranks.eligible_for_formal_ranking]
 if len(candidates)!=2:raise RuntimeError('Both CellHint pilot configurations must pass before selection')
 selected=candidates.sort_values(['mean_unified_overall_score','config_id'],ascending=[False,True]).iloc[0].config_id
 oldsel=json.loads((W/'primary/pilot_selected_configurations.json').read_text())['selected']
 newsel=ranks.sort_values('mean_unified_overall_score',ascending=False).groupby('method').first().config_id.to_dict()
 changes={m:{'frozen':c,'unified_pilot':newsel.get(m)} for m,c in oldsel.items() if m in newsel and c!=newsel[m]}
 selection={'status':'completed','selected_cellhint':selected,'selection_rule':'highest mean unified applicable-metric 40/60 across K=10,20,30; config id breaks exact ties',
            'pilot_ranking_sha256':fh(O/'pilot/unified/unified_scib_configuration_ranking.tsv'),
            'previous_pilot_selected_configurations':oldsel,'all_method_unified_pilot_selected':newsel,'changed_previously_selected_configurations':changes}
 (O/'pilot/CELLHINT_PILOT_SELECTION.json').write_text(json.dumps(selection,indent=2))
 wait(O/'primary/shared_harmonization/SHARED_HARMONIZATION_COMPLETE.json','waiting_primary_shared_harmonization')
 status('integrating_selected_primary_cellhint',selected=selected)
 run([OV,str(C/'integrate_cellhint_primary_selected.py')],O/'logs'/f'primary_{selected}_integrate.log','omicverse')
 status('evaluating_primary_cellhint',selected=selected)
 metric=evaluate('primary',O/'primary/embeddings'/f'cellhint__{selected}.json',selected)
 combine_score('primary',[metric])
 ranking_status=json.loads((O/'primary/unified/UNIFIED_SCIB_RANKING_STATUS.json').read_text())
 status('nine_method_scores_completed' if not changes else 'nine_method_scores_completed_parameter_sensitivity_requires_review',
        nine_method_metric_evaluation_complete=ranking_status['formal_ranking_complete'],changed_previously_selected_configurations=changes,
        parameter_policy='Previously pilot-frozen settings retained; unified pilot re-selection differences explicitly recorded',
        result=str(O/'primary/unified/unified_scib_configuration_ranking.tsv'))
if __name__=='__main__':
 try:main()
 except Exception as e:
  status('failed_requires_attention',error=repr(e));raise
