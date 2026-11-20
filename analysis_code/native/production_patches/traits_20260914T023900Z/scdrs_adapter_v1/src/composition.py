from pathlib import Path
import json,subprocess
import pandas as pd
from scipy.stats import false_discovery_control
from .core import Blocked,atomic_json

def execute(config,run_dir):
    run=Path(run_dir);out=run/'composition_models';out.mkdir(exist_ok=False)
    specs=json.loads(Path(config['composition']['models_file']).read_text())
    if not specs:raise Blocked('No approved composition model specifications')
    counts=pd.read_csv(run/'counts/biological_state_counts.tsv',sep='\t')
    cohort=pd.read_csv(run/'cohort/cohort_design.tsv',sep='\t').drop_duplicates('biological_sample_uid')
    results=[];audit=[]
    for spec in specs:
        mid=spec['model_id'];record={'model_id':mid,'status':'NOT_SELECTED' if spec.get('approved') is not True else 'READY'}
        try:
            if spec.get('approved') is not True:
                record['reason']='Composition model pending approval';audit.append(record);continue
            if not mid.replace('_','').isalnum():raise Blocked('Unsafe model ID')
            if spec['level']=='L1':
                tab=counts[counts.level1_annotation==spec['state']].groupby('biological_sample_uid',as_index=False).agg(n_state=('n_cells','sum'),n_parent=('n_all_qc','first'))
            else:
                tab=counts[counts.level2_annotation==spec['state']][['biological_sample_uid','n_cells','n_parent']].rename(columns={'n_cells':'n_state'})
            tab=tab.merge(cohort,on='biological_sample_uid',validate='one_to_one')
            for key,values in spec.get('filters',{}).items():tab=tab[tab[key].isin(values)]
            tab=tab[tab.inference_eligible.astype(str).str.lower().isin(['true','1'])]
            if tab.empty:raise Blocked('No eligible biological units')
            inp=out/(mid+'.tsv');js=out/(mid+'.json');tab.to_csv(inp,sep='\t',index=False);atomic_json(js,spec)
            cmd=[config['composition']['rscript'],str(Path(__file__).resolve().parents[1]/'scripts/composition_glmm.R'),str(inp),str(js),str(out/mid)]
            with open(out/(mid+'.log'),'w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
            record.update(exit_code=r.returncode,log=str(out/(mid+'.log')))
            if r.returncode:
                reason=(out/(mid+'.log')).read_text(errors='replace').strip()
                known=['Insufficient independent units','Confounded fixed-effect design','Missing covariates','Repeated donor requires','No two-sex study support','Random study effect needs','age_z requires explicitly approved','Constant age','Both optimizers failed','Declared tested terms absent','Unknown donor not an independent person','Duplicate biological unit','Invalid integer denominator']
                record.update(status='BLOCKED' if any(x in reason for x in known) else 'FAILED',reason=reason or f'GLMM exit {r.returncode} without diagnostic',failure_kind='NATIVE_MODEL_NOT_ESTIMABLE' if any(x in reason for x in known) else 'NATIVE_PROGRAM_FAILURE');audit.append(record);continue
            results.append(pd.read_csv(out/mid/'effects.tsv',sep='\t'));record.update(status='COMPLETE')
        except Blocked as e:record.update(status='BLOCKED',reason=str(e),failure_kind='MODEL_GATE')
        except Exception as e:record.update(status='FAILED',reason=str(e),exception_type=type(e).__name__,failure_kind='PROGRAM_EXCEPTION')
        audit.append(record)
    if results:
        d=pd.concat(results,ignore_index=True);d['FDR_BH']=d.groupby('FDR_family').P.transform(lambda x:false_discovery_control(x.to_numpy()));d.to_csv(out/'all_effects.tsv',sep='\t',index=False)
    atomic_json(out/'MODEL_AUDIT.json',audit)
    if not results:pd.DataFrame(columns=['model_id','term','estimate','SE','P','FDR_family','FDR_BH']).to_csv(out/'all_effects.tsv',sep='\t',index=False)
    active=[x['status'] for x in audit if x['status']!='NOT_SELECTED']
    status='FAILED' if 'FAILED' in active else ('COMPLETE' if active and all(x=='COMPLETE' for x in active) else ('PARTIAL' if 'COMPLETE' in active else 'BLOCKED'))
    return {'status':status,'models':audit}
