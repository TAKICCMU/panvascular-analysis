"""Explicit model tasks; context meta-analysis does not invent coefficient SEs."""
from pathlib import Path
import json, subprocess, os, gzip
from scipy import sparse,io
import numpy as np
import pandas as pd
from .core import Blocked,atomic_json
from .cohort import select_contrast

def prepare_de_tasks(config,pb_dir,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    contrasts=pd.read_csv(config['statistical_units']['contrasts_file'],sep='\t',keep_default_na=False)
    tasks=[];audit=[]
    for level in ['L1','L2']:
        base=Path(pb_dir)/level;units=pd.read_csv(base/'units.tsv',sep='\t')
        matrix=sparse.load_npz(base/'counts_units_by_genes.npz')
        row_lookup=pd.Series(range(len(units)),index=units.pseudobulk_id)
        for state,g in units.groupby('cell_state',sort=True):
            g=g[g.n_cells>=config['pseudobulk']['min_cells']].copy()
            for _,spec in contrasts.iterrows():
                identity={'level':level,'cell_state':state,'contrast_id':spec.contrast_id}
                try:
                    m,diag=select_contrast(g,spec)
                    formula=str(spec.formula)
                    cov=[v for v in str(spec.get('covariates','')).split(';') if v and v!='nan']
                    expected='~ '+ ' + '.join((['donor_uid'] if diag['paired'] else [])+cov+['group'])
                    if formula.replace(' ','')!=expected.replace(' ',''):raise Blocked(f'Formula/covariate declaration conflict: expected {expected}')
                    data={**identity,'reference':spec.reference,'case':spec.case,'formula':formula,'meta_context':spec.meta_context,'min_tested_genes':config['pseudobulk']['min_tested_genes'],'rows':json.loads(m.to_json(orient='records')),'diagnostics':diag}
                    task_id=f'de_{len(tasks):06d}';spec_path=out/(task_id+'.json');atomic_json(spec_path,data)
                    inputs=out/(task_id+'_inputs');inputs.mkdir()
                    small=matrix[row_lookup.loc[m.pseudobulk_id].to_numpy()].T
                    with gzip.open(inputs/'counts.mtx.gz','wb') as f:io.mmwrite(f,small,field='integer')
                    m.to_csv(inputs/'units.tsv',sep='\t',index=False)
                    cmd=[config['pseudobulk']['rscript'],str(Path(__file__).resolve().parents[1]/'scripts/pseudobulk_edgeR.R'),str(inputs/'counts.mtx.gz'),str(inputs/'units.tsv'),str(base/'genes.tsv'),str(spec_path),str(out/task_id)]
                    tasks.append({'id':task_id,'argv':cmd,**identity});audit.append({**identity,'status':'READY_FOR_APPROVED_EXECUTION','reason':''})
                except Blocked as e:audit.append({**identity,'status':'NOT_SELECTED' if str(spec.get('approved')).lower() not in ['true','1'] else 'BLOCKED_MODEL','reason':str(e)})
    pd.DataFrame(audit,columns=['level','cell_state','contrast_id','status','reason']).to_csv(out/'MODEL_ELIGIBILITY.tsv',sep='\t',index=False)
    atomic_json(out/'TASKS.json',tasks)
    return tasks

def execute_de(config,run_dir):
    tasks=prepare_de_tasks(config,Path(run_dir)/'pseudobulk',Path(run_dir)/'de')
    report=[]
    for task in tasks:
        log=Path(run_dir)/'de'/(task['id']+'.log');record={**task,'log':str(log)}
        try:
            with open(log,'w') as f:r=subprocess.run(task['argv'],stdout=f,stderr=subprocess.STDOUT,env=os.environ.copy(),check=False)
            record['exit_code']=r.returncode
            if r.returncode:
                reason=log.read_text(errors='replace').strip()
                known=['Non-identifiable design','Insufficient residual degrees','Insufficient expressed genes','Missing covariates silently','Explicit simple case-reference','Not integer pseudobulk counts']
                expected=any(x in reason for x in known)
                record.update(status='BLOCKED' if expected else 'FAILED',reason=reason or f'edgeR exit {r.returncode} without diagnostic',failure_kind='NATIVE_MODEL_NOT_ESTIMABLE' if expected else 'NATIVE_PROGRAM_FAILURE')
            else:record['status']='COMPLETE'
        except Exception as e:record.update(status='FAILED',reason=str(e),exception_type=type(e).__name__,failure_kind='PROGRAM_EXCEPTION')
        report.append(record)
    atomic_json(Path(run_dir)/'de/EXECUTION.json',report)
    eligibility=pd.read_csv(Path(run_dir)/'de/MODEL_ELIGIBILITY.tsv',sep='\t').fillna('')
    models=[dict(x,status='BLOCKED' if x['status']=='BLOCKED_MODEL' else x['status']) for x in eligibility.to_dict('records') if x['status']!='READY_FOR_APPROVED_EXECUTION']+report
    atomic_json(Path(run_dir)/'de/MODEL_AUDIT.json',models)
    active=[x['status'] for x in models if x['status']!='NOT_SELECTED']
    status='FAILED' if 'FAILED' in active else ('COMPLETE' if active and all(x=='COMPLETE' for x in active) else ('PARTIAL' if 'COMPLETE' in active else 'BLOCKED'))
    return {'status':status,'tasks':report,'models':models}

def meta_input(config,run_dir):
    """Preserve historical approximation only when explicitly approved; separate resource-dependent step."""
    out=Path(run_dir)/'meta';out.mkdir(exist_ok=False)
    files=sorted((Path(run_dir)/'de').glob('de_*/DE.tsv'))
    if not files:raise Blocked('No completed DE inputs')
    d=pd.concat([pd.read_csv(p,sep='\t') for p in files],ignore_index=True)
    d.to_csv(out/'all_DE.tsv.gz',sep='\t',index=False)
    if not config['pseudobulk'].get('approve_legacy_ql_se_approximation',False):raise Blocked('Meta SE approximation needs review; exact QL-F is not a Wald coefficient SE')
    if config['pseudobulk'].get('context_synthesis_strategy','one_approved_model_per_context') != 'one_approved_model_per_context':
        raise Blocked('Within-context synthesis strategy not implemented/approved; use one approved joint model per context')
    # One result per state/context; overlapping donor sets or multiple tissue contrasts require an approved synthesis plan.
    keys=['level','cell_state','gene_id','meta_context']
    if d.duplicated(keys).any():raise Blocked('Multiple effects per context require explicit within-context synthesis')
    donors=[]
    for p in files:
        m=pd.read_csv(p.parent/'included_units.tsv',sep='\t');t=pd.read_csv(p,sep='\t',nrows=1)
        donors.extend([(t.level.iloc[0],t.cell_state.iloc[0],t.meta_context.iloc[0],v) for v in set(m.donor_uid)])
    audit=pd.DataFrame(donors,columns=['level','state','context','donor'])
    if (audit.groupby(['level','state','donor']).context.nunique()>1).any():raise Blocked('Cross-context repeated donors: independence assumption violated')
    resource=config['pseudobulk'].get('gwas_priority_table')
    expected=config['pseudobulk'].get('gwas_priority_sha256')
    from .core import sha
    if not resource or not Path(resource).exists() or not expected or sha(resource)!=expected:raise Blocked('Validated GWAS priority table missing')
    master=pd.read_csv(resource,sep='\t')
    if not {'SYMBOL','primary_FDR','primary_Bonferroni'}.issubset(master):raise Blocked('Wrong GWAS priority schema')
    symbol=config['pseudobulk']['gene_symbol_column']
    if symbol not in d:raise Blocked(f'Gene symbol map missing: {symbol}')
    keep=master.loc[master.primary_FDR.astype(str).str.lower().isin(['true','1']) | master.primary_Bonferroni.astype(str).str.lower().isin(['true','1']),'SYMBOL']
    d[d[symbol].isin(keep)].to_csv(out/'priority_meta_input.tsv',sep='\t',index=False)
    parameters = meta_parameters(config)
    atomic_json(out/'EFFECTIVE_MODEL_PARAMETERS.json',parameters)
    cmd=[config['pseudobulk']['rscript'],str(Path(__file__).resolve().parents[1]/'scripts/cross_context_meta.R'),str(out/'priority_meta_input.tsv'),str(out/'meta_results.tsv'),str(out/'EFFECTIVE_MODEL_PARAMETERS.json')]
    with open(out/'meta.log','w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=False)
    if r.returncode:raise Blocked('Meta-analysis failed; see meta.log')
    fitted=pd.read_csv(out/'meta_results.tsv',sep='\t')
    if fitted.empty or not fitted.status.eq('COMPLETE').any():raise Blocked('No gene/state has an estimable multi-context meta; complete NOT_TESTED table retained')
    return {'status':'COMPLETE','argv':cmd,'method':'REML','SE_kind':'APPROXIMATE_FROM_QL_F_NOT_COEFFICIENT_SE','parameters':parameters}


def meta_parameters(config):
    p=config['pseudobulk']
    n=int(p['meta_min_contexts']);strict=int(p['strict_shared_min_contexts']);direction=float(p['direction_consistency'])
    if n<2 or strict<n or not 0.5<=direction<=1:raise Blocked('Invalid explicit meta thresholds')
    if p.get('approve_legacy_ql_se_approximation') is not True:raise Blocked('Legacy QL-F approximation remains unapproved')
    return {'meta_min_contexts':n,'strict_shared_min_contexts':strict,'direction_consistency':direction,'fdr_threshold':float(p.get('meta_fdr_threshold',0.05)),'method':'REML','SE_kind':'APPROXIMATE_FROM_QL_F_NOT_COEFFICIENT_SE','approve_legacy_ql_se_approximation':True,'context_synthesis_strategy':'one_approved_model_per_context','FDR_family':'all tested priority genes within level and state','effect_scale':'log2_fold_change'}


def contrast_feasibility(cohort,contrasts):
    """Read-only upper bound. Never changes approvals or bypasses execution gates."""
    from .core import missing
    rows=[]
    for _,spec in contrasts.iterrows():
        m=cohort.copy()
        for col in ['study','vascular_bed','cohort_disease_context','assay_type']:
            if not missing(spec.get(col)):m=m[m[col].astype(str).eq(str(spec[col]))]
        if not missing(spec.get('study_allowlist','')):m=m[m.study.isin(json.loads(spec.study_allowlist))]
        m=m[m.comparison_group.isin([spec.case,spec.reference])].copy()
        resolved=~m.donor_uid.map(missing)&~m.donor_uid.astype(str).str.contains('PROVISIONAL|UNRESOLVED|PLACEHOLDER',case=False,regex=True)
        m=m[resolved];case=set(m.loc[m.comparison_group.eq(spec.case),'donor_uid']);ref=set(m.loc[m.comparison_group.eq(spec.reference),'donor_uid'])
        paired=str(spec.get('paired')).lower() in ['true','1'];minimum=int(spec.min_independent_per_arm)
        donor_pass=len(case)>=minimum and len(ref)>=minimum and (not paired or len(case&ref)>=minimum)
        rows.append({'contrast_id':spec.contrast_id,'study':spec.study,'meta_context':spec.meta_context,'n_case_donors_recomputed':len(case),'n_reference_donors_recomputed':len(ref),'n_complete_pairs_recomputed':len(case&ref),'both_arms_nonempty':bool(case and ref),'paired':paired,'minimum_per_arm':minimum,'passes_donor_upper_bound':donor_pass,'proposed_grade':'PRIMARY_CANDIDATE_PENDING' if donor_pass else ('EXPLORATORY_ONLY_PENDING' if case and ref else 'NOT_TESTABLE_NO_REFERENCE_OR_CASE'),'formula':spec.formula,'covariates':spec.get('covariates',''),'reference_meaning':spec.get('reference_interpretation','REQUIRES_REVIEW'),'approved':str(spec.approved).lower() in ['true','1'],'execution_status':'BLOCKED_PENDING_DESIGN_AND_STATE_COVERAGE','limitations':'Upper bound before approved technical merges, state min cells, complete-case covariates and matrix rank; not execution eligibility'})
    return pd.DataFrame(rows)
