"""Targeted scopes refit the same approved downstream models; no re-integration."""
from pathlib import Path
import json,hashlib,shutil
import pandas as pd
from scipy import sparse
from .core import Blocked,atomic_json


def scoped_tables(counts,cohort,rule):
    m=cohort.copy()
    if rule['kind']=='assay':m=m[m[rule['column']].isin(rule['values'])]
    elif rule['kind']=='exclude_study':m=m[~m.study.isin(rule['values'])]
    else:raise Blocked('Unknown targeted sensitivity rule')
    if m.empty:raise Blocked('No samples in requested sensitivity')
    t=counts[counts.sample_uid.isin(m.sample_uid)].copy()
    changes=cohort.groupby(['study','vascular_bed','assay_type'],dropna=False).size().rename('before').reset_index()
    after=m.groupby(['study','vascular_bed','assay_type'],dropna=False).size().rename('after').reset_index()
    return t,m,changes.merge(after,how='left').fillna({'after':0})


def complete_biosample_scope(cohort,selected):
    """Refuse partial technical-library inclusion; reused PB units must be exact."""
    chosen=set(selected.sample_uid.astype(str));bios=[]
    for bio,g in cohort.groupby('biological_sample_uid',dropna=False):
        ids=set(g.sample_uid.astype(str));hit=ids&chosen
        if hit and hit!=ids:raise Blocked(f'Sensitivity splits an aggregated biological sample: {bio}')
        if hit:bios.append(bio)
    return bios


def filter_pseudobulk(source,dest,bios):
    dest=Path(dest);dest.mkdir()
    for level in ['L1','L2']:
        src=Path(source)/level;out=dest/level;out.mkdir()
        u=pd.read_csv(src/'units.tsv',sep='\t');keep=u.biological_sample_uid.isin(bios).to_numpy()
        x=sparse.load_npz(src/'counts_units_by_genes.npz')
        if x.shape[0]!=len(u):raise Blocked('Pseudobulk row identity mismatch')
        sparse.save_npz(out/'counts_units_by_genes.npz',x[keep]);u.loc[keep].to_csv(out/'units.tsv',sep='\t',index=False)
        shutil.copy2(src/'genes.tsv',out/'genes.tsv')


def compare_effects(primary,scoped,keys,effect,se=None):
    """Outer join preserves unavailable fits; no selection by significance."""
    if primary.duplicated(keys).any() or scoped.duplicated(keys).any():raise Blocked('Ambiguous sensitivity effect keys')
    cols=keys+[effect]+([se] if se and se in primary and se in scoped else [])
    d=primary[cols].merge(scoped[cols],on=keys,how='outer',suffixes=('_primary','_sensitivity'),indicator=True)
    d['comparison_status']=d['_merge'].map({'both':'COMPARABLE','left_only':'NOT_ESTIMABLE_IN_SCOPE','right_only':'NO_PRIMARY_EFFECT'})
    d['effect_difference']=d[effect+'_sensitivity']-d[effect+'_primary']
    d['same_direction']=pd.NA
    ok=d['_merge'].eq('both');d.loc[ok,'same_direction']=(d.loc[ok,effect+'_sensitivity']*d.loc[ok,effect+'_primary']>0).to_numpy()
    return d.drop(columns='_merge')


def _load_effects(root,module):
    root=Path(root)
    if module=='de':
        paths=sorted((root/'de').glob('de_*/DE.tsv'));keys=['level','cell_state','gene_id','contrast_id'];effect='logFC'
    else:
        p=root/'composition_models/all_effects.tsv';paths=[p] if p.exists() else [];keys=['model_id','term'];effect='estimate'
    tables=[pd.read_csv(p,sep='\t') for p in paths]
    d=pd.concat(tables,ignore_index=True) if tables else pd.DataFrame(columns=keys+[effect])
    if not set(keys+[effect])<=set(d):raise ValueError('Effect schema mismatch for '+module)
    return d,keys,effect,bool(paths and len(d))


def _model_audit(root,module,fallback=None):
    p=Path(root)/('de' if module=='de' else 'composition_models')/'MODEL_AUDIT.json'
    if p.exists():return json.loads(p.read_text())
    return (fallback or {}).get('models',[])


def _annotate_comparison(config,primary,scope,module,result):
    a,keys,effect,has_primary=_load_effects(primary,module);b,_,_,_= _load_effects(scope,module)
    d=compare_effects(a,b,keys,effect)
    idcols=['level','cell_state','contrast_id'] if module=='de' else ['model_id']
    pa=_model_audit(primary,module);sa=_model_audit(scope,module,result)
    def lookup(audit):return {tuple(str(x.get(k,'')) for k in idcols):x for x in audit}
    pl,sl=lookup(pa),lookup(sa)
    if module=='composition' and d.empty:
        specs=json.loads(Path(config['composition']['models_file']).read_text())
        d=pd.DataFrame([{'model_id':x['model_id'],'term':t,'comparison_status':'NO_PRIMARY_OR_SCOPE_EFFECT','estimate_primary':float('nan'),'estimate_sensitivity':float('nan'),'effect_difference':float('nan'),'same_direction':pd.NA} for x in specs if x.get('approved') is True for t in x.get('tested_terms',[])])
    rows=[]
    for row in d.to_dict('records'):
        identity=tuple(str(row.get(k,'')) for k in idcols);pr=pl.get(identity,{});sr=sl.get(identity,{})
        missing_primary=pd.isna(row.get(effect+'_primary'));missing_scope=pd.isna(row.get(effect+'_sensitivity'))
        row['primary_model_status']=pr.get('status','NO_PRIMARY_EFFECT' if missing_primary else 'COMPLETE')
        row['scope_model_status']=sr.get('status',result.get('status','UNKNOWN') if missing_scope else 'COMPLETE')
        row['primary_reason']=pr.get('reason','MISSING_PRIMARY_RESULTS' if not has_primary else 'FEATURE_NOT_TESTED_IN_PRIMARY') if missing_primary else ''
        row['scope_reason']=sr.get('reason',result.get('reason','FEATURE_NOT_TESTED_IN_SCOPE' if row['scope_model_status']=='COMPLETE' else 'MODEL_AUDIT_REQUIRED')) if missing_scope else ''
        row['scope_failure_kind']=sr.get('failure_kind',result.get('failure_kind','')) if missing_scope else ''
        rows.append(row)
    final=pd.DataFrame(rows) if rows else d.assign(primary_model_status=pd.Series(dtype=str),scope_model_status=pd.Series(dtype=str),primary_reason=pd.Series(dtype=str),scope_reason=pd.Series(dtype=str),scope_failure_kind=pd.Series(dtype=str))
    final.to_csv(Path(scope)/(module+'_primary_vs_scope.tsv'),sep='\t',index=False)
    return has_primary,len(final)


def _overall(statuses):
    active=[x for x in statuses if x!='NOT_SELECTED']
    if 'FAILED' in active:return 'FAILED'
    if active and all(x=='COMPLETE' for x in active):return 'COMPLETE'
    if 'COMPLETE' in active or 'PARTIAL' in active:return 'PARTIAL'
    return 'BLOCKED'


def run_same_models(config,primary,scope,bios):
    """Preserve successful effects and exact failure reasons without passing partials."""
    from . import inference,composition,covarnet
    modules=config['sensitivity'].get('model_modules',[])
    if not modules or not set(modules)<= {'de','composition','covarnet'}:raise Blocked('Explicit approved sensitivity model_modules (de/composition/covarnet) required')
    report=[]
    for module in modules:
        result={}
        try:
            if module=='de':
                filter_pseudobulk(Path(primary)/'pseudobulk',Path(scope)/'pseudobulk',bios)
                result=inference.execute_de(config,scope)
            elif module=='composition':
                source=pd.read_csv(Path(primary)/'counts/biological_state_counts.tsv',sep='\t')
                source[source.biological_sample_uid.isin(bios)].to_csv(Path(scope)/'counts/biological_state_counts.tsv',sep='\t',index=False)
                result=composition.execute(config,scope)
            else:
                ontology=Path(config.get('covarnet',{}).get('ontology_path') or Path(primary)/'counts/ontology.tsv')
                shutil.copy2(ontology,Path(scope)/'counts/ontology.tsv');result=covarnet.execute(config,scope)
                if result.get('status')!='COMPLETE':raise Blocked('CoVarNet scope only prepared, not fitted')
                compare_nmf_outputs(primary,scope).to_csv(Path(scope)/(module+'_primary_vs_scope.tsv'),sep='\t',index=False)
        except Blocked as e:result={'status':'BLOCKED','reason':str(e),'failure_kind':'MODEL_GATE'}
        except Exception as e:result={'status':'FAILED','reason':str(e),'exception_type':type(e).__name__,'failure_kind':'PROGRAM_EXCEPTION'}
        record={'module':module,'status':result.get('status','FAILED'),'same_config':True,'executor_receipt':result}
        if module in {'de','composition'}:
            try:
                has_primary,n=_annotate_comparison(config,primary,scope,module,result)
                record.update(comparison_rows=n,comparison_generated=True)
                if not has_primary:
                    record['comparison_status']='BLOCKED_MISSING_PRIMARY'
                    if record['status']!='FAILED':record['status']='BLOCKED'
                    record['reason']='No completed primary effects; scope-only results retained'
            except Exception as e:
                record.update(status='FAILED',comparison_generated=False,reason=str(e),exception_type=type(e).__name__,failure_kind='COMPARISON_PROGRAM_EXCEPTION')
        report.append(record)
    atomic_json(Path(scope)/'MODEL_EXECUTION.json',report)
    return report

def execute(config,run_dir):
    p=Path(run_dir);out=p/'sensitivity';out.mkdir(exist_ok=False)
    cfg=config['sensitivity']
    if cfg.get('run') is not True:raise Blocked('Targeted sensitivity not approved')
    t=pd.read_csv(p/'counts/sample_state_counts.tsv',sep='\t');m=pd.read_csv(p/'cohort/cohort_design.tsv',sep='\t')
    rules=[{'id':'scRNA_only','kind':'assay','column':cfg['assay_column'],'values':cfg['scRNA_only_values']},{'id':'exclude_target_study','kind':'exclude_study','values':cfg['exclude_studies']}]
    seen={};reports=[]
    for r in rules:
        count,meta,changes=scoped_tables(t,m,r);d=out/r['id'];d.mkdir();(d/'counts').mkdir();(d/'cohort').mkdir()
        signature=tuple(sorted(meta.sample_uid.astype(str)));equivalent=seen.get(signature);seen.setdefault(signature,r['id'])
        count.to_csv(d/'counts/sample_state_counts.tsv',sep='\t',index=False);meta.to_csv(d/'cohort/cohort_design.tsv',sep='\t',index=False);changes.to_csv(d/'joint_study_bed_assay_coverage.tsv',sep='\t',index=False)
        receipt={'rule':r,'identical_sample_scope_to':equivalent,'sample_uid_sha256':hashlib.sha256('\n'.join(signature).encode()).hexdigest(),'n_samples':len(meta),'interpretation':'Study removal may remove vascular bed with assay. Same approved formula and support threshold; no L1/L2 changes.'}
        try:
            bios=complete_biosample_scope(m,meta)
            if equivalent:
                original=next(x for x in reports if x['rule']['id']==equivalent)
                receipt.update(status=original['status'],phase='REUSED_IDENTICAL_SCOPE',results_from=equivalent)
                if original['status']!='COMPLETE':receipt['reason']='Equivalent source scope did not complete'
            else:
                fits=run_same_models(config,p,d,bios);receipt.update(status=_overall([x['status'] for x in fits]),models=fits)
        except Blocked as e:receipt.update(status='BLOCKED',reason=str(e))
        atomic_json(d/'SCOPE.json',receipt);reports.append(receipt)
    atomic_json(out/'EXECUTION.json',reports)
    return {'status':_overall([x['status'] for x in reports]),'scopes':reports,'models_refitted':cfg['model_modules']}


def compare_nmf_outputs(primary,scope):
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    def load(root):
        paths=list((Path(root)/'covarnet/stable_NMF').glob('module_node_weights_K*.tsv'))
        # Only primary locked weights; multiple files require an explicit primary manifest.
        if len(paths)!=1:raise Blocked('Expected one fitted primary CoVarNet weight table')
        w=pd.read_csv(paths[0],sep='\t',index_col=0)
        if w.index.duplicated().any():raise Blocked('Duplicate CoVarNet state IDs')
        return w
    a,b=load(primary),load(scope);common=a.index.intersection(b.index)
    if len(common)<2:raise Blocked('Too few shared state IDs for NMF sensitivity comparison')
    A=a.loc[common].to_numpy(float);B=b.loc[common].to_numpy(float)
    sim=A.T@B/np.maximum(np.linalg.norm(A,axis=0)[:,None]*np.linalg.norm(B,axis=0)[None,:],np.finfo(float).tiny)
    ma=np.argmax(sim,axis=1);ri,ci=linear_sum_assignment(-sim);one=dict(zip(ri,ci));rows=[]
    for i,name in enumerate(a.columns):
        coverage=float(a.loc[common,name].sum()/a[name].sum()) if a[name].sum()>0 else 0.
        rows.append({'reference_module':name,'many_to_one_module':b.columns[ma[i]],'many_to_one_cosine':sim[i,ma[i]],'many_to_one_duplicate':int(np.sum(ma==ma[i]))>1,'one_to_one_module':b.columns[one[i]] if i in one else None,'one_to_one_cosine':sim[i,one[i]] if i in one else np.nan,'reference_weight_coverage':coverage,'missing_state_ids':';'.join(a.index.difference(common).astype(str)),'comparison_scope':'SHARED_STATES_ONLY','evidence_status':'DESCRIPTIVE_FULL_COVERAGE' if coverage>=1-1e-12 else 'LIMITED_MISSING_WEIGHT'})
    return pd.DataFrame(rows)
