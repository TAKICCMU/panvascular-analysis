"""Atlas external-method adapters. Pure functions; no import-time data access/jobs."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy import sparse
MISSING = {'', 'NR', 'NA', 'nan', 'None', 'not_available', 'Unknown', 'Not_assessed'}

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024**2), b''): h.update(block)
    return h.hexdigest()

def qualified_states(obs, l1='level1_annotation', l2='level2_annotation'):
    """State identity is parent+child, never merge similarly named states across L1."""
    a=obs[l1].astype(str); b=obs[l2].astype(str)
    return (a+' | '+b).where(~a.isin(MISSING)&~b.isin(MISSING))

def validate_sample_metadata(obs, sample='sample_uid', columns=()):
    if obs.index.has_duplicates: raise ValueError('Duplicate cell IDs')
    if obs[sample].isna().any() or obs[sample].astype(str).isin(MISSING).any():
        raise ValueError('Missing sample identity')
    for c in columns:
        bad=obs.groupby(sample,observed=True)[c].nunique(dropna=False)
        if (bad>1).any(): raise ValueError(f'{c} is not constant within {sample}: {bad[bad>1].index.tolist()}')
    return obs.groupby(sample,observed=True).first()[list(columns)]

def communication_eligibility(obs, biological_sample='biological_sample_uid',
                              level='level1_annotation', min_cells=50):
    """Zero count and not-estimable pair remain explicit; no pooling across people."""
    validate_sample_metadata(obs, biological_sample)
    counts=pd.crosstab(obs[biological_sample],obs[level])
    counts=counts.loc[:,~counts.columns.astype(str).isin(MISSING)]
    rows=[]
    for sample, vals in counts.iterrows():
        for sender in counts:
            for receiver in counts:
                ns,nr=int(vals[sender]),int(vals[receiver])
                rows.append((sample,sender,receiver,ns,nr,ns>=min_cells and nr>=min_cells))
    return pd.DataFrame(rows,columns=[biological_sample,'sender','receiver','sender_cells','receiver_cells','eligible'])

def specificity_from_counts(X, obs, genes, labels, study_key='study', chunk=20000,
                            min_cells=500,min_studies=3,min_study_cells=20):
    """Exact inherited CP1e6/log1p -> group means -> gene specificity; no atlas mutation.
    Study-balanced score averages group-specific study means with >=20 cells.
    Ineligible states are reported, excluded identically from both specificity variants.
    Unmapped output genes are audited and omitted only AFTER full-library normalization.
    Mapped Ensembl features must be unique; unknown mapping is never guessed.
    """
    genes=pd.Index(genes)
    tokens=pd.Series(genes).fillna('').astype(str)
    mapped=tokens.str.fullmatch(r'ENSG\d+').to_numpy()
    if tokens[mapped].duplicated().any():raise ValueError('Duplicate valid Ensembl mapping requires explicit resolution')
    if not mapped.any():raise ValueError('No mapped Ensembl genes')
    gene_audit=pd.DataFrame({'input_feature_position':np.arange(len(genes)),'ensembl_id':tokens,
       'retained_in_specificity_output':mapped,'reason':np.where(mapped,'mapped','missing_or_invalid_Ensembl_not_in_MAGMA_output')})
    if X.shape!=(len(obs),len(genes)): raise ValueError('Matrix/metadata shape mismatch')
    labels=pd.Series(np.asarray(labels),index=obs.index)
    valid=labels.notna()&~labels.astype(str).isin(MISSING)
    grouped=pd.DataFrame({'group':labels,'study':obs[study_key]},index=obs.index).loc[valid]
    if grouped.study.isna().any() or grouped.study.astype(str).isin(MISSING).any():
        raise ValueError('Missing study identity')
    counts=grouped.groupby(['group','study'],observed=True).size().unstack(fill_value=0)
    audit=pd.DataFrame({'n_cells':counts.sum(axis=1),'n_studies':(counts>=min_study_cells).sum(axis=1)})
    audit['eligible']=(audit.n_cells>=min_cells)&(audit.n_studies>=min_studies)
    groups=audit.index[audit.eligible].tolist(); studies=sorted(grouped.study.unique())
    if not groups: raise ValueError('No eligible states under inherited thresholds')
    gi=labels.map(dict(zip(groups,range(len(groups))))).fillna(-1).to_numpy(int)
    si=obs[study_key].map(dict(zip(studies,range(len(studies))))).fillna(-1).to_numpy(int)
    gs=np.zeros((len(groups)*len(studies),len(genes)),np.float64)
    nc=np.zeros(len(gs),np.int64)
    for start in range(0,len(obs),chunk):
        stop=min(start+chunk,len(obs)); keep=gi[start:stop]>=0
        if not keep.any(): continue
        b=sparse.csr_matrix(X[start:stop])[keep].astype(float)
        if not np.isfinite(b.data).all() or (b.data<0).any() or not np.equal(b.data,np.floor(b.data)).all():
            raise ValueError('X must contain nonnegative integer raw counts')
        totals=np.asarray(b.sum(axis=1)).ravel()
        if (totals<=0).any(): raise ValueError('Zero-library cell cannot enter specificity')
        b=sparse.diags(1e6/totals)@b; b.data=np.log1p(b.data)
        code=gi[start:stop][keep]*len(studies)+si[start:stop][keep]
        d=sparse.csr_matrix((np.ones(len(code)),(code,np.arange(len(code)))),shape=(len(gs),len(code)))
        rows=np.unique(code); gs[rows]+=(d@b)[rows].toarray(); nc+=np.bincount(code,minlength=len(gs))
    gs=gs.reshape(len(groups),len(studies),len(genes)); nc=nc.reshape(len(groups),len(studies))
    cell=gs.sum(axis=1)/nc.sum(axis=1)[:,None]
    study=np.stack([(gs[i,nc[i]>=min_study_cells]/nc[i,nc[i]>=min_study_cells,None]).mean(axis=0) for i in range(len(groups))])
    mapping=pd.DataFrame({'variable':[f'S{i:04d}' for i in range(len(groups))],'group':groups})
    def finish(m):
        den=m.sum(axis=0); k=(den>0)&mapped
        out=pd.DataFrame((m[:,k]/den[k]).T,columns=mapping.variable)
        out.insert(0,'Average',m[:,k].mean(axis=0)); out.insert(0,'GENE',genes[k]); out.attrs['gene_mapping_audit']=gene_audit; return out
    return finish(cell),finish(study),mapping,audit.reset_index()

def prepare_scdrs_counts(X,obs,gene_map):
    """Aggregate count columns with same reviewed gene symbol BEFORE scDRS filtering.
    gene_map columns feature_id/symbol must follow exact X feature order.
    """
    if len(gene_map)!=X.shape[1] or gene_map.feature_id.duplicated().any(): raise ValueError('Gene map identity invalid')
    symbols=gene_map['symbol'] if 'symbol' in gene_map else gene_map['gene_symbol']
    if symbols.isna().any() or symbols.astype(str).isin(MISSING).any(): raise ValueError('Unresolved symbol mapping')
    code,names=pd.factorize(symbols,sort=True)
    d=sparse.csr_matrix((np.ones(len(code)),(np.arange(len(code)),code)),shape=(len(code),len(names)))
    out=sparse.csr_matrix(X)@d
    return out,pd.Index(names,name='gene_symbol')

def scdrs_commands(executable,h5ad,cov,geneset,out):
    """Native scoring and native group tests, not mean-score substitute p values."""
    out=Path(out)
    return [[str(executable),'compute-score','--h5ad-file',str(h5ad),'--h5ad-species','human',
             '--gs-file',str(geneset),'--gs-species','human','--cov-file',str(cov),
             '--ctrl-match-opt','mean_var','--weight-opt','vs','--flag-filter-data','True',
             '--flag-raw-count','True','--n-ctrl','1000','--flag-return-ctrl-raw-score','False',
             '--flag-return-ctrl-norm-score','True','--out-folder',str(out/'scores')]]

def magma_command(executable,gene_results,specificity,out):
    return list(map(str,[executable,'--gene-results',gene_results,'--gene-covar',specificity,
                        '--model','direction=pos','condition=Average','--out',out]))

def audit_resources(rows,hash_existing=False):
    out=[]
    for row in rows:
        r=dict(row); p=Path(r['path']) if r.get('path') else None
        r['exists']=bool(p and p.is_file()); r['observed_sha256']=''
        if r['exists'] and hash_existing: r['observed_sha256']=sha256(p)
        r['status']='BLOCKED_MISSING_RESOURCE' if not r['exists'] else 'PRESENT_UNVERIFIED'
        if r['observed_sha256'] and r.get('expected_sha256'):
            r['status']='VERIFIED' if r['observed_sha256']==r['expected_sha256'] else 'BLOCKED_HASH_MISMATCH'
        out.append(r)
    return out

def target_evidence_join(tables):
    """Outer-join distinct evidence layers. Missing evidence is NA, never negative/zero.
    No hard-coded target list, old significance, or new causal tier is manufactured.
    """
    result=None
    for name,frame in tables.items():
        if 'gene_id' not in frame or frame.gene_id.duplicated().any():
            raise ValueError(f'{name}: require one reviewed summary per stable gene_id')
        x=frame.rename(columns={c:f'{name}__{c}' for c in frame if c!='gene_id'})
        result=x if result is None else result.merge(x,on='gene_id',how='outer',validate='one_to_one')
    if result is None: raise ValueError('No evidence tables')
    return result

def render_spatial_adapter(source,atlas_path,output_root,manifest_path,l1_key,l2_key,study_key,
                           approved_forced_genes,atlas_sha256,config_path=None,training_config=None):
    if __package__:
        from .spatial_renderer import render
    else:
        from spatial_renderer import render
    return render(source,atlas_path,output_root,manifest_path,l1_key,l2_key,study_key,
                  approved_forced_genes,atlas_sha256,config_path,training_config)

def execute(config,run_dir):
    """Main-pipeline design-only adapter: writes resource audit, never scientific scores."""
    out=Path(run_dir)/'external_preflight';out.mkdir(parents=True,exist_ok=False)
    inventory=Path(__file__).resolve().parents[1]/'external_resource_inventory.tsv'
    rows=pd.read_csv(inventory,sep='\t',keep_default_na=False).to_dict('records')
    audit=audit_resources(rows,hash_existing=False)
    pd.DataFrame(audit).to_csv(out/'resources.tsv',sep='\t',index=False)
    result={'status':'WAITING_FOR_FINAL_L2_ATLAS','external_models_executed':False,'resource_count':len(audit)}
    (out/'STATUS.json').write_text(json.dumps(result,indent=2));return result
