import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,subprocess,time
import pandas as pd,numpy as np
from scipy.stats import false_discovery_control
B=Path(__file__).resolve().parents[1];W=B.parent
S=W/'pseudobulk_manuscript_strategy_20260916_v1'
E=W/'remaining_downstream_20260916_v1/enrichment'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rec=json.loads((B/'RECEIPT.json').read_text())
for p,h in rec['bound_files'].items():assert sha(p)==h,p
assert json.loads((B/'tests/DL_PASS.json').read_text())['status']=='PASS'
os.environ['R_LIBS_USER']=str(E/'Rlib')
for folder in ['inputs','tables','outputs','logs']:(B/folder).mkdir(exist_ok=False)
frames=[];donors=[]
for p in sorted((S/'results').glob('*/DE.tsv.gz')):
    d=pd.read_csv(p,sep='\t');d['source_model']=p.parent.name;frames.append(d)
    m=pd.read_csv(p.with_name('included_units.tsv'),sep='\t');donors.append(m[['disease_context','donor_uid_v2']].drop_duplicates())
d=pd.concat(frames,ignore_index=True);assert len(frames)==79
dm=pd.concat(donors).drop_duplicates();overlap=dm.groupby('donor_uid_v2').disease_context.nunique();assert overlap.le(1).all(),'Cross-context donor overlap requires covariance handling'
dm.to_csv(B/'tables/donor_context_independence_audit.tsv',sep='\t',index=False)
d.to_csv(B/'inputs/all_manuscript_scope_DE.tsv.gz',sep='\t',index=False)
with (B/'logs/full_gene_meta.log').open('x') as f:subprocess.run(['/usr/bin/Rscript',str(B/'code/meta_DL.R'),str(B),'run'],check=True,stdout=f,stderr=subprocess.STDOUT)
meta=pd.read_csv(B/'tables/full_gene_DL_meta.tsv.gz',sep='\t')
for _,z in meta.groupby(['annotation_level','cell_state']):assert np.allclose(false_discovery_control(z.p.to_numpy(dtype=float)),z.FDR_all_tested_genes,rtol=1e-7,atol=1e-12)
jobs=[];excluded=[]
def ranktable(z,aid,role,stat,pcol,bcol):
    x=pd.DataFrame(dict(gene=z.gene_symbol,statistic=stat,selected_up=(z[pcol]<.05)&(z[bcol]>0),selected_down=(z[pcol]<.05)&(z[bcol]<0)))
    bad=x.gene.isna()|x.gene.eq('')|x.gene.duplicated(False)|~np.isfinite(x.statistic)
    if bad.any():
        ex=x[bad].copy();ex['analysis_id']=aid;excluded.append(ex)
    x=x[~bad];assert len(x)>100 and x.gene.is_unique
    x['analysis_id']=aid;p=B/'inputs'/(aid+'.tsv');x.to_csv(p,sep='\t',index=False)
    jobs.append(dict(id=aid,role=role,annotation_level=str(z.annotation_level.iloc[0]),cell_state=str(z.cell_state.iloc[0]),contexts=sorted(z.disease_context.unique().tolist()) if 'disease_context' in z else 'at_least_two_original_contexts_per_gene',ranks=str(p),ranks_sha256=sha(p)))
for (level,state),z in meta.groupby(['annotation_level','cell_state']):
    aid='PRIMARY_META_DL_'+level+'_'+hashlib.sha256(state.encode()).hexdigest()[:12]
    ranktable(z,aid,'primary_original_pathway_DL',z.z,'FDR_all_tested_genes','beta')
for model,z in d.groupby('source_model'):
    assert np.allclose(z.loc[z.F.lt(0),'PValue'],1)
    ranktable(z,'SUPPORTING_'+model,'supporting_single_context',np.sign(z.logFC)*np.sqrt(z.F.clip(lower=0)),'FDR','logFC')
if excluded:pd.concat(excluded).to_csv(B/'tables/ambiguous_or_invalid_symbol_exclusions.tsv.gz',sep='\t',index=False)
(B/'ANALYSES.json').write_text(json.dumps(jobs,indent=2))
done=[]
for j in jobs:
    out=B/'outputs'/j['id']
    with (B/'logs'/(j['id']+'.log')).open('x') as f:subprocess.run(['/usr/bin/Rscript',str(B/'code/native_pathway.R'),j['ranks'],rec['gene_sets'],str(out)],stdout=f,stderr=subprocess.STDOUT,check=True)
    ora=pd.read_csv(out/'ora_all_tested_before_display_filter.tsv.gz',sep='\t')
    for _,z in ora.groupby(['analysis_id','database','direction']):assert np.allclose(false_discovery_control(z.pval.to_numpy(dtype=float)),z.padj,rtol=1e-7,atol=1e-12)
    g=pd.read_csv(out/'fgsea_all.tsv.gz',sep='\t');bad=g[['NES','pval','padj']].isna().any(axis=1)
    assert g.loc[bad,'pval'].isna().all() and g.loc[bad,'padj'].isna().all()
    assert np.isfinite(g.loc[~bad,['NES','pval','padj']]).all().all()
    for _,z in g[~bad].groupby(['analysis_id','database']):assert np.allclose(false_discovery_control(z.pval.to_numpy(dtype=float)),z.padj,rtol=1e-6,atol=1e-12)
    g[bad].to_csv(out/'native_unestimable_pathways.tsv',sep='\t',index=False)
    item=dict(id=j['id'],role=j['role'],status='COMPLETE_WITH_NATIVE_NA' if bad.any() else 'COMPLETE',unestimable=int(bad.sum()),files={str(p):sha(p) for p in out.iterdir() if p.is_file()})
    (out/'VALIDATED_COMPLETE.json').write_text(json.dumps(item,indent=2));done.append(item)
    (B/'PROGRESS.json').write_text(json.dumps(dict(completed=len(done),expected=len(jobs),last=j['id'],primary_meta_analyses=sum(x['role']=='primary_original_pathway_DL' for x in jobs),supporting_contexts=79),indent=2));print('COMPLETE',j['id'],len(done),flush=True)
(B/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE_MANUSCRIPT_PRIMARY_PSEUDOBULK_ENRICHMENT',results=done,approximate_SE_limitation_retained=True,previous_DE_REML_results_unchanged=True),indent=2))
