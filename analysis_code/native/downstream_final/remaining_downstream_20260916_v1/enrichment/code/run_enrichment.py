import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,subprocess,sys,time
from pathlib import Path
P=Path(__file__).resolve().parents[1]
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
os.environ['R_LIBS_USER']=str(P/'Rlib')
rec=json.loads((P/'PREPARATION.json').read_text())
for p,h in rec['bound_files'].items():assert sha(p)==h,p
assert (P/'tests/PASS.json').exists()
binding={str(p):sha(p) for p in [P/'PREPARATION.json',P/'ANALYSES.json',P/'code/native.R',Path(__file__),P/'tests/PASS.json']}
(P/'ACTUAL_INPUT_RECEIPT.json').write_text(json.dumps(binding,indent=2))
jobs=json.loads((P/'ANALYSES.json').read_text());results=[]
(P/'outputs').mkdir(exist_ok=True);(P/'logs').mkdir(exist_ok=True)
for j in jobs:
    out=P/'outputs'/j['id'];done=out/'VALIDATED_COMPLETE.json'
    if done.exists():
        d=json.loads(done.read_text());assert d['rank_sha256']==sha(j['ranks']) and d['native_sha256']==sha(P/'code/native.R')
        for p,h in d['outputs'].items():assert sha(p)==h
    else:
        assert not out.exists(),'Incomplete output requires independent recovery: '+str(out)
        with (P/'logs'/(j['id']+'.log')).open('x') as log:
            subprocess.run(['/usr/bin/Rscript',str(P/'code/native.R'),j['ranks'],rec['gene_sets'],str(out)],check=True,stdout=log,stderr=subprocess.STDOUT)
        import pandas as pd,numpy as np
        from scipy.stats import false_discovery_control
        ora=pd.read_csv(out/'ora_all_tested_before_display_filter.tsv.gz',sep='\t')
        for _,z in ora.groupby(['analysis_id','database','direction']):assert np.allclose(false_discovery_control(z.pval.to_numpy()),z.padj.to_numpy(),rtol=1e-7,atol=1e-12)
        gsea=pd.read_csv(out/'fgsea_all.tsv.gz',sep='\t')
        assert gsea.NES.notna().all() and gsea.pval.dropna().between(0,1).all()
        done.write_text(json.dumps(dict(status='COMPLETE_VALIDATED',rank_sha256=sha(j['ranks']),native_sha256=sha(P/'code/native.R'),outputs={str(p):sha(p) for p in out.iterdir() if p.is_file()},n_GSEA=len(gsea),n_ORA=len(ora),GSEA_unresolved_pvalues=int(gsea.pval.isna().sum())),indent=2))
    results.append(dict(id=j['id'],complete=str(done)));(P/'PROGRESS.json').write_text(json.dumps(dict(completed=len(results),expected=len(jobs),last=j['id'],updated=time.time()),indent=2));print('COMPLETE',j['id'],len(results),flush=True)
(P/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE_VALIDATED_CONTEXT_AND_GENETIC_ENRICHMENT',results=results,cross_context_meta_enrichment='NOT_RUN_NO_VALIDATED_SE_INPUT'),indent=2))
