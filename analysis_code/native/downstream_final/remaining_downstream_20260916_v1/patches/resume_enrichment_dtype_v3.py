import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,subprocess,time
from pathlib import Path
import pandas as pd,numpy as np
from scipy.stats import false_discovery_control
N=Path(__file__).resolve().parents[1];P=N/'enrichment';O=P/'dtype_resume_patch_v3'
O.mkdir(exist_ok=False);(O/'outputs').mkdir();(O/'logs').mkdir()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rec=json.loads((P/'PREPARATION.json').read_text())
for p,h in rec['bound_files'].items():assert sha(p)==h,p
assert json.loads((P/'tests/PASS.json').read_text())['native_sha256']==sha(P/'code/native.R')
os.environ['R_LIBS_USER']=str(P/'Rlib')
binding={str(p):sha(p) for p in [P/'PREPARATION.json',P/'ANALYSES.json',P/'code/native.R',P/'tests/PASS.json',Path(__file__)]}
binding['prior_native_results']={str(f):sha(f) for root in [P/'outputs',P/'native_NA_review_patch_v2/outputs'] for f in root.rglob('*') if f.is_file()}
(O/'RECEIPT.json').write_text(json.dumps(binding,indent=2))
jobs=json.loads((P/'ANALYSES.json').read_text());results=[]
for j in jobs:
    candidates=[P/'outputs'/j['id'],P/'native_NA_review_patch_v2/outputs'/j['id']]
    existing=[x for x in candidates if (x/'fgsea_all.tsv.gz').exists()]
    assert len(existing)<=1
    prior=existing[0] if existing else candidates[0];reuse=bool(existing)
    out=prior if reuse else O/'outputs'/j['id']
    if reuse:
        assert j['id'] in {x['id'] for x in jobs[:5]}, 'Unexpected prior partial output'
        required=['fgsea_all.tsv.gz','ora_all_tested_before_display_filter.tsv.gz','sessionInfo.txt','gene_set_version.tsv']
        assert all((out/f).exists() for f in required)
    else:
        with (O/'logs'/(j['id']+'.log')).open('x') as log:subprocess.run(['/usr/bin/Rscript',str(P/'code/native.R'),j['ranks'],rec['gene_sets'],str(out)],check=True,stdout=log,stderr=subprocess.STDOUT)
    ora=pd.read_csv(out/'ora_all_tested_before_display_filter.tsv.gz',sep='\t')
    for _,z in ora.groupby(['analysis_id','database','direction']):assert np.allclose(false_discovery_control(z.pval.to_numpy(dtype=float)),z.padj.to_numpy(),rtol=1e-7,atol=1e-12)
    g=pd.read_csv(out/'fgsea_all.tsv.gz',sep='\t')
    assert g.pval.dropna().between(0,1).all() and g.padj.dropna().between(0,1).all()
    invalid=g[['NES','pval','padj']].isna().any(axis=1)
    # Native fgsea can return unestimable pathways; retain their NA and identify them explicitly.
    assert g.loc[invalid,'pval'].isna().all() and g.loc[invalid,'padj'].isna().all()
    assert np.isfinite(g.loc[~invalid,['NES','pval','padj']]).all().all()
    for _,z in g[~invalid].groupby(['analysis_id','database']):assert np.allclose(false_discovery_control(z.pval.to_numpy(dtype=float)),z.padj.to_numpy(),rtol=1e-6,atol=1e-12)
    dest=O/'outputs'/j['id'];dest.mkdir(exist_ok=True)
    g.loc[invalid].to_csv(dest/'UNESTIMABLE_NATIVE_GSEA_PATHWAYS.tsv',sep='\t',index=False)
    d=dict(status='VALIDATED_WITH_NATIVE_NA_PATHWAYS' if invalid.any() else 'COMPLETE_VALIDATED',id=j['id'],reused_native=reuse,n_GSEA=len(g),n_native_NA=int(invalid.sum()),source_path=str(out),files={str(p):sha(p) for p in out.iterdir() if p.is_file()},rank_sha256=sha(j['ranks']),native_sha256=sha(P/'code/native.R'))
    (dest/'VALIDATED_COMPLETE.json').write_text(json.dumps(d,indent=2));results.append(d)
    (O/'PROGRESS.json').write_text(json.dumps(dict(completed=len(results),expected=len(jobs),last=j['id'],native_NA_total=sum(x['n_native_NA'] for x in results)),indent=2));print('COMPLETE',j['id'],len(results),'native_NA',int(invalid.sum()),flush=True)
(O/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE_WITH_NATIVE_UNESTIMABLE_PATHWAYS_RETAINED',analyses=results,cross_context_meta='NOT_RUN'),indent=2))
