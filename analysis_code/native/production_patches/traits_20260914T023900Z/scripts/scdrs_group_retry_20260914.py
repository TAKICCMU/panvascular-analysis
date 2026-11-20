import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['CUDA_VISIBLE_DEVICES']=''
import anndata as ad,numpy as np,subprocess,json,time,scanpy as sc
from pathlib import Path
W=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/traits_20260914T023900Z')
p=W/'smoke/scdrs_synthetic.h5ad';a=ad.read_h5ad(p);c=a.copy();sc.pp.normalize_total(c,target_sum=1e4);sc.pp.log1p(c);sc.pp.pca(c,n_comps=20,random_state=0);sc.pp.neighbors(c,n_neighbors=15,n_pcs=20,random_state=0);a.obsp['connectivities']=c.obsp['connectivities'];a.obsp['distances']=c.obsp['distances'];a.uns['neighbors']=c.uns['neighbors'];new=W/'smoke/scdrs_synthetic_with_graph.h5ad';a.write_h5ad(new,compression='lzf');out=W/'smoke/scdrs_native/group_retry';out.mkdir();score=next((W/'smoke/scdrs_native/scores').glob('*.full_score.gz'))
cmd=['/home/anzhen_vhd/miniconda3/envs/scdrs_env/bin/scdrs','perform-downstream','--h5ad-file',str(new),'--score-file',str(score),'--out-folder',str(out),'--group-analysis','smoke_group','--flag-filter-data','True','--flag-raw-count','True']
with open(W/'logs/scdrs_group_retry.log','w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
import pandas as pd
z={'returncode':r.returncode,'tiny_only':True,'graph_source':'CPU Scanpy neighbors of synthetic normalized counts/PCA20;k15;seed0','changed_native_scdrs':False,'command':cmd}
if r.returncode==0:
 fs=list(out.glob('*.scdrs_group.*'));assert len(fs)==1;d=pd.read_csv(fs[0],sep='\t');assert {'group','assoc_mcp','hetero_mcp'}<=set(d);assert np.isfinite(d[['assoc_mcp','hetero_mcp']].to_numpy()).all();z['schema_pass']=True;z['result_rows']=len(d)
(W/'evidence/SCDRS_GROUP_NATIVE_RETRY.json').write_text(json.dumps(z,indent=2));print(json.dumps(z))
