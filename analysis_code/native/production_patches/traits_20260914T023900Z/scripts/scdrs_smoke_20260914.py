import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['CUDA_VISIBLE_DEVICES']=''
import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp,json,subprocess,time
from pathlib import Path
W=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/traits_20260914T023900Z')
r=np.random.default_rng(20260914);n=400;g=1600
X=r.poisson(r.gamma(2,1,g)[None,:]*r.uniform(.5,1.5,n)[:,None]).astype('float32')
a=ad.AnnData(X=sp.csr_matrix(X),obs=pd.DataFrame({'smoke_group':pd.Categorical(['A']*200+['B']*200)},index=[f'SMOKE_{i}' for i in range(n)]),var=pd.DataFrame(index=[f'SMOKE_G{i}' for i in range(g)]))
a.write_h5ad(W/'smoke/scdrs_synthetic.h5ad',compression='lzf')
cov=pd.DataFrame({'const':np.ones(n)},index=a.obs_names)
for key,x in {'log1p_counts':np.log1p(X.sum(1)),'log1p_genes':np.log1p((X>0).sum(1)),'pct_mt':X[:,:13].sum(1)/X.sum(1)*100}.items():cov[key]=(x-x.mean())/x.std()
cov.to_csv(W/'smoke/scdrs_synthetic.cov',sep='\t')
(W/'smoke/scdrs_synthetic.gs').write_text('TRAIT\tGENESET\nSMOKE_TRAIT\t'+','.join(f'SMOKE_G{i}:1' for i in range(50))+'\n')
exe='/home/anzhen_vhd/miniconda3/envs/scdrs_env/bin/scdrs'
out=W/'smoke/scdrs_native';(out/'scores').mkdir(parents=True,exist_ok=True)
cmd=[exe,'compute-score','--h5ad-file',str(W/'smoke/scdrs_synthetic.h5ad'),'--h5ad-species','human','--gs-file',str(W/'smoke/scdrs_synthetic.gs'),'--gs-species','human','--cov-file',str(W/'smoke/scdrs_synthetic.cov'),'--ctrl-match-opt','mean_var','--weight-opt','vs','--flag-filter-data','True','--flag-raw-count','True','--n-ctrl','1000','--flag-return-ctrl-raw-score','False','--flag-return-ctrl-norm-score','True','--out-folder',str(out/'scores')]
(W/'smoke/scdrs_compute_command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
with open(W/'logs/scdrs_compute_smoke.log','w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
z={'compute_returncode':p.returncode,'n_cells':n,'n_genes':g,'n_ctrl':1000,'synthetic_only':True,'seconds':time.time()-t}
if p.returncode==0:
 fs=list((out/'scores').glob('*.full_score.gz'));assert len(fs)==1
 d=pd.read_csv(fs[0],sep='\t',index_col=0);z['columns']=len(d.columns);z['norm_ctrl_columns']=len([k for k in d if k.startswith('ctrl_norm_score')]);assert z['norm_ctrl_columns']==1000
 z['full_score_gzip_bytes']=fs[0].stat().st_size;z['estimated_full_atlas_gzip_gib_per_scenario']=fs[0].stat().st_size/n*1175115/1024**3;z['estimate_basis']='synthetic gzip row size, preliminary not biological score compressibility guarantee'
 cmd=[exe,'perform-downstream','--h5ad-file',str(W/'smoke/scdrs_synthetic.h5ad'),'--score-file',str(fs[0]),'--out-folder',str(out/'group'),'--group-analysis','smoke_group','--flag-filter-data','True','--flag-raw-count','True']
 (W/'smoke/scdrs_group_command.json').write_text(json.dumps(cmd,indent=2))
 with open(W/'logs/scdrs_group_smoke.log','w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 z['group_returncode']=p.returncode
 if p.returncode==0:
  ff=list((out/'group').glob('*.scdrs_group.*'));assert len(ff)==1;tab=pd.read_csv(ff[0],sep='\t');assert {'group','assoc_mcp','hetero_mcp'}<=set(tab);z['group_schema_pass']=True
(W/'evidence/SCDRS_NATIVE_SMOKE.json').write_text(json.dumps(z,indent=2));print(json.dumps(z),flush=True)
