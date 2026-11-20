import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='2'
os.environ['CUDA_VISIBLE_DEVICES']=''
import sys,json,time,resource
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
W=Path(__file__).resolve().parents[1];sys.path.insert(0,str(W/'code'))
from scdrs_adapter import prepare_scdrs_counts,build_derived,native_sequence,sha
R=W.parents[1];start=time.time();out=W/'tests'/'actual_obs_smoke_v1';out.mkdir(exist_ok=False)
# Bounded exact-dtype unit cases exercise production helper.
gm=pd.DataFrame({'feature_id':['a','b'],'gene_symbol':['B','A'],'ensembl_id':['ENSG1','ENSG2']})
x=sparse.csr_matrix([[2,3],[5,8]],dtype=np.float32);y,n,audit=prepare_scdrs_counts(x,gm,['a','b']);assert y.dtype==np.float32 and np.array_equal(y.toarray(),[[3,2],[8,5]])
dup=gm.copy();dup.gene_symbol='A';dup.ensembl_id='ENSG1';y,_,_=prepare_scdrs_counts(x,dup,['a','b']);assert np.array_equal(y.toarray(),[[5],[13]])
bad=dup.copy();bad.loc[1,'ensembl_id']='ENSG2'
try:prepare_scdrs_counts(x,bad,['a','b']);raise AssertionError('ambiguous duplicate accepted')
except ValueError:pass
big=sparse.csr_matrix([[2**24+1,1]],dtype=np.int64);y,_,_=prepare_scdrs_counts(big,gm,['a','b']);assert y.dtype==np.int64 and int(y[0,1])==2**24+1
atlas=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad';a=ad.read_h5ad(atlas,backed='r')
mapfile=W/'inputs/gene_map_scDRS_READY.tsv';assert sha(mapfile)=='55378df49e46c988dd641a633c041951281af5189447ca952fc805bfad85e24c'
gm=pd.read_csv(mapfile,sep='\t');assert gm.feature_id.astype(str).tolist()==a.var_names.astype(str).tolist()
parent=a.obs.level1_annotation.astype(str);ng=a.uns['neighbors_cpu_n30'];graph=a.obsp[ng['connectivities_key']].tocsr();rng=np.random.default_rng(20260914)
selected=[];selection={}
for label in parent.value_counts().index[:2]:
 eligible=np.flatnonzero(parent.to_numpy()==label);allowed=np.zeros(a.n_obs,bool);allowed[eligible]=True
 chosen=set();queue=[]
 # Start many observed-study seeds then expand real frozen neighbors.
 for _,group in a.obs.iloc[eligible].groupby('study',observed=True):
  pos=a.obs_names.get_indexer(group.index);queue.extend(rng.choice(pos,min(2,len(pos)),replace=False).tolist())
 while len(chosen)<700:
  if not queue:queue.append(int(rng.choice(eligible)))
  q=queue.pop(0)
  if q in chosen:continue
  chosen.add(q)
  neigh=graph.indices[graph.indptr[q]:graph.indptr[q+1]];queue.extend([int(v) for v in neigh if allowed[v] and int(v) not in chosen])
 selected.extend(chosen);selection[label]=len(chosen)
rows=np.array(sorted(selected));counts=a.X[rows,:].tocsr();obs=a.obs.iloc[rows].copy();small=ad.AnnData(X=counts,obs=obs,var=a.var.copy())
small.uns['neighbors_cpu_n30']=ng
for key in [ng['connectivities_key'],ng['distances_key']]:small.obsp[key]=a.obsp[key][rows,:][:,rows].copy()
a.file.close();del a,graph
covpath=R/'production_launches/traits_20260914T023900Z/inputs/scDRS_raw_verified.cov.gz';assert sha(covpath)=='826c7d919686301f548decf35685417d0d41e4b25ad67ca715dbdb4d87305037'
cov=pd.read_csv(covpath,sep='\t',index_col=0).loc[small.obs_names]
report=build_derived(small,gm,cov,out/'derived.h5ad')
b=ad.read_h5ad(out/'derived.h5ad');order=pd.Index(gm.scdrs_gene_symbol).get_indexer(b.var_names)
assert np.array_equal(counts[:,order].toarray(),b.X.toarray());assert b.obs_names.equals(small.obs_names)
for k in ['level1_annotation','level2_annotation']:assert b.obs[k].astype(str).equals(small.obs[k].astype(str))
assert (b.obsp['connectivities']!=small.obsp[ng['connectivities_key']]).nnz==0
manifest=json.loads((R/'enrichment/scdrs_label_independent/genesets/GENESETS_FROZEN.json').read_text())
scenarios=[{'name':s['scenario'],'geneset':s['gs_file'],'sha256':s['gs_sha256']} for s in manifest['scenarios']]
records=native_sequence(out/'derived.h5ad',out/'native_covariates.cov.gz',scenarios,out/'native','/home/anzhen_vhd/miniconda3/envs/scdrs_env/bin/scdrs',report['group_keys'],2)
report.update({'status':'PASS_NATIVE_SHARED_ADAPTER_ACTUAL_OBS','selection':selection,'seconds':time.time()-start,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'n_native_commands':len(records),'graph_note':'Induced subgraph of frozen_global_n30, no recomputed PCA/neighbors. Smoke is compatibility test, not scientific inference. Full production uses complete graph.','unit_tests':4,'atlas_original_shape':[1175115,37336],'actual_obs_smoke_shape':list(small.shape)})
(out/'SMOKE_RESULT.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
