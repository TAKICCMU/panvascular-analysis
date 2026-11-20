import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import fcntl,json,hashlib,time,datetime
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');W=R/'production_launches/traits_20260914T023900Z';L=R/'production_launches/execution_control_20260914/.auxiliary_compute.lock';L.parent.mkdir(exist_ok=True);lock=open(L,'a');print('Waiting auxiliary compute lock',flush=True);fcntl.flock(lock,fcntl.LOCK_EX);print('Acquired auxiliary lock: 2 threads; raw stream begin',flush=True)
import h5py,numpy as np,pandas as pd,scipy.sparse as sp
from anndata.io import read_elem
p=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad';st=p.stat();t=time.time();chunk=10000
with h5py.File(p,'r') as f:
 v=read_elem(f['var']);names=read_elem(f['obs']['_index']);n=len(names);m=len(v);mt=v.gene_symbol.astype('string').str.startswith('MT-',na=False).to_numpy();assert np.array_equal(mt,v.mt.fillna(False).to_numpy());ptr=f['X/indptr'][:];total=np.empty(n);genes=np.empty(n);pct=np.empty(n);samples=[];explicit_zero=0
 X=f['X'];storage={k:int(X[k].id.get_storage_size()) for k in X}
 for start in range(0,n,chunk):
  stop=min(start+chunk,n);lo=int(ptr[start]);hi=int(ptr[stop]);data=X['data'][lo:hi];ix=X['indices'][lo:hi];assert np.isfinite(data).all() and (data>=0).all() and np.equal(data,np.floor(data)).all();b=sp.csr_matrix((data,ix,ptr[start:stop+1]-lo),shape=(stop-start,m));s=np.asarray(b.sum(axis=1,dtype=np.float64)).ravel();assert (s>0).all();total[start:stop]=s;explicit_zero+=int((data==0).sum());b.eliminate_zeros();genes[start:stop]=np.diff(b.indptr);pct[start:stop]=np.asarray(b[:,mt].sum(axis=1,dtype=np.float64)).ravel()/s*100
  samples.append(b[:32].copy())
  if start%100000==0:print(json.dumps({'processed_cells':stop,'elapsed_sec':time.time()-t}),flush=True)
 cov=pd.DataFrame({'const':np.ones(n)},index=pd.Index(names,name='index'));standard={}
 for key,x in {'log1p_counts':np.log1p(total),'log1p_genes':np.log1p(genes),'pct_mt':pct}.items():
  mu=float(np.mean(x));sd=float(np.std(x,ddof=0));standard[key]={'mean':mu,'sd_population_ddof0':sd,'zero_variance':sd==0,'zero_variance_policy':'return exact zero standardized covariate only if actually constant'};cov[key]=(x-mu)/sd if sd else np.zeros(n)
 assert np.isfinite(cov.to_numpy()).all();cov.to_csv(W/'inputs/scDRS_raw_verified.cov.gz',sep='\t',compression='gzip');pd.DataFrame({'total_counts_raw':total,'n_genes_raw':genes,'pct_mt_raw':pct},index=cov.index).to_csv(W/'evidence/raw_covariate_values.tsv.gz',sep='\t',compression='gzip')
 compare={}
 for key,new in [('total_counts',total),('library_size_raw_cpu',total),('n_genes_by_counts',genes),('pct_counts_mt',pct)]:
  old=np.asarray(read_elem(f['obs'][key]),dtype=float);delta=np.abs(old-new);compare[key]={'finite_old':int(np.isfinite(old).sum()),'exact_match':int((delta==0).sum()),'max_abs_difference':float(delta.max()),'within_1e_5':int((delta<=1e-5).sum())}
 fullgraph={k:sum(int(f['obsp'][k][j].id.get_storage_size()) for j in f['obsp'][k]) for k in f['obsp']}
 sample=sp.vstack(samples,format='csr');nnz=sample.nnz
 with h5py.File(W/'smoke/raw_serialization_sample.h5','w') as q:
  for k,arr in [('data',sample.data.astype('float64')),('indices',sample.indices.astype('int64')),('indptr',sample.indptr.astype('int64'))]:q.create_dataset(k,data=arr,compression='lzf')
 with h5py.File(W/'smoke/raw_serialization_sample.h5') as q:compressed={k:q[k].id.get_storage_size() for k in q}
 estimated=(compressed['data']+compressed['indices'])/nnz*int(ptr[-1])+compressed['indptr']/sample.shape[0]*n+sum(fullgraph.values())
 st2=p.stat();assert (st.st_size,st.st_mtime_ns,st.st_ino)==(st2.st_size,st2.st_mtime_ns,st2.st_ino)
 def sha(path):
  h=hashlib.sha256()
  with open(path,'rb') as f:
   for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
  return h.hexdigest()
 report={'status':'PASS_RAW_DERIVED_COVARIATES','n_cells':n,'cell_order_matches_final_obs':True,'atlas_sha256_from_certified_release':'00fc74c1ddaa5e4b32d8664f134e0ec5a0c762f5e9a9fb45b96acbd40804a310','atlas_stat_unchanged':True,'elapsed_seconds':time.time()-t,'mt_gene_count':int(mt.sum()),'explicit_zero_entries_removed_for_detection_count':explicit_zero,'standardization':standard,'metadata_comparison':compare,'covariates':str(W/'inputs/scDRS_raw_verified.cov.gz'),'covariate_sha256':sha(W/'inputs/scDRS_raw_verified.cov.gz'),'scdrs_final_input_not_ready_due_to_gene_map':True,'serialization_sample_cells':sample.shape[0],'serialization_sample_nnz':nnz,'sample_compressed_bytes':compressed,'full_source_X_actual_stored_bytes':storage,'full_frozen_graph_stored_bytes':fullgraph,'estimated_float64_raw_symbols_and_graph_GiB':estimated/1024**3,'estimate_not_a_guarantee':'stratified-by-row 32 cells per 10000; compression, metadata and symbol mapping affect final size','no_full_derived_matrix_written':True}
 (W/'evidence/REAL_COVARIATES_AND_DISK_ESTIMATE.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
fcntl.flock(lock,fcntl.LOCK_UN)
