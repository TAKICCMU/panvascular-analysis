"""Integrate only the pilot-selected CellHint configuration using shared alignment."""
import os
os.environ.update(OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',NUMBA_NUM_THREADS='4',LOKY_MAX_CPU_COUNT='4')
if hasattr(os,'sched_getaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
import json,time,hashlib
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd,cellhint
from scipy import sparse
from prepare_cellhint_primary_shared import fh,nh,R,W,OUT

def main():
 selected=json.loads((OUT/'pilot/CELLHINT_PILOT_SELECTION.json').read_text())['selected_cellhint']
 assert selected in ['meta3','meta5']
 work=OUT/'primary';shared=work/'shared_harmonization';out=work/'embeddings';out.mkdir(parents=True,exist_ok=True)
 cp=out/f'cellhint__{selected}.json'
 if cp.exists():return
 contract=json.loads((shared/'SHARED_HARMONIZATION_COMPLETE.json').read_text())
 assert contract['status']=='completed'
 p=W/'inputs/scib_primary_sampleonly_cap1500_hvg2000.h5ad';a=ad.read_h5ad(p,backed='r');obs=a.obs.copy();a.file.close()
 assert nh(obs.index)==contract['obs_names_sha256']
 pp=W/'omicverse_pca/primary/omicverse_pca.npy';pca=np.load(pp,mmap_mode='r');assert fh(pp)==contract['pca_sha256']
 mapped=pd.read_csv(shared/'cellhint_primary_reannotation.tsv.gz',sep='\t',index_col='cell_id').reindex(obs.index)
 for c in ['reannotation','group']:obs[f'cellhint_harmonized_{c}']=mapped[c].fillna('UNASSIGNED').astype(str)
 data=ad.AnnData(X=np.zeros((len(obs),1),dtype=np.float32),obs=obs);data.obsm['X_pca_omicverse']=np.asarray(pca,dtype=np.float32)
 t=time.time()
 # Preserve the original default n_neighbors=max(15,3*n_batches)=741. This
 # native graph remains broad; only its distance-ranked K subset enters scIB.
 cellhint.integrate(data,batch='sample_uid',cell_type='cellhint_harmonized_group',use_rep='X_pca_omicverse',n_meta_neighbors=int(selected[-1]))
 distkey=data.uns['neighbors'].get('distances_key','distances');connkey=data.uns['neighbors'].get('connectivities_key','connectivities')
 d=data.obsp[distkey].tocsr();d.sum_duplicates()
 sparse.save_npz(out/f'cellhint__{selected}.native_distances.npz',d)
 sparse.save_npz(out/f'cellhint__{selected}.native_connectivities.npz',data.obsp[connkey].tocsr())
 idx=np.empty((len(obs),30),np.int32);vals=np.empty((len(obs),30),np.float32);degrees=[]
 for i in range(len(obs)):
  lo,hi=d.indptr[i:i+2];cols=d.indices[lo:hi];ds=d.data[lo:hi];keep=cols!=i;cols=cols[keep];ds=ds[keep];degrees.append(len(cols))
  if len(cols)<30:raise RuntimeError(f'Native graph row {i} has only {len(cols)} real neighbors; no padding permitted')
  order=np.argsort(ds)[:30];idx[i]=cols[order];vals[i]=ds[order]
 np.savez_compressed(out/f'cellhint__{selected}.graph.npz',indices=idx,distances=vals)
 record={'status':'completed','method':'cellhint','config_id':selected,'output_type':'graph_only','obs_names_sha256':nh(obs.index),'input_shape':[len(obs),2000],
         'parameters':{'n_meta_neighbors':int(selected[-1]),'n_neighbors_native':max(15,obs.sample_uid.nunique()*3),'K_exported':30,'batch_key':'sample_uid','input_label_key':'studywise_marker_level1'},
         'shared_harmonization':contract,'minimum_native_nonself_degree':int(min(degrees)),'elapsed_seconds':time.time()-t,'native_graph_saved':True,'never_self_padded':True}
 cp.write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2),flush=True)
if __name__=='__main__':main()
