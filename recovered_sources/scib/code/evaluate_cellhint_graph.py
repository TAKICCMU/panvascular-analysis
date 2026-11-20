"""Validate genuine CellHint neighbors, then apply the common scIB self convention."""
import os
os.environ.update(OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',NUMBA_NUM_THREADS='4',JAX_PLATFORM_NAME='cpu')
if hasattr(os,'sched_getaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
import argparse,hashlib,json,time
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
import scib_metrics as sm
from scib_metrics.nearest_neighbors import NeighborsResults

def fh(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def nh(v):
 h=hashlib.sha256()
 for x in v:h.update(str(x).encode());h.update(b'\0')
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--metadata',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--stage',required=True);args=ap.parse_args()
 args.output.mkdir(parents=True,exist_ok=True)
 meta=json.loads(args.metadata.read_text());gp=args.metadata.with_name(args.metadata.stem+'.graph.npz')
 signature={'graph_sha256':fh(gp),'metadata_sha256':fh(args.metadata),'input_sha256':fh(args.input),'evaluator':'takeover.graph-self-standard.v1'}
 statusp=args.output/'GRAPH_SCIB_COMPLETE.json'
 if statusp.exists() and json.loads(statusp.read_text()).get('signature')==signature:return
 a=ad.read_h5ad(args.input,backed='r');obs=a.obs.copy();names=a.obs_names.copy();a.file.close()
 assert meta['status']=='completed' and meta['obs_names_sha256']==nh(names)
 z=np.load(gp);inds=np.asarray(z['indices'],dtype=np.int32);dist=np.asarray(z['distances'],dtype=np.float32)
 n=len(names);rows=np.arange(n)[:,None]
 duplicate=np.array([len(r)-len(np.unique(r)) for r in inds])
 audit={'n_cells':n,'input_shape':list(inds.shape),'self_entries':int((inds==rows).sum()),'duplicate_rows':int((duplicate>0).sum()),
        'invalid_indices':int(((inds<0)|(inds>=n)).sum()),'nonfinite_distances':int((~np.isfinite(dist)).sum()),'negative_distances':int((dist<0).sum())}
 audit['passed']=bool(inds.shape==dist.shape and inds.shape==(n,30) and not any(audit[k] for k in ['self_entries','duplicate_rows','invalid_indices','nonfinite_distances','negative_distances']))
 (args.output/'CELLHINT_GRAPH_AUDIT.json').write_text(json.dumps(audit,indent=2))
 if not audit['passed']:raise RuntimeError('CellHint raw export has padding/invalid neighbors; native graph recovery needed')
 # One genuine self entry is required by NeighborsResults, exactly matching the
 # existing embedding PyNNDescent convention. Never pad missing neighbors.
 inds=np.concatenate([rows.astype(np.int32),inds[:,:29]],axis=1)
 dist=np.concatenate([np.zeros((n,1),dtype=np.float32),dist[:,:29]],axis=1)
 labels=obs['studywise_marker_level1'].astype(str).to_numpy();batches=obs['sample_uid'].astype(str).to_numpy()
 t=time.time();records=[]
 for k in [10,20,30]:
  nn=NeighborsResults(indices=inds[:,:k],distances=dist[:,:k])
  row={'stage':args.stage,'method':'cellhint','config_id':meta['config_id'],'output_type':'graph_only','benchmark_track':'label-informed graph',
       'K':k,'n_cells':n,'n_metric_subset_cells':0,'ASW_label':np.nan,'ASW_batch':np.nan,'PCR_batch':np.nan,
       'iLISI_batch':float(sm.ilisi_knn(nn,batches,scale=True)),'cLISI_label':float(sm.clisi_knn(nn,labels,scale=True)),
       'graph_connectivity':float(sm.graph_connectivity(nn,labels)),'embedding_metadata':str(args.metadata),
       'knn_convention':'K includes 1 genuine self and K-1 real neighbors; no padding'}
  for r in [.2,.5,1.,1.5]:
   try:s=sm.nmi_ari_cluster_labels_leiden(nn,labels,optimize_resolution=False,resolution=r,n_jobs=4,seed=0)
   except TypeError as e:
    if 'seed' not in str(e):raise
    s=sm.nmi_ari_cluster_labels_leiden(nn,labels,optimize_resolution=False,resolution=r,n_jobs=4)
   row[f'NMI_resolution_{r:g}']=float(s['nmi']);row[f'ARI_resolution_{r:g}']=float(s['ari'])
  row['NMI_best']=max(row[f'NMI_resolution_{r:g}'] for r in [.2,.5,1.,1.5]);row['ARI_best']=max(row[f'ARI_resolution_{r:g}'] for r in [.2,.5,1.,1.5])
  row['evaluation_seconds']=time.time()-t;records.append(row)
  pd.DataFrame(records).to_csv(args.output/'cellhint_metrics_by_K.partial.tsv',sep='\t',index=False)
  print(json.dumps(row),flush=True)
 d=pd.DataFrame(records);d.to_csv(args.output/'cellhint_metrics_by_K.tsv',sep='\t',index=False)
 statusp.write_text(json.dumps({'status':'completed','signature':signature,'parameters':{'K':[10,20,30],'resolutions':[.2,.5,1.,1.5],'seed':0},'elapsed_seconds':time.time()-t},indent=2))
if __name__=='__main__':main()
