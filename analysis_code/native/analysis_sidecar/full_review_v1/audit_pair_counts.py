import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib
import h5py,numpy as np,pandas as pd
from scipy import sparse
W=Path(__file__).resolve().parent
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
atlas=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad'
cert=json.loads((R/'maintenance_resume_20260914/ATLAS_VERIFIED.json').read_text())
st=atlas.stat()
for field,val in [('size',st.st_size),('inode',st.st_ino),('device',st.st_dev),('mtime_ns',st.st_mtime_ns)]:
 assert cert[field]==val
rows=pd.read_csv(W/'all247_sample_aggregation_audit.tsv',sep='\t',keep_default_na=False)
targets=rows[rows.study.eq('GSE237230')]
def decode(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
def obs_col(h,key):
 x=h['obs'][key]
 if isinstance(x,h5py.Group):
  categories=decode(x['categories'][:]);codes=x['codes'][:];assert (codes>=0).all()
  return categories[codes]
 return decode(x[:])
with h5py.File(atlas,'r') as h:
 samples=obs_col(h,'sample_uid')
 indexkey=h['obs'].attrs['_index'];ids=decode(h['obs'][indexkey][:])
 x=h['X'];assert x.attrs['encoding-type']=='csr_matrix'
 nvar=x.attrs['shape'][1]
 subset={}
 for sample in targets.sample_uid:
  ix=np.flatnonzero(samples==sample);parts=[];barcodes=[]
  for i in ix:
   start,end=x['indptr'][i:i+2]
   parts.append(sparse.csr_matrix((x['data'][start:end],x['indices'][start:end],np.array([0,end-start])),shape=(1,nvar)))
   # Existing cell IDs use final token after | as the source barcode.
   raw=ids[i].split('|')[-1]
   barcodes.append(raw)
  mat=sparse.vstack(parts,format='csr');mat.sort_indices()
  subset[sample]=(ix,barcodes,mat)
 results=[]
 for bio,g in targets.groupby('biological_sample_uid'):
  assert len(g)==2
  a,b=g.sample_uid.tolist();ia,ba,ma=subset[a];ib,bb,mb=subset[b]
  # Record full unmodified barcode tokens; never strip suffixes to force matches.
  assert len(ba)==len(set(ba)) and len(bb)==len(set(bb))
  aa={v:i for i,v in enumerate(ba)};ab={v:i for i,v in enumerate(bb)};common=sorted(set(aa)&set(ab))
  identical=0;corr=[]
  for barcode in common:
   ra=ma[aa[barcode]].toarray().ravel();rb=mb[ab[barcode]].toarray().ravel()
   identical+=int(np.array_equal(ra,rb))
   if np.std(ra)>0 and np.std(rb)>0:corr.append(float(np.corrcoef(ra,rb)[0,1]))
  results.append(dict(biological_sample_uid=bio,sample_a=a,sample_b=b,n_a=len(ba),n_b=len(bb),barcode_example_a=ba[0] if ba else '',barcode_example_b=bb[0] if bb else '',exact_barcode_overlap=len(common),exact_all_gene_count_matches=identical,median_raw_count_pearson=float(np.median(corr)) if corr else None,interpretation='Barcode/count evidence only; overlap is not by itself proof of duplicate cells across independent captures. No rows deleted or aggregated.'))
pd.DataFrame(results).to_csv(W/'GSE237230_retained_counts_pair_audit.tsv',sep='\t',index=False)
(W/'CELL_PAIR_AUDIT.json').write_text(json.dumps(dict(status='BOUNDED_READ_ONLY_AUDIT_COMPLETE',atlas_sha256=cert['observed_sha256'],pairs=results,cell_rows_read=sum(len(x[0]) for x in subset.values()),gene_columns=int(nvar),atlas_unchanged=(atlas.stat().st_mtime_ns==st.st_mtime_ns),models_run=False),indent=2))
print((W/'CELL_PAIR_AUDIT.json').read_text())
