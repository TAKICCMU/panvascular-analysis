"""Original native output checks and BH, isolated from one-thread coordinator."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['MALLOC_CONF']='background_thread:false'
import argparse,json,gzip,re,hashlib
from pathlib import Path
import numpy as np,pandas as pd

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def bh(x):
 x=np.asarray(x,float);o=np.full(len(x),np.nan);ok=np.isfinite(x);p=x[ok];ix=np.argsort(p)
 if len(p):
  z=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];y=np.empty(len(p));y[ix]=np.minimum(z,1);o[ok]=y
 return o

def score(run):
 files=list((run/'scores').glob('*.full_score.gz'))
 if len(files)!=1:raise ValueError('Expected one full native score output')
 header=pd.read_csv(files[0],sep='\t',nrows=0,index_col=0);controls=[c for c in header if str(c).startswith('ctrl_norm_score_')]
 if len(controls)!=1000:raise ValueError('Expected 1000 normalized controls')
 ids=[]
 with gzip.open(files[0],'rt') as fh:
  next(fh)
  for line in fh:
   fields=line.rstrip('\n').split('\t')
   if len(fields)!=len(header.columns)+1:raise ValueError('Truncated score row')
   ids.append(fields[0])
 if len(ids)!=len(set(ids)):raise ValueError('Duplicate scored IDs')
 pd.Series(ids,name='cell_id').to_csv(run/'native_scored_cell_ids.tsv.gz',index=False,sep='\t')
 text=(run/'compute-score.log').read_text();lines=[s for s in text.splitlines() if any(w in s.lower() for w in ['cells','genes','gene set','gene-set','overlap','n_cell','n_gene','shape'])]
 a=re.search(r'--h5ad-file loaded: n_cell=(\d+), n_gene=(\d+)',text);b=re.search(r'Trait=([^,]+), n_gene=(\d+)',text)
 if not lines or not a or not b or int(a.group(1))!=len(ids):raise ValueError('Native filtering/overlap/count evidence absent or inconsistent')
 (run/'NATIVE_FILTER_AND_OVERLAP_LOG.txt').write_text('\n'.join(lines)+'\n')
 d=dict(status='PASS',native_filtered_n_cells=int(a.group(1)),native_filtered_n_genes=int(a.group(2)),native_trait_overlap=int(b.group(2)),n_scored_cells=len(ids),n_control_columns=len(controls),score_columns=list(header.columns[:8]),score_file=str(files[0]),score_bytes=files[0].stat().st_size,score_sha256=sha(files[0]),native_filter_and_overlap_evidence=lines)
 (run/'SCORE_VALIDATION.json').write_text(json.dumps(d,indent=2));return d

def group(run,key):
 folder=run/key;files=list(folder.glob('*.scdrs_group.*'))
 if len(files)!=1:raise ValueError('Native group output missing')
 z=pd.read_csv(files[0],sep='\t')
 if not {'group','assoc_mcp','hetero_mcp'}<=set(z):raise ValueError('Native group schema differs')
 z['formally_eligible']=~z.group.astype(str).str.startswith('Excluded_');z.to_csv(folder/'all_native_groups_with_scope.tsv',sep='\t',index=False)
 kept=z[z.formally_eligible].copy()
 for p in ['assoc_mcp','hetero_mcp']:
  x=kept[p].to_numpy(float)
  if not np.isfinite(x).all() or ((x<0)|(x>1)).any():raise ValueError('Nonfinite/out-of-range eligible P')
  kept[p+'_BH']=bh(x)
 kept.to_csv(folder/'eligible_native_group_BH.tsv',sep='\t',index=False)
 return dict(n_all_groups=len(z),n_eligible_groups=len(kept),nonfinite_assoc=int((~np.isfinite(z.assoc_mcp)).sum()),nonfinite_hetero=int((~np.isfinite(z.hetero_mcp)).sum()))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['score','group']);p.add_argument('--run',required=True);p.add_argument('--key');a=p.parse_args();run=Path(a.run)
 d=score(run) if a.mode=='score' else group(run,a.key)
 if a.mode=='group':(run/(a.key+'_VALIDATION.json')).write_text(json.dumps(d,indent=2))
 print(json.dumps(d))
