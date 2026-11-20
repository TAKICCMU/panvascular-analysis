"""Full sparse spatial content and exact identifier audit; no atlas counts read."""
import os,json,hashlib,time,traceback
from pathlib import Path
import h5py,numpy as np,pandas as pd
from anndata._io.specs import read_elem
from scipy.io import mmread
ROOT=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
JOB=Path(__file__).resolve().parents[1]
ATLAS=ROOT/'panvascular_CPU_theta2_L1L2_final_v1.h5ad'
MAN=Path('/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v3/00_inputs/spatial_sample_manifest_curated.tsv')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def positions(p):
 a=pd.read_csv(p,header=None,nrows=1)
 t=pd.read_csv(p) if str(a.iloc[0,0]).lower() in ('barcode','barcodes') else pd.read_csv(p,header=None)
 t=t.iloc[:,:6];t.columns=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres'];t['barcode']=t.barcode.astype(str)
 return t
with h5py.File(ATLAS,'r') as f:v=read_elem(f['var'])
ids=v.stable_ensembl_id.astype('string');symbols=v.display_gene_symbol.astype('string')
id_ok=ids.notna()&~ids.isin(['','NR','nan'])&~ids.duplicated(keep=False)
sym_ok=symbols.notna()&~symbols.isin(['','NR','nan'])&~symbols.duplicated(keep=False)
idmap={str(x):str(i) for i,x in ids[id_ok].items()};symmap={str(x):str(i) for i,x in symbols[sym_ok].items()}
man=pd.read_csv(MAN,sep='\t');results=[];files={str(MAN):sha(MAN)}
for _,r in man.iterrows():
 out={'section_id':r.section_id,'GSE':r.GSE,'status':'PASS'}
 try:
  for k in ['matrix','barcodes','features','positions']:files[str(r[k])]=sha(r[k])
  m=mmread(r.matrix).tocsr();f=pd.read_csv(r.features,sep='\t',header=None,dtype=str,keep_default_na=False);b=pd.read_csv(r.barcodes,sep='\t',header=None,dtype=str,keep_default_na=False).iloc[:,0];p=positions(r.positions)
  assert m.shape==(len(f),len(b)), 'matrix axes mismatch'
  assert np.isfinite(m.data).all() and (m.data>=0).all() and np.equal(m.data,np.floor(m.data)).all(),'nonraw counts'
  assert b.is_unique and p.barcode.is_unique,'duplicate barcode'
  assert b.isin(p.barcode).all(),'barcodes without exact positions; no suffix fallback'
  sourceid=f.iloc[:,0];symbol=f.iloc[:,1] if f.shape[1]>1 else sourceid
  # Stable IDs dominate. A symbol fallback is allowed only for non-Ensembl IDs,
  # exact unique atlas symbol and consistent source identity; never casefold/suffix.
  rows=[]
  for j,(i,s) in enumerate(zip(sourceid,symbol)):
   gene=idmap.get(i);method='exact_stable_ensembl_id'
   if gene is None and not str(i).startswith('ENSG'):
    gene=symmap.get(s);method='exact_unique_symbol_nonensembl_feature'
   rows.append({'source_row':j,'source_id':i,'source_symbol':s,'atlas_gene':gene or '', 'method':method if gene else 'unmatched'})
  gm=pd.DataFrame(rows)
  hit=gm[gm.atlas_gene!=''];amb=[]
  for gene,t in hit.groupby('atlas_gene'):
   if len(t)>1 and t.source_id.nunique()>1:amb.append(gene)
  assert not amb, 'ambiguous different source IDs map to same atlas gene: '+repr(amb[:20])
  assert len(hit)>1000,'too few exact genes'
  gm.to_csv(JOB/'evidence'/f'{r.GSE}__{str(r.section_id).replace(":", "_")}__gene_map.tsv.gz',sep='\t',index=False)
  in_tissue=p.set_index('barcode').loc[b,'in_tissue'].astype(int).to_numpy()==1
  x=m.T.tocsr()[in_tissue];n=np.asarray(x.sum(1)).ravel();g=np.diff((x>0).tocsr().indptr)
  keep=(n>=(500 if m.shape[0]>10000 else 100))&(g>=(200 if m.shape[0]>10000 else 50));fallback=int(keep.sum())<50
  if fallback:keep=n>0
  out.update(n_spots=len(b),n_genes=len(f),nnz=int(m.nnz),mapped_genes=hit.atlas_gene.nunique(),n_unmatched=int((gm.atlas_gene=='').sum()),n_same_id_duplicate_rows=int(hit.atlas_gene.duplicated().sum()),in_tissue=int(in_tissue.sum()),qc_retained=int(keep.sum()),qc_fallback=fallback,raw_matrix_bytes=int(m.data.nbytes+m.indices.nbytes+m.indptr.nbytes))
  print(json.dumps(out),flush=True)
 except Exception as e:out.update(status='BLOCKED_INPUT',error=repr(e));print(traceback.format_exc(),flush=True)
 results.append(out)
pd.DataFrame(results).to_csv(JOB/'evidence'/'SPATIAL_CONTENT_VALIDATION.tsv',sep='\t',index=False)
(JOB/'evidence'/'SPATIAL_INPUT_HASHES.json').write_text(json.dumps(files,indent=2))
(JOB/'evidence'/'INPUT_VALIDATION.json').write_text(json.dumps({'status':'PASS' if all(x['status']=='PASS' for x in results) else 'PARTIAL_BLOCKED_INPUT','sections':len(results),'passed':sum(x['status']=='PASS' for x in results),'raw_counts_read':True,'exact_barcode_matching':True,'gene_rule':'stable ID preferred; exact unique symbol only for non-Ensembl feature IDs; no casefold/version stripping','atlas_var_sha256':hashlib.sha256(v.to_csv().encode()).hexdigest(),'files_sha256':sha(JOB/'evidence'/'SPATIAL_INPUT_HASHES.json')},indent=2))
