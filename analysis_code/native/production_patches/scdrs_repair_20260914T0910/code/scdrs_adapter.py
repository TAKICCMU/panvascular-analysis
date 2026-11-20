"""Shared smoke/production scDRS adapter; no launch or scientific work on import."""
from pathlib import Path
import hashlib,json,os,subprocess,time,copy,re,gzip
import numpy as np,pandas as pd
from scipy import sparse

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()

def inspect_integer_values(data):
 maximum=0
 for start in range(0,len(data),1000000):
  x=data[start:start+1000000]
  if not np.isfinite(x).all() or (x<0).any() or not np.equal(x,np.floor(x)).all():raise ValueError('Counts must be finite nonnegative integers')
  if len(x):maximum=max(maximum,int(x.max()))
 return maximum

def prepare_scdrs_counts(X,gene_map,feature_names):
 """Reorder unique source symbols without float64 projection. Verified duplicates sum in int64.
 For duplicates, same nonmissing stable Ensembl ID required; no symbol-alias inference.
 Float32 output only after exact integer round-trip check; otherwise preserve int64.
 """
 gm=gene_map.copy()
 if gm.feature_id.astype(str).tolist()!=list(map(str,feature_names)):raise ValueError('Gene-map feature order differs from counts')
 symbol_column='scdrs_gene_symbol' if 'scdrs_gene_symbol' in gm else 'gene_symbol'
 names=gm[symbol_column]
 if names.isna().any() or names.astype(str).str.strip().isin(['','NR','nan','None','not_available']).any():raise ValueError('Missing symbol mapping')
 names=names.astype(str);X=sparse.csr_matrix(X);max_value=inspect_integer_values(X.data)
 codes,symbols=pd.factorize(names,sort=True);sizes=np.bincount(codes)
 duplicates=gm[names.duplicated(keep=False)].copy()
 if len(duplicates):
  for symbol,g in gm.groupby(symbol_column,sort=False):
   if len(g)<2:continue
   if 'ensembl_id' not in g:raise ValueError('Duplicate symbols require same verified gene identity')
   ids=g.ensembl_id.fillna('').astype(str).str.replace(r'\.\d+$','',regex=True)
   if ids.nunique()!=1 or not ids.iloc[0].startswith('ENSG'):raise ValueError('Ambiguous duplicate symbol: '+str(symbol))
  if max_value*int(sizes.max())>np.iinfo(np.int64).max:raise ValueError('Potential integer sum overflow')
  P=sparse.csr_matrix((np.ones(len(codes),dtype=np.int64),(np.arange(len(codes)),codes)),shape=(len(codes),len(symbols)))
  result=X.astype(np.int64)@P
 else:
  order=np.argsort(codes);result=X if np.array_equal(order,np.arange(len(codes))) else X[:,order].tocsr()
 # Do not allocate a float64 expression/projection matrix. Exact check is block bounded.
 f32_safe=True
 for start in range(0,len(result.data),1000000):
  z=result.data[start:start+1000000]
  if not np.array_equal(z,z.astype(np.float32).astype(z.dtype)):
   f32_safe=False;break
 if f32_safe:result=result.astype(np.float32,copy=False)
 elif np.issubdtype(result.dtype,np.integer):result=result.astype(np.int64,copy=False)
 # Already floating counts cannot become inexact here: retain the exact original dtype.
 audit={'input_dtype':str(X.dtype),'output_dtype':str(result.dtype),'input_nnz':int(X.nnz),'output_nnz':int(result.nnz),
  'n_input_features':X.shape[1],'n_output_symbols':len(symbols),'n_duplicate_symbol_features':len(duplicates),
  'float32_integer_roundtrip_exact':f32_safe,'source_symbol_only_features':int(gm.get('ensembl_id',pd.Series('',index=gm.index)).fillna('').astype(str).eq('').sum()),
  'symbol_semantics':'unique source symbols retained; not automatically canonical HGNC; duplicate aggregation only same verified gene'}
 return result,pd.Index(symbols,name='gene_symbol'),audit

def compatible_obs(obs):
 """Lossless values in native-readable ordinary NumPy/categorical encodings; no imputation."""
 out=obs.copy();audit=[]
 for col in out:
  s=out[col];original=str(s.dtype);missing=int(s.isna().sum())
  if pd.api.types.is_numeric_dtype(s.dtype):
   if pd.api.types.is_extension_array_dtype(s.dtype):out[col]=s.to_numpy(dtype=float,na_value=np.nan)
  else:
   # Non-numeric missing values remain NA in a standard categorical encoding.
   vals=s.astype(object);vals=vals.where(pd.notna(vals),None)
   values=[None if v is None else str(v) for v in vals]
   out[col]=pd.Categorical(values)
  audit.append({'column':col,'source_dtype':original,'derived_dtype':str(out[col].dtype),'n_missing_before':missing,'n_missing_after':int(out[col].isna().sum())})
  if missing!=int(out[col].isna().sum()):raise ValueError('Obs missingness changed')
 return out,pd.DataFrame(audit)

def build_derived(adata,gm,cov,output,graph_key='neighbors_cpu_n30',min_cells=500,min_studies=3,min_study_cells=20):
 import anndata as ad
 output=Path(output)
 if output.exists():raise FileExistsError(output)
 output.parent.mkdir(parents=True,exist_ok=True)
 X,names,dtype_audit=prepare_scdrs_counts(adata.X[:],gm,adata.var_names)
 obs,encoding_audit=compatible_obs(adata.obs)
 # Preserve both original formal labels. Qualified state is a separate derived test label.
 l1='level1_annotation';l2='level2_annotation';parent=obs[l1].astype(str);child=obs[l2].astype(str)
 valid=~parent.isin(['Unknown','Not_assessed']) & ~child.isin(['Unknown','Not_assessed'])
 obs['_scdrs_state']=(parent+' | '+child).where(valid,'Excluded_Not_assessed')
 b=ad.AnnData(X=X,obs=obs,var=pd.DataFrame(index=names));group_keys=[]
 for field,label in [(l1,l1),(l2,'_scdrs_state')]:
  support=b.obs.groupby([label,'study'],observed=True).size().unstack(fill_value=0)
  eligible=support.index[(support.sum(axis=1)>=min_cells)&((support>=min_study_cells).sum(axis=1)>=min_studies)]
  eligible=[g for g in eligible if str(g) not in ['Unknown','Not_assessed','Excluded_Not_assessed']]
  key=field+'_scdrs_test';group_keys.append(key)
  b.obs[key]=b.obs[label].astype(str).where(b.obs[label].isin(eligible),'Excluded_Not_eligible')
  support.to_csv(output.parent/f'{field}_study_support.tsv',sep='\t')
  pd.DataFrame({'group':support.index,'n_cells':support.sum(axis=1),'n_studies_ge20':(support>=min_study_cells).sum(axis=1),'eligible':support.index.isin(eligible)}).to_csv(output.parent/f'{field}_eligibility.tsv',sep='\t',index=False)
 ng=adata.uns[graph_key];ck=ng['connectivities_key'];dk=ng['distances_key']
 if int(ng['params']['n_neighbors'])!=30:raise ValueError('Frozen graph must use approved n30')
 for key,new in [(ck,'connectivities'),(dk,'distances')]:
  graph=adata.obsp[key]
  if graph.shape!=(b.n_obs,b.n_obs) or not np.isfinite(graph.data).all():raise ValueError('Frozen graph identity/values invalid')
  b.obsp[new]=graph.copy()
 b.uns['neighbors']={'connectivities_key':'connectivities','distances_key':'distances','params':copy.deepcopy(ng['params'])}
 if not b.obs_names.equals(adata.obs_names) or not cov.index.equals(b.obs_names):raise ValueError('Cell/covariate order mismatch')
 if not np.isfinite(cov.to_numpy(dtype=float)).all():raise ValueError('Invalid covariates')
 cov.to_csv(output.parent/'native_covariates.cov.gz',sep='\t',compression='gzip')
 b.uns['scdrs_adapter']={'graph_policy':'frozen_global_n30','graph_source_uns':graph_key,'annotation_source':'final formal L1/L2 unchanged','group_thresholds':{'min_cells':min_cells,'min_studies':min_studies,'min_study_cells':min_study_cells}}
 b.write_h5ad(output,compression='lzf')
 encoding_audit.to_csv(output.parent/'obs_encoding_audit.tsv',sep='\t',index=False)
 gm.to_csv(output.parent/'gene_mapping_used.tsv.gz',sep='\t',index=False)
 audit=dict(dtype_audit,shape=list(b.shape),n_obs_columns=len(b.obs.columns),graph_policy='frozen_global_n30',connectivities_nnz=int(b.obsp['connectivities'].nnz),distances_nnz=int(b.obsp['distances'].nnz),group_keys=group_keys,
  derived_h5ad_sha256=sha(output),covariates_sha256=sha(output.parent/'native_covariates.cov.gz'),derived_bytes=output.stat().st_size)
 (output.parent/'DERIVED_OBJECT.json').write_text(json.dumps(audit,indent=2))
 return audit

def bh(x):
 x=np.asarray(x,float);o=np.full(len(x),np.nan);ok=np.isfinite(x);p=x[ok];ix=np.argsort(p)
 if len(p):
  z=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];y=np.empty(len(p));y[ix]=np.minimum(z,1);o[ok]=y
 return o

def native_sequence(h5ad,cov,scenarios,output,executable,group_keys,threads,n_ctrl=1000,validate_existing=False,resource_check=None):
 """Same function used by real-obs smoke and full run; all3 traits strictly serial."""
 if threads not in [2,12] or n_ctrl!=1000:raise ValueError('Only approved thread budgets and n_ctrl1000')
 output=Path(output);output.mkdir(parents=True,exist_ok=validate_existing);records=[]
 env=os.environ.copy();env['CUDA_VISIBLE_DEVICES']=''
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]=str(threads)
 def call(cmd,log):
  t=time.time()
  with open(log,'w') as f:r=subprocess.run(list(map(str,cmd)),env=env,stdout=f,stderr=subprocess.STDOUT)
  records.append({'cmd':list(map(str,cmd)),'returncode':r.returncode,'seconds':time.time()-t,'log':str(log)})
  (output/'COMMANDS.json').write_text(json.dumps(records,indent=2))
  if r.returncode:raise RuntimeError(f'Native failure {r.returncode}; {log}')
 for trait in scenarios:
  if resource_check is not None:resource_check()
  (output/'CURRENT_STAGE.json').write_text(json.dumps({'trait':trait['name'],'status':'COMPUTE_SCORE_START','time':time.time()}))
  if sha(trait['geneset'])!=trait['sha256']:raise ValueError('Gene-set SHA differs')
  run=output/trait['name'];run.mkdir(exist_ok=validate_existing);scores=run/'scores';scores.mkdir(exist_ok=validate_existing)
  cmd=[executable,'compute-score','--h5ad-file',h5ad,'--h5ad-species','human','--gs-file',trait['geneset'],'--gs-species','human','--cov-file',cov,'--ctrl-match-opt','mean_var','--weight-opt','vs','--flag-filter-data','True','--flag-raw-count','True','--n-ctrl','1000','--flag-return-ctrl-raw-score','False','--flag-return-ctrl-norm-score','True','--out-folder',scores]
  if not validate_existing:call(cmd,run/'compute-score.log')
  files=list(scores.glob('*.full_score.gz'))
  if len(files)!=1:raise ValueError('Expected one full native score output')
  header=pd.read_csv(files[0],sep='\t',nrows=0,index_col=0)
  controls=[c for c in header if str(c).startswith('ctrl_norm_score_')]
  if len(controls)!=1000:raise ValueError(f'Expected1000 normalized controls, found{len(controls)}')
  score_ids=[]
  with gzip.open(files[0],'rt') as fh:
   next(fh)
   for line in fh:
    fields=line.rstrip('\n').split('\t')
    if len(fields)!=len(header.columns)+1:raise ValueError('Truncated score row')
    score_ids.append(fields[0])
  if len(score_ids)!=len(set(score_ids)):raise ValueError('Duplicate scored cell IDs')
  pd.Series(score_ids,name='cell_id').to_csv(run/'native_scored_cell_ids.tsv.gz',index=False,sep='\t')
  text=(run/'compute-score.log').read_text()
  native_lines=[line for line in text.splitlines() if any(w in line.lower() for w in ['cells','genes','gene set','gene-set','overlap','n_cell','n_gene','shape'])]
  (run/'NATIVE_FILTER_AND_OVERLAP_LOG.txt').write_text('\n'.join(native_lines)+'\n')
  if not native_lines:raise ValueError('Native filtering/overlap evidence absent')
  shape_match=re.search(r'--h5ad-file loaded: n_cell=(\d+), n_gene=(\d+)',text)
  overlap_match=re.search(r'Trait=([^,]+), n_gene=(\d+)',text)
  if not shape_match or not overlap_match:raise ValueError('Native shape/trait overlap absent')
  if int(shape_match.group(1))!=len(score_ids):raise ValueError('Score rows differ from native filtered cell count')
  schema={'native_filtered_n_cells':int(shape_match.group(1)),'native_filtered_n_genes':int(shape_match.group(2)),'native_trait_overlap':int(overlap_match.group(2)),'n_scored_cells':len(score_ids),'native_filter_and_overlap_evidence':native_lines,'trait':trait['name'],'n_control_columns':len(controls),'score_columns':list(header.columns[:8]),'score_bytes':files[0].stat().st_size}
  for key in group_keys:
   group=run/key;group.mkdir(parents=True,exist_ok=validate_existing)
   cmd=[executable,'perform-downstream','--h5ad-file',h5ad,'--score-file',files[0],'--out-folder',group,'--group-analysis',key,'--flag-filter-data','True','--flag-raw-count','True']
   if not validate_existing:call(cmd,run/(key+'.log'))
   f=list(group.glob('*.scdrs_group.*'))
   if len(f)!=1:raise ValueError('Native group output missing')
   z=pd.read_csv(f[0],sep='\t')
   if not {'group','assoc_mcp','hetero_mcp'}<=set(z):raise ValueError('Native group schema differs')
   z['formally_eligible']=~z.group.astype(str).str.startswith('Excluded_')
   z.to_csv(group/'all_native_groups_with_scope.tsv',sep='\t',index=False)
   kept=z[z.formally_eligible].copy()
   for p in ['assoc_mcp','hetero_mcp']:
    values=kept[p].to_numpy(float)
    if not np.isfinite(values).all() or ((values<0)|(values>1)).any():raise ValueError('Nonfinite/out-of-range eligible native '+p)
    kept[p+'_BH']=bh(values)
   kept.to_csv(group/'eligible_native_group_BH.tsv',sep='\t',index=False)
   schema[key]={'n_all_groups':len(z),'n_eligible_groups':len(kept),'nonfinite_assoc':int((~np.isfinite(z.assoc_mcp)).sum()),'nonfinite_hetero':int((~np.isfinite(z.hetero_mcp)).sum())}
  (run/'SCHEMA_VALIDATION.json').write_text(json.dumps(schema,indent=2))
  (output/'CURRENT_STAGE.json').write_text(json.dumps({'trait':trait['name'],'status':'TRAIT_COMPLETE','time':time.time()}))
 return records
