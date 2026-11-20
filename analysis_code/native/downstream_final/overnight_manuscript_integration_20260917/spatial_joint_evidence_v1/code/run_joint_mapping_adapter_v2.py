import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','POLARS_MAX_THREADS','RAYON_NUM_THREADS','ARROW_NUM_THREADS','BLOSC_NTHREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import json,ast,time,hashlib,fcntl,traceback,sys,gc,warnings,math
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--mapping-root',type=Path,default=None)
parser.add_argument('--output-root',type=Path)
parser.add_argument('--dispatcher-admitted',action='store_true')
parser.add_argument('--self-test',action='store_true')
parser.add_argument('--factor-path',type=Path,default=None)
args=parser.parse_args()
BASE=Path(__file__).resolve().parents[1];W=BASE.parent.parent;N=W/'remaining_downstream_20260916_v1';O=args.output_root
def verify_dispatcher_parent():
 expected=str((N/'code/dispatch.py').resolve())
 cmd=[x.decode() for x in Path(f'/proc/{os.getppid()}/cmdline').read_bytes().split(b'\0') if x]
 if expected not in cmd:raise RuntimeError('dispatcher-admitted requires actual direct parent '+expected+'; got '+repr(cmd))
 return {'pid':os.getppid(),'cmdline':cmd,'expected_dispatcher':expected}
if args.self_test:
 try:verify_dispatcher_parent()
 except RuntimeError:print('PASS: non-dispatcher actual parent rejected')
 else:raise AssertionError('Self-test must be launched outside dispatcher')
 print('PASS: parser and admission checks loaded without heavy imports');sys.exit(0)
if O is None:raise RuntimeError('--output-root required')
parent_admission=verify_dispatcher_parent() if args.dispatcher_admitted else None
if O.exists() and list(O.iterdir()):raise RuntimeError('Output directory must be new/empty; will not overwrite another branch')
O.mkdir(exist_ok=True,parents=True)

def dump(name,x): (O/name).write_text(json.dumps(x,indent=2,ensure_ascii=False,default=str))
def actual_threads():
 total=0
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   c=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode()
   if str(W.parent) not in c or any(x in c for x in ['python3 -','bash -c','sh -c','sshd:','scp ']):continue
   s=dict(x.split(':',1) for x in (p/'status').read_text().splitlines() if ':' in x)
   if s['State'].strip().startswith(('T','Z')):continue
   total+=int(s['Threads'])
  except (OSError,UnicodeError,KeyError):pass
 return total
lock=None
while lock is None and not args.dispatcher_admitted:
 n=actual_threads()
 if n+1<=40:
  for i in range(4):
   h=open(W/'locks'/f'cpu_heavy_{i}.lock','a')
   try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB);lock=h;break
   except BlockingIOError:h.close()
 if lock is None:
  dump('STATUS.json',{'status':'WAITING_RESOURCES','actual_threads':n,'pid':os.getpid()});time.sleep(10)
dump('STATUS.json',{'status':'RUNNING','pid':os.getpid(),'lock':lock.name if lock else 'owned_by_verified_dispatcher_parent','dispatcher_admission':parent_admission,'threads':1,'started':time.time()})
import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
from scipy.io import mmread
from scipy.stats import norm,spearmanr
import statsmodels.api as sm
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
threadpool_limits(1)
K_NEIGHBORS=6
src=BASE/'code/legacy_spatial.py'; tree=ast.parse(src.read_text()); nodes=[]
for node in tree.body:
 if isinstance(node,ast.FunctionDef):nodes.append(node)
 if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PATHOLOGY_SIGNATURES' for t in node.targets):nodes.append(node)
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(src),'exec'),globals())
# Explicit-query nearest neighbours removes actual self; inherited kneighbors(None) already excludes self and then wrongly drops one nearest neighbour.
def neighbour_lag(values,coords,k=6):
 n=len(coords); ids=NearestNeighbors(n_neighbors=min(k+1,n)).fit(coords).kneighbors(coords,return_distance=False)
 return np.vstack([np.mean(values[np.array([j for j in row if j!=i][:k])],axis=0) for i,row in enumerate(ids)])
hashes={}
def sha(p):
 p=Path(p)
 if str(p) not in hashes:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(8<<20),b''):h.update(b)
  hashes[str(p)]=h.hexdigest()
 return hashes[str(p)]
def save(df,name):df.to_csv(O/name,sep='\t',index=False)
try:
 manifestpath=Path('/home/anzhen_vhd/scRNAseq/codex/spatial_validation_v3/00_inputs/spatial_sample_manifest_curated.tsv')
 ds_path=W/'spatial_mapping_final_review_20260916_v1/tables/dataset_summary.tsv'
 targetpath=BASE.parent/'candidate_prioritization/Tier1_spatial_interface.tsv';wtpath=W/'covarnet/stable_NMF/module_node_weights_K14.tsv';nodepath=W/'covarnet/CoVarNet_node_parent_map.tsv';factorpath=W/'spatial_gpu_resume_v2/outputs/reference_factors.tsv';tierpath=BASE.parent/'covarnet_final_tiers/14_final_module_evidence_and_tiers.tsv'
 if args.factor_path is not None:factorpath=args.factor_path
 for p in [src,Path(__file__),manifestpath,ds_path,targetpath,wtpath,nodepath,factorpath,tierpath]:sha(p)
 if args.mapping_root is not None:
  newroot=args.mapping_root.parent
  vpath=newroot/'validation_v2/VALIDATION.json';validation=json.loads(vpath.read_text())
  assert validation.get('status')=='PASS' and validation.get('reference_retrained') is True, 'New mapping integrity validation is not PASS'
  assert validation.get('n_targets')==29 and validation.get('n_datasets')==7 and validation.get('n_sections')==49
  assert json.loads((newroot/'RUN_STATUS.json').read_text()).get('status')=='COMPLETE', 'Production incomplete'
  receiptpath=newroot/'configs/PRODUCTION_RECEIPT.json';receipt=json.loads(receiptpath.read_text())
  for bp,bh in receipt['bound_files'].items():assert sha(bp)==bh,'New mapping input changed: '+bp
  rcpath=args.mapping_root/'reference_COMPLETE.json';rc=json.loads(rcpath.read_text())
  assert rc['reference_retrained'] is True and rc['target_excluded_before_model_setup'] is True
  assert sha(args.mapping_root/'reference_signatures.tsv.gz')==rc['signature_sha']
  assert sha(receiptpath)==rc['receipt_sha']
  assert sha(newroot/'code/spatial_target_excluded.py')==rc['code_sha']
  cnew=json.loads((newroot/'configs/config.json').read_text())
  assert set(cnew['target_ids'])==set(pd.read_csv(targetpath,sep='\t').gene_id.astype(str)), 'New excluded target identities differ'
  sha(vpath);sha(rcpath);sha(newroot/'validation_v2/file_hashes_and_bytes.tsv')
 man=pd.read_csv(manifestpath,sep='\t');ds=pd.read_csv(ds_path,sep='\t');targets=pd.read_csv(targetpath,sep='\t');factors=pd.read_csv(factorpath,sep='\t');fmap=dict(zip(factors.reference_state,factors.reference_group)); nm=pd.read_csv(nodepath,sep='\t');nm['label']=nm.parent_lineage+' | '+nm.level2_annotation
 weights=pd.read_csv(wtpath,sep='\t').set_index('cell_state');weights.index=weights.index.map(nm.set_index('cell_state').label);assert weights.index.notna().all();weights=weights[['CM09']]
 assert len(targets)==29 and targets[['gene_symbol','cell_state']].drop_duplicates().shape[0]==29
 signatures=PATHOLOGY_SIGNATURES;genes=sorted(set(targets.gene_symbol)|set(sum(signatures.values(),[])))
 rows=[];detection=[];audit=[];module_rows=[];spot_parts=[];coverage=[];match=[];rng=np.random.default_rng(20260810)
 for _,d in ds.iterrows():
  path=Path(d.posterior_path)
  if args.mapping_root is not None:
   path=args.mapping_root/f'{d.GSE}_posterior.h5ad'
   provenance=args.mapping_root/f'{d.GSE}_COMPLETE.json'
   assert provenance.exists(), 'New mapping requires its own per-dataset provenance'
   sha(provenance);prov=json.loads(provenance.read_text());assert sha(path)==prov['posterior_sha'];assert prov['reference_sha']==rc['signature_sha'] and prov['receipt_sha']==rc['receipt_sha'];assert prov['excluded_targets_present']==[]
  else:assert sha(path)==d.posterior_sha256
  print('LOAD',d.GSE,flush=True);a=ad.read_h5ad(path);assert a.n_obs==d.spots
  ab={}
  for post,key in [('q05','q05_cell_abundance_w_sf'),('mean','means_cell_abundance_w_sf')]:
   z=a.obsm[key].copy();z.columns=[fmap.get(strip_cell2location_prefix(x),strip_cell2location_prefix(x)) for x in z.columns];assert z.shape[1]==128 and z.columns.is_unique and z.index.equals(a.obs_names)
   assert set(weights.index)<=set(z.columns);ab[post]=z
  for _,mr in man[man.GSE==d.GSE].iterrows():
   if str(mr.include_mapping).lower() not in ['true','1','yes']:continue
   section=str(mr.section_id);idx=np.flatnonzero(a.obs.section_id.astype(str).to_numpy()==section);assert len(idx)>0;obs=a.obs.iloc[idx].copy()
   assert obs.biological_unit_id.astype(str).eq(str(mr.biological_unit_id)).all()
   gmfile=W/'spatial_gpu_resume_v2/evidence'/f'{d.GSE}__{section.replace(":","_")}__gene_map.tsv.gz'; assert gmfile.exists();sha(gmfile)
   for p in [mr.matrix,mr.barcodes,mr.features,mr.positions]:sha(p)
   gm=pd.read_csv(gmfile,sep='\t',keep_default_na=False);mx=mmread(mr.matrix).tocsr().T;bars=pd.read_csv(mr.barcodes,sep='\t',header=None,dtype=str).iloc[:,0];assert len(bars)==mx.shape[0]
   assert np.isfinite(mx.data).all() and (mx.data>=0).all() and np.max(np.abs(mx.data-np.rint(mx.data)))<1e-5
   rawids=pd.Index(section+'::'+bars.astype(str));sel=rawids.get_indexer(obs.index);assert (sel>=0).all() and len(set(sel))==len(sel)
   lib=np.asarray(mx.sum(axis=1)).ravel()[sel];assert np.allclose(lib,obs.total_counts_original.to_numpy(),rtol=0,atol=0)
   Y=np.zeros((len(sel),len(genes)));present={}
   for j,g in enumerate(genes):
    columns=gm.loc[gm.atlas_gene==g,'source_row'].to_numpy(int);present[g]=len(columns)>0
    if len(columns):Y[:,j]=np.asarray(mx[sel][:,columns].sum(axis=1)).ravel()
    else:Y[:,j]=np.nan
   expr=pd.DataFrame(np.log1p(Y*1e4/np.maximum(lib,1)[:,None]),index=obs.index,columns=genes)
   fallback=bool(obs.QC_fallback.astype(bool).any());blocks=spatial_blocks(obs);assert {'array_row','array_col'}<=set(obs) and len(set(blocks))>=2
   audit.append(dict(GSE=d.GSE,section_id=section,biological_unit_id=mr.biological_unit_id,metadata_status=mr.metadata_status,n_spots=len(idx),QC_fallback=fallback,n_spatial_blocks=len(set(blocks)),raw_integer_counts=True,original_full_library_verified=True,min_library=float(min(lib)),max_library=float(max(lib)),raw_genes=mx.shape[1],posterior_genes=a.n_vars,n_targets_present=sum(present[g] for g in targets.gene_symbol)))
   meta=dict(GSE=d.GSE,section_id=section,biological_unit_id=str(mr.biological_unit_id),metadata_status=mr.metadata_status,QC_fallback=fallback)
   # Signature genes scored within section. Each target is removed from its own regional program to avoid direct mathematical self-inclusion.
   sigscore={}
   for s,gs in signatures.items():
    valid=[g for g in gs if present[g]];sigscore[s]=np.mean(np.column_stack([zscore(expr[g]) for g in valid]),axis=1) if len(valid)>=3 else np.full(len(idx),np.nan)
    coverage.append({**meta,'signature':s,'n_valid':len(valid),'genes':';'.join(valid)})
   for post in ['q05','mean']:
    abundance=ab[post].iloc[idx];lineage=lineage_abundance(abundance);scores=module_scores(abundance,weights)
    for _,t in targets.iterrows():
     gene=t.gene_symbol;ct=t.cell_state;signature={'Fibroblast':'fibrous_cap_core','Endothelial cell':'endothelial_activation','Macrophage':'lipid_foam'}[ct]
     gs=[g for g in signatures[signature] if present[g] and g!=gene];pat=np.mean(np.column_stack([zscore(expr[g]) for g in gs]),axis=1) if len(gs)>=3 else np.full(len(idx),np.nan)
     y=Y[:,genes.index(gene)];dm={**meta,'posterior':post,'gene':gene,'cell_context':ct,'pathology':signature,'expected_direction':int(np.sign(t.meta_logFC)),'target_removed_from_program':gene in signatures[signature],'n_program_genes':len(gs),'n_spots':len(idx),'n_positive_spots':int(np.sum(y>0)),'gene_present_in_source':present[gene]}
     eligible=present[gene] and len(gs)>=3 and len(idx)>=100 and np.sum(y>0)>=20 and ct in lineage
     rr=fit_spatial_count_model(y,lib,lineage[ct].to_numpy(),pat,blocks) if eligible else []
     good=rr and all(np.isfinite(v['beta']) and np.isfinite(v['se']) and v['se']>0 for v in rr)
     dm['model_status']='FIT_OK' if good else ('FIT_FAILED' if eligible else 'INPUT_GATE_NOT_MET');detection.append(dm)
     if good:
      for v in rr:rows.append({**dm,**v})
    co=coordinates(a,idx);logv=np.log1p(abundance.to_numpy());corr=matrix_correlations(logv);cross=spatial_cross_matrix(logv,co);ordered=weights.CM09.sort_values(ascending=False);top=[x for x in ordered.index if ordered[x]>0][:10];inds=[abundance.columns.get_loc(x) for x in top]
    within=upper_mean(corr,inds);lag=upper_mean(cross,inds);nullw=[];nulll=[]
    for k in range(1000):
     st=sample_matched_null(top,list(abundance.columns),rng)
     if len(st)==len(top):
      ni=[abundance.columns.get_loc(x) for x in st];nullw.append(upper_mean(corr,ni));nulll.append(upper_mean(cross,ni))
    mod={**meta,'posterior':post,'module':'CM09','n_top_nodes':len(top),'top_nodes':';'.join(top),'within_spot_coherence':within,'spatial_lag_coherence':lag,'within_spot_empirical_p':(1+np.sum(np.array(nullw)>=within))/(1+len(nullw)),'spatial_lag_empirical_p':(1+np.sum(np.array(nulll)>=lag))/(1+len(nulll)),'n_null':len(nullw),'mapped_weight_fraction':1.,'moran_I':moran_i(scores.CM09.to_numpy(),co)}
    for s,vals in sigscore.items():mod['rho__'+s]=float(spearmanr(scores.CM09,vals).statistic) if np.isfinite(vals).all() else np.nan
    module_rows.append(mod)
    if post=='q05':
     spdf=expr.add_prefix('expr__');spdf['CM09']=scores.CM09;spdf['section_id']=section;spdf['biological_unit_id']=mr.biological_unit_id;spdf['library_size_original']=lib;spdf['spot_id']=obs.index;spdf['x']=co[:,0];spdf['y']=co[:,1]
     for s,val in sigscore.items():spdf['program__'+s]=val
     for ct in targets.cell_state.unique():spdf['abundance__'+ct]=lineage[ct].to_numpy()
     spot_parts.append(spdf)
   print('DONE_SECTION',section,len(idx),flush=True);save(pd.DataFrame(audit),'section_input_audit.tsv');save(pd.DataFrame(detection),'target_model_roster.tsv');save(pd.DataFrame(rows),'target_section_models.tsv');save(pd.DataFrame(module_rows),'CM09_section_projection.tsv');dump('SOURCE_SHA256.json',hashes)
   del mx;gc.collect()
  del a;gc.collect()
 raw=pd.DataFrame(rows);unitrows=[];metarows=[]
 for sensitivity in ['all_sections','exclude_QC_fallback']:
  rr=raw if sensitivity=='all_sections' else raw[~raw.QC_fallback]
  for key,x in rr.groupby(['posterior','gene','cell_context','pathology','term','GSE','biological_unit_id']):
   b,se=fixed_effect(x.beta,x.se);unitrows.append(dict(zip(['posterior','gene','cell_context','pathology','term','GSE','biological_unit_id'],key))|{'sensitivity':sensitivity,'beta':b,'se':se,'n_sections':len(x),'serial_section_covariance_unavailable':len(x)>1})
 unit=pd.DataFrame(unitrows)
 for key,x in unit.groupby(['posterior','gene','cell_context','pathology','term','sensitivity']):
  z=random_effects(x.beta,x.se);z.update(dict(zip(['posterior','gene','cell_context','pathology','term','sensitivity'],key)));z['n_datasets']=x.GSE.nunique();z['CI_low']=z['beta']-1.96*z['se'];z['CI_high']=z['beta']+1.96*z['se'];metarows.append(z)
 meta=pd.DataFrame(metarows);meta['FDR_BH_all_29_pairs_x_3_terms']=np.nan
 # All 87 prespecified pair-term tests are family denominator, including NA/insufficient tests conservatively as p=1.
 for key,x in meta.groupby(['posterior','sensitivity']):
  p=np.r_[x.p.to_numpy(),np.ones(max(0,87-len(x)))];meta.loc[x.index,'FDR_BH_all_29_pairs_x_3_terms']=bh_fdr(p)[:len(x)]
 save(unit,'target_biological_unit_models.tsv');save(meta,'target_meta_descriptive.tsv');save(pd.DataFrame(coverage),'program_gene_coverage.tsv');pd.concat(spot_parts).to_csv(O/'spot_joint_source.tsv.gz',sep='\t',index=False)
 modules=pd.DataFrame(module_rows);munit=modules.groupby(['posterior','GSE','biological_unit_id'],as_index=False).median(numeric_only=True);save(munit,'CM09_biological_unit_descriptive.tsv')
 save(modules.groupby(['posterior','GSE'],as_index=False).median(numeric_only=True),'CM09_dataset_descriptive.tsv')
 save(weights.reset_index(),'CM09_projection_weights.tsv');dump('SOURCE_SHA256.json',hashes)
 summary={'status':'TECHNICAL_COMPLETE_DESCRIPTIVE_NOT_INDEPENDENT_VALIDATION','mapping_root_override':str(args.mapping_root),'datasets':len(ds),'sections':len(audit),'spots':sum(x['n_spots'] for x in audit),'biological_units':pd.DataFrame(audit).biological_unit_id.nunique(),'QC_fallback_sections':sum(x['QC_fallback'] for x in audit),'target_pairs':29,'target_section_posterior_models_ok':sum(x['model_status']=='FIT_OK' for x in detection),'target_model_rows':len(raw),'target_meta_rows':len(meta),'reference_target_excluded':args.mapping_root is not None,'threads':1,'ended':time.time()};dump('COMPLETE.json',summary);dump('STATUS.json',summary)
except Exception as e:
 dump('STATUS.json',{'status':'FAILED','error':repr(e),'traceback':traceback.format_exc()});dump('SOURCE_SHA256.json',hashes);raise
