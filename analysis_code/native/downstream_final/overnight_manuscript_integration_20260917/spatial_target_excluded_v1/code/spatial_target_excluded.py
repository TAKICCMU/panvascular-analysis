"""Separate authorized descriptive cell2location runner. Never imports CPU bootstrap."""
from pathlib import Path
import os,json,hashlib,sys,traceback,time,gc,ast,inspect
JOB=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(name,value):
 p=JOB/name;p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(value,indent=2,default=str))
def guard():
 r=json.loads((JOB/'configs/PRODUCTION_RECEIPT.json').read_text());c=json.loads((JOB/'configs/config.json').read_text())
 assert r['scope']==['target_excluded_reference','target_excluded_mapping']
 assert c['reference_mode']=='target_excluded_both_stages' and len(c['target_ids'])==29
 for p,h in r['bound_files'].items():
  assert Path(p).is_file() and sha(p)==h, 'Bound input changed '+p
 assert os.environ.get('CUDA_VISIBLE_DEVICES')==c['gpu_uuid']
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
  assert os.environ[k]=='4'
 return c

def execute(stage):
 c=guard()
 import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp,scanpy as sc,torch,scvi
 from scipy.io import mmread
 try:
  from scipy.io._fast_matrix_market import _fmm_core
 except ModuleNotFoundError:
  pass
 from threadpoolctl import threadpool_limits
 threadpool_limits(limits=4)
 from cell2location.models import RegressionModel,Cell2location
 from cell2location.utils.filtering import filter_genes

 torch.set_num_threads(4);torch.set_num_interop_threads(1);scvi.settings.dl_num_workers=0
 torch.set_float32_matmul_precision('high');scvi.settings.seed=c['seed'];np.random.seed(c['seed']);torch.manual_seed(c['seed'])
 assert torch.cuda.is_available(), 'CUDA required; no CPU fallback'
 torch.cuda.set_device(0);test=torch.arange(128,device='cuda').float();assert float((test+1).sum().cpu())==8256
 api={n:str(inspect.signature(m)) for n,m in [('RegressionModel.train',RegressionModel.train),('RegressionModel.sample_posterior',RegressionModel.sample_posterior),('Cell2location.train',Cell2location.train),('Cell2location.sample_posterior',Cell2location.sample_posterior)]}
 assert 'return_samples' in inspect.signature(Cell2location.sample_posterior).parameters
 dump('evidence/CUDA_API_ACTUAL.json',{'api':api,'device':torch.cuda.get_device_name(0),'tensor_device':str(test.device),'torch':torch.__version__,'pid':os.getpid(),'visible':os.environ['CUDA_VISIBLE_DEVICES'],'threads':torch.get_num_threads()})
 out=JOB/('smoke' if stage=='smoke' else 'outputs');out.mkdir(exist_ok=True)
 def history(model,stem):
  rows=[]
  for metric,values in model.history.items():
   vals=np.asarray(values).ravel()
   rows += [{'metric':metric,'epoch':i,'value':float(x)} for i,x in enumerate(vals)]
  assert rows,'Missing training history'
  pd.DataFrame(rows).to_csv(out/f'{stem}_history.tsv',sep='\t',index=False)
  assert np.isfinite([r['value'] for r in rows]).all(), 'Nonfinite training history; full values saved'
  dump(str(out.relative_to(JOB)/f'{stem}_convergence_note.json'),{'classification':'technical diagnostics only; scientific convergence not certified','epochs':len(rows),'no_arbitrary_ELBO_threshold':True,'requires_review':['ELBO trajectory','reconstruction QC','prior/posterior scale']})
  return rows
 def train(model,stem,epochs,batch):
  model.train(max_epochs=epochs,batch_size=batch,train_size=1,accelerator='gpu',device=[0])
  tensors=list(model.module.parameters())+list(model.module.buffers());devices=sorted({str(t.device) for t in tensors});assert devices and all(d.startswith('cuda') for d in devices),'Model tensors not on CUDA: '+str(devices);device=devices[0]
  hh=history(model,stem);model.save(str(out/f'{stem}_model'),overwrite=False)
  dump(str(out.relative_to(JOB)/f'{stem}_device.json'),{'device':device,'pid':os.getpid(),'epochs':epochs,'batch_size':batch,'gpu_peak_bytes':torch.cuda.max_memory_allocated(),'history_rows':len(hh)})
 def posterior(model,a,kind,stem):
  # Native sampling produces actual means/std/quantiles and no persisted draws.
  kw={'num_samples':8 if stage=='smoke' else c['posterior_samples'],'batch_size':min(2500,a.n_obs) if kind=='reference' else a.n_obs,'accelerator':'gpu','device':[0],'return_samples':False}
  a=model.export_posterior(a,sample_kwargs=kw,use_quantiles=False)
  assert 'posterior_samples' not in model.samples,'Unexpected raw posterior draws'
  factor_names=list(map(str,a.uns['mod']['factor_names']))
  vals=a.varm['means_per_cluster_mu_fg'] if kind=='reference' else a.obsm['means_cell_abundance_w_sf']
  vv=np.asarray(vals);expected=(a.n_vars if kind=='reference' else a.n_obs,len(factor_names))
  assert vv.shape==expected,'Posterior axis mismatch'
  assert np.isfinite(vv).all() and (vv>=0).all(),'Invalid negative/nonfinite posterior mean'
  rowsum=vv.sum(1);colsum=vv.sum(0);rowlabels=a.var_names if kind=='reference' else a.obs_names
  pd.DataFrame({'axis_id':rowlabels,'total_mean':rowsum,'zero_total':rowsum==0}).to_csv(out/f'{stem}_row_total_audit.tsv.gz',sep='\t',index=False)
  pd.DataFrame({'factor':factor_names,'total_mean':colsum,'zero_factor':colsum==0}).to_csv(out/f'{stem}_factor_total_audit.tsv',sep='\t',index=False)
  assert vv.sum()>0,'Completely zero posterior means'
  dump(str(out.relative_to(JOB)/f'{stem}_degeneracy_audit.json'),{'expected_shape':expected,'actual_shape':vv.shape,'zero_rows':int((rowsum==0).sum()),'zero_factors':int((colsum==0).sum()),'overall_positive':True,'total_mean':float(vv.sum()),'interpretation':'Technical checks only; zero factors/spots explicitly flagged for scientific review'})
  for name,d in model.samples.items():
   if isinstance(d,dict):
    for k,v in d.items():assert np.isfinite(v).all(), 'Nonfinite posterior '+name+':'+k
  # Native defaults export q05/q95; q50 is not claimed as available or as mean.
  dump(str(out.relative_to(JOB)/f'{stem}_posterior_policy.json'),{'num_samples':kw['num_samples'],'return_samples':False,'summaries':['means','stds','q05','q95'],'q50':'not requested; native default summaries retained','use_quantiles':False,'saved_raw_draws':False})
  import matplotlib;matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  model.plot_QC(use_n_obs=min(1000,a.n_obs));plt.savefig(out/f'{stem}_reconstruction_QC.png',dpi=160,bbox_inches='tight');plt.close('all')
  for k in list(a.layers):del a.layers[k]
  a.write_h5ad(out/f'{stem}_posterior.h5ad',compression='gzip')
  return a
 if stage=='smoke':
  rng=np.random.default_rng(c['seed']);x=rng.poisson(np.vstack([np.r_[np.repeat(6.,20),np.repeat(1.,20)],np.r_[np.repeat(1.,20),np.repeat(6.,20)]])[np.arange(80)%2])
  a=ad.AnnData(sp.csr_matrix(x.astype('float32')),obs=pd.DataFrame({'state':['A','B']*40,'batch':['s1']*40+['s2']*40},index=[f'c{i}' for i in range(80)]),var=pd.DataFrame(index=[f'g{i}' for i in range(40)]))
  RegressionModel.setup_anndata(a,batch_key='batch',labels_key='state');m=RegressionModel(a);train(m,'tiny_reference',2,40);a=posterior(m,a,'reference','tiny_reference');sig=a.varm['means_per_cluster_mu_fg'].copy();sig.columns=a.uns['mod']['factor_names'];del m;gc.collect();torch.cuda.empty_cache()
  b=ad.AnnData(sp.csr_matrix(rng.poisson(np.tile(sig.to_numpy().sum(1)*2,(24,1))).astype('float32')),obs=pd.DataFrame({'sample_id':['p1']*12+['p2']*12},index=[f'spot{i}' for i in range(24)]),var=pd.DataFrame(index=sig.index))
  Cell2location.setup_anndata(b,batch_key='sample_id');m=Cell2location(b,cell_state_df=sig,N_cells_per_location=12,detection_alpha=20);train(m,'tiny_mapping',2,None);b=posterior(m,b,'mapping','tiny_mapping');assert b.obsm['means_cell_abundance_w_sf'].shape==(24,2)
  dump('smoke/PASS.json',{'status':'PASS','native_reference':True,'native_mapping':True,'CUDA_tensor_and_model':True,'full_analysis':False,'code_sha':sha(__file__)})
  return
 assert json.loads((JOB/'smoke/PASS.json').read_text())['code_sha']==sha(__file__),'Smoke code identity mismatch'
 from exclusion import excluded_mask, assert_excluded
 excluded=set(c['excluded_feature_tokens'])
 inp=pd.read_csv(JOB/'evidence/SPATIAL_CONTENT_VALIDATION.tsv',sep='\t');man=pd.read_csv(c['spatial_manifest'],sep='\t')
 old=Path(c['reuse_completed_root'])
 reference_done=out/'reference_COMPLETE.json'
 if reference_done.exists():
  z=json.loads(reference_done.read_text());assert z['code_sha']==sha(__file__) and z['receipt_sha']==sha(JOB/'configs/PRODUCTION_RECEIPT.json')
  assert z['signature_sha']==sha(out/'reference_signatures.tsv.gz')
  sig=pd.read_csv(out/'reference_signatures.tsv.gz',sep='\t',index_col=0)
 else:
  # New model from raw counts, with exactly the original cells, fixed labels and retained non-target genes.
  assert not (out/'reference_model').exists(), 'Partial reference model: preserve and review before restart'
  ref=ad.read_h5ad(c['reference_counts'])
  mask=excluded_mask(ref.var,excluded)
  pd.DataFrame({'gene':ref.var_names,'excluded':mask}).to_csv(out/'reference_exclusion_audit.tsv',sep='\t',index=False)
  ref=ref[:,~mask].copy();assert_excluded(ref.var,excluded)
  for k in list(ref.uns):
   if k.startswith('_scvi'):del ref.uns[k]
  RegressionModel.setup_anndata(ref,batch_key='reference_study',labels_key='reference_state')
  ref.write_h5ad(out/'reference_counts.h5ad',compression='gzip')
  m=RegressionModel(ref);train(m,'reference',250,2500);ref=posterior(m,ref,'reference','reference')
  sig=ref.varm['means_per_cluster_mu_fg'].copy();sig.columns=ref.uns['mod']['factor_names'];sig.to_csv(out/'reference_signatures.tsv.gz',sep='\t')
  assert len(sig.columns)==128 and not set(map(str.upper,sig.index)) & excluded
  del m,ref;gc.collect();torch.cuda.empty_cache()
  dump('outputs/reference_COMPLETE.json',{'signature_sha':sha(out/'reference_signatures.tsv.gz'),'receipt_sha':sha(JOB/'configs/PRODUCTION_RECEIPT.json'),'code_sha':sha(__file__),'reference_retrained':True,'target_excluded_before_model_setup':True})
 factors=pd.read_csv(old/'outputs/reference_factors.tsv',sep='\t');factors.to_csv(out/'reference_factors.tsv',sep='\t',index=False)
 nuisance=set(factors.loc[factors.reference_group.str.startswith('Other | '),'reference_state'])
 assert not set(map(str.upper,sig.index)) & excluded
 for gse,rows in man[man.include_mapping.astype(str).str.lower().isin(['true','1','yes'])].groupby('GSE',sort=True):
  completed=out/f'{gse}_COMPLETE.json'
  if completed.exists():
   z=json.loads(completed.read_text());assert z['code_sha']==sha(__file__) and z['receipt_sha']==sha(JOB/'configs/PRODUCTION_RECEIPT.json') and z['reference_sha']==sha(out/'reference_signatures.tsv.gz') and z['posterior_sha']==sha(out/f'{gse}_posterior.h5ad');continue
  assert not (out/f'{gse}_model').exists(), 'Partial mapping model: preserve and review before restart '+gse
  if not inp[inp.GSE==gse].status.eq('PASS').all():
   dump('outputs/'+gse+'_BLOCKED.json',{'reason':'spatial exact content validation failed'});continue
  try:
   sections=[]
   for _,row in rows.iterrows():
    gm=pd.read_csv(JOB/'evidence'/f'{gse}__{str(row.section_id).replace(":", "_")}__gene_map.tsv.gz',sep='\t',keep_default_na=False)
    mx=mmread(row.matrix).tocsr().T;bar=pd.read_csv(row.barcodes,sep='\t',header=None,dtype=str).iloc[:,0];first=pd.read_csv(row.positions,header=None,nrows=1)
    pos=pd.read_csv(row.positions) if str(first.iloc[0,0]).lower() in ['barcode','barcodes'] else pd.read_csv(row.positions,header=None)
    pos=pos.iloc[:,:6];pos.columns=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres'];pos=pos.set_index('barcode').loc[bar]
    tissue=pos.in_tissue.astype(int).to_numpy()==1;origcount=np.asarray(mx.sum(1)).ravel();ngene=np.asarray((mx>0).sum(1)).ravel();qc=(origcount>=(500 if mx.shape[1]>10000 else 100))&(ngene>=(200 if mx.shape[1]>10000 else 50))&tissue
    fallback=qc.sum()<50
    if fallback:qc=(origcount>0)&tissue
    # Full source-feature UMI includes targets and unmapped features for independent outcomes.
    rawobs=pos.iloc[np.flatnonzero(qc)].copy();rawobs.index=pd.Index([str(row.section_id)+'::'+str(b) for b in bar[qc]])
    assert gm.source_row.tolist()==list(range(mx.shape[1])), 'Gene map does not cover source matrix in order'
    rawvar=gm.copy().set_index('source_row');rawvar.index=rawvar.index.astype(str)
    raw=ad.AnnData(mx[qc].astype('float32'),obs=rawobs,var=rawvar)
    raw.uns['scope']='outcome expression only; never input to target-excluded model'
    rawdir=out/'full_umi';rawdir.mkdir(exist_ok=True)
    rawpath=rawdir/(str(row.section_id).replace(':','_')+'.h5ad')
    raw.write_h5ad(rawpath,compression='gzip');del raw
    hit=gm[gm.atlas_gene!=''];genes=pd.Index(hit.atlas_gene.unique());cols=genes.get_indexer(hit.atlas_gene);agg=sp.csr_matrix((np.ones(len(hit)),(hit.source_row.to_numpy(),cols)),shape=(mx.shape[1],len(genes)))
    sx=(mx[qc]@agg).astype('float32');so=pos.iloc[np.flatnonzero(qc)].copy();so.index=pd.Index([str(row.section_id)+'::'+str(b) for b in bar[qc]])
    for k in ['GSE','sample_id','section_id','biological_unit_id','metadata_status','condition','tissue','dataset_role']:so[k]=str(row[k])
    so['QC_fallback']=fallback;so['total_counts_original']=origcount[qc];so['n_genes_original']=ngene[qc];section=ad.AnnData(sx,obs=so,var=pd.DataFrame(index=genes));section.obsm['spatial']=so[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();sections.append(section)
   spatial=ad.concat(sections,join='outer',fill_value=0,merge='same');spatial.write_h5ad(out/f'{gse}_full_atlas_gene_UMI.h5ad',compression='gzip');common=spatial.var_names.intersection(sig.index);common=common[[str(x).upper() not in excluded for x in common]];assert not set(map(str.upper,common)) & excluded;assert len(common)>=1000, 'Insufficient common genes'
   spatial=spatial[:,common].copy();spatial.write_h5ad(out/f'{gse}_counts.h5ad',compression='gzip');Cell2location.setup_anndata(spatial,batch_key='sample_id');m=Cell2location(spatial,cell_state_df=sig.loc[common],N_cells_per_location=12,detection_alpha=20);train(m,gse,30000,None);spatial=posterior(m,spatial,'mapping',gse)
   for stat in ['means','stds','q05','q95']:
    val=spatial.obsm[stat+'_cell_abundance_w_sf'];named=[x for x in val.columns if not any(n in str(x) for n in nuisance)];other=[x for x in val.columns if x not in named]
    val[named].to_csv(out/f'{gse}_{stat}_named.tsv.gz',sep='\t');val[other].to_csv(out/f'{gse}_{stat}_nuisance.tsv.gz',sep='\t')
   dump('outputs/'+gse+'_provenance.json',{'atlas_sha':c['atlas_sha256'],'receipt_sha':sha(JOB/'configs/PRODUCTION_RECEIPT.json'),'reference_sha':sha(out/'reference_signatures.tsv.gz'),'posterior_sha':sha(out/f'{gse}_posterior.h5ad'),'n_spots':spatial.n_obs,'n_genes':spatial.n_vars,'QC_fallback_sections':spatial.obs.loc[spatial.obs.QC_fallback,'section_id'].unique().tolist(),'scope':'target-excluded reference and mapping conditional on frozen atlas labels; full UMI reserved for outcomes','status':'TECHNICAL_COMPLETE_REQUIRES_SCIENTIFIC_REVIEW'})
   dump('outputs/'+gse+'_COMPLETE.json',{'code_sha':sha(__file__),'receipt_sha':sha(JOB/'configs/PRODUCTION_RECEIPT.json'),'reference_sha':sha(out/'reference_signatures.tsv.gz'),'posterior_sha':sha(out/f'{gse}_posterior.h5ad'),'n_genes':spatial.n_vars,'excluded_targets_present':[]})
   del m,spatial,sections;gc.collect();torch.cuda.empty_cache()
  except Exception as e:
   dump('outputs/'+gse+'_FAILED.json',{'error':repr(e),'traceback':traceback.format_exc()});print(traceback.format_exc(),flush=True)
   if isinstance(e,torch.cuda.OutOfMemoryError):torch.cuda.empty_cache()
 dump('RUN_STATUS.json',{'status':'COMPLETE' if not list(out.glob('*_FAILED.json')) and not list(out.glob('*_BLOCKED.json')) and len(list(out.glob('GSE*_COMPLETE.json')))==7 else 'PARTIAL_FAILED','scientific_review_required':True,'completed_at':time.time()})
if __name__=='__main__':
 try:execute(sys.argv[1])
 except Exception as e:
  dump('RUN_STATUS.json',{'status':'FAILED','error':repr(e),'traceback':traceback.format_exc(),'stage':sys.argv[1]});raise
