from pathlib import Path
import os,time,json,hashlib,shutil,inspect,importlib,importlib.metadata as im
os.environ['CUDA_VISIBLE_DEVICES']=''
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='16'
import numpy as np,pandas as pd,anndata as ad
R=Path(__file__).resolve().parents[1]; O=R.parent/'reliable_L1_v2_20260909'; BASE=R.parents[2]; W=BASE/'02_scib_benchmark/studywise_sampleonly_v1'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def nh(names):return hashlib.sha256(b''.join(str(s).encode()+b'\0' for s in names)).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
def status(s):dump(R/'STATUS.json',dict(stage=s,pid=os.getpid(),time=time.time()))
for d in ['inputs','embeddings','metrics','logs','audit','figures','tables','cache']: (R/d).mkdir(exist_ok=True)
status('VALIDATING_FROZEN_INPUT')
checks={}
try:
 prep=json.loads((O/'inputs/PREPARE_COMPLETE.json').read_text()); ev=json.loads((O/'metrics/evaluation_manifest.json').read_text()); reuse=json.loads((O/'audit/REUSED_EMBEDDINGS.json').read_text()); hm=next(x for x in reuse if x['method']=='harmony'); pm=next(x for x in reuse if x['method']=='pca'); pmeta=json.loads(Path(pm['metadata']['parameters']['pca_metadata']).read_text())
 source=O/'inputs/primary.h5ad'; original=Path(pmeta['fit_input'])
 for label,p,expected in [('frozen_input',source,prep['input_sha256']),('original_PCA_input',original,pmeta['fit_input_sha256']),('PCA',O/'embeddings/pca.npy',ev['pca_sha256']),('GPU_Harmony',O/'embeddings/harmony.npy',hm['sha256'])]:
  actual=sha(p);checks[label]={'path':str(p),'sha256':actual,'expected':expected,'passed':actual==expected};assert actual==expected,label
 a=ad.read_h5ad(source,backed='r');old=ad.read_h5ad(original,backed='r');obs=a.obs.copy()
 assert a.n_obs==315949==prep['n_cells'];assert a.n_vars==2000;assert obs.index.is_unique;assert nh(obs.index)==prep['cell_ids_sha256']==hm['metadata']['obs_names_sha256'];assert obs.index.equals(old.obs_names);assert a.var_names.equals(old.var_names);assert nh(a.var_names)==pmeta['var_names_sha256']
 for col in ['sample_uid','studywise_marker_level1']:assert np.array_equal(obs[col].astype(str),old.obs[col].astype(str)),col
 labs=pd.read_csv(O/'labels/primary_supervision_labels.tsv.gz',sep='\t',index_col='cell_id');assert obs.index.equals(labs.index)
 for col in ['reliable_L1','supervision_label','trusted']:assert np.array_equal(obs[col].astype(str),labs[col].astype(str)),col
 mask=obs.trusted.to_numpy(bool);selected=obs.loc[mask,['reliable_L1','trusted','sample_uid']].copy();selected['reliable_L1']=selected.reliable_L1.astype(str);selected['sample_uid']=selected.sample_uid.astype(str)
 membership=hashlib.sha256(json.dumps({'all_cells':list(obs.index.astype(str)),'trusted':mask.tolist(),'labels':selected.reliable_L1.tolist(),'batches':selected.sample_uid.tolist()},sort_keys=True).encode()).hexdigest()
 assert membership==ev['membership_sha256'];saved=pd.read_csv(O/'metrics/evaluated_cells.tsv.gz',sep='\t',index_col=0);assert saved.index.equals(selected.index)
 for c in selected:assert np.array_equal(selected[c].astype(str),saved[c].astype(str)),c
 assert int(mask.sum())==182968;assert obs.sample_uid.nunique()==247
 checks['identity_and_metadata']={'passed':True,'n_primary':len(obs),'n_evaluated':int(mask.sum()),'n_samples':247,'n_HVG':a.n_vars,'n_PCs':30,'membership_sha256':membership,'cell_ids_sha256':nh(obs.index),'HVG_order_sha256':nh(a.var_names),'evaluation_labels':'reliable_L1','integration_label':'not used','batch':'sample_uid'}
 a.file.close();old.file.close()
 # Freeze audit of all historical benchmark files touched/read. Check again at completion.
 hist={str(p):sha(p) for d in ['metrics','code','labels','embeddings'] for p in (O/d).iterdir() if p.is_file() and p.suffix not in ['.lock']}
 hist[str(source)]=sha(source);dump(R/'audit/historical_hashes.json',hist)
 for name,target in [('inputs/primary.h5ad',source),('embeddings/pca.npy',O/'embeddings/pca.npy'),('embeddings/harmony_gpu.npy',O/'embeddings/harmony.npy')]:
  dst=R/name
  if not dst.exists():dst.symlink_to(target)
 for f in ['evaluation_manifest.json','raw_scib_metrics.tsv','official_scores.tsv','final_ranking.tsv']:shutil.copy2(O/'metrics'/f,R/'audit'/('original_'+f))
 obs.to_csv(R/'tables/frozen_primary_metadata.tsv.gz',sep='\t',index_label='cell_id');a=None
 dump(R/'audit/INPUT_VALIDATION.json',{'passed':True,'checks':checks})
 (R/'benchmark_input_validation_report.md').write_text('# 冻结 benchmark 输入核验\n\nPASS：315,949 个 Primary 细胞、247 个 sample_uid、2000 个 HVG 及其顺序、30 维原 PCA、182,968 个可信细胞评价集合均与历史冻结记录一致。\n\n逐细胞比对 cell_id、sample_uid、原标签、reliable_L1、supervision_label、trusted 与独立保存的标签和评价清单；文件 SHA256 与历史 manifest 一致。表达源为冻结 primary.h5ad；本轮仅引用，不从全量 atlas 重建。\n\n整合 seed=0，评分 seed=0；官方指标专用邻居数 15/50/90。主分数不做跨方法 min-max，biology/batch 权重 0.6/0.4。\n\n旧 HVG/PCA 原样保留，包括此前报告的 38 个基因差异；本轮不实施 clip 更新。\n\n详见 audit/INPUT_VALIDATION.json 和 historical_hashes.json。\n')
except Exception as e:
 dump(R/'audit/INPUT_VALIDATION.json',{'passed':False,'checks':checks,'error':repr(e)});(R/'benchmark_input_validation_report.md').write_text('核验失败，已停止，未运行 Harmony。\n'+repr(e));status('BLOCKED_INPUT_MISMATCH');raise
status('CPU_HARMONY')
import omicverse as ov,torch
ov.settings.cpu_init();torch.set_num_threads(16)
m=importlib.import_module('omicverse.external.harmony.harmony');cpu=importlib.import_module('omicverse.external.harmony._harmony_cpu')
assert im.version('omicverse')==hm['metadata']['environment']['packages']['omicverse']
for mod in [m,cpu]:shutil.copy2(inspect.getfile(mod),R/'audit'/Path(inspect.getfile(mod)).name)
shutil.copy2(BASE/'code/07_run_pca_harmony_scanorama.py',R/'audit/original_integration_runner.py')
params=dict(theta=1.,lamb=None,sigma=.1,nclust=None,tau=0,block_size=.05,max_iter_harmony=20,max_iter_kmeans=20,epsilon_cluster=1e-5,epsilon_harmony=1e-4,alpha=.2,verbose=True,random_state=0,device='cpu')
# Historical runner specified theta, max_iter, random_state, device, verbose;
# assert every formerly implicit setting still equals installed default.
for k in ['lamb','sigma','nclust','tau','block_size','max_iter_kmeans','epsilon_cluster','epsilon_harmony','alpha']:assert params[k]==inspect.signature(m.run_harmony).parameters[k].default,k
config={'parameters':params,'implicit_defaults_provenance':'historical runner + same installed OmicVerse version; historic source hash not archived, cannot prove bitwise historical source identity','versions':{k:im.version(k) for k in ['omicverse','numpy','scipy','torch','anndata']},'threads':16,'cuda_visible_devices':os.environ['CUDA_VISIBLE_DEVICES'],'GPU_Harmony_used':False,'source':inspect.getfile(m),'implementation':'OmicVerse HarmonyCPU NumPy','historical_GPU_metadata':hm}
dump(R/'audit/harmony_config.json',config)
if not (R/'embeddings/harmony_cpu.npy').exists():
 z=np.load(R/'embeddings/pca.npy');assert z.shape==(315949,30);obs.sample_uid=obs.sample_uid.astype(str).astype('category');start=time.time();res=m.run_harmony(z,obs,'sample_uid',**params);assert isinstance(res,cpu.HarmonyCPU)
 out=np.asarray(res.result(),dtype=np.float32);assert out.shape==z.shape and np.isfinite(out).all();np.save(R/'embeddings/harmony_cpu.npy',out)
 obj={k:np.asarray(getattr(res,k)).tolist() for k in ['objective_harmony','objective_kmeans'] if hasattr(res,k)}
 dump(R/'audit/harmony_runtime.json',dict(seconds=time.time()-start,backend=type(res).__name__,embedding_sha256=sha(R/'embeddings/harmony_cpu.npy'),shape=list(out.shape),objectives=obj))
status('CPU_HARMONY_COMPLETE')
