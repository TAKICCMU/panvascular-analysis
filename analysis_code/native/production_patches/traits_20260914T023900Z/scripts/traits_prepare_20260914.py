import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['CUDA_VISIBLE_DEVICES']=''
import sys,json,hashlib,datetime,subprocess,time
from pathlib import Path
import numpy as np,pandas as pd,h5py
from anndata.io import read_elem
ROOT=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
CODE=ROOT/'downstream_CPU_L1L2_v2_1_release_compat_20260914'
WORK=ROOT/'production_launches'/'traits_20260914T023900Z'
for d in ['inputs','evidence','smoke','logs','scripts']:(WORK/d).mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
source=ROOT/'panvascular_CPU_theta2_L1L2_final_v1.h5ad'
with h5py.File(source) as f:
 v=read_elem(f['var']);sample_obs={k:read_elem(f['obs'][k]) for k in ['_index','total_counts','n_genes_by_counts','pct_counts_mt','library_size_raw_cpu']}
 schema={k:{'shape':list(f['X'][k].shape),'dtype':str(f['X'][k].dtype),'compression':f['X'][k].compression} for k in f['X']}
 gm=pd.DataFrame({'feature_id':v.index,'gene_symbol':v.gene_symbol.astype('string').to_numpy(),'ensembl_id':v.stable_ensembl_id.astype('string').to_numpy(),'display_symbol':v.display_gene_symbol.astype('string').to_numpy(),'gene_mapping_status':v.gene_mapping_status.astype('string').to_numpy()})
 gm['valid_ensembl']=gm.ensembl_id.str.fullmatch(r'ENSG\d+',na=False)
 assert not gm.loc[gm.valid_ensembl,'ensembl_id'].duplicated().any()
 assert not gm.loc[gm.gene_symbol.notna(),'gene_symbol'].duplicated().any()
 prop=pd.read_csv(CODE/'configs/gene_map_PROPOSED.tsv',sep='\t',keep_default_na=False)
 eq=gm.feature_id.tolist()==prop.feature_id.tolist() and gm.gene_symbol.fillna('').tolist()==prop.gene_symbol.tolist() and gm.ensembl_id.fillna('').tolist()==prop.ensembl_id.tolist()
 assert eq
 gm.to_csv(WORK/'inputs/gene_map_MAGMA_READY.tsv',sep='\t',index=False)
 gm.loc[~gm.valid_ensembl|gm.gene_symbol.isna()].to_csv(WORK/'evidence/unresolved_gene_mapping.tsv',sep='\t',index=False)
 v.loc[v.mt.fillna(False),['gene_symbol','stable_ensembl_id']].to_csv(WORK/'evidence/mitochondrial_gene_panel.tsv',sep='\t')
 rows=json.loads((CODE/'evidence/CPU_TRAIT_INPUT_CHECK.json').read_text())
 for r in rows:
  if r.get('path'):r['observed_sha256']=sha(r['path']);assert r['observed_sha256']==r['sha256']
 (WORK/'evidence/TRAIT_RESOURCES_VERIFIED.json').write_text(json.dumps(rows,indent=2))
 summary={'atlas':str(source),'atlas_stat':{'size':source.stat().st_size,'mtime_ns':source.stat().st_mtime_ns,'inode':source.stat().st_ino},'n_features':len(v),'valid_ensembl':int(gm.valid_ensembl.sum()),'unresolved_symbol':int(gm.gene_symbol.isna().sum()),'PROPOSED_exact_match':eq,'raw_schema':schema,'MAGMA_mapping_status':'READY_missing_ENSG_audited_after_full_library_normalization','scDRS_mapping_status':'BLOCKED_INPUT_735_symbol_fallback_without_canonical_symbol','raw_scan_performed':False,'covariates_prepared':False,'mt_symbols_match_var_mt':bool(np.array_equal(v.mt.fillna(False).to_numpy(),v.gene_symbol.astype('string').str.startswith('MT-',na=False).to_numpy()))}
 for k,arr in sample_obs.items():
  if k!='_index':
   a=np.asarray(arr,dtype=float);summary[k]={'n':len(a),'finite':int(np.isfinite(a).sum()),'min':float(np.nanmin(a)),'max':float(np.nanmax(a))}
 (WORK/'evidence/GENE_MAP_AUDIT.json').write_text(json.dumps(summary,indent=2))
 # Native MAGMA smoke: actual approved gene-results, synthetic covariates only, no atlas-specific findings.
 rng=np.random.default_rng(20260914)
 genes=[]
 with open(rows[0]['path']) as g:
  for line in g:
   if not line.startswith('#') and line.strip():genes.append(line.split()[0])
 pd.DataFrame({'GENE':genes,'Average':rng.uniform(size=len(genes)),'SMOKE_A':rng.uniform(size=len(genes)),'SMOKE_B':rng.uniform(size=len(genes))}).to_csv(WORK/'smoke/magma_covariates.tsv',sep='\t',index=False)
 cmd=['/home/anzhen_vhd/scRNAseq/magma','--gene-results',rows[0]['path'],'--gene-covar',str(WORK/'smoke/magma_covariates.tsv'),'--model','direction=pos','condition=Average','--out',str(WORK/'smoke/magma_native')]
 (WORK/'smoke/magma_command.json').write_text(json.dumps(cmd,indent=2))
 t=time.time()
 with open(WORK/'logs/magma_smoke.log','w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 result={'returncode':r.returncode,'seconds':time.time()-t,'is_scientific_result':False,'synthetic_covariates':True}
 if r.returncode==0:
  out=pd.read_csv(WORK/'smoke/magma_native.gsa.out',sep=r'\s+',comment='#');assert {'VARIABLE','P'}<=set(out);assert out.P.notna().all();result['schema_pass']=True;result['rows']=len(out)
 (WORK/'evidence/MAGMA_NATIVE_SMOKE.json').write_text(json.dumps(result,indent=2));print(str(WORK),json.dumps(result),flush=True)
