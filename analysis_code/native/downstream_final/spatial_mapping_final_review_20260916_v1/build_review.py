import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import hashlib,json,base64,html,shutil,zipfile,datetime
import numpy as np,pandas as pd,h5py
from anndata.io import read_elem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
W=Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1')
B=W/'spatial_mapping_final_review_20260916_v1';CUR=W/'spatial_gpu_resume_v2';OLD=W.parent/'spatial_gpu_data_v1/shared_gpu_r3'
for d in ['tables','figures','evidence']:(B/d).mkdir(parents=True,exist_ok=True)
assert not (B/'VALIDATION.json').exists()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(8<<20),b''):h.update(x)
 return h.hexdigest()
c=json.loads((CUR/'configs/config.json').read_text());run=json.loads((CUR/'RUN_STATUS.json').read_text());assert run['status']=='COMPLETE'
manifest=pd.read_csv(c['spatial_manifest'],sep='\t');sig=pd.read_csv(OLD/'outputs/reference_signatures.tsv.gz',sep='\t',index_col=0);assert sig.shape==(20428,128)
refhash=sha(OLD/'outputs/reference_signatures.tsv.gz');assert sha(CUR/'outputs/reference_signatures.tsv.gz')==refhash
summary=[];sections=[];posterior_checks=[];hashes={};elbolast=[];figs=[]
for gse in ['GSE241346','GSE243179','GSE248608','GSE274572','GSE293427','GSE306288','GSE314851']:
 root=OLD if gse in ['GSE241346','GSE243179','GSE248608'] else CUR;out=root/'outputs';prov=json.loads((out/f'{gse}_provenance.json').read_text())
 posterior=out/f'{gse}_posterior.h5ad';ph=sha(posterior);assert ph==prov['posterior_sha'];assert prov['atlas_sha']==c['atlas_sha256'] and prov['reference_sha']==refhash
 assert sha(root/'configs/PRODUCTION_RECEIPT.json')==prov['receipt_sha']
 hashes[str(posterior)]=ph
 with h5py.File(posterior,'r') as f:
  obs=read_elem(f['obs']);assert obs.index.is_unique and len(obs)==prov['n_spots'];assert obs.in_tissue.astype(int).eq(1).all()
  vals={k:read_elem(f['obsm'][k+'_cell_abundance_w_sf']) for k in ['means','stds','q05','q95']}
  for stat,a in vals.items():
   assert a.shape==(len(obs),128) and a.index.equals(obs.index)
   factors=[x.split('cell_abundance_w_sf_',1)[1] for x in a.columns];assert factors==list(sig.columns)
   v=a.to_numpy();assert np.isfinite(v).all() and (v>=0).all()
   posterior_checks.append(dict(GSE=gse,summary=stat,shape=str(a.shape),finite_nonnegative=True,cell_order_match=True,factor_order_match=True))
  means=vals['means'].to_numpy();rowtotal=means.sum(1);assert (rowtotal>0).all() and (means.sum(0)>0).all()
  assert (vals['q05'].to_numpy()<=vals['q95'].to_numpy()).all()
  obs['posterior_mean_total_abundance']=rowtotal
  for sid,v in obs.groupby('section_id',observed=True):
   sections.append(dict(GSE=gse,section_id=sid,sample_id=';'.join(v.sample_id.astype(str).unique()),biological_unit_id=';'.join(v.biological_unit_id.astype(str).unique()),retained_spots=len(v),QC_fallback=bool(v.QC_fallback.any()),median_raw_UMI=float(v.total_counts_original.median()),median_detected_genes=float(v.n_genes_original.median()),median_mean_abundance=float(v.posterior_mean_total_abundance.median()),p05_mean_abundance=float(v.posterior_mean_total_abundance.quantile(.05)),p95_mean_abundance=float(v.posterior_mean_total_abundance.quantile(.95))))
 observed=set(obs.section_id.astype(str));expected=set(manifest.loc[manifest.GSE.eq(gse),'section_id'].astype(str));assert observed==expected,(gse,observed^expected)
 hist=pd.read_csv(out/f'{gse}_history.tsv',sep='\t');h=hist[hist.metric.eq('elbo_train')];assert len(h)==30000 and np.isfinite(h.value).all()
 last=h.tail(10).copy();last['GSE']=gse;elbolast.append(last)
 device=json.loads((out/f'{gse}_device.json').read_text());assert device['device']=='cuda:0'
 fallback=sorted(obs.loc[obs.QC_fallback.astype(bool),'section_id'].astype(str).unique());assert fallback==sorted(prov['QC_fallback_sections'])
 summary.append(dict(GSE=gse,source='reused_previous' if root==OLD else 'new_completed',sections=len(observed),spots=len(obs),genes=prov['n_genes'],factors=128,epochs=len(h),posterior_draws=1000,device=device['device'],median_mean_abundance=float(np.median(rowtotal)),p05_mean_abundance=float(np.quantile(rowtotal,.05)),p95_mean_abundance=float(np.quantile(rowtotal,.95)),prior_N_cells_per_location=12,QC_fallback_sections=len(fallback),status='TECHNICAL_COMPLETE_SCIENTIFIC_REVIEW_PENDING',posterior_path=str(posterior),posterior_sha256=ph))
 fig,axs=plt.subplots(1,2,figsize=(10,3));axs[0].plot(h.epoch,h.value,lw=.7);axs[0].set(xlabel='Epoch',ylabel='Training ELBO loss',title=gse)
 axs[1].hist(rowtotal,bins=60,color='#3a7da6');axs[1].axvline(12,color='red',ls='--',label='Prior expected 12');axs[1].set(xlabel='Posterior mean total abundance',ylabel='Spots');axs[1].legend();fig.tight_layout();p=B/'figures'/f'{gse}_diagnostics.png';fig.savefig(p,dpi=160);plt.close(fig);figs.append(p)
 for suffix in ['provenance.json','device.json','posterior_policy.json','convergence_note.json','degeneracy_audit.json','history.tsv','factor_total_audit.tsv','row_total_audit.tsv.gz']:
  p=out/f'{gse}_{suffix}'
  if p.exists():shutil.copy2(p,B/'evidence'/p.name);hashes[str(p)]=sha(p)
 p=out/f'{gse}_reconstruction_QC.png';shutil.copy2(p,B/'figures'/p.name);figs.append(B/'figures'/p.name)
 print('VALIDATED',gse,len(obs),flush=True)
s=pd.DataFrame(summary);sec=pd.DataFrame(sections);checks=pd.DataFrame(posterior_checks)
s.to_csv(B/'tables/dataset_summary.tsv',sep='\t',index=False);sec.to_csv(B/'tables/section_QC.tsv',sep='\t',index=False);checks.to_csv(B/'tables/posterior_checks.tsv',sep='\t',index=False);pd.concat(elbolast).to_csv(B/'tables/last10_ELBO.tsv',sep='\t',index=False)
fallback=sec[sec.QC_fallback];fallback.to_csv(B/'tables/QC_fallback_sections.tsv',sep='\t',index=False)
assert len(sec)==49
text=f'''# 空间 reference/mapping 完成审核

7/7数据集、{len(sec)}切片、{int(s.spots.sum())}个保留spots已生成描述性mapping输出。3组核验复用，4组新计算；reference不重复训练。所有模型均CUDA，mapping各30000epochs，后验1000抽样，仅存means/stds/q05/q95，不声称存有q50或全部抽样。

本包只读取输出并绘图，没有重新拟合、改变QC、重新训练或开展空间关联/靶基因检验。冻结atlas保持不变。

## 技术检查
七个posterior SHA与原provenance一致；atlas/reference身份一致，receipt哈希一致；各输出spot顺序、128factor顺序、形状、非负/有限数值、q05≤q95检查通过；没有全零spot或factor。49切片与批准manifest精确匹配。ELBO全曲线及末10轮保留，重建QC图来自原生运行。

## 需要科学审核的两类问题
1. {len(fallback)}个切片触发原代码QC fallback：严格QC剩余少于50 spots后，保留组织内非零counts spots。见QC_fallback_sections.tsv。这不是所有spot通过原严格QC，且未自动删除或重新计算。
2. N_cells_per_location=12为先验期望，不是硬上限。后验总丰度及原始UMI/检出基因分布逐切片列出，需结合组织图像及spot大小检查尺度。完整训练轮数、有限ELBO和数值有效不等于科学收敛或准确细胞身份已获证明。

本轮状态：TECHNICAL_COMPLETE_PENDING_SPATIAL_SCIENTIFIC_REVIEW。未进行CKB空间关联、细胞来源验证或独立生物学验证。
'''
(B/'REPORT_CN.md').write_text(text)
body='<html><meta charset="utf-8"><style>body{font:16px sans-serif;max-width:1250px;margin:25px auto;line-height:1.6}table{border-collapse:collapse;font-size:13px}td,th{padding:5px;border:1px solid #ddd}pre{white-space:pre-wrap}img{max-width:100%}</style><h1>七组空间mapping技术审核</h1><pre>'+html.escape(text)+'</pre>'
body+='<h2>数据集汇总</h2>'+s.drop(columns=['posterior_path','posterior_sha256']).to_html(index=False,float_format=lambda v:f'{v:.2f}')
body+='<h2>触发放宽QC的切片</h2>'+fallback.to_html(index=False,float_format=lambda v:f'{v:.2f}')
body+='<h2>全部切片QC</h2>'+sec.to_html(index=False,float_format=lambda v:f'{v:.2f}')
for p in figs:body+='<h3>'+p.stem+'</h3><img src="data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()+'">'
(B/'REPORT_CN.html').write_text(body)
val=dict(status='ALL_SEVEN_TECHNICAL_OUTPUTS_VALIDATED_SCIENTIFIC_REVIEW_PENDING',datasets=7,sections=49,spots=int(s.spots.sum()),posterior_summaries_checked=28,hashes_verified=hashes,QC_fallback_sections=fallback.section_id.tolist(),analysis_rerun=False,association_tests_run=False,checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(B/'VALIDATION.json').write_text(json.dumps(val,ensure_ascii=False,indent=2)+'\n')
(B/'FILES_SHA256.json').write_text(json.dumps({str(p.relative_to(B)):sha(p) for p in B.rglob('*') if p.is_file()},indent=2)+'\n')
with zipfile.ZipFile(B/'spatial_mapping_review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in B.rglob('*'):
  if p.is_file() and p.name!='spatial_mapping_review.zip':z.write(p,p.relative_to(B))
print(s.to_string(index=False));print('FALLBACK',fallback.to_dict('records'));print('ZIP_SHA',sha(B/'spatial_mapping_review.zip'))
