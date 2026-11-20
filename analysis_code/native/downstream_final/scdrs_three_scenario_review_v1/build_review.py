import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,gzip,itertools,html,base64,zipfile
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/scdrs_storage_v1')
OUT=ROOT.parent/'expanded_parallel_20260916_v1/scdrs_three_scenario_review_v1';OUT.mkdir(exist_ok=True)
assert not (OUT/'VALIDATION.json').exists(),'No overwrite'
S={'primary_r0.1':ROOT/'native/primary_r0.1','primary_r0.6':ROOT/'native/primary_r0.6','leave_one_CAD_r0.1':ROOT/'third_score_completion_v3'}
KEYS={'L1':'level1_annotation_scdrs_test','L2':'level2_annotation_scdrs_test'}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8<<20),b''):h.update(x)
 return h.hexdigest()
def ids_sha(p):
 h=hashlib.sha256()
 with gzip.open(p,'rb') as f:
  for x in iter(lambda:f.read(8<<20),b''):h.update(x)
 return h.hexdigest()
def bh(v):
 v=np.asarray(v);order=np.argsort(v);adj=np.minimum.accumulate((v[order]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1].clip(0,1);z=np.empty(len(v));z[order]=adj;return z
frames={};manifest={};summary=[];id_hashes={}
for scenario,p in S.items():
 schema=json.loads((p/'SCHEMA_VALIDATION.json').read_text());assert schema['status']=='PASS' and schema['n_scored_cells']==1175115 and schema['n_control_columns']==1000
 manifest[str(p/'SCHEMA_VALIDATION.json')]=sha(p/'SCHEMA_VALIDATION.json');id_hashes[scenario]=ids_sha(p/'native_scored_cell_ids.tsv.gz')
 for level,key in KEYS.items():
  f=p/key/'eligible_native_group_BH.tsv';d=pd.read_csv(f,sep='\t');assert len(d)==(15 if level=='L1' else 103) and d.group.is_unique
  assert d.n_ctrl.eq(1000).all()
  assert np.isfinite(d[['assoc_mcp','assoc_mcz','hetero_mcp','hetero_mcz']].values).all()
  for stem in ['assoc_mcp','hetero_mcp']:
   assert d[stem].between(0,1).all() and np.allclose(d[stem+'_BH'],bh(d[stem]),atol=1e-7,rtol=1e-6)
  d=d.sort_values('group').reset_index(drop=True);frames[(scenario,level)]=d;manifest[str(f)]=sha(f)
  summary.append({'scenario':scenario,'level':level,'eligible_groups':len(d),'associated_groups_FDR05':int(d.assoc_mcp_BH.lt(.05).sum()),'heterogeneity_groups_FDR05':int(d.hetero_mcp_BH.lt(.05).sum()),'scored_cells':schema['n_scored_cells'],'native_gene_overlap':schema['native_trait_overlap']})
assert len(set(id_hashes.values()))==1,'Different scored cell identity/order'
metrics=[];overlaps=[];wideframes={}
for level in KEYS:
 base=frames[('primary_r0.1',level)]
 for scenario in S:
  d=frames[(scenario,level)];assert d.group.tolist()==base.group.tolist() and np.array_equal(d.n_cell,base.n_cell)
 wide=base[['group','n_cell']].copy()
 for scenario in S:
  d=frames[(scenario,level)]
  for c in ['assoc_mcp','assoc_mcz','assoc_mcp_BH','hetero_mcp','hetero_mcz','hetero_mcp_BH']:wide[scenario+'__'+c]=d[c]
 wide.to_csv(OUT/(level+'_all_scenario_metrics.tsv'),sep='\t',index=False);wideframes[level]=wide
 sets={s:set(frames[(s,level)].loc[frames[(s,level)].assoc_mcp_BH.lt(.05),'group']) for s in S}
 for a,b in itertools.combinations(S,2):
  x,y=frames[(a,level)],frames[(b,level)];metrics.append({'level':level,'scenario_a':a,'scenario_b':b,'assoc_mcz_spearman':float(x.assoc_mcz.rank().corr(y.assoc_mcz.rank())),'FDR05_overlap':len(sets[a]&sets[b]),'FDR05_union':len(sets[a]|sets[b])})
 for g in base.group:overlaps.append({'level':level,'group':g,**{s:g in sets[s] for s in S},'significant_in_all_three':all(g in sets[s] for s in S)})
pd.DataFrame(summary).to_csv(OUT/'scenario_summary.tsv',sep='\t',index=False);pd.DataFrame(metrics).to_csv(OUT/'scenario_pairwise_comparison.tsv',sep='\t',index=False);pd.DataFrame(overlaps).to_csv(OUT/'group_significance_overlap.tsv',sep='\t',index=False)
# Native association statistics only; no changes to tests, covariates or scores.
for level in KEYS:
 wide=wideframes[level];z=np.column_stack([wide[s+'__assoc_mcz'] for s in S]);sig=np.column_stack([wide[s+'__assoc_mcp_BH']<.05 for s in S])
 fig,ax=plt.subplots(figsize=(10,max(5,len(wide)*.22)));lim=max(abs(np.nanmin(z)),abs(np.nanmax(z)))
 im=ax.imshow(z,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto');ax.set_yticks(range(len(wide)));ax.set_yticklabels(wide.group,fontsize=7);ax.set_xticks(range(3));ax.set_xticklabels(list(S),rotation=20,ha='right');ax.set_title(level+' scDRS native association z; * BH FDR < 0.05')
 for i,j in zip(*np.where(sig)):ax.text(j,i,'*',ha='center',va='center',fontsize=8)
 fig.colorbar(im,ax=ax,label='assoc_mcz',shrink=.5);fig.tight_layout();fig.savefig(OUT/(level+'_scenario_heatmap.png'),dpi=160);plt.close(fig)
report='''# scDRS三个遗传输入场景：完成与对照报告

第三场景 leave_one_CAD_r0.1 的原生L1/L2分组检验已完成，原验证器通过。未重新计算任何单细胞score。
三场景均对相同的1,175,115细胞评分；本轮还独立核对了评分cell_id的完整顺序哈希、15个L1/103个L2资格组、各组细胞数、1000个control、P值范围和组内BH结果。

三个场景均分别包含L1和L2，第三场景不是“L2场景”。第三场景改变遗传输入（leave-one-CAD），不是将疾病样本和健康样本分成第三组。

## 结果与解释边界
- primary_r0.1：L1 5组、L2 31组association BH FDR<0.05。
- primary_r0.6：L1 6组、L2 31组association BH FDR<0.05。
- leave_one_CAD_r0.1：L1 0组、L2 6组association BH FDR<0.05。
- 三场景共同显著的L2仅有NRG1 matrix-associated endothelial cell和OMD matrix-associated endothelial cell；没有共同显著L1。
- 说明部分遗传关联结果对场景选择敏感，不能概括为所有场景均稳定。显著组数减少本身不能证明某个疾病或基因具有因果作用。
- 三个场景共享细胞和部分遗传信号，不是三次独立验证。关联检验与异质性检验分别保留，不混为同一个显著性结论。
- scDRS的多基因遗传关联结果与CKB单基因疾病差异表达是不同证据；本表不是CKB特异的遗传检验。
- 临床分组及供者修订不会自动改变这些固定单细胞score；如日后进行供者/疾病层面的关联，需要使用新sidecar另行建模。

## 可追溯性
原生结果和旧验证文件只读。本轮保存全部原生组指标、场景对照表、展示图、输入SHA256和核验信息。全量5GB score文件不复制到审核包。本轮没有重读并重新哈希完整score数值文件，沿用原生验证器已完成的score验证；评分cell_id顺序则在本轮重新全量核对。
'''
# Check the report's overlap statement against actual data rather than hard-code conclusions.
shared=[x['group'] for x in overlaps if x['level']=='L2' and x['significant_in_all_three']]
assert shared==['Endothelial cell | NRG1 matrix-associated endothelial cell','Endothelial cell | OMD matrix-associated endothelial cell'],shared
(OUT/'REPORT_CN.md').write_text(report)
body='<html lang="zh"><meta charset="utf-8"><style>body{font:16px sans-serif;margin:30px}pre{white-space:pre-wrap}td,th{border:1px solid #ccc;padding:5px}table{border-collapse:collapse;font-size:12px}img{max-width:100%}</style><pre>'+html.escape(report)+'</pre>'+pd.DataFrame(summary).to_html(index=False)+pd.DataFrame(metrics).to_html(index=False)
for level in KEYS:
 img=base64.b64encode((OUT/(level+'_scenario_heatmap.png')).read_bytes()).decode();body+='<h2>'+level+'</h2><img src="data:image/png;base64,'+img+'">'+wideframes[level].to_html(index=False,float_format=lambda v:f'{v:.4g}')
(OUT/'REPORT_CN.html').write_text(body+'</html>')
(OUT/'INPUT_SHA256.json').write_text(json.dumps(manifest,indent=2));(OUT/'VALIDATION.json').write_text(json.dumps({'status':'PASS_THREE_SCENARIO_NATIVE_GROUP_COMPARISON','scored_cell_order_sha256':id_hashes,'identical_groups_and_ncells':True,'rechecked_BH':True,'retrained_or_rescored':False,'full_score_numeric_hash_recomputed':False,'same_all_three_L2':shared},indent=2))
(OUT/'FILES_SHA256.json').write_text(json.dumps({p.name:sha(p) for p in OUT.iterdir() if p.is_file()},indent=2))
with zipfile.ZipFile(OUT/'scdrs_three_scenario_review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in OUT.iterdir():
  if p.is_file() and p.suffix!='.zip':z.write(p,p.name)
print(pd.DataFrame(summary).to_string(index=False));print(pd.DataFrame(metrics).to_string(index=False));print(shared)
