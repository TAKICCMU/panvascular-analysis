import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,base64,html,zipfile,shutil
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path(__file__).resolve().parents[1];W=B.parent

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for a in iter(lambda:f.read(8<<20),b''):h.update(a)
 return h.hexdigest()
def save(name,x): (B/name).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n')
r=pd.read_csv(B/'tables/model_roster.tsv',sep='\t',keep_default_na=False);cks=[];status=[];hashes={};los=[];scopechecks=[]
for row in r[r.status.eq('READY')].itertuples():
 out=B/'results'/row.model_id;st=json.loads((out/'MODEL.json').read_text());status.append(dict(model_id=row.model_id,level=row.level,context=row.context,state=row.cell_state,**{k:v for k,v in st.items() if k!='LOS'}))
 if st['status']!='COMPLETE':continue
 de=pd.read_csv(out/'DE.tsv.gz',sep='\t');assert de.gene_id.is_unique and np.isfinite(de.PValue).all() and de.PValue.between(0,1).all();assert np.allclose(de.FDR,false_discovery_control(de.PValue.to_numpy()),rtol=1e-10,atol=1e-13)
 hashes[str(out/'DE.tsv.gz')]=sha(out/'DE.tsv.gz');cc=de[de.gene_symbol.eq('CKB')].copy()
 if len(cc):cks.append(cc)
 inc=pd.read_csv(out/'included_units.tsv',sep='\t');assert inc.n_cells.ge(20).all() and not inc.study_canonical.str.contains('HCA').any()
 for study,v in inc.groupby('study_canonical'):
  scopechecks.append(dict(model_id=row.model_id,study=study,n_case=int(v.condition.eq('disease').sum()),n_reference=int(v.condition.eq('healthy').sum())))
 for p in out.glob('LOS_*.tsv.gz'):
  z=pd.read_csv(p,sep='\t');los.append(z[z.gene_symbol.eq('CKB')]);hashes[str(p)]=sha(p)
new=pd.concat(cks,ignore_index=True);new.to_csv(B/'tables/CKB_new_original_strategy_context.tsv',sep='\t',index=False)
meta=pd.concat([pd.read_csv(p,sep='\t') for p in (B/'meta').glob('*.tsv')],ignore_index=True)
for _,z in meta.groupby(['annotation_level','cell_state']):assert z.gene_id.is_unique and np.allclose(z.target_set_meta_FDR,false_discovery_control(z.meta_P.to_numpy()),rtol=1e-10,atol=1e-13)
meta.to_csv(B/'tables/all_priority_meta.tsv.gz',sep='\t',index=False)
mk=meta[meta.gene_symbol.eq('CKB')].copy();mk.to_csv(B/'tables/CKB_new_original_strategy_meta.tsv',sep='\t',index=False)
st=pd.DataFrame(status);st.to_csv(B/'tables/model_status.tsv',sep='\t',index=False);pd.DataFrame(scopechecks).to_csv(B/'tables/post_design_study_counts.tsv',sep='\t',index=False)
if los:pd.concat(los,ignore_index=True).to_csv(B/'tables/CKB_leave_one_study.tsv',sep='\t',index=False)
old=pd.read_csv(B/'evidence/old_CKB_all_disease_context_edgeR_results.tsv',sep='\t');oldmeta=pd.read_csv(B/'evidence/old_CKB_GWAS_prioritised_cross_disease_meta_analysis.tsv',sep='\t')
old['analysis']='original_old_atlas';oldmeta['analysis']='original_old_atlas';new['analysis']='original_strategy_corrected_atlas';mk['analysis']='original_strategy_corrected_atlas'
old.annotation_level=old.annotation_level.replace({'Level1':'L1','Level2':'L2'});oldmeta.annotation_level=oldmeta.annotation_level.replace({'Level1':'L1','Level2':'L2'})
comp=pd.concat([oldmeta,mk],ignore_index=True);comp.to_csv(B/'tables/CKB_old_vs_new_meta.tsv',sep='\t',index=False)
current=pd.read_csv(W/'pseudobulk_3backgrounds_20260916_v1/tables/CKB_review_results.tsv',sep='\t');current['analysis']='expanded_3background_current';current=current.rename(columns={'level':'annotation_level','context':'disease_context','logFC':'logFC'})
ctx=pd.concat([old,new,current],ignore_index=True);ctx.to_csv(B/'tables/CKB_all_versions_context.tsv',sep='\t',index=False)
# Display only L1 estimates; all L2 and leave-study rows remain downloadable.
plot=comp[comp.annotation_level.eq('L1') & comp.cell_state.isin(['Fibroblast','Macrophage'])].copy();plot['name']=plot.cell_state+' / '+plot.analysis
fig,ax=plt.subplots(figsize=(9,3.6));y=np.arange(len(plot));colors=['#888888' if x=='original_old_atlas' else '#2166ac' for x in plot.analysis]
for i,(_,a) in enumerate(plot.iterrows()):ax.errorbar(a.meta_logFC,i,xerr=[[a.meta_logFC-a.meta_CI_low],[a.meta_CI_high-a.meta_logFC]],fmt='o',color=colors[i],capsize=3)
ax.set_yticks(y,plot.name);ax.axvline(0,color='black',lw=.8);ax.set_xlabel('CKB pooled log2 fold change (disease / reference), legacy REML 95% CI');fig.tight_layout();fig.savefig(B/'tables/CKB_meta_comparison.png',dpi=180);fig.savefig(B/'tables/CKB_meta_comparison.pdf');plt.close(fig)
summary=pd.read_csv(B/'tables/included_study_summary.tsv',sep='\t');scope=pd.read_csv(B/'tables/sample_scope_audit.tsv',sep='\t');hca=pd.read_csv(B/'tables/HCA_original_scope.tsv',sep='\t');block=st.get('block_omitted',pd.Series(dtype=bool)).fillna(False).astype(bool).sum()
main=comp[comp.annotation_level.eq('L1') & comp.cell_state.isin(['Fibroblast','Macrophage'])][['analysis','cell_state','n_disease_contexts','meta_logFC','meta_CI_low','meta_CI_high','meta_P','target_set_meta_FDR','I2']]
oldc=old[old.annotation_level.eq('L1') & old.cell_state.isin(['Fibroblast','Macrophage'])][['cell_state','disease_context','logFC','PValue','FDR']]
newc=new[new.annotation_level.eq('L1') & new.cell_state.isin(['Fibroblast','Macrophage'])][['cell_state','disease_context','n_case','n_control','logFC','PValue','FDR','design_formula']]
text=f'''# 原文 pseudobulk 策略复现与 CKB 核查

本版使用当前冻结CPU L1/L2与正确counts，复现最后投稿包的主分析范围和原始R函数，不恢复旧counts/旧注释/已知错误供者信息。既有三背景分析、HCA信息和冻结图谱均未修改。

## 原文为什么显著
原文表4中的CKB是跨背景REML合并效应，并在MAGMA预选靶基因集合内进行BH校正；不是单疾病背景全基因FDR。旧成纤维细胞两个背景及旧巨噬细胞三个背景的单背景CKB均未达到全基因FDR<0.05，但效应同向，合并后分别报告约0.00345和0.03946。二者检验问题与校正集合不同，不能直接拿此前单背景FDR判断“显著性丢失”。数据去重、细胞注释、供者配对和病例对照范围也改变，不能未经单因素实验把全部差异归因于HCA。

## 原文HCA分组
服务器原始sample_disease_map.tsv中HCA为healthy，disease_context空白、include_primary=FALSE；原始建图代码和实际模型审核表一致。HCA参加图谱描述，但未作为本次六研究主pseudobulk的病例或对照。此次在当前数据中排除全部{len(hca)}个HCA技术记录，未将其分配到三个主背景。原版HCA标本合并不能当作本轮已确认的新标本关系。

## 原文门槛与分组
六canonical studies、53技术记录（旧表52生物样本记录，GSE210248两capture合并一个供者组织单位）。Atherosclerosis: GSE159677+GSE248610（来源GSE248609）；Aortic_dissection: GSE213740+GSE222318+GSE254132；Pulmonary_arterial_hypertension: GSE228644（含GSE210248、GSE228643）。此旧定义包含PHPF病例，保留原分组并记录，不将其重新解释为同一种病因。动脉瘤、瓣膜病、其他外部健康图谱不进入原主分析。

raw counts按供者×组织×condition×cell_state求和，聚合后≥20细胞；每cell_state/context先保留至少2病例和2对照单位的研究，再要求背景总计≥5病例和≥5对照单位。原代码实际门槛统计的是聚合单位，不是原05a注释中泛称的3供者门槛。校正后的供者数单独输出，不把同供者两capture当两人。GSE248609六标本保留真实四供者，旧版六供者错误不恢复。

edgeR4.0.16，filterByExpr(min.count=10)，至少100表达基因，TMM，robust estimateDisp/QLF。原prepare_design和run_one_de从归档源码直接提取、未经改写：配对donor×tissue block，非配对study block，sex/age完整度≥90%且满秩才纳入，年龄≥5不同值。原策略没有额外tissue/assay项。原函数存在无法纳入block时舍弃它的行为，本轮检测到{int(block)}个完成模型遗漏多水平block；结果需据此解释，不静默修复再称原文复现。

至少2背景才计算REML。严格保留原文SE近似abs(logFC)/sqrt(QLF F)、metafor默认z检验与每level/state内MAGMA预选基因BH。该SE是近似，不能称edgeR直接提供的精确标准误，原meta不确定性处理仍需方法学审核。保留原始2199个MAGMA候选集合以固定复现的检验范围；每state实际可合并数量和BH分母随表达/背景覆盖变化。L2遵循“已登记资格映射”原则使用当前103标签，不强设旧71或把旧标签改名冒充新版。

## 实际完成
原文策略DE完成{int(st.status.eq('COMPLETE').sum())}个模型；不可估计{int(st.status.ne('COMPLETE').sum())}；未过输入门槛{int(r.status.ne('READY').sum())}个组合。全基因DE及原文leave-one-study重新拟合，后者只保留MAGMA候选行。CKB有{len(new)}个背景结果、{len(mk)}个可合并meta结果。未满足表达门槛的CKB不记为零效应或不显著。

## 审核界限
这是原文策略在正确新图谱上的复现；不是对原错误图谱的逐数字重现，也不是认定原meta近似最优。原文与新版L2不能假设一对一。未运行HVG/PCA/Harmony、重注释、CoVarNet或空间关联。结果不等于CKB因果证据。
'''
(B/'REPORT_CN.md').write_text(text)
body='<html><head><meta charset="utf-8"><style>body{font:16px sans-serif;max-width:1400px;margin:30px auto;line-height:1.6}table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #ddd;padding:6px}pre{white-space:pre-wrap}img{max-width:100%}</style></head><body><h1>原文策略复现：CKB审核</h1>'
for heading,df in [('旧版与新版的同口径meta比较',main),('旧版各背景原始结果',oldc),('新版原文策略各背景结果',newc),('实际纳入研究与分组',summary),('全部CKB meta（含L2）',mk)]:body+='<h2>'+heading+'</h2>'+df.to_html(index=False,float_format=lambda v:f'{v:.5g}')
body+='<img src="data:image/png;base64,'+base64.b64encode((B/'tables/CKB_meta_comparison.png').read_bytes()).decode()+'"><h2>方法、区别与局限</h2><pre>'+html.escape(text)+'</pre></body></html>'; (B/'REPORT_CN.html').write_text(body)
validation=dict(status='COMPLETE_VALIDATED',DE_complete=int(st.status.eq('COMPLETE').sum()),not_estimable=int(st.status.ne('COMPLETE').sum()),ineligible=int(r.status.ne('READY').sum()),CKB_contexts=len(new),CKB_meta=len(mk),no_HCA_in_models=True,all_gene_BH_rechecked=True,target_set_BH_rechecked=True,block_omitted=int(block),result_sha256=hashes,tests=json.loads((B/'tests/PASS.json').read_text()))
save('VALIDATION.json',validation);save('COMPLETE.json',{k:v for k,v in validation.items() if k not in ['result_sha256','tests']})
with zipfile.ZipFile(B/'CKB_manuscript_strategy_review.zip','w',zipfile.ZIP_DEFLATED) as f:
 for p in [B/'REPORT_CN.html',B/'REPORT_CN.md',B/'VALIDATION.json',B/'COMPLETE.json',B/'RECEIPT.json']+list((B/'tables').glob('*'))+list((B/'evidence').glob('*'))+list((B/'code').glob('*'))+list((B/'tests').glob('*')):
  if p.is_file():f.write(p,p.relative_to(B))
print(main.to_string(index=False));print(json.dumps({k:v for k,v in validation.items() if k not in ['result_sha256','tests']}))
