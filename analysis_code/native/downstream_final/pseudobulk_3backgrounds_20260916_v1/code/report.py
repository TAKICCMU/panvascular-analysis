import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,csv,hashlib,html,base64
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path(__file__).resolve().parents[1]
rows=pd.read_csv(B/'tables/model_roster.tsv',sep='\t',keep_default_na=False);records=[];models=[]
for r in rows.to_dict('records'):
 p=B/'results'/r['model_id'];f=p/'MODEL.json'
 if not f.exists():models.append({**r,'actual_status':'PENDING' if r['status']=='READY_FOR_DESIGN_CHECK' else r['status']});continue
 d=json.loads(f.read_text());models.append({**r,'actual_status':d['status']})
 if d['status']!='COMPLETE':continue
 a=pd.read_csv(p/'DE.tsv',sep='\t');assert a.gene_id.is_unique and np.isfinite(a.PValue).all() and a.PValue.between(0,1).all()
 # Independent BH check, without relying on the R output's own calculation.
 pv=a.PValue.to_numpy();order=np.argsort(pv);qsort=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1].clip(0,1);q=np.empty(len(pv));q[order]=qsort;assert np.allclose(a.FDR,q,rtol=1e-6,atol=1e-12)
 ck=a[a.gene_symbol.eq('CKB')];assert len(ck)<=1
 record={**{k:r[k] for k in ['model_id','level','context','cell_state']},'n_case_donors':d['n_case_donors'],'n_reference_donors':d['n_reference_donors'],'n_units':d['n_units'],'n_studies':d['n_studies'],'status':'TESTED' if len(ck) else 'FILTERED_LOW_EXPRESSION','covariates':json.dumps(d['covariates'],ensure_ascii=False)}
 if len(ck):record.update(ck.iloc[0].to_dict());record['fold_change']=2**record['logFC']
 records.append(record)
mods=pd.DataFrame(models);mods.to_csv(B/'tables/MODEL_STATUS.tsv',sep='\t',index=False)
res=pd.DataFrame(records)
if not len(res):raise RuntimeError('No native models completed yet')
res.to_csv(B/'tables/CKB_review_results.tsv',sep='\t',index=False)
figs=B/'figures';figs.mkdir(exist_ok=True)
z=res[(res.level=='L1')&res.status.eq('TESTED')].copy();order=sorted(z.cell_state.unique());palette={'G1':'#b74c49','G2':'#4176a6','G3':'#438668'}
fig,ax=plt.subplots(figsize=(10,max(4,len(order)*.43+1.5)))
for i,k in enumerate(['G1','G2','G3']):
 s=z[z.context==k];ys=[order.index(v)+(i-1)*.2 for v in s.cell_state];ax.scatter(s.logFC,ys,s=48,color=palette[k],label=k,alpha=.85)
 for (_,r),y in zip(s.iterrows(),ys):
  if r.FDR<.05:ax.text(r.logFC,y,' *',va='center',fontsize=12)
ax.set_yticks(range(len(order)));ax.set_yticklabels(order);ax.axvline(0,c='grey',lw=.8);ax.set_xlabel('CKB adjusted log2 fold change (disease / healthy)');ax.legend();ax.set_title('L1 CKB: * within-model all-gene FDR < 0.05');fig.tight_layout();fig.savefig(figs/'CKB_L1_logFC.png',dpi=220);fig.savefig(figs/'CKB_L1_logFC.pdf');plt.close(fig)
cols=[c for c in ['context','cell_state','logFC','fold_change','PValue','FDR','n_case_donors','n_reference_donors','status'] if c in res]
report='# 三大疾病背景 pseudobulk：CKB阶段结果\n\n状态：'+('所有计划主模型已结束。' if (B/'COMPLETE.json').exists() else '部分结果，剩余模型仍在队列或运行。')+'\n\n'
report+='G1=动脉粥样硬化及钙化性主动脉瓣病；G2=主动脉夹层及动脉瘤；G3=肺高压。logFC>0表示疾病组更高，<0表示更低。FDR为每个细胞群/背景内全部受检基因BH；不是只对CKB或显著基因校正。未通过表达过滤不代表无变化。\n\n'
report+='## 方法和解释边界\n\n原raw计数按供者×组织×条件汇总；20细胞、5独立供者/臂门槛。edgeR 4.0.16、TMM、稳健QL。保留study/重复供者、组织及assay调整空间，缺失人口学协变量按原90%完整度规则处理并逐模型记录。只有冗余设计列被移除，完全混杂的疾病效应不强行检验。\n\n外部健康对照包含在候选输入与可用拟合中，但只有健康或只有病例的研究在study调整后不能单独识别疾病效应。study_condition_information.tsv提供各study对条件设计的支持；不能把新增健康细胞数说成等量新增疾病效应证据。保留儿童健康瓣膜，年龄混杂和不同组织的共同效应假设仍是限制。跨状态共享供者，因此本轮各状态/背景的显著性不自动等同独立重复验证。\n\n'
report+='已验证 '+str(len(records))+' 个模型的完整DE表、唯一基因及BH一致性。CKB当前结果见表格。未运行新的跨背景meta；不使用abs(logFC)/sqrt(F)伪装精确标准误。\n'
(B/'CKB_REVIEW_CN.md').write_text(report)
image=base64.b64encode((figs/'CKB_L1_logFC.png').read_bytes()).decode()
body='<html lang="zh"><meta charset="utf-8"><style>body{font:16px sans-serif;margin:30px}pre{white-space:pre-wrap}td,th{border:1px solid #ddd;padding:5px}table{border-collapse:collapse;font-size:13px}</style><pre>'+html.escape(report)+'</pre><img style="max-width:100%" src="data:image/png;base64,'+image+'"><h2>L1</h2>'+res[res.level.eq('L1')][cols].to_html(index=False,float_format=lambda v:f'{v:.4g}')+'<h2>L2</h2>'+res[res.level.eq('L2')][cols].to_html(index=False,float_format=lambda v:f'{v:.4g}')+'</html>'
(B/'CKB_REVIEW_CN.html').write_text(body)
(B/'RESULT_VALIDATION.json').write_text(json.dumps({'status':'PASS_COMPLETED_MODELS','validated_models':len(records),'model_status_counts':mods.actual_status.value_counts().to_dict(),'completed_all':(B/'COMPLETE.json').exists()},indent=2))
print(res[res.level.eq('L1')][cols].to_string(index=False))
