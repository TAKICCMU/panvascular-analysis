import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,csv,hashlib
O=Path(__file__).resolve().parents[1]
def read(n):return list(csv.DictReader((O/n).open(),delimiter='\t'))
c=json.loads((O/'COMPLETE.json').read_text());a=read('section_input_audit.tsv');r=read('target_model_roster.tsv');m=read('target_meta_descriptive.tsv');mod=read('CM09_dataset_descriptive.tsv')
assert len(a)==49 and sum(int(x['n_spots']) for x in a)==45348 and len({x['GSE'] for x in a})==7
assert all(x['raw_integer_counts']=='True' and x['original_full_library_verified']=='True' for x in a)
assert sum(x['QC_fallback']=='True' for x in a)==2
assert len(r)==49*29*2 and len(m)==29*3*2*2
assert len({(x['gene'],x['cell_context']) for x in r})==29
import math
assert all(math.isfinite(float(x['beta'])) and math.isfinite(float(x['se'])) and float(x['se'])>0 for x in read('target_section_models.tsv'))
assert len(mod)==14
stats={}
for term in ['cell_abundance','pathology','interaction']:
 t=[x for x in m if x['posterior']=='q05' and x['sensitivity']=='all_sections' and x['term']==term]
 stats[term]={'tested_pairs':len(t),'FDR05':sum(float(x['FDR_BH_all_29_pairs_x_3_terms'])<.05 for x in t)}
selected=[x for x in m if x['gene'] in ['CKB','ZNF589','COL4A1','ARHGEF12'] and x['posterior']=='q05' and x['sensitivity']=='all_sections']
md=['# 新候选与CM09空间联合结果（描述性）','**状态：已完成并验证7研究、49切片、45,348 spots、29个gene-state对。该结果不是独立空间验证。**','原始counts与全库offset逐切片通过，全部posterior SHA匹配既有review；128 factors到现行标签精确映射；CM09所有NMF权重节点可映射。','生物单位数：'+str(c['biological_units'])+'（沿manifest，不等于全部已确证独立供者）。2个QC fallback切片单独标记，并已生成排除后的全部敏感性meta。',f"目标切片模型：{c['target_section_posterior_models_ok']} / {len(r)} 可估计，其余按低检测/输入门槛或拟合失败保留，不填零效应。共{c['target_model_rows']}条term估计，{len(m)}条跨单位meta（含2posterior×2QC范围）。",'## 候选联合关联','q05、全部切片：'+json.dumps(stats,ensure_ascii=False)+'。上述显著数属于固定87检验家族的描述性统计，不等于候选独立验证通过数。','|gene|cell state|term|beta|95% CI|FDR|units|datasets|','|---|---|---|---:|---|---:|---:|---:|']
for x in selected:md.append('|'+ '|'.join([x['gene'],x['cell_context'],x['term'],f"{float(x['beta']):.4g}",f"{float(x['CI_low']):.4g} to {float(x['CI_high']):.4g}",f"{float(x['FDR_BH_all_29_pairs_x_3_terms']):.4g}",x['k'],x['n_datasets']])+'|')
md+=['## CM09空间投影','14组研究×posterior、98切片×posterior记录均已输出。CM09_dataset_descriptive.tsv为切片指标在研究内的中位数；其中经验P的中位数不作为合并P，不据此宣称研究层面显著。模块区域程序关系仅作描述，无新增定向病理假设。','|study|posterior|within-spot coherence median|neighbour-lag coherence median|','|---|---|---:|---:|']
for x in mod:md.append(f"|{x['GSE']}|{x['posterior']}|{float(x['within_spot_coherence']):.4f}|{float(x['spatial_lag_coherence']):.4f}|")
md+=['## 改稿边界与方法变更','现有cell2location参考未排除29目标，因此丰度与表达可能循环；禁止写成独立验证，禁止target_excluded=true。真正目标排除mapping须在独立输出目录复算；adapter只读新provenance，不自动授予独立性。','继承原稿Poisson+空间block稳健SE与序列切片固定合并/DL汇总；同供者序列切片残差协方差未显式建模、study层共性也未建随机层，故本版汇总是描述性关联。已知组织/标本混杂、metadata_status不确定性须保留。','明确改变：每目标从自身区域signature移除，避免直接自包含；区域基因z-score在切片内计算；邻居显式查询再排自身，修复旧代码重复排自身导致少1邻居。CoVarNet分数也按切片CLR/pseudocount计算；这些不是逐字无改动复现。','候选来自本轮29对，CM09来自本轮唯一高置信分级；KCNK3非本轮Tier1未作为预设主要空间目标，没有新空间阴性结论。','完整方法、命令见METHODS_CN.md；可复算代码code/run_joint.py，原稿快照code/legacy_spatial.py，逐文件来源哈希SOURCE_SHA256.json；所有参数与真实状态见COMPLETE.json。']
md += ['## Fallback与CKB敏感性核查','两个fallback分别为GSE243179::GSM7780296_Sample_1（2615 spots）和GSE306288::GSM9759911_Patient_23V（1512 spots），后者不是512。','CKB Fib丰度项q05全部切片beta=-0.136687、FDR=0.0449355；排除fallback beta=-0.155295、FDR=0.0286694。mean分别FDR=0.110930/0.069136；I2约95%，说明显著性依赖posterior摘要，且单位效应高度异质。不能将q05单个FDR作为稳健的独立验证。']
text='\n\n'.join(md)
import re
text=re.sub(r'(?<=\|)\n\n(?=\|)','\n',text)
(O/'REPORT_CN.md').write_text(text);(O/'VALIDATION.json').write_text(json.dumps({'status':'PASS','datasets':7,'sections':49,'spots':45348,'target_pairs':29,'model_roster_rows':len(r),'meta_rows':len(m),'descriptive_FDR_counts':stats,'not_independent_validation':True},indent=2))
files={str(p.relative_to(O)):hashlib.sha256(p.read_bytes()).hexdigest() for p in O.rglob('*') if p.is_file() and p.name not in ['OUTPUT_SHA256.json','run.log']};(O/'OUTPUT_SHA256.json').write_text(json.dumps(files,indent=2))
print(json.dumps(c));print(stats)
