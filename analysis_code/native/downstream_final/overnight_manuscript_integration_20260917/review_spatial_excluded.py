import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[k]='1'
from pathlib import Path
import json, hashlib, html
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parent
P=O/'spatial_target_excluded_joint_v1'
V=O/'spatial_target_excluded_joint_review_v1'
V.mkdir(exist_ok=False)
def read(root,name):return pd.read_csv(root/name,sep='\t')
def save(x,name):x.to_csv(V/name,sep='\t',index=False)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
c=json.loads((P/'COMPLETE.json').read_text())
iv=json.loads((O/'spatial_target_excluded_v1/validation_v2/VALIDATION.json').read_text())
assert iv['status']=='PASS' and c['reference_target_excluded'] is True
a=read(P,'section_input_audit.tsv');r=read(P,'target_model_roster.tsv');m=read(P,'target_meta_descriptive.tsv')
b=read(P,'target_biological_unit_models.tsv');s=read(P,'target_section_models.tsv')
assert len(a)==49 and a.n_spots.sum()==45348 and a.GSE.nunique()==7
assert a.raw_integer_counts.all() and a.original_full_library_verified.all() and a.QC_fallback.sum()==2
keys=['posterior','gene','cell_context','pathology','term','sensitivity']
assert len(r)==2842 and len(m)==348 and not m.duplicated(keys).any()
assert len(r[['gene','cell_context']].drop_duplicates())==29
assert not r.duplicated(['section_id','posterior','gene','cell_context']).any()
assert r.model_status.eq('FIT_OK').sum()==c['target_section_posterior_models_ok']
assert len(s)==3*r.model_status.eq('FIT_OK').sum()==c['target_model_rows']
assert np.isfinite(s[['beta','se','p']].to_numpy()).all() and (s.se>0).all() and s.p.between(0,1).all()
assert not s.duplicated(['section_id','posterior','gene','cell_context','term']).any()
assert np.allclose(s.p,2*norm.sf(abs(s.beta/s.se)),rtol=1e-7,atol=1e-13)
check=[]
for row in m.itertuples():
    mask=np.ones(len(b),dtype=bool)
    for k in keys:mask &= b[k].eq(getattr(row,k)).to_numpy()
    x=b.loc[mask];assert len(x)==row.k and x.GSE.nunique()==row.n_datasets
    w=1/x.se.to_numpy()**2;y=x.beta.to_numpy();fixed=np.sum(w*y)/w.sum()
    Q=np.sum(w*(y-fixed)**2);den=w.sum()-(w*w).sum()/w.sum()
    tau=max(0,(Q-(len(x)-1))/den) if len(x)>1 and den>0 else 0
    wr=1/(x.se.to_numpy()**2+tau);beta=np.sum(wr*y)/wr.sum();se=np.sqrt(1/wr.sum())
    vals=[beta,se,2*norm.sf(abs(beta/se)),tau,max(0,(Q-(len(x)-1))/Q)*100 if len(x)>1 and Q>0 else 0]
    assert np.allclose(vals,[row.beta,row.se,row.p,row.tau2,row.I2],rtol=1e-7,atol=1e-11)
    assert np.allclose([row.CI_low,row.CI_high],[beta-1.96*se,beta+1.96*se],rtol=1e-7,atol=1e-11)
    check.append(dict(zip(keys,[getattr(row,k) for k in keys]))|{'meta_recomputed':'PASS'})
for _,x in m.groupby(['posterior','sensitivity']):
    assert len(x)==87
    assert np.allclose(multipletests(x.p.to_numpy(),method='fdr_bh')[1],x.FDR_BH_all_29_pairs_x_3_terms,rtol=1e-7,atol=1e-13)
save(pd.DataFrame(check),'meta_recomputation.tsv')
old=read(O/'spatial_joint_evidence_v1','target_meta_descriptive.tsv')
comp=old.merge(m,on=keys,suffixes=('_old','_excluded'),validate='one_to_one',how='outer',indicator=True)
assert comp._merge.eq('both').all()
comp['beta_difference']=comp.beta_excluded-comp.beta_old
comp['sign_changed']=np.sign(comp.beta_excluded)!=np.sign(comp.beta_old)
comp['FDR05_old']=comp.FDR_BH_all_29_pairs_x_3_terms_old<.05
comp['FDR05_excluded']=comp.FDR_BH_all_29_pairs_x_3_terms_excluded<.05
save(comp,'all_348_old_vs_target_excluded.tsv')
ckb=comp[(comp.gene=='CKB')&(comp.term=='cell_abundance')].copy();assert len(ckb)==4
save(ckb,'CKB_four_settings.tsv')
rr=read(O/'spatial_joint_evidence_v1','target_model_roster.tsv')
change=rr.merge(r,on=['section_id','posterior','gene','cell_context'],suffixes=('_old','_excluded'),validate='one_to_one')
save(change[change.model_status_old!=change.model_status_excluded],'model_eligibility_changes.tsv')
counts=m.assign(FDR05=m.FDR_BH_all_29_pairs_x_3_terms<.05).groupby(['posterior','sensitivity','term'],as_index=False).agg(tested=('gene','size'),FDR05=('FDR05','sum'))
save(counts,'FDR_counts_all_settings.tsv')
source=json.loads((P/'SOURCE_SHA256.json').read_text())
for i,(path,h) in enumerate(source.items()):
    assert sha(path)==h,'Changed source: '+path
    if i%50==0:print('HASHES',i,len(source),flush=True)
mod=read(P,'CM09_dataset_descriptive.tsv');assert len(mod)==14
text=['# 目标基因排除后的空间联合分析核验与新旧对照',
'计算与统计一致性核验通过：7研究、49切片、45,348 spots、29个预设基因–细胞组合。reference与mapping均排除目标基因并独立重训；未改冻结图谱及临床信息。',
f"{c['target_section_posterior_models_ok']}/2842切片–后验–目标模型可估，{len(s)}条term估计、348条meta。旧版可估2046个，新版2044个，具体变化见model_eligibility_changes.tsv。失败/门槛不足仍保留，不填零。",
'独立复核了完整模型清单、唯一性、有限正标准误、双侧P值、全部348条DL汇总、95%CI及4个固定87项检验家族的BH。来源文件实时SHA全部匹配；这些是完整性与统计复算检查，不证明模型收敛或因果关系。',
'## CKB与成纤维细胞丰度关联',
'|后验|切片范围|旧beta|新beta|旧FDR|新FDR|新I²|',
'|---|---|---:|---:|---:|---:|---:|']
for x in ckb.itertuples():text.append(f'|{x.posterior}|{x.sensitivity}|{x.beta_old:.6f}|{x.beta_excluded:.6f}|{x.FDR_BH_all_29_pairs_x_3_terms_old:.6g}|{x.FDR_BH_all_29_pairs_x_3_terms_excluded:.6g}|{x.I2_excluded:.2f}%|')
text+=['方向在四个设置下均为负。目标排除后，q05两种切片范围以及mean去fallback达到FDR<0.05；mean全部切片未达门槛。I²仍约95%。因此可写“存在方向一致但显著性依赖分析设定、且高度异质的区域关联”，不能写所有敏感性均稳健。',
'## 全部候选与CM09',counts.to_markdown(index=False),
'上述计数是每个固定87检验家族中的显著term数，不是独立因果基因个数；不能将不同后验/QC设置相加。完整348行对照包含效应、FDR、方向变化，不选择性仅汇报显著项。',
'CM09仍使用当前唯一高置信模块，全部权重节点映射；研究表为切片指标的研究内中位数。经验P中位数不是合并P，共享解卷积可提高表观一致性，不等于独立组织学验证。',
'## 正文P38–P40修改边界',
'P38可更新为128 factors、7研究49切片45348 spots，reference250轮、各mapping30000轮、后验1000次；本次29目标双阶段排除版本已完成。',
'P39应使用本表新的CKB Fib效应，保留四设置与异质性；不沿用旧Mac/SPP1效应，也不把未纳入预设Tier1空间清单的KCNK3写成新阴性结果。',
'P40改用新版CM09，不能继承旧模块编号和病理相关显著结论。',
'目标排除降低了目标表达参与解卷积所致的直接循环性，但标签仍来自同一图谱；切片内及研究内相关、组织/疾病混杂和未确认供者关系仍存在。序列切片采用固定合并，未显式拟合其残差协方差；随后跨生物单位DL也未额外建立研究随机层。因此本报告仍是区域关联证据，不是细胞自主作用或因果验证。',
'两QC fallback切片与旧版本一致，完整counts库大小为offset，目标从自身signature中移除。未删除任何切片以追求显著性；保留all_sections和exclude_QC_fallback两套结果。']
(V/'REPORT_CN.md').write_text('\n\n'.join(text))
(V/'review.html').write_text('<!doctype html><meta charset="utf-8"><title>空间目标排除核验</title><style>body{max-width:1100px;margin:2em auto;font:16px sans-serif}pre{white-space:pre-wrap}td,th{padding:6px;border:1px solid #ddd}</style><h1>空间目标排除核验</h1><pre>'+html.escape('\n\n'.join(text))+'</pre>'+ckb.to_html(index=False))
(V/'VALIDATION.json').write_text(json.dumps({'status':'PASS','datasets':7,'sections':49,'spots':45348,'models_ok':int(c['target_section_posterior_models_ok']),'meta_rows':348,'BH_family_size':87,'all_meta_recomputed':True,'source_hashes_verified':len(source),'not_causal_or_independent_validation':True},indent=2))
files={str(f.relative_to(V)):sha(f) for f in V.iterdir() if f.is_file()}
files['review_script']=sha(__file__)
(V/'OUTPUT_SHA256.json').write_text(json.dumps(files,indent=2))
print('REVIEW_PASS',V,flush=True)
