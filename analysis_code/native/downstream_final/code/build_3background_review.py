import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import pandas as pd,json,hashlib,html,shutil,zipfile
from datetime import datetime,timezone
D=Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1')
W=D/'expanded_parallel_20260916_v1'
OUT=W/'disease_background_3groups_review_v1'
OUT.mkdir(exist_ok=False)
source=D/'covarnet_pseudobulk_launch_review_20260915/updated_scope_audit_v2/sample_sidecar_updated.tsv'
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
atlas=R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad';before=atlas.stat()
m=pd.read_csv(source,sep='\t',keep_default_na=False)
assert len(m)==247 and m.sample_uid.is_unique and m.n_cells.sum()==1175115
original=m.copy(deep=True)
groups={
 'G1':['GSE131778','GSE155512','GSE159677','GSE184073','GSE224273','GSE247238','GSE248610','GSE253904','GSE255895','GSE220774','GSE273980','PRJNA562645','GSE213128'],
 'G2':['GSE155468','GSE166676','GSE224587','GSE237230','GSE213740','GSE222318','GSE224559','GSE254132'],
 'G3':['GSE228644','GSE271435'],
 'REFERENCE_ONLY':['GSE198750','GSE266682','HCA-2245bca0-6563-4e88-ab26-f2b8f60383a7']}
labels={'G1':'动脉粥样硬化及钙化性主动脉瓣病','G2':'主动脉夹层及动脉瘤','G3':'肺高压','REFERENCE_ONLY':'通用健康参照：暂不指定疾病背景'}
lookup={study:k for k,studies in groups.items() for study in studies}
assert len(lookup)==26 and set(m.study)==set(lookup)
m['proposed_background_id']=m.study.map(lookup)
m['proposed_background_CN']=m.proposed_background_id.map(labels)
m['diagnosis_subtype_preserved']=m.cohort_disease_context
m['clinical_group_CN']=m.group_binary.map({'disease':'疾病组','healthy':'健康/参照组'})
m['context_review_status']='PENDING_USER_REVIEW'
m['pseudobulk_run_authorized']=False
m['review_notes']='归属仅用于本轮审核；是否进入模型还需按原门槛检查每个细胞群。'
m.loc[m.study.eq('GSE248610'),'review_notes']='糖尿病大血管病的钙化/未钙化血管段；暂列G1延续原稿范围，是否符合本轮G1请审核；未钙化3例保持用户已确认健康组。'
m.loc[m.study.eq('PRJNA562645'),'review_notes']='均列G1并保持现有健康/疾病分组；狭窄与年龄研究分支保留区别。'
ageing=m.study.eq('PRJNA562645')&m.cohort_disease_context.str.contains('ageing',regex=False)
m.loc[ageing,'review_notes']='年龄研究的非CAVD健康瓣膜；6/6/11岁3例及46/57/47岁3例。建议先作为描述性参照，不自动与狭窄研究对照合并；待审核。'
m.loc[m.study.eq('GSE213128'),'review_notes']='健康主动脉瓣，按组织对应G1参照；该study无病例，不能按原规则直接替其他study提供主对照。'
m.loc[m.study.eq('GSE255895'),'review_notes']='冠心病供者的无斑块胸廓内动脉，全部保持健康/参照；该study无病例，分类归G1不等于具备原规则的疾病比较资格。'
m.loc[m.study.eq('GSE220774'),'review_notes']='归G1，健康3记录仅1名供者；不能当成3名对照。此前排除仅针对CoVarNet，未扩大到本轮分类。'
m.loc[m.study.eq('GSE224587'),'review_notes']='归G2；GSM7048230/32/34保持健康，其余瘤样本保持疾病；同供者配对信息保留。'
m.loc[m.study.eq('GSE155468'),'review_notes']='归G2；2个非动脉瘤患病对照保持用户指定健康/参照，不重改分组。'
m.loc[m.study.eq('GSE228644'),'review_notes']='PAH与肺纤维化相关肺高压统一归G3，但保留具体诊断；同一健康供者的2个capture仍只对应1人。'
m.loc[m.study.eq('GSE271435'),'review_notes']='CTEPH归G3；只有病例，不能仅因合并大类而借用其他study对照。'
m.loc[m.proposed_background_id.eq('REFERENCE_ONLY'),'review_notes']='没有对应疾病病例背景的健康动脉；建议保留通用参照，不根据血管部位强行指定G1/G2/G3，不当作第四类疾病。'
assert m[original.columns].equals(original)
assert m.loc[m.sample_uid.str.contains('GSM7048230|GSM7048232|GSM7048234'), 'group_binary'].eq('healthy').all()
assert m.loc[m.study.eq('GSE255895'),'group_binary'].eq('healthy').all()
m.to_csv(OUT/'sample_background_proposal_247.tsv',sep='\t',index=False)
friendly=m[['study','sample_uid','donor_uid','proposed_background_CN','clinical_group_CN','diagnosis_subtype_preserved','tissue','age','sex','n_cells','review_notes']].copy()
friendly.columns=['研究','样本','供者','拟疾病大类','疾病或健康组','保留的具体诊断','组织','年龄','性别','细胞数','说明']
friendly.to_csv(OUT/'逐样本审核表_247.csv',index=False,encoding='utf-8-sig')
study=[]
for s,z in m.groupby('study',sort=True):
 d=z[z.group_binary.eq('disease')];h=z[z.group_binary.eq('healthy')]
 study.append({'研究':s,'拟疾病大类':z.proposed_background_CN.iloc[0],'疾病样本数':len(d),'健康参照样本数':len(h),'疾病供者数':d.donor_uid.nunique(),'健康参照供者数':h.donor_uid.nunique(),'供者总数（不重复）':z.donor_uid.nunique(),'两组共有供者数':len(set(d.donor_uid)&set(h.donor_uid)),'细胞数':int(z.n_cells.sum()),'具体诊断':'; '.join(sorted(z.cohort_disease_context.unique())),'归属不等于模型纳入':True})
study=pd.DataFrame(study);study.to_csv(OUT/'研究层面审核表_26.csv',index=False,encoding='utf-8-sig')
summ=[]
for k,z in m.groupby('proposed_background_id',sort=True):
 d=z[z.group_binary.eq('disease')];h=z[z.group_binary.eq('healthy')]
 summ.append({'大类':labels[k],'疾病样本':len(d),'健康参照样本':len(h),'全部样本':len(z),'研究数':z.study.nunique(),'去重供者':z.donor_uid.nunique(),'细胞数':int(z.n_cells.sum())})
summary=pd.DataFrame(summ);summary.to_csv(OUT/'三大背景样本汇总.csv',index=False,encoding='utf-8-sig')
issues=pd.DataFrame([
 {'事项':'糖尿病血管钙化6样本','涉及':'GSE248609（canonical study GSE248610）','建议':'延续原稿范围列入第一组；保留糖尿病大血管病诊断，不把它改写成单纯粥样硬化。','待审核':'确认第一组包含这个血管钙化分支。'},
 {'事项':'瓣膜年龄研究6个健康参照','涉及':'PRJNA562645_2，SRX20865347–SRX20865352','建议':'疾病背景可列第一组；作为描述性参照保留，先不自动充当狭窄主比较的健康组，其中3人为6、6、11岁。','待审核':'确认这6例在主疾病比较中的使用范围，不能靠扩大分组掩盖年龄或研究分支差异。'},
 {'事项':'无明确疾病对应的健康动脉52样本','涉及':'GSE198750、GSE266682、HCA项目','建议':'保留为通用健康参照，疾病背景暂不指定，不强行按血管部位归类。','待审核':'确认这批样本只做图谱/表达描述，暂不作为跨study疾病对照。'}
])
issues.to_csv(OUT/'需要审核的3项.csv',index=False,encoding='utf-8-sig')
def mdtable(df):
 rows=[[str(v).replace('|','/').replace('\n',' ') for v in r] for r in [df.columns.tolist(),*df.values.tolist()]]
 return '\n'.join(['| '+' | '.join(rows[0])+' |','| '+' | '.join(['---']*len(rows[0]))+' |']+['| '+' | '.join(r)+' |' for r in rows[1:]])
report='''# 三大疾病背景划分：待用户审核

仅整理分组建议，没有启动新的pseudobulk、meta或CKB统计；没有改冻结obs、样本健康/疾病状态或供者身份。

按用户最新指令统一为三个研究用大类：
1. 动脉粥样硬化及钙化性主动脉瓣病：动脉粥样硬化/冠状动脉相关血管病与主动脉瓣狭窄放在同一研究大类；原稿糖尿病血管钙化分支暂列此组供审核。不是把所有诊断都重新命名为“粥样硬化”。
2. 主动脉夹层及动脉瘤：夹层、胸主动脉瘤、腹主动脉瘤合并为一个背景。
3. 肺高压：PAH、肺纤维化相关肺高压、CTEPH合并为一个背景，具体诊断仍单独保留。

同一大类不再按study或血管部位拆成不同背景；study、原始研究分支、具体诊断、组织和供者仍用于审核及后续模型，不能直接不加调整地混池。

## 数据范围
'''+mdtable(summary)+'''

以上为metadata层面样本数，不是每个L1/L2的最终pseudobulk数量。同一个人可能有疾病和参照两个标本，不能把两组人数相加当总供者。52个通用健康参照不构成第四类疾病。

## 需要审核的事项
'''+mdtable(issues)+'''

## 已明确保留的决定
- GSE248609未钙化段、GSE155468的2个非动脉瘤患病对照、GSE255895全部12样本，仍按用户确认作为健康/参照。
- GSE224587近端GSM7048230/32/34为健康，同供者瘤标本为疾病，保持配对。
- GSE220774和GSM7025944/45的排除只属于CoVarNet；没有自动扩大为本轮疾病分类排除。GSE220774的3个健康记录对应1名供者。
- GSE255895虽然都归第一组健康，但原稿主比较仍要求study内有病例和参照，不能拿它给不相关的病例-only study凑对照。CTEPH-only等研究同理。
- CKB继续列为重点关注基因；本轮不根据CKB方向、P值或表达结果改变疾病分组。

## 原门槛与本轮新分类的区别
原05a脚本主分析限定6个study、3个背景；这是旧纳入范围，不应在用户新三大背景决定后仍静默沿用。新范围须本次审核后重新确定。
原05b脚本在每个细胞群中先要求每个study至少2病例和2参照，再要求疾病背景合计至少5病例和5参照pseudobulk单位。原始计数按供者×组织×条件汇总，每个单位/细胞群至少20细胞；filterByExpr min.count=10。此前“每个study各5供者”的实现过严，旧25个结果保留为单研究结果，不当成新三大背景分析。
正式建模时仍需确认配对/重复标本关系、年龄/性别完整度和模型可识别性；分类归属并不自动代表通过统计门槛。

## 状态
PAUSED_FOR_DISEASE_BACKGROUND_3GROUP_REVIEW

等待用户审核三大背景的逐样本归属和上述3项，不自动启动pseudobulk。其他不依赖本轮背景的结构/描述/遗传与mapping任务继续。
'''
(OUT/'REVIEW_CN.md').write_text(report)
body='<html lang="zh"><meta charset="utf-8"><title>三大疾病背景审核</title><style>body{font:16px sans-serif;margin:30px}table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #ddd;padding:6px}th{position:sticky;top:0;background:#edf3f8}.scroll{overflow:auto}pre{white-space:pre-wrap}</style><h1>三大疾病背景审核（尚未运行pseudobulk）</h1>'
body+='<pre>'+html.escape(report)+'</pre><h2>26个研究</h2><div class="scroll">'+study.to_html(index=False)+'</div><h2>247个样本</h2><div class="scroll">'+friendly.to_html(index=False)+'</div></html>'
(OUT/'REVIEW_CN.html').write_text(body)
after=atlas.stat();assert (before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns)==(after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
validation={'status':'PASS_METADATA_REVIEW_ONLY','approval_status':'PAUSED_FOR_DISEASE_BACKGROUND_3GROUP_REVIEW','sample_rows':len(m),'unique_samples':m.sample_uid.nunique(),'cells':int(m.n_cells.sum()),'studies':m.study.nunique(),'donors':m.donor_uid.nunique(),'three_background_assigned_samples':int(m.proposed_background_id.ne('REFERENCE_ONLY').sum()),'reference_only_samples':int(m.proposed_background_id.eq('REFERENCE_ONLY').sum()),'all_previous_fields_unchanged':m[original.columns].equals(original),'atlas_stat_unchanged':True,'pseudobulk_started':False,'source':str(source),'source_sha256':sha(source),'created_utc':datetime.now(timezone.utc).isoformat()}
(OUT/'VALIDATION.json').write_text(json.dumps(validation,indent=2,ensure_ascii=False))
(OUT/'FILES_SHA256.json').write_text(json.dumps({str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file()},indent=2))
with zipfile.ZipFile(W/'disease_background_3groups_review_v1.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in OUT.iterdir():
  if p.is_file():z.write(p,p.name)
print(summary.to_string(index=False));print(json.dumps(validation,ensure_ascii=False))
