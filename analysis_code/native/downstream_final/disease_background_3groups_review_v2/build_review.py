from pathlib import Path
import csv,json,hashlib,html,collections,sys,zipfile
from datetime import datetime,timezone
base=Path(sys.argv[1]) if len(sys.argv)>1 else Path('.')
src=base/'disease_background_3groups_review_v1/sample_background_proposal_247.tsv'
out=base/'disease_background_3groups_review_v2';out.mkdir(exist_ok=True)
assert not (out/'VALIDATION.json').exists(), 'Do not overwrite completed release'
r=list(csv.DictReader(src.open(),delimiter='\t'));old=[x.copy() for x in r]
labels={'G1':'动脉粥样硬化及钙化性主动脉瓣病','G2':'主动脉夹层及动脉瘤','G3':'肺高压','EXCLUDED_UNMATCHED':'排除主分析：无相应部位病例'}
origcols=list(r[0])[:list(r[0]).index('proposed_background_id')]
match={'carotid artery':'G1','coronary artery':'G1','aorta':'G2','ascending aorta':'G2','descending aorta':'G2','aortic arch':'G2','abdominal aorta':'G2','pulmonary artery':'G3'}
refs=[]
for x in r:
 isref=x['proposed_background_id']=='REFERENCE_ONLY'
 x['previous_background_id']=x['proposed_background_id']
 if isref:
  k=match.get(x['tissue'],'EXCLUDED_UNMATCHED')
  assert x['tissue'] in match or x['tissue_detail'] in ['brachial artery','iliac artery']
  x['proposed_background_id']=k
 x['proposed_background_CN']=labels[x['proposed_background_id']]
 x['main_analysis_scope_candidate']=str(x['proposed_background_id']!='EXCLUDED_UNMATCHED').lower()
 x['context_review_status']='USER_RULE_APPLIED_AWAITING_ROW_MAPPING_REVIEW'
 x['pseudobulk_run_authorized']='False'
 x['grouping_decision_source']='User approved three backgrounds and specific control rules; current row mapping pending final review'
 x['review_notes']='按三大疾病背景归类；保留具体诊断、组织、study及供者。分组范围不等于每个细胞群已满足统计门槛。'
 if isref:
  x['review_notes']='按现有组织字段匹配疾病背景；供者和具体血管部位保留，不将多个标本当成多个供者。' if x['main_analysis_scope_candidate']=='true' else '肱动脉或髂动脉没有对应部位病例：按用户规则排除主分析，仍保留图谱及描述数据。'
 if x['study']=='GSE248610':
  x['review_notes']='用户确认：钙化3样本为G1疾病组，未钙化3样本为健康组；共4位患者，保留真实配对关系。'
 if x['study']=='PRJNA562645' and 'ageing' in x['cohort_disease_context']:
  x['review_notes']='用户明确纳入主要病例对照健康组，包括年龄6、6、11、46、57、47岁；不以儿童或年龄研究分支为由自动排除，保留年龄及研究分支用于模型核验。'
 if x['study']=='GSE220774':x['review_notes']+=' GSE220774排除仅限CoVarNet，不扩展到本轮pseudobulk；健康3条记录仅1名供者。'
 if isref:refs.append(x)
assert all(all(a[k]==b[k] for k in origcols) for a,b in zip(old,r))
assert len(r)==247 and len({x['sample_uid'] for x in r})==247 and sum(int(x['n_cells']) for x in r)==1175115
z=[x for x in r if x['study']=='GSE248610'];assert collections.Counter(x['group_binary'] for x in z)=={'disease':3,'healthy':3};assert len({x['donor_uid'] for x in z})==4
age=[x for x in r if x['study']=='PRJNA562645' and 'ageing' in x['cohort_disease_context']];assert len(age)==6 and all(x['group_binary']=='healthy' and x['main_analysis_scope_candidate']=='true' for x in age)
assert collections.Counter(x['proposed_background_id'] for x in refs)=={'G1':17,'G2':23,'G3':6,'EXCLUDED_UNMATCHED':6}
def write(name,rows,sep=','):
 with (out/name).open('w',encoding='utf-8-sig' if sep==',' else 'utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter=sep);w.writeheader();w.writerows(rows)
write('sample_background_sidecar_v2.tsv',r,'\t')
def simple(x):return dict(zip(['研究','样本','供者','疾病背景','疾病或健康','原始组织','具体组织','原始诊断','主分析候选','说明'],[x[k] for k in ['study','sample_uid','donor_uid','proposed_background_CN','group_binary','tissue','tissue_detail','source_diagnosis','main_analysis_scope_candidate','review_notes']]))
write('逐样本分组表_247.csv',[simple(x) for x in r]);write('健康动脉52样本归位表.csv',[simple(x) for x in refs])
summ=[]
for k in ['G1','G2','G3','EXCLUDED_UNMATCHED']:
 z=[x for x in r if x['proposed_background_id']==k]
 summ.append({'疾病背景':labels[k],'疾病样本':sum(x['group_binary']=='disease' for x in z),'健康样本':sum(x['group_binary']=='healthy' for x in z),'样本总数':len(z),'去重供者':len({x['donor_uid'] for x in z}),'细胞数':sum(int(x['n_cells']) for x in z)})
write('三大背景更新汇总.csv',summ)
def table(rows):
 keys=list(rows[0]);return '| '+' | '.join(keys)+' |\n| '+' | '.join(['---']*len(keys))+' |\n'+'\n'.join('| '+' | '.join(str(x[k]).replace('|','/') for k in keys)+' |' for x in rows)
report='''# 疾病背景分组 v2：落实用户决定后的逐样本审核

已落实用户明确决定：
- GSE248609（汇总study为GSE248610）：3个钙化标本归第一类疾病组，3个未钙化标本归健康组。6标本来自4人，未增加供者。
- PRJNA562645_2：6个健康瓣膜全部纳入主要病例对照候选范围，包括3个儿童样本，不自动排除。
- 52个通用健康动脉按组织位置归位，无法匹配部位者排除主分析。

## 更新后的样本范围
'''+table(summ)+'''

上述是样本记录数，非每个细胞群最终合格的pseudobulk单位数。供者可跨组织和疾病背景，不能把各类供者数相加当独立总人数。

## 52个健康动脉的归位
- 17个冠状动脉/颈动脉标本 → 动脉粥样硬化背景健康组。
- 23个主动脉标本（升、降、弓、腹及其他主动脉） → 夹层/动脉瘤背景健康组。具体主动脉部位仍保留，不将主动脉组织当瓣膜组织。
- 6个肺动脉标本 → 肺高压背景健康组。肺动脉口径信息仍保留。
- 3个肱动脉标本（GSE266682）和3个髂动脉标本（HCA）无相应部位病例 → 排除主分析，仍保留图谱和描述数据。

## 建模衔接（未执行）
用户已明确批准纳入额外健康参照和6例健康瓣膜，因此旧脚本“study内必须同时有足量病例/对照”的筛选不得静默把这些新增参照全部排除。后续需在这个已更新的候选范围核验study、组织、年龄、性别与疾病变量是否可识别；无法识别的具体对比须单列报告，不能静默删去study/供者项来得到P值。
同供者不同标本保留配对/重复测量关系，健康组为本研究对照定义，不等同于供者完全没有基础病。儿童瓣膜年龄保留，不通过重命名或忽略年龄混入成人。
此处仅更新外部分析sidecar和审核材料，不更改冻结h5ad、L1/L2、临床原始字段、CoVarNet范围或已运行代码。旧分组v1保留。新pseudobulk、疾病通讯和CKB推断分析仍未启动。

状态：PAUSED_FOR_DISEASE_BACKGROUND_3GROUP_REVIEW
请审核本次52例归位明细，再执行新的pseudobulk。
'''
(out/'REVIEW_CN.md').write_text(report)
body='<html lang="zh"><meta charset="utf-8"><title>疾病背景分组v2</title><style>body{font:16px sans-serif;margin:32px}pre{white-space:pre-wrap}td,th{padding:6px;border:1px solid #ccc}table{border-collapse:collapse;font-size:13px}</style><pre>'+html.escape(report)+'</pre><h2>52个健康动脉逐条归位</h2><table>'
fr=[simple(x) for x in refs]
body+='<tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in fr[0])+'</tr>'
for x in fr:body+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in x.values())+'</tr>'
(out/'REVIEW_CN.html').write_text(body+'</table></html>')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
v={'created_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(src),'rows':247,'cells':1175115,'source_metadata_fields_unchanged':True,'donor_ids_unchanged':True,'main_scope_candidate_samples':241,'main_excluded_unmatched_samples':6,'reference52_mapping':dict(collections.Counter(x['proposed_background_id'] for x in refs)),'gse248610_samples':6,'gse248610_donors':4,'six_age_branch_valves_in_main_scope':True,'scientific_analysis_started':False,'status':'PAUSED_FOR_DISEASE_BACKGROUND_3GROUP_REVIEW'}
(out/'VALIDATION.json').write_text(json.dumps(v,indent=2,ensure_ascii=False))
(out/'FILES_SHA256.json').write_text(json.dumps({p.name:sha(p) for p in out.iterdir() if p.is_file() and p.suffix!='.zip' and p.name!='FILES_SHA256.json'},indent=2))
with zipfile.ZipFile(out/'disease_background_3groups_review_v2.zip','w',zipfile.ZIP_DEFLATED) as zz:
 for p in out.iterdir():
  if p.is_file() and p.suffix!='.zip':zz.write(p,p.name)
print(json.dumps(v,ensure_ascii=False));print(table(summ))
