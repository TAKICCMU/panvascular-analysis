from pathlib import Path
import json,hashlib,datetime,html,zipfile
import pandas as pd
from audit_scope import R,D,W,PHASE,OLD,read,sha

def main():
    scope=json.loads((W/'SCOPE_AUDIT.json').read_text())
    node=json.loads((W/'NODE_AUDIT.json').read_text())
    m=read(W/'sample_sidecar_updated.tsv')
    rev=read(W/'L2_node_eligibility_review.tsv')
    comm=read(W/'communication_old_vs_new_247.tsv')
    assert len(rev)==116 and rev.level2_state_id.is_unique
    assert len(comm)==247 and comm.sample_uid.is_unique
    assert not m.group_review_required.astype(str).str.lower().eq('true').any()
    assert m.loc[m.sample_uid.isin(['GSE224587__GSM7048230','GSE224587__GSM7048232','GSE224587__GSM7048234']),'group_binary'].eq('healthy').all()
    for donor,center,proximal in [('AAA1','GSM7048229','GSM7048230'),('AAA2','GSM7048231','GSM7048232'),('AAA3','GSM7048233','GSM7048234')]:
        a=m.set_index('sample_uid').loc[['GSE224587__'+center,'GSE224587__'+proximal]]
        assert a.donor_uid.nunique()==1 and a.donor_uid.iloc[0].endswith('::'+donor)
        assert a.group_binary.tolist()==['disease','healthy']
    findings=read(W/'grouping_findings.tsv')
    eligible=comm.new_separate_sample_eligible.astype(str).str.lower().eq('true')
    olddone=comm.old_completed.astype(str).str.lower().eq('true')
    assert (olddone & ~eligible).sum()==0
    assert int(eligible.sum())==229 and int((eligible&~olddone).sum())==26
    assert node['profiles']==223 and node['nodes_audited']==116 and not node['NMF_run']
    labels={'RECOMMEND_CORE':'建议核心节点','EXCLUDE_LOW_REPRESENTATION':'数量或覆盖不足，不纳入核心','SENSITIVITY_ONLY_BROAD_PARENT':'宽泛父类，建议仅作敏感性补充','HOLD_LIMITED_PROGRAM_REPLICATION':'表达重复性证据不足，暂缓核心纳入'}
    rev['审核建议']=rev.recommendation.map(labels)
    simple=rev[['level1_annotation','level2_annotation','n_cells','n_samples','n_studies','n_resolved_donors','审核建议','replicated_positive_genes']].rename(columns={'level1_annotation':'L1','level2_annotation':'L2节点','n_cells':'细胞数','n_samples':'样本数','n_studies':'研究数','n_resolved_donors':'供者数','replicated_positive_genes':'跨研究支持基因'})
    simple.to_csv(W/'L2网络节点审核_简表.csv',index=False,encoding='utf-8-sig')
    counts=node['recommendations']
    hold=rev[rev.recommendation.eq('HOLD_LIMITED_PROGRAM_REPLICATION')].level2_annotation.tolist()
    report=f'''# CoVarNet新范围审核与通讯样本追溯

本轮只做配置更新、读取既有结果及描述性核验。冻结h5ad、L1/L2及既有通讯结果未改；没有运行CoVarNet、敏感性模型、疾病DE或新的通讯模型。

## 1. 分组

已按用户确认将GSM7048230、GSM7048232、GSM7048234三个近端样本归为健康分析组；与各自瘤中心保持AAA1、AAA2、AAA3供者配对。三个近端属于AAA患者自身参照，不代表健康供者。

排除GSE220774的15个样本，以及GSM7025944/GSM7025945。拟纳入230个样本、1,100,684个细胞：123个疾病组样本、107个健康/参照组样本。最终临床来源的供者字段均有记录；旧记忆中GSE213128缺供者的情况已被GSA Individual映射补齐，本轮读取到HRI274922/23/24/26/27，不从样本名推断供者。既有单供者内年龄/性别字段未发现冲突；这不等于独立基因型身份验证。

### 明显需要区分的“健康”含义

|研究/样本|现有分组|来源实际含义|审核结论|
|---|---|---|---|
|GSE248609（归入GSE248610）：GSM7917903/04/05|健康|糖尿病大血管病患者的未钙化胫前动脉段|未钙化不等于无疾病；建议标记为疾病内参照，不能称一般健康|
|GSE155468：GSM4704932/33|健康|原始GEO明确写Non-aneurysmal diseased control|非动脉瘤的患病对照；不能当健康供者|
|GSE255895：12个样本|健康|冠心病CABG患者的无斑块胸廓内动脉|可作非病变组织参照，但不是健康人|
|GSE159677近端3个样本、GSE224587近端3个样本|健康|病人自身近端/邻近参照|保留配对，不能解释为疾病供者对健康供者|

其他移植/心脏疾病背景参照也逐样本列在grouping_findings.tsv。所有Non-Diseased归入健康分析组，保留原二元分组；建议使用“健康/非病变参照”并保留reference类型与供者原诊断。前两项仍需明确是否纳入统一健康对照解释。疾病关联分析需要处理study、血管床、assay和donor，不能将跨组织所有健康样本作为可交换对照。网络节点资格本轮没有使用疾病关联显著性来筛选。

原始GEO来源已存public_sources，包括原始gzip、文本和SHA256。GSE255895本身是snRNA；assay类型保留，不能把技术差异解释成疾病差异。

## 2. L2网络节点资格核验

继续使用原门槛：每profile至少100个总细胞；每节点至少500个细胞、10个有该节点的profile、3个研究。按用户指定保留独立样本，本轮profile=sample_uid，同供者关系另存donor_block。7个小样本不足100细胞，共524细胞，仅在网络输入准备中不纳入，原始atlas不删细胞。核验输入为223个profile、1,100,160个细胞。

全部116个命名L2逐项核验完成：

|审核建议|节点数|
|---|---:|
|建议核心节点|{counts.get('RECOMMEND_CORE',0)}|
|数量/覆盖不足，不进入核心|{counts.get('EXCLUDE_LOW_REPRESENTATION',0)}|
|宽泛父类，建议仅作敏感性补充|{counts.get('SENSITIVITY_ONLY_BROAD_PARENT',0)}|
|表达重复性证据不足，暂缓核心纳入|{counts.get('HOLD_LIMITED_PROGRAM_REPLICATION',0)}|

暂缓节点：{'; '.join(hold)}。

频率过关共{node['frequency_pass']}个，不能直接都叫“稳定状态”。表达核验复用已完成sample-level原始counts伪bulk：每份样本内状态对剩余同L1比较全基因文库CPM，多份同供者样本先取中位数，再在研究间描述表达方向。保留67,266条样本-基因核验记录（实际条数以NODE_AUDIT.json为准）。基因集合来自既有within-parent top marker及锚定基因，排除核糖体/线粒体家族；具体规则见NODE_REVIEW_RULES.json。

表达方向审核采用项目内明确的保守建议规则：至少两个基因分别在至少三个研究中方向为正，且至少三个研究支持多基因方向。它不是CoVarNet官方资格标准，不是独立验证，也不是要求每个marker在所有细胞出现。没有借用疾病差异P值、未运行显著性检验。主来源集中度连续报告，没有按任意集中比例自动删节点。父类剩余群不能只因名称宽泛就修改L2，其标签和分母全部保留。

当前核心候选各profile均有非零总数。候选counts、以全部atlas细胞为分母的比例矩阵已保存，仅供审核。NMF模块稳定性、留一研究、供者成块bootstrap及原有分母/父谱系敏感性尚未运行，不能把本轮描述性审核称为已经通过这些验证。

**运行前适配点**：旧src/covarnet.py默认按study×donor×tissue×condition自动合并profile，与用户现在保留独立样本的选择不一致。后续批准运行时必须在独立版本使用显式profile键和donor block，不能靠改旧aggregate_approved标志绕过；正在运行的代码不改。

## 3. 为什么通讯是203/203

实际生产prepare.py与冻结输入逐项对上：

247条样本记录 → 41条因旧aggregate_approved=False暂缓 → 206个各自独立的生物样本 → 3个不足通讯细胞类型门槛 → 203个进入LIANA和CellChat并完成。

41条包含GSE131778的7条、GSE210248的2条、GSE247238/GSE294291的22条、HCA的2条及GSE237230的8条。它们恰好是用户最新明确取样关系的41条。旧203个中不存在把多个sample简单池化成一个生物样本的情况；未进入的41条当时被直接暂缓。

原通讯门槛是：排除Unknown/Not_assessed，每个样本至少2类L1，每类至少50细胞；对合格类型每样本每类最多抽200细胞用于通讯。203/203表示原计划的203个都完成，不代表全部247个样本都被分析，也不是全部细胞不经抽样进入通讯。

原3个不达标：GSE224273 GSM7018579、GSM7018582只有T细胞达到50；GSE247238 GSM7885778（CAR8）没有L1达到50。

按新确认关系保留每个sample，复用原50细胞门槛：旧203个均仍合格；原暂缓41个中26个新合格，15个仍不足。全通讯范围为229个可分析样本、18个不足；可新增26个，旧203个成员关系不变可复用。新增分布：GSE131778 7、GSE210248 2、GSE237230 7、GSE247238/GSE294291 8、HCA 2。

本次GSE220774和两个内/外膜样本排除仅适用于CoVarNet，不自动扩展到通讯。如果未来明确让通讯也沿用这一范围，则旧203减17再加26=212。本轮只生成新增资格表，没有自动补跑、删除或重算已有通讯；分组更名本身不改变样本内L1描述通讯分数。

## 4. 文件

- sample_sidecar_updated.tsv / 研究_样本_供者_分组.csv：更新分组、原诊断、样本/供者关系。
- grouping_findings.tsv：逐样本异常及解释。
- L2_node_eligibility_review.tsv / L2网络节点审核_简表.csv：全116节点审核及建议。
- node_marker_sample_evidence.tsv / node_marker_study_replication.tsv：表达证据。
- communication_old_vs_new_247.tsv：全247样本通讯旧新版对照。
- communication_newly_eligible.tsv / communication_still_ineligible.tsv：26个可新增及18个不足样本。
- code/、日志、规则与来源哈希：可复核过程。

状态：AUDIT_COMPLETE_WAITING_FOR_USER_REVIEW。未启动CoVarNet、正式敏感性或通讯补跑。
'''
    report=report.replace('67,266条',f"{node['evidence_rows']:,}条")
    (W/'REVIEW_CN.md').write_text(report)
    h='<meta charset="utf-8"><style>body{font-family:system-ui;margin:24px}table{border-collapse:collapse}td,th{padding:7px;border:1px solid #ddd}th{background:#eef3fa;position:sticky;top:0}pre{white-space:pre-wrap;font-family:inherit;line-height:1.6}</style><h1>CoVarNet与通讯核验</h1><pre>'+html.escape(report)+'</pre><h2>全116个L2节点</h2>'+simple.to_html(index=False,escape=True)
    (W/'离线审核.html').write_text(h)
    validation={'status':'PASS_AUDIT_INTEGRITY_NOT_MODEL_APPROVAL','sample_keys_unique':True,'paired_three_center_proximal_checked':True,'excluded_17_checked':True,'published_counts_match_updated_sample_totals':True,'pseudobulk_unit_counts_match_final_state_ledger':True,'all_116_nodes_checked':True,'old_203_reconciled':True,'new_26_eligible':True,'frozen_atlas_opened_for_expression':False,'new_models_run':False,'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (W/'VALIDATION.json').write_text(json.dumps(validation,indent=2))
    sources=json.loads((W/'SOURCE_SHA256.json').read_text())
    pb=R/'pseudobulk/descriptive_final_atlas_v1/outputs'
    for p in [pb/'genes.tsv',pb/'sample_L1/units.tsv.gz',pb/'sample_L2/units.tsv.gz',pb/'sample_L1/counts_units_by_genes.npz',pb/'sample_L2/counts_units_by_genes.npz',OLD/'L2_top20_within_parent_descriptive_markers.tsv',OLD/'L2_anchor_within_parent_evidence.tsv']:
        sources[str(p)]=sha(p)
    (W/'SOURCE_SHA256.json').write_text(json.dumps(sources,indent=2))
    registry=R/'maintenance_resume_20260914/RESUME_REGISTRY.json'
    before=registry.read_bytes();d=json.loads(before)
    (W/'RESUME_REGISTRY_before_audit.json').write_bytes(before)
    d['pending_full_design_review'].update(latest_covarnet_group_config=str(W/'config.json'),latest_node_and_communication_audit=str(W/'REVIEW_CN.md'),status='AUDIT_COMPLETE_WAITING_FOR_USER_REVIEW',models_started=False)
    assert registry.read_bytes()==before
    temp=registry.with_name('RESUME_REGISTRY.audit_review.tmp');temp.write_text(json.dumps(d,indent=2,ensure_ascii=False));temp.replace(registry)
    files={str(p.relative_to(W)):sha(p) for p in W.rglob('*') if p.is_file() and '__pycache__' not in str(p) and p.name!='FILES_SHA256.json'}
    (W/'FILES_SHA256.json').write_text(json.dumps(files,indent=2,ensure_ascii=False))
    z=W.parent/'updated_scope_audit_v2_review.zip'
    with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
        for name in list(files)+['FILES_SHA256.json']:f.write(W/name,arcname='updated_scope_audit_v2/'+name)
    print(json.dumps({'review':str(W/'REVIEW_CN.md'),'zip':str(z),'zip_sha256':sha(z),'node_decisions':node['recommendations']},ensure_ascii=False))

if __name__=='__main__':main()
