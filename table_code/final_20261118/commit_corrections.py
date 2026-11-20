from pathlib import Path
import hashlib,json,os,shutil
from datetime import datetime
from zipfile import ZipFile

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
O=P/'图件编辑资料/补充表核查_20261118';S=O/'staging';B=O/'backup_before_correction'
m=json.loads((O/'correction_input_hashes.json').read_text())
v=json.loads((O/'correction_verification.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert v['status']=='VERIFIED_FOR_IN_PLACE_REPLACEMENT'
for f,h in m['affected'].items():
    assert sha(P/f)==h,('Original changed since verification',f)
    assert sha(B/f)==h,('Backup mismatch',f)
    assert sha(S/f)==v['staged_sha256'][f],('Staged file mismatch',f)
for f,h in m['protected'].items():assert sha(P/f)==h,('Protected file mismatch',f)
replaced=[]
for f in m['affected']:
    dest=P/f;tmp=dest.with_name('.'+dest.name+'.correction.tmp')
    assert not tmp.exists()
    shutil.copy2(S/f,tmp)
    assert sha(tmp)==v['staged_sha256'][f]
    os.replace(tmp,dest)
    assert sha(dest)==v['staged_sha256'][f]
    if dest.suffix in ('.docx','.xlsx'):
        with ZipFile(dest) as z:assert z.testzip() is None
    replaced.append({'file':f,'before_sha256':m['affected'][f],'after_sha256':sha(dest)})
for f,h in m['protected'].items():assert sha(P/f)==h,('Protected file changed',f)
receipt={'status':'COMPLETED_IN_PLACE','completed_at':datetime.now().astimezone().isoformat(),'replaced':replaced,
 'protected_files_unchanged':len(m['protected']),'lightweight_backup':str(B),'science_recomputed':False,
 'visual_review':'Reviewed changed manuscript, legend and supplement pages and representative edited workbook ranges; not a claim of visual review of every data row.'}
(O/'correction_delivery_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
note='''# 本轮修正完成

已直接覆盖当前稿件目录中的工作簿、Manuscript、独立 Figure Legends、Supplementary Information 及对应 PDF；Figure 目录中的图注副本同步更新。未新建投稿包。

- **表8**：192个数值和该工作表的原XML完全保留。目录、字段字典、Methods、Results和图注已明确：表8先在组织与供者标识内平均，再跨供者标识平均；Figure 1e按样本平均。
- **表59**：从5项扩充为18项已有校正家族定义，补入细胞比例、CoVarNet和空间模块等原有边界，包含28项和84项家族。未重算任何P值或FDR。
- **表60**：从13条扩充为63条软件、资源和随机种子记录；原13条保留。新增内容依据原环境、参数和执行记录，按流程区分版本，不以当前安装环境替代原环境。
- **编号与排名**：修正Table 20→当前Table 22、marker Figure 1g→当前Figure 1h。benchmark Source Data按当前8种方法更新rank为1–8；三类得分和行顺序原样保留。
- **数字格式**：表17、18的40个基因组位置或计数由数字文本改为等值整数，采用千位分隔显示；未转换基因ID、版本号或记录标识。
- **表24**：数据行高度统一为60 pt。两列超长组织清单仅取消自动换行，完整字符串仍在原单元格，选中单元格可查看完整内容。

检查覆盖1,536,904个原有非空单元格（不含本次扩充的表59、60），表8全部值与表24全部数据保持一致。保留原筛选范围、冻结窗格及引文域；191个登记保护文件的哈希未变。实际查看了修改页和工作簿代表区域，未宣称逐行目视全部数据。

原文件轻量备份、来源记录和逐项修改保存在本目录。没有重新运行分析，没有改变效应、CI、P/FDR、空间分组或图件。
'''
(O/'修正完成说明.md').write_text(note)
# Remove only this run's redundant full export and automatic inspection dump.
for name in ('artifact_edited.xlsx.inspect.ndjson','artifact_edited.xlsx'):
    p=S/name
    if p.exists():p.unlink()
# The accepted files now exist at their formal locations; retain only lightweight QA materials here.
for f in m['affected']:
    p=S/f
    if p.exists() and sha(p)==v['staged_sha256'][f]:p.unlink()
print(json.dumps({'status':receipt['status'],'files_replaced':len(replaced),'protected_files_unchanged':len(m['protected']),'report':str(O/'修正完成说明.md')},ensure_ascii=False))
