from pathlib import Path
import json,hashlib,zipfile
W=Path(__file__).resolve().parents[1]
command=f"/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B {W}/code/production_runner_r2.py --receipt {W}/approvals/scdrs_APPROVAL_RECEIPT_r3.json --execute"
status=json.loads((W/'STAGE_STATUS.json').read_text());resource=json.loads((W/'RESOURCE_STATUS.json').read_text())
report=f"""# scDRS 正式启动回执

已于2026-09-14启动PID301481。当前实际状态为 `{status['status']}`，不是仅READY。原生全量评分尚未完成，不报告全量遗传结论。

完整命令：
```bash
{command}
```

- 当前生产日志：`{W}/production_r2.log`
- 阶段：`{W}/STAGE_STATUS.json`
- 动态资源：`{W}/RESOURCE_STATUS.json`
- 派生完成后：`{W}/outputs/derived/DERIVED_OBJECT.json`
- 每场景与命令：`{W}/outputs/native/CURRENT_STAGE.json` 与 `COMMANDS.json`
- 原生输出：`{W}/outputs/native/{{scenario}}/scores/`，各L1/L2 group目录
- 最终退出：`{W}/EXIT_STATUS.json`

首次资源门控project39线程预算通过（已SIGSTOP的旧进程单独审计，不恢复或终止）；磁盘和内存通过。每个原生场景前重新门控。三场景严格串行，n_ctrl1000，native CPU12，控制器1线程。

审阅包包含新代码/资源修订、旧adapter来源/diff、实际小测日志与schema、config、各版本receipt、映射审核。不包含counts、派生h5ad、全量cov、控制分数大表。科学模型保持原生，frozen_global_n30明确保留。完整小测细节见`evidence/SCDRS_ADAPTER_REVIEW_CN.md`。
"""
(W/'SCDRS_START_REPORT_CN.md').write_text(report)
(W/'evidence/START_STAGE_SNAPSHOT.json').write_text(json.dumps(status,indent=2));(W/'evidence/START_RESOURCE_SNAPSHOT.json').write_text(json.dumps(resource,indent=2))
files=[]
for dirname in ['code','approvals','configs','audit_scripts','evidence']:
 files += [p for p in (W/dirname).rglob('*') if p.is_file() and '__pycache__' not in str(p) and p.suffix!='.pyc']
files += [p for p in (W/'tests').glob('*') if p.is_file() and p.suffix in ['.py','.log','.json']]
for folder in ['actual_obs_smoke_v1','actual_l2_obs_smoke_v2','actual_l2_obs_smoke_v3']:
 base=W/'tests'/folder
 files += [p for p in base.glob('*') if p.is_file() and p.suffix in ['.json','.tsv']]
 files += [p for p in (base/'native').rglob('*') if p.is_file() and p.suffix in ['.json','.log','.txt','.tsv'] and '/scores/' not in str(p)]
files += [p for p in (W/'inputs').glob('*') if p.is_file() and (p.suffix=='.json' or p.name=='gene_map_scDRS_READY.tsv')]
files += [W/'SCDRS_START_REPORT_CN.md',W/'production_r2.log']
manifest={str(p.relative_to(W)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(files))}
(W/'REVIEW_BUNDLE_SHA256.json').write_text(json.dumps(manifest,indent=2));files.append(W/'REVIEW_BUNDLE_SHA256.json')
z=W/'scDRS_repair_review_bundle_20260914.zip';assert not z.exists()
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
 for p in sorted(set(files)):f.write(p,str(p.relative_to(W)))
print(json.dumps({'zip':str(z),'bytes':z.stat().st_size,'sha256':hashlib.sha256(z.read_bytes()).hexdigest(),'files':len(set(files))},indent=2))
