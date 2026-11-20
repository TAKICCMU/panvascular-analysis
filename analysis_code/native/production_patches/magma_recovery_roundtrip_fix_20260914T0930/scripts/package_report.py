from pathlib import Path
import json,datetime,hashlib,zipfile
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');W=R/'production_launches/magma_recovery_roundtrip_fix_20260914T0930';C=W/'magma_adapter_v3';pid=298454;manifest=json.loads((W/'LAUNCH_MANIFEST.json').read_text());runs=list((W/'outputs/runs').glob('*')) if (W/'outputs/runs').is_dir() else [];run=str(runs[0]) if len(runs)==1 else 'Pending native authorization/preflight';phase='native_authorization_preflight';snapshot={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'shell_pid':298450,'dispatcher_pid':pid,'command':Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0',b' ').decode() if Path(f'/proc/{pid}/cmdline').exists() else 'process exited; inspect status','run_dir_at_snapshot':run,'scientific_completed':False,'resource_gate':json.loads((W/'RESOURCE_GATE.json').read_text())};(W/'ACTUAL_LAUNCH.json').write_text(json.dumps(snapshot,indent=2))
text=f'''# MAGMA r3修复与实际启动说明

最新工作目录：`{W}`。此次只续跑MAGMA及其cohort读取依赖，不重跑五个描述模块，不执行scDRS、临床分析或CoVarNet。

实际启动：shell PID298450，轻控制dispatcher PID298454。启动时间2026-09-14 09:35:07+08；gate通过时disk free74.61GiB，MemAvailable643.19GiB。科学worker12线程；magma父worker为1线程，控制入口为1线程，避免空闲父进程也保留12线程池。新run在本次报告快照为：`{run}`。进度须以真实RUN_STATUS和RECOVERY_PROGRESS为准，PID存在不代表科学计算完成。

完整启动命令：
```bash
nohup bash '{W}/start_magma.sh' > '{W}/logs/magma_dispatcher.log' 2>&1 < /dev/null &
```
实际Python命令：
```bash
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B '{W}/scripts/run_magma_light.py' --config '{W}/configs/magma_config_APPROVED.yaml' --execute --modules cohort magma
```

代码：`{C}`。code SHA256：`{manifest['code_sha256']}`。
配置：`{W}/configs/magma_config_APPROVED.yaml`。
生产回执：`{W}/approvals/magma_approval.json`，绑定代码、实际配置、atlas、来源依赖、复用证书及测试。
日志：`{W}/logs/magma_dispatcher.log`及新run/magma.log、magma/runtime.log、magma/commands.jsonl。
科学产物：`{W}/outputs/runs/`。

修复内容：
1. MAGMA官方manual v1.10第21和24页明确，condition=Average会在每个MODEL中输出一条Average条件行；此前过滤之前按VARIABLE做一对一merge错误。保留所有raw .gsa.out，明确分离条件行，校验(VARIABLE,MODEL)唯一、每目标单模型、每模型单目标及一个Average；未知变量、缺失或重复均拒绝。目标表仍用one_to_one合并，BH仅在同层/场景/规格的目标变量内部计算。
2. 复用表的TSV不保存DataFrame轴名。上轮因columns.name None vs variable报错，该失败未证明科学数值改变。r3仅规范化CSV不存储的轴名/attrs，保留行列身份、顺序、dtype、shape和浮点序列化容差验证。全量旧表的数值等价仍由正式运行重新核验，不预先宣称已完成。

本轮实际测试：20/20 PASS。10项解析回归覆盖真实15模型gsa、Average分离、BH对照、重复目标/模型、未知变量、缺失目标/条件和非有限P。另10项完整闭环测试使用冻结counts中的7500细胞×37336基因（15 L1各500细胞，保留至少3 study且每study20的原门槛），共15,164,153非零值，实际specificity→TSV→read_csv→run_recovery走通。数值、列身份/顺序、行顺序和dtype改动必须失败。完整12路径测试的11次native调用采用已保存真实gsa fixture，不是新科学富集结果；测试L2仅为接口fixture，不写入atlas。

复用与新计算：认证保留旧L1 cell-weighted specificity及primary_r0.1_L1_cell_weighted原始gsa/log；对旧config/code/atlas/preflight、native command及每个产物SHA逐项绑定。L1 study-balanced此前未落盘，须由冻结raw补算；L2两种规格仍需计算。目标是3场景×2层×2规格，共12个family，其中1个native输出复用，11次新的MAGMA调用。原失败版本所有日志、配置、结果保留不覆盖。

科学进度查看：
```bash
/home/anzhen_vhd/miniconda3/envs/omicverse/bin/python -B '{R}/production_launches/execution_control_20260914/check_traits_repair.py'
find '{W}/outputs/runs' -name RUN_STATUS.json -exec cat {{}} \\\;
find '{W}/outputs/runs' -name RECOVERY_PROGRESS.json -exec cat {{}} \\\;
```
'''
(W/'MAGMA_R3_REPAIR_START_REPORT_CN.md').write_text(text,encoding='utf-8')
files=[]
for folder in ['scripts','evidence','configs','approvals']:
 files.extend(p for p in (W/folder).rglob('*') if p.is_file() and p.stat().st_size<5*1024**2)
files.extend([W/'MAGMA_R3_REPAIR_START_REPORT_CN.md',W/'ACTUAL_LAUNCH.json',W/'LAUNCH_MANIFEST.json',W/'RESOURCE_GATE.json',W/'start_magma.sh',C/'src/magma_recovery.py',C/'scripts/external_run_module.py',C/'scripts/run_downstream.py'])
z=W/'MAGMA_r3_repair_start_review.zip'
with zipfile.ZipFile(z,'w',compression=zipfile.ZIP_DEFLATED) as archive:
 for p in files:archive.write(p,str(p.relative_to(W)))
(W/'MAGMA_r3_repair_start_review.zip.sha256').write_text(hashlib.sha256(z.read_bytes()).hexdigest()+'  '+z.name+'\n');print(str(z),z.stat().st_size)
