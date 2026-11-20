from pathlib import Path
import json,csv,hashlib,datetime,subprocess,shlex
ROOT=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');BASE=ROOT/'production_launches';A=BASE/'execution_control_20260914';CPU=BASE/'phase1_/20260913T183720_920345Z';TRAIT=BASE/'traits_20260914T023900Z';GPU=BASE/'spatial_gpu_20260913T183753Z';RUN=CPU/'outputs/runs/20260913T183912_636388_077099775370'
MR=BASE/'magma_category_fix_20260914T0301'
def j(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
subprocess.run(['python3',str(A/'check_progress.py')],stdout=subprocess.DEVNULL,check=True)
snap=j(A/'CURRENT_PROGRESS.json');status=j(RUN/'RUN_STATUS.json');cpur=j(CPU/'approvals/cpu_approval.json');tr=j(MR/'approvals/magma_approval.json');gr=j(GPU/'configs/PRODUCTION_RECEIPT.json')
rows=[]
gpu_codes={str(p.relative_to(GPU)):sha(p) for p in sorted((GPU/'code').rglob('*.py'))}
gpu_code_digest=hashlib.sha256('\n'.join(k+':'+v for k,v in gpu_codes.items()).encode()).hexdigest()
(A/'GPU_CODE_IDENTITY.json').write_text(json.dumps({'sha256':gpu_code_digest,'encoding':'sorted relative path:sha256 joined by newline UTF8 without final newline','files':gpu_codes},indent=2))
def entry(job,role,modules,config,receipt,pid,command,gpu,cpu,log,state,science):
 return {'job_id':job.name,'role':role,'modules':modules,'config':str(config),'config_sha':sha(config),'receipt':str(receipt),'receipt_sha':sha(receipt),'code_sha':j(receipt).get('code_sha256',gpu_code_digest),'atlas_sha':cpur['atlas_sha256'],'pid':pid,'command':command,'gpu_uuid':gpu,'cpu_budget':cpu,'log':str(log),'output':str(job/'outputs'),'snapshot_utc':snap['utc'],'status':state,'science_started':science}
cpucmd=shlex.join(j(CPU/'LAUNCH_MANIFEST.json')['command']);tc=shlex.join(j(MR/'LAUNCH_MANIFEST.json')['command'])
rows.append(entry(CPU,'CPU_PHASE1','cohort,counts,markers,figures,dotplot',CPU/'configs/cpu_config_APPROVED.yaml',CPU/'approvals/cpu_approval.json','3060/3063',cpucmd,'NONE','2 workers x12; dispatcher12 reserved',CPU/'logs/cpu_dispatcher.log','COMPLETE',True))
rows.append(entry(MR,'MAGMA','cohort,magma',MR/'configs/magma_config_APPROVED.yaml',MR/'approvals/magma_approval.json','52263/52267',tc,'NONE','1 worker12 + dispatcher12',MR/'logs/magma_dispatcher.log','RUNNING (inspect current module status)',True))
rows.append(entry(GPU,'RESOURCE_WAITER','native_smoke,descriptive_reference,descriptive_mapping',GPU/'configs/config.json',GPU/'configs/PRODUCTION_RECEIPT.json','40593','/usr/bin/python3 -u '+str(GPU/'code/resource_waiter.py'),'GPU-ee506680-afc9-11fd-3ed0-99b220571975','waiter idle; GPU4 only with auxiliary lock',GPU/'logs/resource_waiter.log',j(GPU/'RUN_STATUS.json')['status'],False))
with (A/'LAUNCHED_JOBS.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
report=f'''# 首批真实执行与并行队列记录

快照时间：{snap['utc']}。这是本轮真实运行记录，不是READY配置准备。

## 已实际执行

五模块均为COMPLETE，CPU_EXIT_STATUS exit_code=0，完成时间{j(CPU/'CPU_EXIT_STATUS.json')['completed_at']}。模块真实耗时：

| 模块 | 状态 | 秒 |
|---|---|---:|
'''
for m in ['cohort','counts','markers','figures','dotplot']:report+=f"| {m} | {status[m]['status']} | {status[m].get('seconds',0):.2f} |\n"
report+=f'''
启动shell PID3060，dispatcher PID3063（现已结束）；marker worker PID5534（已结束）。完整命令及按内容绑定的CPU配置/回执保存在下列目录。
全量计数复核：1,175,115细胞、247样本；L2 Unknown118,220、Not_assessed2,832均包含在全QC分母。parent缺失=NA、parent存在而state未检出=0。206样本沿用既有生物学聚合批准，41条未批样本保留在sample统计，独立清单BIOLOGICAL_AGGREGATION_EXCLUDED_SAMPLES.tsv。inference_eligible仍0，未改任何审批。全部marker为描述性均值/检出率；global-rest与within-parent-rest分表，未运行患者级DE。

MAGMA当前独立修复版：shell PID52263，dispatcher PID52267，已通过资源门禁并进入真实执行入口。原r2 PID42599/42603在真实Categorical编码处失败并退出code3；所有失败日志保留。仅在新magma_adapter_v1/src/external_modules.py中将临时身份编码Series转object后映射；6项回归和真实百万细胞category编码测试PASS，与object输入的四份特异性输出一致。新code digest ac18f657d6197cf3229f64e9578bb14c5f85b846b4917d7f621cb404c783a006；基因范围/阈值/分母/MAGMA模型未变。三个场景限定primary_r0.1、primary_r0.6、leave_one_CAD_r0.1。实际native小测试PASS。当前不是MAGMA科学结果完成；具体worker进展见outputs/runs/*/RUN_STATUS.json。原排队PID34069/34071在尚未启动科学计算时被本任务安全停止，并保留原配置、回执和MAGMA_QUEUE_R1_SUPERSEDED.json；r2仅将磁盘门槛改为包含GPU16GiB预约的40GiB，不修改活跃baseline。

GPU：资源等待进程PID40593已提交；{j(GPU/'RUN_STATUS.json')['status']}。49/49 raw空间矩阵内容、barcode和实际gene-map检查PASS，独立门禁5/5和266-file绑定PASS。CUDA张量、tiny原生reference/mapping及真实训练均尚未执行。GPU由既有PID2760458占用，未抢占。无其他compute PID、利用率<10%、空闲显存>16GiB连续3次60秒检查且余盘>46GiB后，获取共享auxiliary锁，再自动tiny→真实reference→各GSE顺序mapping；tiny失败只阻断GPU。不是49个独立模型。新后验只保存真实means/std/q05/q95，不持久保存全部1000draws。

## 未启动分支和确切原因

scDRS全量评分未启动、未签发生产回执：
- 735个SYMBOL_FALLBACK缺少已核实canonical symbol/Ensembl；未使用显示名称猜测。MAGMA可按既有映射规则处理36,601个唯一Ensembl，并保留全raw库大小分母，所以不被此项拖延。
- 全量raw-derived covariates已完成、读回PASS：const、标准化log1p_counts/log1p_genes/pct_mt（13个明确MT基因）。811细胞旧n_genes_by_counts与raw差1–2，协变量使用raw实算，逐细胞差异保留；未写回metadata。总counts与metadata完全一致。
- 原生compute-score(n_ctrl=1000)与修复输出目录/合成邻居图后的perform-downstream小测试PASS；独立clone最小修复无损LZF、复用冻结n30图与mkdir。完整adapter/最终派生obs对旧scDRS环境尚未跑通全流程，不能声称正式评分通过。
- 实测外推raw+graph15.98GiB、三套full control scores15.48GiB，共31.46GiB另加临时余量；加GPU16GiB与20GiB运行底线，目前不能同时无条件预约。未减少n_ctrl或删除历史结果。

未批准的临床DE、临床composition推断、meta、CoVarNet/NMF、communication、coloc、targets和空间关联检验均未启动。

## 当前资源及统一预算

快照load={snap['load_average']}；可用内存{snap['mem_available_kib']/1024**2:.2f}GiB；项目盘空闲{snap['disk_free_gib']:.2f}GiB。
GPU：{snap['gpu']}。
GPU进程：{snap['gpu_processes']}。
可用affinity60核、无cgroup CPU quota上限；CPU启动时共享负载较高，采用2×12而非2×18，不改运行中配置。按最保守dispatcher12+worker24+aux4=40；phase1结束后MAGMA dispatcher12+worker12+GPU4=28。GPU准备与traits准备通过共享排他flock使用aux4。资源不足只排队本项目新任务，不终止其他任务。
GPU预约16GiB；MAGMA4GiB；未预约磁盘底线20GiB。原phase1预约10GiB，现任务已完成并保留实际产物；scDRS未预约。分支启动时重新检查，监控进程PID41817每60秒记录本次进程及后代线程/RSS、GPU和余盘，日志RESOURCE_PROGRESS.jsonl；监控不是科学计算。

## 完整执行命令（已经提交；不要重复运行）

CPU科学入口：
```bash
{cpucmd}
```
实际后台提交：
```bash
nohup bash {CPU}/start_cpu.sh > {CPU}/logs/cpu_dispatcher.log 2>&1 < /dev/null &
```
MAGMA科学入口：
```bash
{tc}
```
实际后台提交：
```bash
nohup bash {MR}/start_magma.sh > {MR}/logs/magma_dispatcher.log 2>&1 < /dev/null &
```
GPU资源队列：
```bash
nohup /usr/bin/python3 -u {GPU}/code/resource_waiter.py > {GPU}/logs/resource_waiter.log 2>&1 < /dev/null &
```
GPU通过等待门禁后分别以cell2location_GPU Python执行code/spatial_gpu.py smoke和production，显式CUDA_VISIBLE_DEVICES为上面UUID，不调用CPU隐藏CUDA的bootstrap。

## 文件地址、配置和回执
'''
for row in rows:
 report+=f"\n### {row['role']}\n\n- 配置：`{row['config']}`\n- config SHA256：`{row['config_sha']}`\n- 回执：`{row['receipt']}`\n- receipt SHA256：`{row['receipt_sha']}`\n- 日志：`{row['log']}`\n- 输出：`{row['output']}`\n"
report+=f'''
CPU实际结果run：`{RUN}`。
完整marker大表留在该run/markers/；审核ZIP只提取显示用top30（含并列）和完整表SHA/地址，不重算统计、不替换全表。UMAP PNG/PDF、counts/cohort、DotPlot各有独立小于35MiB的审核ZIP。

正式atlas路径：`{ROOT}/panvascular_CPU_theta2_L1L2_final_v1.h5ad`。
认证SHA256保持`{cpur['atlas_sha256']}`；原五模块CPU baseline code digest保持`{cpur['code_sha256']}`。母图谱、L1/L2、临床和全局坐标没有写入。新模型仅使用派生文件；只在production_launches下写本轮产物。

## 下一次查看进度

本地终端可直接执行：
```bash
ssh -o BatchMode=yes codex-shpc-gpu 'python3 {A}/check_progress.py'
```
更多日志：
```bash
ssh -o BatchMode=yes codex-shpc-gpu 'tail -n 40 {MR}/logs/magma_dispatcher.log'
ssh -o BatchMode=yes codex-shpc-gpu 'cat {GPU}/RUN_STATUS.json; tail -n 20 {GPU}/logs/resource_waiter.log'
```
本报告是上述UTC时点快照；将来真实训练PID会登记GPU/smoke_PID.json、production_PID.json，MAGMA模块状态以其run目录记录为准。
'''
(A/'EXECUTION_REPORT_CN.md').write_text(report)
(A/'EXECUTION_REPORT_SNAPSHOT.json').write_text(json.dumps({'snapshot':snap,'jobs':rows,'phase1_modules':status,'no_unapproved_inference':True},indent=2))
print(str(A/'EXECUTION_REPORT_CN.md'))
