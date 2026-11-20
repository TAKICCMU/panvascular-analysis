from pathlib import Path
import json,csv,hashlib,time
J=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/production_launches/spatial_gpu_20260913T183753Z')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
r=json.loads((J/'evidence/PREPARE_RECEIPT_RESULT.json').read_text());status=json.loads((J/'RUN_STATUS.json').read_text());rows=list(csv.DictReader(open(J/'evidence/SPATIAL_CONTENT_VALIDATION.tsv'),delimiter='\t'));fb=[x['section_id'] for x in rows if x['qc_fallback']=='True'];cmd=f'/usr/bin/python3 -u {J}/code/resource_waiter.py'
launch={'job_id':J.name,'role':'RESOURCE_WAITER','modules':'descriptive_reference;descriptive_mapping;native_smoke_first','config_sha':r['config_sha256'],'code_sha':sha(J/'code/resource_waiter.py'),'atlas_sha':'00fc74c1ddaa5e4b32d8664f134e0ec5a0c762f5e9a9fb45b96acbd40804a310','PID':40593,'command':cmd,'GPU_UUID':'GPU-ee506680-afc9-11fd-3ed0-99b220571975','CPU_budget':'waiter1 idle; nativeGPU4 only with auxiliary flock','log':str(J/'logs/resource_waiter.log'),'started_utc':json.loads((J/'logs/resource_waiter.log').read_text().splitlines()[0])['time_utc'],'science_started':False,'config':str(J/'configs/config.json'),'receipt':str(J/'configs/PRODUCTION_RECEIPT.json'),'output':str(J/'outputs')}
with open(J/'LAUNCHED_JOBS.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(launch),delimiter='\t');w.writeheader();w.writerow(launch)
(J/'GPU_BRANCH_HANDOFF.md').write_text(f'''# GPU描述性空间分支执行记录

状态：WAITING_RESOURCE；资源等待PID40593已真实提交。未将排队PID视为科学分析完成。

## 已实际完成
- 49/49空间matrix完整稀疏读取与内容验证；非负有限整数、shape、barcode唯一性及与position精确对应通过。
- 采用确证stable Ensembl ID优先；只有非Ensembl来源feature才允许exact唯一symbol对应。未全大写、未make_unique、未barcode去后缀。49份实际gene map与全部空间输入SHA已冻结。
- {len(fb)}个切片触发旧QC规则的低counts fallback，见下；保留为描述性/QC-review，不自动认定优质。
- 5/5独立门禁测试通过；266文件真实内容绑定guard通过；compile通过。
- 已安装API实际inspect通过，scvi sample_posterior支持return_samples=False。设置实际均值、std、q05、q95汇总，不保存全1000 draws，不将median冒充mean。
- 复用本轮CPU真实preflight的母atlas SHA和stat认证，未再次全文件扫描，也未修改atlas/metadata/L1/L2。

## 未实际执行
实际CUDA张量、tiny原生RegressionModel和Cell2location、真实reference训练、mapping均未启动。当前其他GPU计算PID2760458约90%占用。

等待条件：无其他GPU compute PID、GPU util<10%、freeGPU>16GiB连续3次（60秒间隔），并且free disk>46GiB。获得全项目auxiliary排他flock后再次检查资源；变化则释放锁回到等待循环。
通过后自动：native tiny reference+mapping → native技术检查 → 真实reference → 各GSE sequential mapping。首次tiny失败仅阻断GPU分支，不开始生产。无需CPU绘图完成后才能占用GPU。

## 资源与存储
- 最多1个GPU训练进程，GPU_UUID={launch['GPU_UUID']}。
- CPU4线程及dataloader0worker，独立CUDA_VISIBLE_DEVICES，无CPUbootstrap。
- sparse reference上界来自真实X indptr和每factor最多1500细胞的最大nnz上界，非全量dense。保守预计13.482GiB，预约16GiB；另留CPUphase1 10GiB+绝对floor20GiB。
- 真实posterior summaries保存在服务器，模型和必要raw子集留存；不上传大模型。

## 方法边界及明确修正
116个命名L2是候选reference；min40、max1500/state、max20study/state、seed20260814。Other|parent作为独立nuisance组，非新L2。实际最多128候选reference groups、最多167718入选细胞是提取前上限，不是训练已经完成。
旧采样的补足步骤可能从前20study以外补细胞；新适配明确只在所选最多20study范围补足，以符合此次授权上限，记录此差异。没有病例/组织驱动筛选。
审阅实际继承的02_spatial_v4_cell2location_GPU.py无显式ribo/MT/HB reference排除，故未静默新增排除，filter_genes仍5/0.03/1.12。
reference250epochs、mapping30000epochs、Ncells12、alpha20。模型维度/非负/finite/总量positive检查，zero-factor/zero-spot完整导出；ELBO和重构QC保留，technical complete不代表科学收敛已审核。
只做descriptive；forced_genes=[]、excluded_targets=[]。不得作为相同target的独立空间验证，不运行临床/空间关联推断。

## 关键地址
配置：{launch['config']}
回执：{launch['receipt']}
日志：{launch['log']}
输出：{launch['output']}
真实输入审核：{J}/evidence/INPUT_VALIDATION.json
完整内容表：{J}/evidence/SPATIAL_CONTENT_VALIDATION.tsv
磁盘估计：{J}/evidence/DISK_RESERVATION.json

## 启动及查看
实际启动：nohup {cmd} > {launch['log']} 2>&1 < /dev/null &
查看：cat {J}/RUN_STATUS.json; tail -n 20 {launch['log']}; ps -p 40593 -o pid,etime,pcpu,rss,args
未来真实模型PID写入{J}/smoke_PID.json及production_PID.json；真实日志smoke.log及production.log。

## QC fallback切片
'''+'\n'.join('- '+x for x in fb)+'\n')
print(json.dumps(launch,indent=2))
