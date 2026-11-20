# L2通讯样本间并行调整

服务器：codex-shpc-gpu。

运行目录：`/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1/expanded_parallel_20260916_v1/remaining_downstream_20260916_v1/L2_parallel_v1/`。

调整时已有64/212个样本完成，1个原生样本正在计算；剩余147个样本分为3条互斥队列，每条49个。原调度控制进程1907255暂时SIGSTOP，其正在执行的原生worker 2105597继续计算。三个新增队列与旧worker最多同时处理4个样本；旧worker结束后保持3个样本并行。

三条队列分别绑定CPU [4,8]、[9,10]、[13,14]，避免原包装脚本把所有worker绑定到同两个CPU。原LIANA、CellChat脚本、抽样、随机种子、置换1000次、bootstrap 100次均不变；每个原生步骤继续原单线程设置。全项目仍受40线程和4个重任务槽位约束。

所有已完成结果复用，待运行目录必须为空，否则停止核查。各队列仅写其独有run_id目录和独立进度文件，不竞争全局索引。新增代码及计划、旧原生代码、原输入回执以SHA256绑定；64个真实完成结果核验、互斥分区与覆盖检查、PID身份保护检查通过。

三个队列成功后，恢复任务核验全部212份COMPLETE及输出哈希，再SIGCONT原控制进程，让原逻辑生成ALL_VALID_RESULTS_INDEX.json。L2疾病差异分析沿用原依赖，在完整索引生成后启动。

监测时读取lane_0/1/2_PROGRESS.json及ACTIVE.json，并统计原L2/outputs下COMPLETE.json；原L2/PROGRESS.json在原控制器恢复前暂不更新。不要因为1907255处于T状态而重启它；这是有意暂停。出现队列失败时先核查相应样本日志和部分输出，不强行恢复原控制器，以免重复执行或争用目录。

调整前资源：当前可用CPU集合60个逻辑核；1分钟负载43.46；内存768 GiB、可用745 GiB；/data可用864 GiB；RTX5090利用率8%、显存769/32607 MiB。此原生通讯流程使用CPU，不因GPU空闲切换算法。
