# scIB 后台接管状态

九方法未完成，八方法已有阶段排名；后台依赖控制器已就绪运行。

## 进程归属

- 原主线 CellHint pilot meta3：PID `2403745`，输入114,833-cell pilot。
- 原主线 CellHint pilot meta5：PID `2403746`，同一pilot。
- 本接管 primary shared harmonization：PID `2508566`，输入315,949-cell primary，固定4核；已通过native-vs-cache最终relation/reannotation/阈值完全等价smoke。
- 本接管 continuation controller：PID `2518463`。原控制器2512021仅在等待阶段被本分支替换，以加入上游进程退出检测；未中断任何分析worker。

## 预计耗时

15:18快照原pilot meta3=210/247、meta5=208/247。按近期速度，pilot剩余约2–5小时或更久，**不是完成保证**；后段样本复杂度和随后741-neighbor native图构建可进一步延长。Primary共享harmonization并行运行，目前在247sample重排序阶段，暂无可靠总ETA。

## 控制器已接上的链

自然等待两个pilot完成 → 冻结输出副本 → 拒绝self补位/重复/无效邻居 → 适配标准self+K-1邻居约定 → 对k10/20/30计算CellHint全部适用scIB指标 → pilot统一40/60选meta3/5 → 等待共享primary harmonization → 所选primary原生integrate → scIB评估 → 九方法适用指标/共同图指标/历史分组敏感性排名。

原先8方法的pilot冻结配置保留。若加入CellHint后统一pilot重缩放改变此前最优配置，控制器将明示 `changed_previously_selected_configurations`，不能悄悄当作完全相同的参数筛选策略。此前8方法重缩放暂未改变任何最优配置。

## 监控与恢复

远程根：`/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/scib`。

- 看 `TAKEOVER_STATUS.json` 和 `logs/controller.log`。
- 上游PID退出却没有所需产物时，控制器写 `failed_requires_attention`，不会永远把失败当等待。
- 任一graph审计或指标失败时，异常留在日志，正式分数不生成；禁止自己补零/补self/skipna绕过。
- 共享primary异常查看 `logs/cellhint_primary_shared.log`；完成检查 `primary/shared_harmonization/SHARED_HARMONIZATION_COMPLETE.json`。
- 最终表：`primary/unified/unified_scib_configuration_ranking.tsv`；正式完整性：同目录 `UNIFIED_SCIB_RANKING_STATUS.json`。
- 检查控制器最后phase及parameter-selection change字段后，才据此用于最终方法选择。
- 阶段表和旧稿对照本地已生成，最终新表须再从远端同步。

所有新写入在隔离目录。原输入、PCA、旧embedding/metric只读；23份旧embedding metric cache已验证签名、sample-only实际抽样、cache表值一致并补冻结当前内容hash。


## 既有 Harmony 收敛核查（2026-09-08）

Primary 315,949-cell theta1 的原日志明确在第15轮达到 OmicVerse 原生收敛判据（阈值0.0001，cap20，sample_uid，seed0）；pilot114,833-cell theta1在第17轮达到。因此当前8方法的所选Harmony配置并非cap20截断。此处是原生日志收敛报告，既有cache未保存objective完整轨迹，不能替代绝对改变量审计；full117万细胞需31轮与这两个subset不同，不应外推。未为此补跑或修改冻结cache。完整逐配置证据见 audit/harmony_convergence_existing/EXISTING_HARMONY_CONVERGENCE_AUDIT.json。
