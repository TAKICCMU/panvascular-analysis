# 单细胞 integration benchmark 修订版结果

综合表现最好：**scPoli**；biology 最好：**scPoli**；batch removal 最好：**Harmony**。
无监督整合方法中，综合得分最高为 **Harmony**（PCA 作为未整合基线另列）。

| 排名 | 方法 | Biology | Batch | Overall |
|---:|---|---:|---:|---:|
| 1 | scPoli † | 0.8319 | 0.4407 | 0.6755 |
| 2 | scANVI † | 0.7849 | 0.4412 | 0.6474 |
| 3 | Harmony | 0.7031 | 0.5345 | 0.6357 |
| 4 | scVI | 0.7282 | 0.4540 | 0.6185 |
| 5 | SCALEX | 0.6340 | 0.5291 | 0.5920 |
| 6 | PCA | 0.7371 | 0.3633 | 0.5876 |
| 7 | LIGER | 0.5962 | 0.4836 | 0.5511 |
| 8 | Scanorama | 0.6614 | 0.2630 | 0.5020 |

所有主指标使用共同冻结可信 L1 子集：182,968/315,949 cells（57.91%）。Unknown 不作为真实细胞类型计分；各方法完整 embedding 仍保留全部 Primary cells。该排名描述可信标签子集，不能外推为全部 Unknown 细胞整合质量。

主分析遵循 scib-metrics 0.5.9 官方指标及默认邻居设置：graph connectivity 15、kBET 输入图 50、LISI 90；未执行可选 K=10/20/30，因此不能声称 K 排名稳定。未逐方法优化 K。

独立算术核验通过：5 个 biology 指标等权均值、5 个 batch 指标等权均值，Overall=0.6×Biology+0.4×Batch；八方法指标完整，无缺失分母变化。主排名使用原生指标范围。
跨方法 min-max 辅表会改变排名，应以原生分数和逐项指标解释结果。

† scANVI 和 supervised scPoli 使用同一可信 L1 监督，L1 biology 分数属于标签知情的评价，不能宣称完全独立生物学验证。

CellHint 为 graph-only 输出，未纳入本次统一 embedding-based 八方法排名；此项排除范围需在论文方法及图注明确披露。

本次汇总检查未发现缺失指标静默平均或人为 min-max 主排名等计分错误；这里只核验输出表算术与评价覆盖，细胞一致性仍依赖评分 runner 的成员 hash 校验，不能替代上游数据身份、marker 注释与训练实现审计。保留监督标签循环评价、可信子集覆盖和 graph-only 方法排除的适用范围说明即可。合理的方法路线差异不要求继续无限追加 sensitivity analysis。

原始指标及所有数值见 metrics/raw_scib_metrics.tsv；本报告不生成或推断 p 值。
