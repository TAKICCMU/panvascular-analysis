# CPU-Harmony 与 GPU-Harmony 补充 benchmark

状态：PAUSED_FOR_CPU_HARMONY_BENCHMARK_REVIEW

## 输入与计算

复用同一冻结 Primary：315,949 个细胞、247 个 sample_uid、2000 个历史 HVG、同一 30 维 PCA。所有指标仍只评价同一 182,968 个可信 L1 细胞；未使用全量 atlas 或新 L1/L2 标签。验证通过。

CPU Harmony 实际后端：HarmonyCPU；耗时 1351.9 秒。theta=1、lambda 默认 1、max_iter_harmony=20、seed=0；其余参数详见 audit/harmony_config.json。CPU 与 GPU 实现的数值运算及停止轮次可能不同；未更改停止规则。历史代码未单独保存当时 Harmony 源码哈希，因此只能核验调用脚本、包版本及当前默认参数，不能声称历史源码逐字节相同。

邻居搜索全部使用 CPU PyNNDescent，沿用原 8 jobs、seed=0、官方 15/50/90 指标邻域。原 NMI/ARI 为 KMeans，不使用 Leiden，因此本轮没有新增 Leiden。评分仍用原 scib-metrics 0.5.9 与原兼容环境；JAX 指标计算可用 GPU，与 CPU Harmony/CPU 邻居计算分开。

## 原始指标与排名

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>backend</th>
      <th>biology_score</th>
      <th>batch_score</th>
      <th>overall_score</th>
      <th>rank</th>
    </tr>
    <tr>
      <th>Embedding</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>scpoli</th>
      <td>historical frozen; see provenance</td>
      <td>0.831930</td>
      <td>0.440744</td>
      <td>0.675456</td>
      <td>1</td>
    </tr>
    <tr>
      <th>scanvi</th>
      <td>historical frozen; see provenance</td>
      <td>0.784886</td>
      <td>0.441221</td>
      <td>0.647420</td>
      <td>2</td>
    </tr>
    <tr>
      <th>Harmony_GPU</th>
      <td>OmicVerse PyTorch CUDA</td>
      <td>0.703117</td>
      <td>0.534464</td>
      <td>0.635656</td>
      <td>3</td>
    </tr>
    <tr>
      <th>Harmony_CPU</th>
      <td>OmicVerse NumPy CPU</td>
      <td>0.702881</td>
      <td>0.533952</td>
      <td>0.635309</td>
      <td>4</td>
    </tr>
    <tr>
      <th>scvi</th>
      <td>historical frozen; see provenance</td>
      <td>0.728168</td>
      <td>0.454049</td>
      <td>0.618521</td>
      <td>5</td>
    </tr>
    <tr>
      <th>scalex</th>
      <td>historical frozen; see provenance</td>
      <td>0.634028</td>
      <td>0.529064</td>
      <td>0.592042</td>
      <td>6</td>
    </tr>
    <tr>
      <th>pca</th>
      <td>historical frozen; see provenance</td>
      <td>0.737065</td>
      <td>0.363317</td>
      <td>0.587566</td>
      <td>7</td>
    </tr>
    <tr>
      <th>liger</th>
      <td>historical frozen; see provenance</td>
      <td>0.596185</td>
      <td>0.483552</td>
      <td>0.551132</td>
      <td>8</td>
    </tr>
    <tr>
      <th>scanorama</th>
      <td>historical frozen; see provenance</td>
      <td>0.661424</td>
      <td>0.262971</td>
      <td>0.502043</td>
      <td>9</td>
    </tr>
  </tbody>
</table>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Harmony_GPU</th>
      <th>Harmony_CPU</th>
      <th>CPU_minus_GPU</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>isolated_labels</th>
      <td>0.531380</td>
      <td>0.529758</td>
      <td>-0.001622</td>
    </tr>
    <tr>
      <th>nmi_ari_cluster_labels_kmeans_nmi</th>
      <td>0.754707</td>
      <td>0.754896</td>
      <td>0.000188</td>
    </tr>
    <tr>
      <th>nmi_ari_cluster_labels_kmeans_ari</th>
      <td>0.609617</td>
      <td>0.609867</td>
      <td>0.000249</td>
    </tr>
    <tr>
      <th>silhouette_label</th>
      <td>0.619879</td>
      <td>0.619882</td>
      <td>0.000003</td>
    </tr>
    <tr>
      <th>clisi_knn</th>
      <td>1.000000</td>
      <td>1.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>bras</th>
      <td>0.786350</td>
      <td>0.786162</td>
      <td>-0.000188</td>
    </tr>
    <tr>
      <th>ilisi_knn</th>
      <td>0.029438</td>
      <td>0.029419</td>
      <td>-0.000018</td>
    </tr>
    <tr>
      <th>kbet_per_label</th>
      <td>0.606955</td>
      <td>0.607554</td>
      <td>0.000599</td>
    </tr>
    <tr>
      <th>graph_connectivity</th>
      <td>0.839826</td>
      <td>0.836832</td>
      <td>-0.002993</td>
    </tr>
    <tr>
      <th>pcr_comparison</th>
      <td>0.409750</td>
      <td>0.409794</td>
      <td>0.000044</td>
    </tr>
    <tr>
      <th>Bio conservation</th>
      <td>0.703117</td>
      <td>0.702881</td>
      <td>-0.000236</td>
    </tr>
    <tr>
      <th>Batch correction</th>
      <td>0.534464</td>
      <td>0.533952</td>
      <td>-0.000511</td>
    </tr>
    <tr>
      <th>Total</th>
      <td>0.635656</td>
      <td>0.635309</td>
      <td>-0.000346</td>
    </tr>
  </tbody>
</table>

综合最好：scpoli；biology 最好：scpoli；batch 最好：Harmony_GPU。CPU Harmony 排名 4，GPU Harmony 排名 3。其他模型没有重训或重评分，所有旧指标逐值保留；排名变化仅来自插入 CPU 这一行。

## 图像与解释

CPU 和 GPU 对照 UMAP 均为相同 182,968 个评价细胞，在各自官方 15-NN 子图上使用相同 CPU UMAP 设置（min_dist=0.3、500 epochs、seed=0）新生成，仅用于可视化。GPU 方法仅复用旧 embedding/邻居，没有重新整合。两个 UMAP 独立布局，不能用坐标直接位移评价整合好坏。图按可靠 L1、study 和 sample_uid 着色。

Embedding 原坐标 RMSE=0.0316891，仅作数值差异诊断。

## 解释边界

主分数使用官方原生指标尺度与 60% biology / 40% batch，不进行跨方法 min-max。未额外搜索 K 或 Harmony 参数。可信标签仅覆盖 57.91% 的 Primary，本结果不能代表所有 Unknown 的表现；scANVI/scPoli 使用同一 L1 监督，其 L1 biology 并非完全独立验证。历史 benchmark HVG 的 38 个基因变化及其依赖更新问题仍保留，本轮为公平比较完全复用历史输入，不自动更新 HVG 或排名。CPU atlas 的 max_iter=100 与本 benchmark 的 20 不同，因此本补充分数不等于对全量 CPU atlas 参数的直接验证。

历史文件哈希保持不变；没有修改 L1、L2、全量 atlas 或默认指针。
