from pathlib import Path
import json,time,hashlib,base64,zipfile,html
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];O=R.parent/'reliable_L1_v2_20260909'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
val=json.loads((R/'audit/INPUT_VALIDATION.json').read_text());assert val['passed']
old=pd.read_csv(O/'metrics/official_scores.tsv',sep='\t',index_col=0)
new=pd.read_csv(R/'metrics/official_scores.tsv',sep='\t',index_col=0)
assert list(new.index)==['harmony_cpu'];assert list(old.columns)==list(new.columns)
combined=pd.concat([old,new]);assert np.isfinite(combined.to_numpy(float)).all()
bio=['isolated_labels','nmi_ari_cluster_labels_kmeans_nmi','nmi_ari_cluster_labels_kmeans_ari','silhouette_label','clisi_knn'];batch=['bras','ilisi_knn','kbet_per_label','graph_connectivity','pcr_comparison']
assert np.allclose(combined[bio].mean(axis=1),combined['Bio conservation']);assert np.allclose(combined[batch].mean(axis=1),combined['Batch correction']);assert np.allclose(.6*combined['Bio conservation']+.4*combined['Batch correction'],combined.Total)
combined.index=combined.index.map(lambda s:'Harmony_GPU' if s=='harmony' else 'Harmony_CPU' if s=='harmony_cpu' else s)
combined['rank']=combined.Total.rank(ascending=False,method='min').astype(int);combined=combined.sort_values('rank');combined.to_csv(R/'benchmark_all_metrics_full_precision.tsv',sep='\t',index_label='method')
comp=combined[['Bio conservation','Batch correction','Total','rank']].rename(columns={'Bio conservation':'biology_score','Batch correction':'batch_score','Total':'overall_score'})
comp.insert(0,'backend',[('OmicVerse NumPy CPU' if m=='Harmony_CPU' else 'OmicVerse PyTorch CUDA' if m=='Harmony_GPU' else 'historical frozen; see provenance') for m in comp.index]);comp.to_csv(R/'benchmark_CPU_Harmony_comparison.tsv',sep='\t',index_label='method')
rawold=pd.read_csv(O/'metrics/raw_scib_metrics.tsv',sep='\t',index_col=0);rawnew=pd.read_csv(R/'metrics/raw_scib_metrics.tsv',sep='\t',index_col=0);pd.concat([rawold,rawnew]).to_csv(R/'tables/raw_scib_metrics_all_methods.tsv',sep='\t')
delta=pd.DataFrame({'Harmony_GPU':combined.loc['Harmony_GPU',bio+batch+['Bio conservation','Batch correction','Total']], 'Harmony_CPU':combined.loc['Harmony_CPU',bio+batch+['Bio conservation','Batch correction','Total']]});delta['CPU_minus_GPU']=delta.Harmony_CPU-delta.Harmony_GPU;delta.to_csv(R/'tables/Harmony_CPU_vs_GPU_metric_differences.tsv',sep='\t')
x=np.load(R/'embeddings/harmony_cpu.npy');y=np.load(R/'embeddings/harmony_gpu.npy');assert x.shape==y.shape==(315949,30);diff=x.astype(float)-y;emb={'shape':list(x.shape),'raw_coordinate_RMSE':float(np.sqrt(np.mean(diff**2))),'mean_cell_L2_difference':float(np.linalg.norm(diff,axis=1).mean()),'note':'Same original PCA coordinates; embedding distance is diagnostic, not an integration-quality metric'};dump(R/'audit/embedding_comparison.json',emb)
runtime=json.loads((R/'audit/harmony_runtime.json').read_text())
original=json.loads((R/'audit/historical_hashes.json').read_text());changed=[p for p,h in original.items() if sha(p)!=h];dump(R/'audit/HISTORICAL_IMMUTABILITY.json',{'passed':not changed,'n_files':len(original),'changed':changed});assert not changed
for meth in ['harmony_cpu','harmony_gpu']:assert (R/f'audit/{meth}_visualization.json').exists()
text=f'''# CPU-Harmony 与 GPU-Harmony 补充 benchmark\n\n状态：PAUSED_FOR_CPU_HARMONY_BENCHMARK_REVIEW\n\n## 输入与计算\n\n复用同一冻结 Primary：315,949 个细胞、247 个 sample_uid、2000 个历史 HVG、同一 30 维 PCA。所有指标仍只评价同一 182,968 个可信 L1 细胞；未使用全量 atlas 或新 L1/L2 标签。验证通过。\n\nCPU Harmony 实际后端：{runtime['backend']}；耗时 {runtime['seconds']:.1f} 秒。theta=1、lambda 默认 1、max_iter_harmony=20、seed=0；其余参数详见 audit/harmony_config.json。CPU 与 GPU 实现的数值运算及停止轮次可能不同；未更改停止规则。历史代码未单独保存当时 Harmony 源码哈希，因此只能核验调用脚本、包版本及当前默认参数，不能声称历史源码逐字节相同。\n\n邻居搜索全部使用 CPU PyNNDescent，沿用原 8 jobs、seed=0、官方 15/50/90 指标邻域。原 NMI/ARI 为 KMeans，不使用 Leiden，因此本轮没有新增 Leiden。评分仍用原 scib-metrics 0.5.9 与原兼容环境；JAX 指标计算可用 GPU，与 CPU Harmony/CPU 邻居计算分开。\n\n## 原始指标与排名\n\n{comp.to_html()}\n\n{delta.to_html()}\n\n综合最好：{comp.overall_score.idxmax()}；biology 最好：{comp.biology_score.idxmax()}；batch 最好：{comp.batch_score.idxmax()}。CPU Harmony 排名 {int(comp.loc['Harmony_CPU','rank'])}，GPU Harmony 排名 {int(comp.loc['Harmony_GPU','rank'])}。其他模型没有重训或重评分，所有旧指标逐值保留；排名变化仅来自插入 CPU 这一行。\n\n## 图像与解释\n\nCPU 和 GPU 对照 UMAP 均为相同 182,968 个评价细胞，在各自官方 15-NN 子图上使用相同 CPU UMAP 设置（min_dist=0.3、500 epochs、seed=0）新生成，仅用于可视化。GPU 方法仅复用旧 embedding/邻居，没有重新整合。两个 UMAP 独立布局，不能用坐标直接位移评价整合好坏。图按可靠 L1、study 和 sample_uid 着色。\n\nEmbedding 原坐标 RMSE={emb['raw_coordinate_RMSE']:.6g}，仅作数值差异诊断。\n\n## 解释边界\n\n主分数使用官方原生指标尺度与 60% biology / 40% batch，不进行跨方法 min-max。未额外搜索 K 或 Harmony 参数。可信标签仅覆盖 57.91% 的 Primary，本结果不能代表所有 Unknown 的表现；scANVI/scPoli 使用同一 L1 监督，其 L1 biology 并非完全独立验证。历史 benchmark HVG 的 38 个基因变化及其依赖更新问题仍保留，本轮为公平比较完全复用历史输入，不自动更新 HVG 或排名。CPU atlas 的 max_iter=100 与本 benchmark 的 20 不同，因此本补充分数不等于对全量 CPU atlas 参数的直接验证。\n\n历史文件哈希保持不变；没有修改 L1、L2、全量 atlas 或默认指针。\n'''
(R/'Harmony_CPU_vs_GPU_report.md').write_text(text)
(R/'README.md').write_text('# CPU-Harmony 补充分支\n\n入口：Harmony_CPU_vs_GPU_report.md；离线图集 review.html；全指标 benchmark_all_metrics_full_precision.tsv。输入均为只读历史路径链接。代码与适配差异位于 code/，全部审计位于 audit/。仅新增 CPU Harmony；原 benchmark 不重命名、不覆盖。\n\n指标邻居搜索与 UMAP 为 CPU；Leiden 未执行，原指标使用 KMeans。评分 JAX 保持历史 GPU 后端以保持一致。\n')
parts=['<html><meta charset="utf-8"><title>CPU Harmony 补充审核</title><body><h1>CPU Harmony 补充 benchmark</h1>',comp.to_html(),delta.to_html(),'<pre style="white-space:pre-wrap">'+html.escape(text)+'</pre>']
for p in sorted((R/'figures').glob('*.png')):parts+=['<h2>'+p.name+'</h2>','<img style="max-width:100%" src="data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()+'">']
parts+=['</body></html>'];(R/'review.html').write_text('\n'.join(parts))
dump(R/'STATUS.json',{'stage':'PAUSED_FOR_CPU_HARMONY_BENCHMARK_REVIEW','completed_unix':time.time(),'historical_unchanged':True,'L1_L2_modified':False})
files=[p for p in R.rglob('*') if p.is_file() and not p.is_symlink() and 'cache' not in p.relative_to(R).parts and p.name not in ['SHA256SUMS.txt','review_package.zip'] and p.suffix not in ['.lock']]
(R/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(R))+'\n' for p in sorted(files)))
with zipfile.ZipFile(R/'review_package.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for p in files+[R/'SHA256SUMS.txt']:
  if p.suffix not in ['.npy','.npz'] and p.stat().st_size<100_000_000:z.write(p,p.relative_to(R))
print(comp.to_string(),flush=True)
