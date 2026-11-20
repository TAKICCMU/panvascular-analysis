import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import fcntl,json,hashlib,datetime
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');W=R/'production_launches/traits_20260914T023900Z';lock=open(R/'production_launches/execution_control_20260914/.auxiliary_compute.lock','a');fcntl.flock(lock,fcntl.LOCK_EX)
import h5py,numpy as np,pandas as pd
from anndata.io import read_elem
raw=pd.read_csv(W/'evidence/raw_covariate_values.tsv.gz',sep='\t',index_col=0)
with h5py.File(R/'panvascular_CPU_theta2_L1L2_final_v1.h5ad') as f:old=np.asarray(read_elem(f['obs/n_genes_by_counts']),dtype=float);ids=read_elem(f['obs/_index'])
assert raw.index.astype(str).tolist()==list(ids);mask=old!=raw.n_genes_raw.to_numpy();pd.DataFrame({'cell_id':raw.index[mask],'old_n_genes_by_counts':old[mask],'n_genes_raw':raw.n_genes_raw.to_numpy()[mask]}).to_csv(W/'evidence/n_genes_metadata_raw_discrepancy.tsv',sep='\t',index=False)
cov=pd.read_csv(W/'inputs/scDRS_raw_verified.cov.gz',sep='\t',index_col=0);assert cov.index.equals(raw.index);assert list(cov)==['const','log1p_counts','log1p_genes','pct_mt'];assert np.isfinite(cov.to_numpy()).all();assert (cov.const==1).all();assert np.max(np.abs(cov.iloc[:,1:].mean().to_numpy()))<1e-12;assert np.allclose(cov.iloc[:,1:].std(ddof=0),1,atol=1e-12)
f=W/'evidence/SCDRS_ADAPTER_STATUS.json';a=json.loads(f.read_text());a['blockers']=[x for x in a['blockers'] if not x.startswith('Exact-cell') and not x.startswith('Need full')];a['covariates']='PASS_REAL_RAW_ALL_CELLS';a['disk_estimate_ready']=True;f.write_text(json.dumps(a,indent=2));r=json.loads((W/'evidence/REAL_COVARIATES_AND_DISK_ESTIMATE.json').read_text());s=json.loads((W/'evidence/SCDRS_NATIVE_SMOKE.json').read_text());estimate=r['estimated_float64_raw_symbols_and_graph_GiB']+3*s['estimated_full_atlas_gzip_gib_per_scenario'];(W/'evidence/COVARIATE_READBACK_TEST.json').write_text(json.dumps({'status':'PASS','n_cells':len(cov),'exact_cell_order':True,'finite':True,'constant_column_exact_one':True,'centered_mean_within_1e_12':True,'population_sd_unit_within_1e_12':True,'old_n_genes_discrepancies':int(mask.sum()),'raw_values_used_not_metadata':True},indent=2))
(W/'TRAIT_PREPARATION_REPORT_CN.md').write_text(f'''# 三个遗传场景：输入准备与实际队列状态

MAGMA原生烟雾测试通过，已建立production回执和串行资源门禁队列；并非全量结果已完成。以MAGMA_LAUNCH_MANIFEST_r2.json和MAGMA_QUEUE_STATUS.json为实际入口。场景仅primary_r0.1、primary_r0.6、leave_one_CAD_r0.1，既有GWAS基因统计SHA全部匹配。

最终var共37,336行，与既有PROPOSED顺序及symbol/Ensembl字段完全一致。36,601个有效唯一Ensembl可用于MAGMA；735个SYMBOL_FALLBACK没有canonical gene_symbol和Ensembl。MAGMA继承原算法：完整raw库大小先归一化，未映射特征不输出到MAGMA gene covariates。735不是根据结果丢弃的细胞或新筛选规则。scDRS原adapter要求每个输入特征有已核实symbol；本轮未把显示名称自动升级为canonical映射，也未猜测Ensembl。因此scDRS正式分支仍BLOCKED_INPUT。历史scDRS曾用var_names strip/upper，这只是旧实现事实，不足以在本轮无歧义映射要求下自行认证735行。

全量raw读取验证1,175,115细胞，covariates实算const、log1p_counts、log1p_genes和pct_mt；明确13个MT基因；非有限/负数/非整数检查通过。均值及population SD(ddof=0)见REAL_COVARIATES_AND_DISK_ESTIMATE.json，未用fillna(0)，未加入临床变量。raw总counts与现metadata逐细胞完全一致。811细胞的既有n_genes_by_counts与raw实测相差1–2，逐细胞差异已导出；covariates使用真实raw结果。pctMT差异不超过1.70e-6，为数值精度量级；不因此修改metadata。压缩cov已重新读回，ID顺序、有限值、中心化及单位SD测试均通过。

scDRS原生compute-score已在400×1600合成counts上实际运行，n_ctrl=1000，mean_var/vs及native raw/filter flags，保留1000个normalized control列。原生perform-downstream首次暴露缺少输出目录和connectivities；保留失败日志。新tiny对象加CPU计算的合成图并显式创建目录后，原生assoc_mcp和hetero_mcp均有有限输出。tiny结果不属于科学分析。

独立scdrs_adapter_v1仅补新派生H5AD的LZF无损压缩、复制最终冻结n30 graph到native键和输出目录创建；CPU活跃baseline未改。正式图使用neighbors_cpu_n30_connectivities及distances、完全相同cell order，不重算neighbors。该clone已语法检查；尚未做全量adapter运行，也没有scDRS生产回执。最终全obs对旧scDRS环境的完整派生schema兼容性尚需在映射解决后验证。

磁盘实测抽样3776细胞/9,870,033非零值：float64新raw-symbol对象及冻结graph约15.98GiB。合成native full-score gzip外推每场景5.16GiB，三个场景约15.48GiB，合计{estimate:.2f}GiB，另需metadata、临时文件和安全余量；估算不是精确承诺。在20GiB floor和16GiB GPU预约下，当前磁盘不足以同时安全启动scDRS。不会降低n_ctrl、关闭control输出或删除历史文件。

MAGMA r2 gate保留20GiB floor+16GiB GPU预约+4GiB MAGMA=40GiB free门槛，且须phase1精确退出标记有效、相关进程已退出、RAM>=96GiB；之后自动运行cohort+magma，1worker×12，dispatcher12。现有L1/L2、原H5AD及CPU源码只读。
''',encoding='utf-8')
fcntl.flock(lock,fcntl.LOCK_UN);print('covariate readback/report PASS',flush=True)
