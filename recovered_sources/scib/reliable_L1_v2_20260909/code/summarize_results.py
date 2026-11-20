#!/usr/bin/env python3
"""Validate final official metrics and produce concise publication-review artifacts."""
import argparse
import csv
import json
import math
from pathlib import Path

METHODS = ('pca', 'harmony', 'scvi', 'scanorama', 'liger', 'scalex', 'scanvi', 'scpoli')
BIO = ('isolated_labels', 'nmi_ari_cluster_labels_kmeans_nmi', 'nmi_ari_cluster_labels_kmeans_ari', 'silhouette_label', 'clisi_knn')
BATCH = ('bras', 'ilisi_knn', 'kbet_per_label', 'graph_connectivity', 'pcr_comparison')
SUPERVISED = {'scanvi', 'scpoli'}
PRETTY = {'pca': 'PCA', 'harmony': 'Harmony', 'scvi': 'scVI', 'scanorama': 'Scanorama', 'liger': 'LIGER', 'scalex': 'SCALEX', 'scanvi': 'scANVI †', 'scpoli': 'scPoli †'}


def read_table(path):
    with path.open() as f:
        records = list(csv.DictReader(f, delimiter='\t'))
    if not records:
        raise ValueError(f'Empty table: {path}')
    idcol = list(records[0])[0]
    result = {}
    for record in records:
        name = record.pop(idcol)
        if name in result:
            raise ValueError(f'Duplicate method {name}: {path}')
        result[name] = {k: float(v) for k, v in record.items()}
    if set(result) != set(METHODS):
        raise ValueError(f'Expected exactly eight methods: {path}')
    for method, row in result.items():
        if any(not math.isfinite(v) for v in row.values()):
            raise ValueError(f'Nonfinite value: {path}/{method}')
    return result


def rank(values):
    return {m: 1 + sum(v > values[m] for v in values.values()) for m in values}


def winners(values):
    peak = max(values.values())
    return [m for m, value in values.items() if math.isclose(value, peak, rel_tol=0, abs_tol=1e-12)]


def check_close(a, b, label):
    if not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-9):
        raise ValueError(f'Independent arithmetic mismatch {label}: {a} != {b}')


def audit(raw, official, scaled, final):
    for method in METHODS:
        if set(raw[method]) != set(BIO + BATCH):
            raise ValueError(f'Expected ten raw metrics: {method}')
        # ARI can legitimately be negative; do not clip or substitute.
        for k, value in raw[method].items():
            if not (-1 <= value <= 1 + 1e-6):
                raise ValueError(f'Unexpected native metric range: {method}/{k}={value}')
        biology = sum(raw[method][k] for k in BIO) / len(BIO)
        batch = sum(raw[method][k] for k in BATCH) / len(BATCH)
        expected = {'Bio conservation': biology, 'Batch correction': batch, 'Total': 0.6 * biology + 0.4 * batch}
        for k, value in raw[method].items():
            check_close(value, official[method][k], f'{method}/{k}')
        for k, value in expected.items():
            check_close(value, official[method][k], f'official/{method}/{k}')
            check_close(value, final[method][k], f'final/{method}/{k}')
    native_ranks = rank({m: official[m]['Total'] for m in METHODS})
    for m in METHODS:
        check_close(native_ranks[m], final[m]['rank'], f'rank/{m}')
    # Independently reproduce sklearn per-column minmax, including constant cols.
    for m in METHODS:
        normalized = {}
        for k in BIO + BATCH:
            lo = min(raw[x][k] for x in METHODS)
            hi = max(raw[x][k] for x in METHODS)
            normalized[k] = (raw[m][k] - lo) / (hi - lo) if hi > lo else 0.0
            check_close(normalized[k], scaled[m][k], f'minmax/{m}/{k}')
        bio = sum(normalized[k] for k in BIO) / len(BIO)
        batch = sum(normalized[k] for k in BATCH) / len(BATCH)
        for k, value in {'Bio conservation': bio, 'Batch correction': batch, 'Total': 0.6 * bio + 0.4 * batch}.items():
            check_close(value, scaled[m][k], f'minmax/{m}/{k}')
    scaled_ranks = rank({m: scaled[m]['Total'] for m in METHODS})
    return {'passed': True, 'n_methods': 8, 'n_metrics': 10, 'native_ranks': native_ranks, 'minmax_ranks': scaled_ranks, 'minmax_changes_rank': native_ranks != scaled_ranks, 'winners': {key: winners({m: official[m][key] for m in METHODS}) for key in ('Bio conservation', 'Batch correction', 'Total')}, 'unsupervised_integration_winner': winners({m: official[m]['Total'] for m in METHODS if m not in SUPERVISED and m != 'pca'})}


def make_figures(raw, official, target):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    order = sorted(METHODS, key=lambda m: official[m]['Total'], reverse=True)
    values = np.array([[raw[m][k] for k in BIO + BATCH] for m in order])
    # Common native 0-1 scale; any negative ARI forces a visibly negative lower
    # endpoint. No per-column z-score/minmax or clipping of legitimate values.
    low = min(0.0, float(values.min()))
    fig, ax = plt.subplots(figsize=(14, 6))
    im = ax.imshow(values, vmin=low, vmax=1, cmap='viridis', aspect='auto')
    for i in range(len(order)):
        for j in range(values.shape[1]):
            relative = (values[i, j] - low) / (1 - low)
            ax.text(j, i, f'{values[i,j]:.3f}', ha='center', va='center', fontsize=9, color='black' if relative > .6 else 'white')
    ax.set_xticks(range(10), ['Isolated ASW', 'KMeans NMI', 'KMeans ARI', 'Cell-type ASW', 'cLISI', 'BRAS', 'iLISI', 'kBET', 'Graph connectivity', 'PCR'], rotation=35, ha='right')
    ax.set_yticks(range(8), [PRETTY[m] for m in order])
    ax.axvline(4.5, color='white', linewidth=2)
    ax.set_title('Official scIB metrics: native values, higher is better')
    fig.colorbar(im, ax=ax, label='Native metric value (no across-method minmax)')
    fig.text(.01, .01, '† L1-informed supervision; L1 biology is not an independent validation.', fontsize=9)
    fig.tight_layout(rect=(0,.04,1,1))
    for ext in ('png','pdf'):
        fig.savefig(target / f'raw_metrics_heatmap.{ext}', dpi=220)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(order)); width = .25
    for offset, key, color in [(-1, 'Bio conservation', '#4878D0'), (0, 'Batch correction', '#EE854A'), (1, 'Total', '#6ACC64')]:
        data = [official[m][key] for m in order]
        bars = ax.bar(x + offset*width, data, width, label=key, color=color)
        for bar, value, m in zip(bars, data, order):
            if m in SUPERVISED:
                bar.set_hatch('//')
            ax.annotate(f'{value:.3f}', (bar.get_x()+bar.get_width()/2, value), xytext=(0,3), textcoords='offset points', ha='center', fontsize=7)
    ax.set_xticks(x, [PRETTY[m] for m in order]); ax.set_ylabel('Native aggregate score')
    ax.set_ylim(min(0, min(official[m][k] for m in order for k in ('Bio conservation','Batch correction','Total')))-.02,1.06)
    ax.legend(); ax.set_title('Biology / Batch / Overall (60% Biology + 40% Batch)')
    fig.text(.01,.01,'† Hatched: supervised methods evaluated against the same L1 label system.',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,1))
    for ext in ('png','pdf'):
        fig.savefig(target / f'aggregate_scores.{ext}', dpi=220)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--no-plots', action='store_true', help='Arithmetic/markdown only; explicitly marks figures not generated')
    args = p.parse_args(); root = args.root
    target = root / 'reports'; target.mkdir(exist_ok=True)
    metrics = root / 'metrics'
    try:
        raw, official, scaled, final = [read_table(metrics / name) for name in ('raw_scib_metrics.tsv','official_scores.tsv','supplementary_minmax_scores.tsv','final_ranking.tsv')]
        qa = audit(raw, official, scaled, final)
        manifest = json.loads((metrics / 'evaluation_manifest.json').read_text())
        n, total = manifest['n_evaluated'], manifest['n_primary']
        if n != 182968 or total != 315949:
            raise ValueError(f'Unexpected frozen evaluation coverage: {n}/{total}; update documented freeze before release')
        qa['coverage'] = {'trusted': n, 'primary': total, 'fraction': n/total}
        qa['figures_generated'] = not args.no_plots
        if not args.no_plots:
            make_figures(raw, official, target)
        names = lambda seq: '、'.join(PRETTY[m].replace(' †','') for m in seq)
        order = sorted(METHODS, key=lambda m: final[m]['rank'])
        lines = ['# 单细胞 integration benchmark 修订版结果', '',
                 f"综合表现最好：**{names(qa['winners']['Total'])}**；biology 最好：**{names(qa['winners']['Bio conservation'])}**；batch removal 最好：**{names(qa['winners']['Batch correction'])}**。",
                 f"无监督整合方法中，综合得分最高为 **{names(qa['unsupervised_integration_winner'])}**（PCA 作为未整合基线另列）。", '',
                 '| 排名 | 方法 | Biology | Batch | Overall |', '|---:|---|---:|---:|---:|']
        lines += [f"| {int(final[m]['rank'])} | {PRETTY[m]} | {official[m]['Bio conservation']:.4f} | {official[m]['Batch correction']:.4f} | {official[m]['Total']:.4f} |" for m in order]
        lines += ['', f'所有主指标使用共同冻结可信 L1 子集：{n:,}/{total:,} cells（{100*n/total:.2f}%）。Unknown 不作为真实细胞类型计分；各方法完整 embedding 仍保留全部 Primary cells。该排名描述可信标签子集，不能外推为全部 Unknown 细胞整合质量。', '',
                  '主分析遵循 scib-metrics 0.5.9 官方指标及默认邻居设置：graph connectivity 15、kBET 输入图 50、LISI 90；未执行可选 K=10/20/30，因此不能声称 K 排名稳定。未逐方法优化 K。', '',
                  '独立算术核验通过：5 个 biology 指标等权均值、5 个 batch 指标等权均值，Overall=0.6×Biology+0.4×Batch；八方法指标完整，无缺失分母变化。主排名使用原生指标范围。',
                  ('跨方法 min-max 辅表会改变排名，应以原生分数和逐项指标解释结果。' if qa['minmax_changes_rank'] else '跨方法 min-max 辅表与主分析排名一致。'), '',
                  '† scANVI 和 supervised scPoli 使用同一可信 L1 监督，L1 biology 分数属于标签知情的评价，不能宣称完全独立生物学验证。', '',
                  'CellHint 为 graph-only 输出，未纳入本次统一 embedding-based 八方法排名；此项排除范围需在论文方法及图注明确披露。', '',
                  '本次汇总检查未发现缺失指标静默平均或人为 min-max 主排名等计分错误；这里只核验输出表算术与评价覆盖，细胞一致性仍依赖评分 runner 的成员 hash 校验，不能替代上游数据身份、marker 注释与训练实现审计。保留监督标签循环评价、可信子集覆盖和 graph-only 方法排除的适用范围说明即可。合理的方法路线差异不要求继续无限追加 sensitivity analysis。', '',
                  '原始指标及所有数值见 metrics/raw_scib_metrics.tsv；本报告不生成或推断 p 值。']
        if args.no_plots:
            lines += ['', '本次使用 --no-plots，图形未生成。']
        (target / 'RESULTS_CN.md').write_text('\n'.join(lines)+'\n')
        (target / 'QA.json').write_text(json.dumps(qa, ensure_ascii=False, indent=2)+'\n')
        print(target / 'RESULTS_CN.md')
    except Exception as exc:
        (target / 'QA.json').write_text(json.dumps({'passed':False,'error':str(exc)},ensure_ascii=False,indent=2)+'\n')
        raise


if __name__ == '__main__':
    main()
