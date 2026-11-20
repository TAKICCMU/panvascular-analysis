#!/usr/bin/env python3
"""Resumable scib-metrics 0.5.9 official benchmark; never loads expression X.

Run one method per process sequentially, then --aggregate. Use --method all only
when all eight embeddings exist. Main scores use native metric scaling; optional
across-method minmax is supplementary. Unknown cells remain in saved embeddings,
but all primary metrics use the same frozen trusted-cell evaluation population.
"""
import argparse
import hashlib
import importlib.metadata
import json
import logging
import math
import os
from pathlib import Path
import sys
import time

METHODS = ('harmony_cpu',)
BIO = ('isolated_labels', 'nmi_ari_cluster_labels_kmeans_nmi',
       'nmi_ari_cluster_labels_kmeans_ari', 'silhouette_label', 'clisi_knn')
BATCH = ('bras', 'ilisi_knn', 'kbet_per_label', 'graph_connectivity', 'pcr_comparison')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def save_json(path, obj):
    tmp = path.with_name(path.name + f'.{os.getpid()}.{time.time_ns()}.tmp')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n')
    tmp.replace(path)


def validate_scores(rows):
    """No per-method missing-metric denominator or silent NaN averaging."""
    for method in METHODS:
        if method not in rows:
            raise ValueError(f'Missing method: {method}')
        if set(rows[method]) != set(BIO + BATCH):
            raise ValueError(f'Metric set mismatch: {method}')
        for metric, value in rows[method].items():
            if value is None or not math.isfinite(value):
                raise ValueError(f'Nonfinite metric blocks aggregate: {method}/{metric}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--method', choices=(*METHODS, 'all'), default='all')
    p.add_argument('--aggregate', action='store_true')
    p.add_argument('--prepare-only', action='store_true', help='CPU neighbor preparation only; never mark metrics complete')
    p.add_argument('--batch-key', default='sample_uid')
    p.add_argument('--seed', type=int, default=20260909)
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--allow-cpu', action='store_true', help='Explicitly permit CPU metrics if JAX GPU is unavailable')
    args = p.parse_args()
    if args.prepare_only and args.aggregate:
        p.error('--prepare-only and --aggregate are mutually exclusive')
    if args.self_test:
        good = {m: {k: 0.5 for k in BIO + BATCH} for m in METHODS}
        validate_scores(good)
        good['pca']['bras'] = None
        try:
            validate_scores(good)
        except ValueError:
            print('PASS: complete scores accepted; missing metric blocks aggregate')
            return
        raise AssertionError('Missing metric accepted')
    # Set before any scientific imports. No full-VRAM preallocation; concurrent
    # GPU training should still be avoided because ASW creates temporary arrays.
    if args.prepare_only or args.aggregate:
        os.environ['JAX_PLATFORMS'] = 'cpu'
    elif not args.allow_cpu:
        os.environ['JAX_PLATFORMS'] = 'cuda'
    os.environ.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')
    os.environ.setdefault('XLA_PYTHON_CLIENT_MEM_FRACTION', '0.30')
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        os.environ.setdefault(key, str(args.jobs))
    # Initialize JAX before Benchmarker/Scanpy can load incompatible conda cuBLASLt.
    import jax
    jax.config.update('jax_compilation_cache_dir', str(args.root / 'cache' / 'jax'))
    initial_devices = jax.devices()
    import numpy as np
    import pandas as pd
    import anndata as ad
    from scib_metrics.benchmark import Benchmarker
    from scib_metrics.nearest_neighbors import NeighborsResults, pynndescent
    version = importlib.metadata.version('scib-metrics')
    if version != '0.5.9':
        raise RuntimeError(f'Expected scib-metrics 0.5.9, got {version}')
    root = args.root.resolve()
    out = root / 'metrics'
    out.mkdir(exist_ok=True)
    (root / 'logs').mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler(), logging.FileHandler(root / 'logs' / f'evaluate_{args.method}_{int(time.time())}.log')])
    versions = {x: importlib.metadata.version(x) for x in ('scib-metrics', 'anndata', 'numpy', 'pandas', 'scipy', 'jax', 'jaxlib', 'scanpy', 'pynndescent', 'umap-learn', 'scikit-learn')}
    logging.info('Versions=%s devices=%s', versions, jax.devices())
    if not args.aggregate and not args.prepare_only and not args.allow_cpu and not any(d.platform == 'gpu' for d in jax.devices()):
        raise RuntimeError('JAX GPU unavailable; fix runtime or explicitly pass --allow-cpu; no silent CPU workload')
    save_json(root / 'logs' / f'backend_{args.method}_{time.time_ns()}.json', {'devices': [str(d) for d in jax.devices()], 'prepare_only': args.prepare_only, 'aggregate_only': args.aggregate, 'versions': versions, 'pid': os.getpid()})
    np.random.seed(args.seed)
    source = ad.read_h5ad(root / 'inputs' / 'primary.h5ad', backed='r')
    obs = source.obs.copy()
    source.file.close()
    if not obs.index.is_unique:
        raise ValueError('Duplicate cell identities')
    required = ['reliable_L1', 'trusted', args.batch_key]
    if obs[required].isna().any().any():
        raise ValueError('Missing labels, trust flags, or batches')
    if not set(obs.trusted.unique()).issubset({True, False, 0, 1}):
        raise ValueError('trusted must be bool/0/1, not strings')
    mask = obs.trusted.to_numpy(dtype=bool)
    if (obs.loc[mask, 'reliable_L1'].astype(str) == 'Unknown').any():
        raise ValueError('Unknown marked trusted')
    selected = obs.loc[mask, required].copy()
    selected['reliable_L1'] = selected.reliable_L1.astype(str)
    selected[args.batch_key] = selected[args.batch_key].astype(str)
    if selected.reliable_L1.nunique() < 2 or selected[args.batch_key].nunique() < 2:
        raise ValueError('Need at least two labels and batches')
    membership = digest({'all_cells': list(obs.index.astype(str)), 'trusted': mask.tolist(), 'labels': selected.reliable_L1.tolist(), 'batches': selected[args.batch_key].tolist()})
    prefile = root / 'embeddings' / 'pca.npy'
    pre = np.load(prefile, mmap_mode='r')
    if len(pre) != len(obs) or not np.isfinite(pre).all():
        raise ValueError('PCA invalid or wrong cell count')
    config = {'membership_sha256': membership, 'pca_sha256': sha(prefile), 'versions': versions, 'seed': args.seed, 'neighbors': [15, 50, 90], 'neighbor_random_state': 0, 'batch_key': args.batch_key, 'n_primary': len(obs), 'n_evaluated': int(mask.sum()), 'coverage': float(mask.mean()), 'unknown_policy': 'All metrics on identical trusted subset; Unknown never a type', 'k_sensitivity': 'Not executed; official metric-specific defaults retained', 'min_max_scale_main': False, 'weights': {'biology': 0.6, 'batch': 0.4}}
    historical = json.loads((root / 'audit' / 'original_evaluation_manifest.json').read_text())
    if config != historical:
        raise RuntimeError('Frozen evaluation configuration/version mismatch; stopping')
    runhash = digest(config)
    manifest = out / 'evaluation_manifest.json'
    if manifest.exists() and json.loads(manifest.read_text()) != config:
        raise RuntimeError('Existing evaluation configuration differs; use new directory')
    save_json(manifest, config)
    celltmp = out / f'evaluated_cells.{os.getpid()}.tsv.gz'
    selected.to_csv(celltmp, sep='\t')
    celltmp.replace(out / 'evaluated_cells.tsv.gz')
    support = selected.groupby('reliable_L1', observed=True).agg(n_cells=(args.batch_key, 'size'), n_batches=(args.batch_key, 'nunique'))
    supporttmp = out / f'label_batch_support.{os.getpid()}.tsv'
    support.to_csv(supporttmp, sep='\t')
    supporttmp.replace(out / 'label_batch_support.tsv')
    locks = []  # Hold method locks until process exit, including prepare-only jobs.
    if not args.aggregate:
        import fcntl
        for method in METHODS if args.method == 'all' else [args.method]:
            lock = open(out / f'{method}.lock', 'a')
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                lock.close()
                raise RuntimeError(f'Method already running: {method}')
            locks.append(lock)
            embfile = root / 'embeddings' / f'{method}.npy'
            provenance = {'config_sha256': runhash, 'embedding_sha256': sha(embfile)}
            cache = out / f'{method}.json'
            if cache.exists():
                if json.loads(cache.read_text())['provenance'] != provenance:
                    raise RuntimeError(f'Stale metric cache: {method}')
                logging.info('Reusing completed metrics: %s', method)
                continue
            emb = np.load(embfile, mmap_mode='r')
            if emb.ndim != 2 or len(emb) != len(obs) or not np.isfinite(emb).all():
                raise ValueError(f'Invalid embedding: {method}')
            tiny = ad.AnnData(obs=selected.copy())
            tiny.obsm['X_pre'] = np.asarray(pre[mask], dtype=np.float32)
            tiny.obsm[method] = np.asarray(emb[mask], dtype=np.float32)
            bm = Benchmarker(tiny, batch_key=args.batch_key, label_key='reliable_L1', embedding_obsm_keys=[method], pre_integrated_embedding_obsm_key='X_pre', n_jobs=args.jobs)
            graphfile = out / f'{method}.neighbors.npz'
            graphmeta = out / f'{method}.neighbors.json'
            def neighbor_computer(x, k):
                if graphfile.exists() and graphmeta.exists() and json.loads(graphmeta.read_text()) == provenance:
                    with np.load(graphfile) as z:
                        return NeighborsResults(indices=z['indices'], distances=z['distances'])
                result = pynndescent(x, n_neighbors=k, random_state=0, n_jobs=args.jobs)
                np.savez(graphfile, indices=result.indices, distances=result.distances)
                save_json(graphmeta, provenance)
                return result
            start = time.time()
            bm.prepare(neighbor_computer=neighbor_computer)
            if args.prepare_only:
                logging.info('Prepared official neighbors only: %s; no metric completion cache written', method)
                del bm, tiny
                import gc
                gc.collect()
                continue
            try:
                bm.benchmark()
            except Exception:
                bm._results.to_csv(out / f"{method}.partial_official_metrics.tsv", sep="\t")
                logging.exception("Official metric failure; partial values preserved, no aggregate")
                raise
            table = bm.get_results(min_max_scale=False, clean_names=False)
            row = {key: float(table.loc[method, key]) for key in BIO + BATCH}
            if any(not math.isfinite(v) for v in row.values()):
                save_json(out / f'{method}.failure.json', {'provenance': provenance, 'scores': {k: v if math.isfinite(v) else None for k, v in row.items()}, 'reason': 'Nonfinite metric; aggregate prohibited, investigate structural NA explicitly'})
                raise RuntimeError(f'Nonfinite official metric for {method}; no ranking produced')
            save_json(cache, {'provenance': provenance, 'metrics': row, 'seconds': time.time() - start})
            del bm, tiny
            import gc
            gc.collect()
    if not args.prepare_only and (args.aggregate or args.method == 'all'):
        rows = {}
        for method in METHODS:
            record = json.loads((out / f'{method}.json').read_text())
            expected = {'config_sha256': runhash, 'embedding_sha256': sha(root / 'embeddings' / f'{method}.npy')}
            if record['provenance'] != expected:
                raise RuntimeError(f'Stale cache: {method}')
            rows[method] = record['metrics']
        validate_scores(rows)
        raw = pd.DataFrame.from_dict(rows, orient='index')
        raw.to_csv(out / 'raw_scib_metrics.tsv', sep='\t')
        # Let the installed official aggregator calculate both versions. This
        # adapter only supplies validated cached outputs, never changes metrics.
        bm = Benchmarker(ad.AnnData(), batch_key=args.batch_key, label_key='reliable_L1', embedding_obsm_keys=list(METHODS))
        bm._results = raw.T
        bm._results['Metric Type'] = ['Bio conservation' if k in BIO else 'Batch correction' for k in bm._results.index]
        for scale, filename in [(False, 'official_scores.tsv'), (True, 'supplementary_minmax_scores.tsv')]:
            result = bm.get_results(min_max_scale=scale, clean_names=False).drop(index='Metric Type').astype(float)
            result.sort_values('Total', ascending=False).to_csv(out / filename, sep='\t')
            if not scale:
                rank = result[['Bio conservation', 'Batch correction', 'Total']].copy()
                rank['rank'] = rank.Total.rank(method='min', ascending=False).astype(int)
                rank.sort_values('rank').to_csv(out / 'final_ranking.tsv', sep='\t')
        logging.info('Complete: all eight methods, identical metric denominator, native scaling main')


if __name__ == '__main__':
    main()
