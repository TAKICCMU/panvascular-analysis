# Reliable L1 integration benchmark v2 — 2026-09-09

This new version preserves all historical files. It retains the frozen 315,949-cell Primary set (247 sample_uid batches, seed 0, at most 1,500 cells per sample), its 2,000 HVGs, and validated unsupervised embeddings. No new pilot search is performed. Historical configuration selection used a nested pilot, so this comparison is not an independent held-out generalization test.

## Labels

The existing 14 canonical broad-lineage marker programs and integration-independent scores are retained. A trusted assignment now requires >=3 detected markers, >=40% program marker detection, and top/second score >=2. All other cells become Unknown and use the literal __unlabeled__ for supervised models. The retained score is the IDF-weighted mean of full-library-normalized log1p marker expression, multiplied by sqrt(detected program markers / program size). IDF=log((N+1)/(marker-positive cells+1))+1 is estimated once across the full atlas and normalized within each program; identical weights are applied to each study_file. Score ratios are heuristics, not calibrated probabilities. Exact marker/weight tables and original scoring code are in labels/source_reference. Thresholds are frozen before this version's scoring, and are not optimized to maximize annotation rate or benchmark scores. No L2 is generated. Rare programs are displayed for review but are not automatically added to the historical 14-class vocabulary.

Full atlas: 1,175,115 cells, 512,399 Unknown (43.60%). Primary: 315,949 cells, 182,968 trusted and 132,981 Unknown (42.09%). See labels/*summary.tsv, per-study counts, complete supervision labels, marker means/detection, and marker_support_panel*.pdf. These are conservative marker-supported provisional identities, not author-curated ground truth. Fib/SMC/FMC overlap remains a limitation; uncertain cells are not forced into a class.

## Integration

PCA, Harmony, scVI, Scanorama, LIGER and SCALEX are reused via read-only source symlinks with source cell-order and embedding SHA verification. They are not retrained. scANVI loads the frozen scVI model, checks latent reproducibility on 1,000 deterministic cells, and trains with trusted labels, 50 maximum epochs and early stopping. scPoli uses the same supervision labels, latent10/condition embedding10, NB, eta5, minibatch512 with 8 data-loader workers (fixed before scoring for GPU throughput; historical default was128), and the official tutorial's 50 total / 40 pretraining epochs (80% unsupervised pretraining); this replaces historical 60/40 without a parameter search. Unknown prototypes do not participate in supervised training. All embeddings retain every Primary cell.

Historical CellHint is graph-only and label-informed. Its old graph uses obsolete supervision, so it is not reused or ranked against the eight embeddings with a different metric denominator. Its historical output is preserved. This run does not claim to be a completed new nine-method comparison.

## Official metrics

Use installed scib-metrics 0.5.9 Benchmarker defaults. Its metrics use different neighborhood sizes: graph connectivity 15, kBET 50 and LISI 90; the latter corresponds to perplexity30. Every method receives identical per-metric parameters. Do not describe this as a single K=30 benchmark. Default KMeans NMI/ARI is used instead of selecting each method's best Leiden resolution. Optional K10/20/30 sensitivity is not part of this primary run: preserving official metric defaults takes priority, and no rank-stability claim will be made without such a run.

All primary metrics use the same 182,968 trusted cells, including label-dependent batch metrics. Unknown is never treated as a cell type. All methods were integrated on the identical full Primary population; evaluation covers 57.91%, so scores cannot establish biological preservation within Unknown. Biology: isolated labels, KMeans NMI/ARI, label ASW, cLISI. Batch: BRAS, iLISI, kBET per label, graph connectivity, PCR comparison against the matching unintegrated PCA.

Save all raw metrics. Main aggregation uses official get_results(min_max_scale=False), category means and 0.6 biology + 0.4 batch; across-method minmax output is supplementary. Missing metrics must be explicitly investigated; no method-specific denominator or silently omitted kBET is allowed. Label-informed L1 scores are not independent validation of scANVI/scPoli.

## Compute and state

Preparation/scANVI use omicverse; scPoli uses its existing compatible environment; metrics use existing scib_GPU 0.5.9 because omicverse lacks scib-metrics. No environment is overwritten. JAX GPU requires LD_LIBRARY_PATH unset and initialization before importing Benchmarker/Scanpy, which otherwise loads an incompatible conda cuBLASLt; retain version checks and record actual backend. GPU-heavy tasks share GPU_HEAVY_STAGE.lock and use bounded memory, CPU thread caps and intermediate embedding/neighbor caches. Raw million-cell counts are not reread; label review uses the existing small marker cache.

Completed: all eight methods and all ten official metrics passed aggregation and independent arithmetic QA. Label identity, saved-model supervision, embedding integrity, and successful GPU runtime checks also passed. See reports/RESULTS_CN.md, reports/QA.json and audit/*QA.json. Overall and biology winner: supervised scPoli; batch and unsupervised overall winner: Harmony. Optional K sensitivity was not run. Main scores use native metric values; supplementary minmax swaps PCA/SCALEX only, leaving the top four unchanged. The previously paused stromal/mural reintegration and its renderer remain paused.

Official references: https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.benchmark.Benchmarker.html ; https://docs.scvi-tools.org/en/latest/api/reference/scvi.model.SCANVI.html ; https://docs.scarches.org/en/latest/scpoli_surgery_pipeline.html
