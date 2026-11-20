"""Independent arithmetic and completeness audit; no integration/model fitting."""
from pathlib import Path
import csv, math, hashlib, json, statistics, sys
B=Path(sys.argv[1]); P=B/'primary'; U=P/'unified'
def rows(p):
    with p.open() as f: return list(csv.DictReader(f,delimiter='\t'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def finite(v):
    try: return math.isfinite(float(v))
    except (ValueError,TypeError): return False
raw=rows(P/'combined_raw_metrics_by_K.tsv'); final=rows(U/'unified_scib_configuration_ranking.tsv')
status=json.loads((U/'UNIFIED_SCIB_RANKING_STATUS.json').read_text())
assert sha(P/'combined_raw_metrics_by_K.tsv')==status['input_sha256']
methods={'pca','harmony','scanorama','liger','scvi','scanvi','scpoli','scalex','cellhint'}
assert len(raw)==27 and len(final)==9 and {r['method'] for r in raw}==methods
assert len({(r['method'],r['config_id'],r['K']) for r in raw})==27
assert {int(r['n_cells']) for r in raw}=={315949}
batch=['PCR_batch','iLISI_batch','ASW_batch','graph_connectivity']
bio=['ASW_label','NMI_best','ARI_best','cLISI_label']; cb=['iLISI_batch','graph_connectivity']; cy=['NMI_best','ARI_best','cLISI_label']
calculated=[]
for k in [10,20,30]:
    part=[r for r in raw if int(r['K'])==k]; assert len(part)==9
    bounds={m:[float(r[m]) for r in part if finite(r[m])] for m in batch+bio}
    for r in part:
        required=cb+cy if r['method']=='cellhint' else batch+bio
        assert all(finite(r[m]) for m in required)
        if r['method']=='cellhint': assert all(not finite(r[m]) for m in ['PCR_batch','ASW_batch','ASW_label'])
        scaled={}
        for m in batch+bio:
            if not finite(r[m]): continue
            lo,hi=min(bounds[m]),max(bounds[m])
            scaled[m]=.5 if math.isclose(lo,hi,rel_tol=1e-5,abs_tol=1e-8) else (float(r[m])-lo)/(hi-lo)
        for metric in ['NMI','ARI']:
            vals=[float(r[f'{metric}_resolution_{res}']) for res in ['0.2','0.5','1','1.5']]
            assert math.isclose(float(r[metric+'_best']),max(vals),abs_tol=1e-12)
        mean=lambda ms:statistics.mean(scaled[m] for m in ms if m in scaled)
        calculated.append({'method':r['method'],'K':k,'unified':.4*mean(batch)+.6*mean(bio),'common_graph':.4*mean(cb)+.6*mean(cy)})
for r in final:
    assert r['eligible_for_formal_ranking']=='True' and r['K_values_present']=='10;20;30'
    for metric in ['unified','common_graph']:
        value=statistics.mean(x[metric] for x in calculated if x['method']==r['method'])
        assert math.isclose(value,float(r[f'mean_{metric}_overall_score']),abs_tol=1e-12)
per_k=[]
for k in [10,20,30]:
    h=next(x for x in calculated if x['method']=='harmony' and x['K']==k)
    s=next(x for x in calculated if x['method']=='scvi' and x['K']==k)
    per_k.append({'K':k,'Harmony_unified':h['unified'],'scVI_unified':s['unified'],'Harmony_minus_scVI':h['unified']-s['unified'],'Harmony_common_graph':h['common_graph'],'scVI_common_graph':s['common_graph']})
qa={'status':'PASS_independent_arithmetic_and_completeness','n_methods':9,'n_K_rows':27,'n_graph_cells':315949,
    'Leiden_4_resolution_best_values_verified':True,'same_within_K_scaling_verified':True,'equal_K_weight_verified':True,
    'graph_only_structural_missingness_verified':True,'source_raw_sha256':sha(P/'combined_raw_metrics_by_K.tsv'),
    'ranking_sha256':sha(U/'unified_scib_configuration_ranking.tsv'),'code_sha256':sha(Path(__file__)),
    'Harmony_vs_scVI_by_K':per_k,'scope':'Arithmetic audit of frozen raw metrics; no new integration, donor bootstrap or significance claim.',
    'files':{p.name:sha(p) for p in sorted(U.glob('*')) if p.is_file()}}
O=B/'audit';O.mkdir(exist_ok=True);(O/'ROOT_PRIMARY_NINE_METHOD_RANKING_QA.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa,indent=2))
