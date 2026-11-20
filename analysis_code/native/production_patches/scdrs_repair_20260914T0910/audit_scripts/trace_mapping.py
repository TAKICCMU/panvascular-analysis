import csv,gzip,re,collections,json
from pathlib import Path
p=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907')
with gzip.open(p/'01_corrected_counts/audit/full_gene_identifier_mapping.tsv.gz','rt') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
target={r['feature_name'] for r in rows if r['gene_mapping_status']=='symbol_fallback'}
print('targets',len(target),flush=True)
allhits={}
for fn in ['/home/anzhen_vhd/Reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf','/home/anzhen_vhd/software/scE2G/resources/genome_annotations/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz']:
    hits=collections.defaultdict(set);op=gzip.open if fn.endswith('gz') else open
    with op(fn,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            a=line.split('\t')
            if len(a)<9 or a[2]!='gene':continue
            attrs=dict(re.findall(r'(\w+) "([^"]*)"',a[8]))
            if attrs.get('gene_name') in target:hits[attrs['gene_name']].add(attrs.get('gene_id'))
    allhits[fn]={k:sorted(v) for k,v in hits.items()}
    print(fn,'matched',len(hits),'multiple',{k:sorted(v) for k,v in hits.items() if len(v)>1},flush=True)
union=set().union(*(set(x) for x in allhits.values()))
print('unmatched',len(target-union),sorted(target-union))
out=p/'06_reanalysis_takeover_20260908/production_launches/scdrs_repair_20260914T0910/inputs'
out.mkdir(parents=True,exist_ok=True)
(out/'local_gtf_exact_symbol_hits.json').write_text(json.dumps(allhits,indent=2))
