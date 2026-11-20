import csv,gzip,json,hashlib
from pathlib import Path
D=Path(__file__).parent.resolve();R=D/'coloc_existing17';(R/'00_metadata').mkdir(parents=True,exist_ok=True)
def read(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
t={x['gene_symbol']:x for x in read(D/'Tier1_spatial_interface.tsv')};m=read(D.parent/'Tier1_full_colocalisation_v1/00_metadata/qtl_download_manifest.tsv');old=read(D.parent/'GWAS_eQTL_colocalisation/00_metadata/qtl_download_manifest.tsv')
for x in old:
 if x['gene']!='CKB':continue
 x['context']=x['tissue'];x['source']=x['study'];x['n_qtl_rows']=len(read(x['output']));m.append(x)
m=[x for x in m if x['gene'] in t]
for x in m:x['nominated_lineage']=t[x['gene']]['cell_state']
keys=sorted(set().union(*(x.keys() for x in m)))
with open(R/'00_metadata/qtl_download_manifest.tsv','w') as f:w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(m)
source=D.parent/'Tier1_full_colocalisation_v1/scripts/02_run_full_coloc.R';s=source.read_text().replace('ROOT <- "/Users/taki/Desktop/Codex任务1/Tier1_full_colocalisation_v1"','ROOT <- '+json.dumps(str(R),ensure_ascii=False));s=s.replace('library(data.table)','library(data.table)\n  setDTthreads(1)');s=s.split('# Plot gene-level maxima')[0]+'\nwriteLines(capture.output(sessionInfo()), file.path(ROOT, "00_metadata", "R_sessionInfo.txt"))\n';(D/'rerun_existing17_coloc.R').write_text(s)
print('manifest rows',len(m),'genes',len(set(x['gene'] for x in m)))
