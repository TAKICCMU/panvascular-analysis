import csv,gzip,json,hashlib
from pathlib import Path
D=Path(__file__).parent;root=D.parent/'Tier1_full_colocalisation_v1'
def read(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rs):
 with open(D/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(1024**2),b''):h.update(c)
 return h.hexdigest()
t=read(D/'Tier1_spatial_interface.tsv');m=read(root/'00_metadata/qtl_download_manifest.tsv');m+=read(D.parent/'GWAS_eQTL_colocalisation/00_metadata/qtl_download_manifest.tsv');c=read(root/'02_results/Tier1_full_coloc_abf_results.tsv.gz');rows=[];hr=[]
for x in m:
 if 'context' not in x:
  x['context']=x['tissue'];x['source']=x['study'];x['n_qtl_rows']=len(read(x['output'])) if x['output'] and Path(x['output']).exists() else 0
for r in t:
 g=r['gene_symbol'];qm=[x for x in m if x['gene']==g];ready=[x for x in qm if x['output'] and Path(x['output']).exists() and float(x.get('n_qtl_rows') or 0)>=100]; cr=[x for x in c if x['gene']==g]
 rows.append({'gene_symbol':g,'cell_state':r['cell_state'],'n_manifest_contexts':len(qm),'n_existing_ge100_qtl_contexts':len(ready),'n_existing_coloc_rows':len(cr),'QTL_status':'EXISTING_ARCHIVE_REUSE_REQUIRES_GWAS_PROVENANCE_REVIEW' if ready else 'NEW_COMPLETE_CIS_QTL_REQUIRED','max_historical_ABF_H4':max((float(x['PP.H4']) for x in cr),default='NA'),'regulatory_interpretation':'Historical coloc is gene/tissue evidence, not proof of current nominated cell state; not used to assign Tier'})
 for x in ready:hr.append({'gene':g,'context':x['context'],'path':x['output'],'sha256':sha(x['output']),'source_url':x['source_url'],'source_build':x['source_build'],'n_qtl_rows':x['n_qtl_rows']})
write('Tier1_eQTL_resource_audit.tsv',rows);write('available_QTL_SHA256.tsv',hr)
write('Tier1_historical_coloc_subset.tsv',[x for x in c if x['gene'] in {r['gene_symbol'] for r in t}])
print('QTL gap',[x['gene_symbol'] for x in rows if x['n_existing_ge100_qtl_contexts']==0]);print('available',sum(x['n_existing_ge100_qtl_contexts']>0 for x in rows));print('GWAS exist',[(str(p),p.exists()) for p in [Path('/Users/taki/Desktop/FUMA/full_bifactor_G_FUMA_all.tsv.gz'),Path('/Users/taki/Desktop/FUMA/LOT_CAD/leave_one_CAD_bifactor_G_FUMA_all.tsv.gz')]])
