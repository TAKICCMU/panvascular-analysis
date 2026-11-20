import importlib.util,csv,json,urllib.request,time,gzip,subprocess
from pathlib import Path
D=Path(__file__).parent;Q=D/'new_qtl';Q.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('legacy',D.parent/'Tier1_full_colocalisation_v1/scripts/01_fetch_full_qtl.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.OUT=Q;m.REQUEST_TIMEOUT_SECONDS=45;m.MAX_ATTEMPTS=1
r=list(csv.DictReader(open(D/'Tier1_eQTL_resource_audit.tsv'),delimiter='\t'));ts={x['gene_symbol']:x for x in csv.DictReader(open(D/'Tier1_spatial_interface.tsv'),delimiter='\t')};rows=[];coords=[]
for x in r:
 if x['n_existing_ge100_qtl_contexts']!='0':continue
 t=ts[x['gene_symbol']];gid=t['gene_id'];cache=Q/(gid+'.ensembl.json')
 try:
  if not cache.exists():
   req=urllib.request.Request('https://rest.ensembl.org/lookup/id/'+gid+'?content-type=application/json',headers={'Accept':'application/json'});cache.write_bytes(urllib.request.urlopen(req,timeout=20).read());time.sleep(.5)
  e=json.loads(cache.read_text());assert e['assembly_name']=='GRCh38'
  gene=dict(gene=x['gene_symbol'],gene_id=gid,nominated_lineage=x['cell_state'],chromosome=str(e['seq_region_name']),start_grch38=str(e['start']),end_grch38=str(e['end']),start_grch37='',end_grch37='');coords.append(gene)
 except Exception as ex:
  rows.append(dict(gene=x['gene_symbol'],status='COORDINATE_FAILED',error=str(ex)));continue
 # Original seven-context GTEx panel, serial, limited to complete cis regions.
 for context,(qid,n) in m.GTEX.items():
  out=m.fetch_one(gene,context,'GTEx_eQTL_Catalogue',m.SUMSTATS_BASE+'/QTS000015/'+qid+'/'+qid+'.all.tsv.gz',n,'GRCh38');rows.append(out);print(x['gene_symbol'],context,out['status'],out.get('n_qtl_rows'),flush=True)
  (D/'NEW_QTL_PROGRESS.json').write_text(json.dumps(rows,indent=2));time.sleep(2.5)
  # Halt external requests after three consecutive service failures.
  if len(rows)>=3 and all(z['status']=='TABIX_FAILED' for z in rows[-3:]):
   (D/'NEW_QTL_BLOCKED.json').write_text(json.dumps({'reason':'three_consecutive_TABIX_FAILED','rows':rows},indent=2));raise SystemExit(2)
if rows:
 keys=sorted(set().union(*(x.keys() for x in rows)))
 with open(D/'new_qtl_manifest.tsv','w') as f:w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(rows)
(D/'NEW_QTL_COMPLETE.json').write_text(json.dumps({'rows':len(rows),'coordinates':coords},indent=2))
