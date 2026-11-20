import importlib.util,csv,json,urllib.request,time,sys,hashlib,shutil
from pathlib import Path
V=Path(__file__).parent.resolve();D=V.parent;Q=V/'new_qtl';Q.mkdir(exist_ok=True)
EXPECTED='a872fc169066e37de51bc5292c10d030094d87ec676bdb6bd6cf2319afe3711b'
assert hashlib.sha256((D/'Tier1_spatial_interface.tsv').read_bytes()).hexdigest()==EXPECTED
spec=importlib.util.spec_from_file_location('legacy',D/'original_fetch_full_qtl.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.OUT=Q;m.REQUEST_TIMEOUT_SECONDS=360;m.MAX_ATTEMPTS=1
old=json.loads((D/'NEW_QTL_PROGRESS.json').read_text());known={(x['gene'],x['context']):x for x in old if x['status'] in ['DOWNLOADED','REUSED','NO_TARGET_GENE_QTL']}
for k,x in known.items():
 p=Path(x['output']);new=Q/p.name
 if not new.exists():shutil.copyfile(p,new)
 done=Path(str(p)+'.done');nd=Path(str(new)+'.done')
 if done.exists() and not nd.exists():shutil.copyfile(done,nd)
r=list(csv.DictReader(open(D/'Tier1_eQTL_resource_audit.tsv'),delimiter='\t'));ts={x['gene_symbol']:x for x in csv.DictReader(open(D/'Tier1_spatial_interface.tsv'),delimiter='\t')}
def gene_record(g):
 t=ts[g];gid=t['gene_id'];cache=Q/(gid+'.ensembl.json');orig=D/'new_qtl'/cache.name
 if not cache.exists() and orig.exists():shutil.copyfile(orig,cache)
 if not cache.exists():
  req=urllib.request.Request('https://rest.ensembl.org/lookup/id/'+gid+'?content-type=application/json',headers={'Accept':'application/json'});cache.write_bytes(urllib.request.urlopen(req,timeout=30).read());time.sleep(1)
 e=json.loads(cache.read_text());assert e['assembly_name']=='GRCh38'
 return dict(gene=g,gene_id=gid,nominated_lineage=t['cell_state'],chromosome=str(e['seq_region_name']),start_grch38=str(e['start']),end_grch38=str(e['end']),start_grch37='',end_grch37='')
def fetch(g,ctx):
 qid,n=m.GTEX[ctx];return m.fetch_one(gene_record(g),ctx,'GTEx_eQTL_Catalogue',m.SUMSTATS_BASE+'/QTS000015/'+qid+'/'+qid+'.all.tsv.gz',n,'GRCh38')
if '--pilot' in sys.argv:
 started=time.time();x=fetch('SKI','Artery_Aorta');x['elapsed_seconds']=time.time()-started;x['timeout_seconds']=360;(V/'PILOT_RESULT.json').write_text(json.dumps(x,indent=2));print(json.dumps(x,indent=2),flush=True);raise SystemExit(0 if x['status'] in ['DOWNLOADED','REUSED','NO_TARGET_GENE_QTL'] else 2)
rows=[];coords=[]
for x in r:
 if x['n_existing_ge100_qtl_contexts']!='0':continue
 g=x['gene_symbol'];coords.append(gene_record(g))
 for ctx in m.GTEX:
  if (g,ctx) in known:
   out=dict(known[(g,ctx)]);out['output']=str(Q/Path(out['output']).name);out['status']='REUSED_SUCCESSFUL_V1' if int(out['n_qtl_rows'])>0 else 'NO_TARGET_GENE_QTL';out['source_version']='v1_success_preserved'
  else:out=fetch(g,ctx);time.sleep(2.5)
  rows.append(out);(V/'NEW_QTL_PROGRESS.json').write_text(json.dumps(rows,indent=2));print(g,ctx,out['status'],out.get('n_qtl_rows'),flush=True)
  if len(rows)>=3 and all(z['status']=='TABIX_FAILED' for z in rows[-3:]):
   (V/'NEW_QTL_BLOCKED.json').write_text(json.dumps({'reason':'three_consecutive_real_service_failures_at_360s','rows':rows},indent=2));raise SystemExit(2)
keys=sorted(set().union(*(x.keys() for x in rows)))
with open(V/'new_qtl_manifest.tsv','w') as f:w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(rows)
assert len(rows)==84
(V/'NEW_QTL_COMPLETE.json').write_text(json.dumps({'rows':len(rows),'coordinates':coords,'timeout_seconds':360,'reused_v1_success':len(known)},indent=2))
