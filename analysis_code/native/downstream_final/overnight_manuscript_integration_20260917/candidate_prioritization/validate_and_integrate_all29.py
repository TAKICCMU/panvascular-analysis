import csv,gzip,json,hashlib
from pathlib import Path
D=Path(__file__).parent.resolve()
def read(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs):
 with open(D/p,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
new=read(D/'coloc_new12/02_results/Tier1_full_coloc_abf_results.tsv.gz');existing=read(D/'coloc_existing17/02_results/Tier1_full_coloc_abf_results.tsv.gz');allr=existing+new
for x in allr:
 assert float(x['n_overlap'])>=100
 assert abs(sum(float(x['PP.H'+str(i)]) for i in range(5))-1)<1e-8
assert len({(x['scenario'],x['gene'],x['context']) for x in allr})==len(allr)
write('all29_complete_cis_coloc_results.tsv',allr)
cov=read(D/'coloc_existing17/02_results/Tier1_gene_level_colocalisation_coverage.tsv')+read(D/'coloc_new12/02_results/Tier1_gene_level_colocalisation_coverage.tsv');write('all29_complete_cis_coloc_coverage.tsv',cov)
t=read(D/'Tier1_integrated_genetic_expression_regulatory_evidence.tsv')
for x in t:
 for s in ['primary','leave_one_CAD']:
  cc=next(c for c in cov if c['gene']==x['gene_symbol'] and c['scenario']==s)
  for k in ['n_contexts_tested','n_strong','maximum_PP.H4','best_context','best_classification']:x[s+'_recomputed_'+k]=cc[k]
 x['QTL_status']='NEW_AND_EXISTING_COMPLETE_CIS_ABF_ANALYZED_SEE_COVERAGE';x['regulatory_interpretation']='Original ABF method with fixed sdY=1 and prior sensitivity; no new SuSiE or matched LD; tissue maxima are descriptive, not multiple-testing adjusted.'
write('Tier1_integrated_evidence_FINAL.tsv',t)
man=read(D/'new_qtl_manifest.tsv');nq=[x for x in man if x.get('output') and Path(x['output']).exists()];write('new_QTL_SHA256.tsv',[dict(gene=x['gene'],context=x['context'],source_url=x['source_url'],source_build=x['source_build'],path=x['output'],n_qtl_rows=x['n_qtl_rows'],sha256=hashlib.sha256(Path(x['output']).read_bytes()).hexdigest()) for x in nq])
s={'status':('COMPLETE_VALIDATED' if all(x['status'] in ['DOWNLOADED','REUSED','NO_TARGET_GENE_QTL'] for x in man) else 'ANALYSIS_COMPLETE_WITH_QTL_RESOURCE_GAPS'),'new_QTL_manifest_rows':len(man),'new_QTL_statuses':{k:sum(x['status']==k for x in man) for k in {x['status'] for x in man}},'new_ABF_results':len(new),'existing_ABF_results':len(existing),'all_ABF_results':len(allr),'coverage_rows':len(cov),'genes':len(t),'no_new_SuSiE':True,'strong_primary_genes':[x['gene'] for x in cov if x['scenario']=='primary' and int(x['n_strong'])>0],'strong_leave_CAD_genes':[x['gene'] for x in cov if x['scenario']=='leave_one_CAD' and int(x['n_strong'])>0]};(D/'ALL29_COMPLETE.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
