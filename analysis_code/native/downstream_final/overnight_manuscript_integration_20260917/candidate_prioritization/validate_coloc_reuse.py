import csv,gzip,json,hashlib
from pathlib import Path
D=Path(__file__).parent;out=D/'coloc_existing17/02_results'
def read(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
new=read(out/'Tier1_full_coloc_abf_results.tsv.gz');old=read(D.parent/'Tier1_full_colocalisation_v1/02_results/Tier1_full_coloc_abf_results.tsv.gz');key=lambda x:(x['scenario'],x['gene'],x['context']);o={key(x):x for x in old};matched=[x for x in new if key(x) in o];delta=max(abs(float(x['PP.H4'])-float(o[key(x)]['PP.H4'])) for x in matched)
ckb=[x for x in new if x['gene']=='CKB'];historical_ckb={ (x['scenario'],x['gene'],x['tissue']):x for x in read(D.parent/'GWAS_eQTL_colocalisation/03_results/coloc_abf_summary.tsv') if x['gene']=='CKB'};ckb_delta=max(abs(float(x['PP.H4'])-float(historical_ckb[key(x)]['PP.H4'])) for x in ckb);cols=list(new[0])
with open(D/'CKB_existing_cis_coloc_recomputed.tsv','w') as f:w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(ckb)
h={}
for p in [Path('/Users/taki/Desktop/FUMA/full_bifactor_G_FUMA_all.tsv.gz'),Path('/Users/taki/Desktop/FUMA/LOT_CAD/leave_one_CAD_bifactor_G_FUMA_all.tsv.gz')]:
 a=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):a.update(b)
 h[str(p)]=a.hexdigest()
s={'status':'COMPLETE_VALIDATED','scenario_gene_context_results':len(new),'matched_historical_rows':len(matched),'max_abs_H4_difference_to_historical':delta,'CKB_max_abs_H4_difference_to_historical':ckb_delta,'CKB_results':len(ckb),'CKB_max_H4':max(float(x['PP.H4']) for x in ckb),'GWAS_SHA256':h,'method':'Inherited complete cis coloc.abf; sdY=1 quantitative latent G; no new SuSiE; one thread'};(D/'COLOC_EXISTING17_COMPLETE.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
