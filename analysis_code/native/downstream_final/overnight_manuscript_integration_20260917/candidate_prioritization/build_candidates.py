import csv,gzip,json,hashlib
from pathlib import Path
from collections import Counter
D=Path(__file__).parent
old=D/'OLD_GWAS_pseudobulk_integrated_candidates.tsv.gz'
def read(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,cols=None):
 if not cols:cols=list(rows[0])
 with open(D/p,'w') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
g={r['SYMBOL']:r for r in read(D/'original_GWAS_master.tsv')}
T1='Tier1_GWAS_Bonferroni_GWAS-aligned_DE_LOT-CAD';T2='Tier2_GWAS_FDR_plus_multidisease_DE';T0='GWAS_prioritised_without_shared_DE_support'
def classify(r):
 f=lambda k:float(r[k]); b=lambda k:str(r[k]).lower()=='true'
 tier=T0
 if b('primary_FDR') and f('target_set_meta_FDR')<.05 and f('n_disease_contexts')>=2 and f('direction_consistency')>=.75:tier=T2
 if b('primary_Bonferroni') and f('target_set_meta_FDR')<.05 and f('n_GWAS_aligned_contexts')>=2 and f('direction_consistency')>=.75 and b('retained_leave_one_CAD_nominal'):tier=T1
 return tier
rows=read(D/'all_priority_meta.tsv.gz')
for r in rows:
 r.update(g[r['gene_symbol'].strip().upper()]);r['Tier']=classify(r);r['strict_three_context_shared']=float(r['n_disease_contexts'])>=3 and float(r['direction_consistency'])>=.75
rows.sort(key=lambda r:([T1,T2,T0].index(r['Tier']),float(r['target_set_meta_FDR']),float(r['meta_P'])))
write('GWAS_pseudobulk_integrated_candidates_NEW.tsv',rows)
sel=[r for r in rows if r['Tier']!=T0];t1=[r for r in rows if r['Tier']==T1]
write('prioritised_target_candidates_NEW.tsv',sel,list(rows[0]));write('Tier1_spatial_interface.tsv',t1,list(rows[0]))
olds=read(old);assert all(classify(r)==r['Tier'] for r in olds),'Old rule mismatch'
key=lambda r:({'Level1':'L1','Level2':'L2'}.get(r['annotation_level'],r['annotation_level']),r['cell_state'],r['gene_symbol'])
a={key(r):r for r in olds};b={key(r):r for r in rows};cmp=[]
for k in sorted(set(a)|set(b)):
 o=a.get(k,{});n=b.get(k,{})
 cmp.append(dict(zip(['annotation_level','cell_state','gene_symbol'],k),old_Tier=o.get('Tier','NO_MATCHED_GENE_STATE_META_ROW'),new_Tier=n.get('Tier','NO_MATCHED_GENE_STATE_META_ROW'),old_meta_logFC=o.get('meta_logFC',''),new_meta_logFC=n.get('meta_logFC',''),old_target_FDR=o.get('target_set_meta_FDR',''),new_target_FDR=n.get('target_set_meta_FDR','')))
write('all_old_vs_new.tsv',cmp);write('old_or_new_Tier1_changes.tsv',[r for r in cmp if r['old_Tier']==T1 or r['new_Tier']==T1]);write('CKB_position.tsv',[r for r in rows if r['gene_symbol']=='CKB'])
s={'new_pairs':dict(Counter(r['Tier'] for r in rows)),'old_pairs':dict(Counter(r['Tier'] for r in olds)),'new_Tier1_genes':sorted(set(r['gene_symbol'] for r in t1)),'old_Tier1_genes':sorted(set(r['gene_symbol'] for r in olds if r['Tier']==T1)),'old_rules_reproduced_all_rows':True,'new_meta_rows':len(rows),'new_gwas_genes':len(set(r['gene_symbol'] for r in rows))}
(D/'SUMMARY.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
