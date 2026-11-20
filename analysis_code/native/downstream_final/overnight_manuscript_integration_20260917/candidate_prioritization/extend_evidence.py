exec(open(__file__.replace('extend_evidence.py','audit_qtl.py')).read().split('t=read')[0])
fuma=Path('/Users/taki/Desktop/FUMA/FUMA结果');t=read(D/'Tier1_spatial_interface.tsv');q={r['gene_symbol']:r for r in read(D/'Tier1_eQTL_resource_audit.tsv')};fs={s:read(fuma/f/'eqtl.txt') for s,f in [('primary','FUMA_LDr20.1Q_bifactor'),('leave_CAD','FUMA_LOTCAD_LDr20.1Q_bifactor')]}
for r in t:
 r.update(q[r['gene_symbol']])
 for s,a in fs.items():
  hits=[x for x in a if x['symbol']==r['gene_symbol'] and x['eqtlMapFilt']=='1']
  r[s+'_FUMA_eqtl_mapping_rows']=len(hits);r[s+'_FUMA_eqtl_tissues']=';'.join(sorted(set(x['tissue'] for x in hits)))
write('Tier1_integrated_genetic_expression_regulatory_evidence.tsv',t)
