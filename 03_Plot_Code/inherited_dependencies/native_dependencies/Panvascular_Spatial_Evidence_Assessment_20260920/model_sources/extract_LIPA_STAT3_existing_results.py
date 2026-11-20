"""Read-only result inventory. No model fitting, meta-analysis or significance recomputation."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,itertools
O=Path(__file__).resolve().parents[1];S=O/'model_sources';T=O/'tables';T.mkdir(exist_ok=True);genes=['LIPA','STAT3']
read=lambda x:pd.read_csv(S/x,sep='\t',float_precision='round_trip')
m=read('target_meta_descriptive.tsv');m=m[m.gene.isin(genes)].copy();m['Figure5c_row']=m.gene.map({'LIPA':14,'STAT3':23});m['is_Figure5c_primary']=m.posterior.eq('q05')&m.sensitivity.eq('all_sections')&m.term.eq('cell_abundance');m.to_csv(T/'模型LIPA_STAT3_全部24项meta.tsv',sep='\t',index=False,float_format='%.17g')
m[m.is_Figure5c_primary].to_csv(T/'模型LIPA_STAT3_Figure5c主分析精确值.tsv',sep='\t',index=False,float_format='%.17g')
r=read('target_model_roster.tsv');r=r[r.gene.isin(genes)].copy();a=read('section_input_audit.tsv');ss=read('target_section_models.tsv');ss=ss[ss.gene.isin(genes)];u=read('target_biological_unit_models.tsv');u=u[u.gene.isin(genes)]
r.to_csv(T/'模型LIPA_STAT3_全部196切片后验资格记录.tsv',sep='\t',index=False)
rr=r.merge(pd.DataFrame({'term':['cell_abundance','pathology','interaction']}),how='cross');keys=['GSE','section_id','biological_unit_id','posterior','gene','cell_context','pathology','term'];rr=rr.merge(ss[keys+['beta','se','z','p','model']],on=keys,how='left',validate='one_to_one');rr.to_csv(T/'模型LIPA_STAT3_全部588切片term含不可估.tsv',sep='\t',index=False,float_format='%.17g')
a.to_csv(T/'模型全部49切片35记录单位元数据.tsv',sep='\t',index=False)
# Preserve every unit and sensitivity setting, including units with no eligible sections.
units=a[['GSE','biological_unit_id','metadata_status']].drop_duplicates();assert len(units)==35
configs=pd.DataFrame(list(itertools.product(genes,['q05','mean'],['all_sections','exclude_QC_fallback'],['cell_abundance','pathology','interaction'])),columns=['gene','posterior','sensitivity','term']);uc=units.merge(configs,how='cross');uc=uc.merge(u,on=['GSE','biological_unit_id','gene','posterior','sensitivity','term'],how='left',validate='one_to_one');uc['available_original_unit_estimate']=uc.beta.notna();uc.to_csv(T/'模型LIPA_STAT3_全部840记录单位term含不可估.tsv',sep='\t',index=False,float_format='%.17g')
# Study summaries describe the FULL distribution of existing estimates; no new pooled beta/P.
rows=[]
for gene,post,sens,term,gse in itertools.product(genes,['q05','mean'],['all_sections','exclude_QC_fallback'],['cell_abundance','pathology','interaction'],sorted(a.GSE.unique())):
 ro=r[(r.gene==gene)&(r.posterior==post)&(r.GSE==gse)];sec=ss[(ss.gene==gene)&(ss.posterior==post)&(ss.term==term)&(ss.GSE==gse)]
 if sens=='exclude_QC_fallback':ro=ro[~ro.QC_fallback];sec=sec[~sec.QC_fallback]
 un=u[(u.gene==gene)&(u.posterior==post)&(u.sensitivity==sens)&(u.term==term)&(u.GSE==gse)]
 z=dict(gene=gene,posterior=post,sensitivity=sens,term=term,GSE=gse,n_sections_in_scope=len(ro),n_fit_sections=len(sec),n_gate_not_met=int(ro.model_status.eq('INPUT_GATE_NOT_MET').sum()),n_failed=int(ro.model_status.eq('FIT_FAILED').sum()),n_unit_estimates=len(un),study_pooled_beta='NOT_COMPUTED_IN_ORIGINAL_PIPELINE',study_pooled_P='NOT_AVAILABLE')
 for level,df in [('section',sec),('unit',un)]:
  z.update({f'{level}_beta_min':df.beta.min(),f'{level}_beta_q25':df.beta.quantile(.25),f'{level}_beta_median_descriptive':df.beta.median(),f'{level}_beta_q75':df.beta.quantile(.75),f'{level}_beta_max':df.beta.max(),f'{level}_positive_count':int(df.beta.gt(0).sum()),f'{level}_negative_count':int(df.beta.lt(0).sum()),f'{level}_zero_count':int(df.beta.eq(0).sum())})
 rows.append(z)
pd.DataFrame(rows).to_csv(T/'模型LIPA_STAT3_全部168研究分布_非研究meta.tsv',sep='\t',index=False,float_format='%.17g')
# Verification against the adopted submission source table, not inference from rendered rounding.
legacy=Path('/Users/taki/Desktop/投稿包_最新版_20260919/Source_Data_Figures/Figure5/target_meta_descriptive.tsv');lm=pd.read_csv(legacy,sep='\t',float_precision='round_trip');cols=['posterior','gene','cell_context','pathology','term','sensitivity'];q=m.merge(lm,on=cols,suffixes=('_copied_server','_adopted_figure'),validate='one_to_one');fields=['beta','se','p','CI_low','CI_high','FDR_BH_all_29_pairs_x_3_terms'];assert all(np.array_equal(q[c+'_copied_server'],q[c+'_adopted_figure']) for c in fields)
verification={'adopted_source':str(legacy),'all_24_rows_exact_numeric_match':True,'matched_fields':fields,'model_roster_rows':len(r),'section_term_rows_including_ineligible':len(rr),'unit_term_rows_including_unavailable':len(uc),'study_distribution_rows':len(rows),'new_models_or_P_tests':False}
(O/'notes/model_inventory_checks.json').write_text(json.dumps(verification,indent=2))
# Local preservation hashes and server-bound executable hashes.
hashes=json.loads((S/'SOURCE_SHA256.json').read_text());manifest=[]
for p in S.iterdir():
 if not p.is_file():continue
 h=hashlib.sha256(p.read_bytes()).hexdigest();matches=[k for k,v in hashes.items() if Path(k).name==p.name and v==h];manifest.append(dict(local_file=str(p),sha256=h,bound_server_path=';'.join(matches)))
pd.DataFrame(manifest).to_csv(T/'模型原始来源文件sha256.tsv',sep='\t',index=False)
print(verification)
