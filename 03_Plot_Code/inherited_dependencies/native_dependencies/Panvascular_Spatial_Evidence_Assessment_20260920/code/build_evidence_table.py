"""Join existing results and raw-count descriptive audit. No hypothesis tests or refits."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path(__file__).resolve().parents[1];J=Path('/Users/taki/Desktop/Codex任务1/overnight_manuscript_integration_20260917/spatial_target_excluded_joint_v1');C=Path('/Users/taki/Documents/Panvascular_Figure_Archive_20260919/Panvascular_Nature_Reconstruction_20260918/SourceData/spatial/candidate29.tsv')
c=pd.read_csv(C,sep='\t')[['gene_id','gene_symbol','cell_state']].rename(columns={'gene_symbol':'gene','cell_state':'cell_context'});assert len(c)==29 and c.gene.is_unique
r=pd.read_csv(O/'tables/existing_raw_gene_section_descriptives.tsv',sep='\t');assert len(r)==29*49
s=pd.read_csv(J/'section_input_audit.tsv',sep='\t');m=pd.read_csv(O/'code/native_sources/sample_manifest.tsv',sep='\t');assert s.section_id.is_unique and m.section_id.is_unique
# Audit metadata_status is the one bound to the actual fitted result; retain curated status separately.
d=r.merge(c,on='gene',validate='many_to_one').merge(s[['section_id','biological_unit_id','metadata_status','QC_fallback','n_spatial_blocks']],on='section_id',validate='many_to_one')
d=d.merge(m[['section_id','sample_id','metadata_status','condition','tissue','metadata_note']].rename(columns={'metadata_status':'curated_metadata_status'}),on='section_id',validate='many_to_one')
d['donor_identity_status']=np.where(d.metadata_status.str.contains('unique_library_identifier',na=False),'未提供_文库身份不等于已确认供者','按既有记录单位映射_保留原metadata_status')
d['gene_id_scope']='candidate Ensembl ID; source raw artifact uses symbol; native ID separately retained'
roster=pd.read_csv(J/'target_model_roster.tsv',sep='\t').query("posterior=='q05'");assert not roster.duplicated(['section_id','gene']).any()
rc=['section_id','gene','pathology','target_removed_from_program','n_program_genes','gene_present_in_source','model_status'];d=d.merge(roster[rc].rename(columns={'pathology':'model_expression_program','model_status':'current_model_status'}),on=['section_id','gene'],validate='one_to_one')
models=pd.read_csv(J/'target_section_models.tsv',sep='\t').query("posterior=='q05'")
for term in ['cell_abundance','pathology','interaction']:
 q=models[models.term.eq(term)][['section_id','gene','beta','se','z','p','model']];q=q.rename(columns={k:'section_'+term+'_'+k for k in ['beta','se','z','p','model']});d=d.merge(q,on=['section_id','gene'],how='left',validate='one_to_one')
 d['section_'+term+'_CI_status']='未提供_原逐切片表无CI_本轮未推导';d['section_'+term+'_FDR_status']='未提供_原逐切片表只提供名义P'
d['section_association_status']=np.where(d.current_model_status.eq('FIT_OK'),'已有逐切片模型_原P未更改','未检验或不可估_'+d.current_model_status.fillna('未提供'))
meta=pd.read_csv(J/'target_meta_descriptive.tsv',sep='\t').query("posterior=='q05' and sensitivity=='all_sections' and term=='cell_abundance'")
meta=meta[['gene','beta','se','p','CI_low','CI_high','FDR_BH_all_29_pairs_x_3_terms','I2','k','n_datasets']].rename(columns={k:'overall_record_unit_meta_'+k for k in ['beta','se','p','CI_low','CI_high','FDR_BH_all_29_pairs_x_3_terms','I2','k','n_datasets']});d=d.merge(meta,on='gene',how='left',validate='many_to_one');d['overall_meta_scope']='总体记录单位meta_不是该切片显著性_不是研究层meta'
nat=pd.read_csv(O/'tables/native_candidate_structure_long.tsv',sep='\t');nat=nat.drop(columns=['GSE']);assert not nat.duplicated(['gene','section_id']).any();d=d.merge(nat,on=['gene','section_id'],how='left',validate='one_to_one')
d['native_structure_status']=np.where(d.native_I.notna(),'已有作者原生表达Moran_切片集合/变换/邻接未完整提供','未提供_本轮未运行空间结构检验')
his=pd.read_csv(O/'tables/historical_candidate_residual_structure_join.tsv',sep='\t');hc=['gene','section_id']+[k for k in his.columns if k.startswith('historical_')];d=d.merge(his[hc],on=['gene','section_id'],how='left',validate='one_to_one')
d['historical_residual_status']=np.where(d.historical_nb_residual_moran_I.notna(),'已有历史残差诊断_不是当前候选表达自相关检验','未提供_无匹配历史结果')
p=pd.read_csv(O/'code/native_sources/pathology_availability.tsv',sep='\t').rename(columns={'GSE':'study'});d=d.merge(p[['study','public_spot_pathology_table','public_nonspot_label']],on='study',validate='many_to_one');d['histology_region_enrichment_status']='未提供_独立逐spot组织分区与对应区域富集检验未找到';d['model_program_is_histology_label']=False
# Native gene-symbol count matches are only descriptive consistency; do not infer barcode identity.
d['native_current_target_total_UMI_difference']=d.target_total_UMI-d.native_gene_counts_under_tissue;d['native_current_positive_spots_difference']=d.n_positive-d.native_detected_barcodes
d['current_display_section']=((d.gene.eq('LIPA')&d.section_id.eq('GSE248608::GSM7917895_CA2'))|(d.gene.eq('STAT3')&d.section_id.eq('GSE314851::GSM9414344_LCA')))
order={g:i for i,g in enumerate(c.gene)};d['_order']=d.gene.map(order);d=d.sort_values(['_order','study','section_id']).drop(columns='_order')
front=['gene','gene_id','cell_context','study','section_id','sample_id','biological_unit_id','metadata_status','donor_identity_status','condition','tissue','n_spots','gene_available','n_positive','positive_fraction','target_total_UMI','raw_count_value_status','current_model_status','model_expression_program','section_cell_abundance_beta','section_cell_abundance_se','section_cell_abundance_p','section_cell_abundance_CI_status','section_cell_abundance_FDR_status','native_structure_status','native_I','native_P_original','native_FDR_original','native_spot_universe','native_transform','native_adjacency','historical_residual_status','historical_nb_residual_moran_I','historical_nb_residual_moran_permutation_p','historical_nb_residual_moran_FDR','histology_region_enrichment_status']
d=d[front+[k for k in d.columns if k not in front]];d.to_csv(O/'tables/gene_section_spatial_evidence.tsv',sep='\t',index=False,na_rep='未提供')
d[d.gene.isin(['LIPA','STAT3'])].to_csv(O/'tables/LIPA_STAT3_全部98个gene_section证据.tsv',sep='\t',index=False,na_rep='未提供')
d[d.current_display_section].to_csv(O/'tables/LIPA_STAT3_当前展示切片证据.tsv',sep='\t',index=False,na_rep='未提供')
summary=[]
for g,q in d.groupby('gene',sort=False):
 rr=meta[meta.gene.eq(g)].iloc[0];summary.append({'gene':g,'cell_context':q.cell_context.iloc[0],'n_sections':len(q),'n_gene_measured_sections':int(q.gene_available.sum()),'n_sections_with_detected_expression':int((q.n_positive>0).sum()),'n_current_association_fitted_sections':int(q.current_model_status.eq('FIT_OK').sum()),'native_expression_Moran_results':int(q.native_I.notna().sum()),'historical_residual_Moran_results':int(q.historical_nb_residual_moran_I.notna().sum()),'histology_confirmed_region_enrichment':'未提供',**{k:rr[k] for k in rr.index if k!='gene'}})
pd.DataFrame(summary).to_csv(O/'tables/candidate_evidence_level_summary.tsv',sep='\t',index=False,na_rep='未提供')
checks={'gene_section_rows':len(d),'unique_genes':d.gene.nunique(),'sections':d.section_id.nunique(),'studies':d.study.nunique(),'recorded_units':d.biological_unit_id.nunique(),'rows_with_source_gene_missing':int((~d.gene_available).sum()),'native_rows':int(d.native_I.notna().sum()),'historical_residual_rows':int(d.historical_nb_residual_moran_I.notna().sum()),'new_tests_or_P_adjustments':False,'max_existing_normalized_value_abs_error':float(d.normalised_max_abs_difference.max()),'candidates_source':str(C),'candidates_sha256':hashlib.sha256(C.read_bytes()).hexdigest()};(O/'notes/evidence_table_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2));print(checks)
