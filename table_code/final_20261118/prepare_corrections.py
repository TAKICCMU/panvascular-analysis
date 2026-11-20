from pathlib import Path
import json, hashlib, shutil, re
import openpyxl

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
O=P/'图件编辑资料/补充表核查_20261118'
R=Path('/Users/taki/Desktop/投稿包_最新版_20260919/图件编辑资料/全流程代码审核_20260922/Panvascular_code_review')
W=openpyxl.load_workbook(P/'Supplementary_Tables.xlsx',read_only=False,data_only=False)
S=O/'staging'; S.mkdir(exist_ok=True)
B=O/'backup_before_correction'; B.mkdir(exist_ok=True)
changed_files=['Supplementary_Tables.xlsx','Manuscript.docx','Manuscript.pdf','Figure_Legends.docx','Figure_Legends.pdf','Supplementary_Information.docx','Supplementary_Information.pdf','Figure/Figure_Legends.docx','Figure/Figure_Legends.md','Source_Data_Figures/ExtendedData1/benchmark_publication_CPU_only.tsv']
hashes={}
for f in changed_files:
    p=P/f
    if not p.exists():continue
    hashes[f]=hashlib.sha256(p.read_bytes()).hexdigest()
    q=B/f;q.parent.mkdir(parents=True,exist_ok=True)
    if q.exists():assert hashlib.sha256(q.read_bytes()).hexdigest()==hashes[f],f
    else:shutil.copy2(p,q)
protected={str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for directory in ['Figure','Figures','Extended_Data','Source_Data_Figures','Supporting_Visuals'] for p in (P/directory).rglob('*') if p.is_file() and str(p.relative_to(P)) not in changed_files}
(O/'correction_input_hashes.json').write_text(json.dumps({'affected':hashes,'protected':protected},ensure_ascii=False,indent=2))

changes=[]
def cell(sheet,address,value):
    old=W[sheet][address].value
    if old!=value:changes.append(dict(sheet=sheet,cell=address,before=old,after=value))

t8='L1 fractions averaged first across samples within each tissue and recorded donor identifier, then across donor identifiers. Figure 1e instead averages across samples.'
cell('Information','C9',t8)
cell('Information','D9','One tissue category; columns are mean L1 fractions across donor identifiers (0–1).')
for row in range(113,130):
    name=W['Data dictionary'].cell(row,2).value
    if row==113:
        description='Tissue category used for the donor-averaged composition summary.';unit='tissue category'
    else:
        description=f'Mean fraction of {name} cells. Per-sample fractions use all atlas cells in that sample as denominator; fractions are averaged within each tissue and recorded donor identifier, then across donor identifiers. Figure 1e uses sample averages instead.';unit='fraction (0–1)'
    cell('Data dictionary',f'C{row}',description);cell('Data dictionary',f'D{row}',unit)
cell('Data dictionary','C83',W['Data dictionary']['C83'].value.replace('Figure 1g','Figure 1h'))

families=[
['MAGMA cell enrichment','GWAS scenario × annotation level × weighting','BH','15 L1 or 103 L2 eligible categories. The Average covariate is excluded from the testing family.'],
['scDRS association and heterogeneity','GWAS scenario × annotation level × test type','BH','Association and heterogeneity P values are corrected separately across the 15 L1 or 103 L2 eligible categories.'],
['Pseudobulk candidate meta-analysis','Cell group × tested genes in the GWAS-priority set','BH','Target-set family used for candidate selection. This is separate from the full transcriptome pathway ranking.'],
['Spatial candidate associations','Posterior summary × section/QC setting','BH','29 candidate pairs × 3 model terms = 87 tests per setting. Nonestimable terms contribute P=1 without an imputed effect.'],
['Pseudobulk GSEA','Analysis ID × pathway database','Native BH within each fgsea call','All testable pathways are included before display selection. Supplementary Table 22 is the index and significant-result view.'],
['Pseudobulk differential expression','Individual disease contrast × annotation level × cell group','BH','All genes passing the original expression filter within the fitted model.'],
['MAGMA gene significance','GWAS scenario × genes tested in that scenario','Bonferroni; BH where reported','The native gene-test universe defines the denominator; the separate 17,885-symbol pathway rank is not substituted for this family.'],
['MAGMA GSEA','Genetic ranking × pathway database','Native BH within each fgsea call','All gene sets meeting the original size limits in that database; displayed pathways do not define a new family.'],
['MAGMA overrepresentation','Selected gene set × pathway database','BH after overlap ≥3 filter','The retained implementation filters to overlap ≥3 before correction within each database; Supplementary Table 21.'],
['Pseudobulk overrepresentation','Analysis ID × pathway database × effect direction','BH before display overlap filter','All nonempty gene sets evaluable against the original background. The ≥3 overlap rule is applied only to the display/index; Supplementary Table 22.'],
['CoVarNet phenotype associations','Phenotype variable','BH','Module association tests are corrected within each phenotype variable. The paired condition comparison has its own family.'],
['CoVarNet matched null','Null statistic','BH','All 14 modules for each metric, using the saved empirical P values; Supplementary Table 33.'],
['CoVarNet coherence within studies','All evaluable module–node-pair meta-analysis tests','BH','The final coherence table is corrected jointly across module and node-pair combinations.'],
['CoVarNet disease associations','Pseudocount setting × result layer','BH','Within-study/tissue/background tests and cross-study meta-analysis tests are corrected separately, each across all stored module comparisons.'],
['CoVarNet coherence within tissue','All evaluable module × tissue × node-pair meta-analysis tests','BH','Joint correction across the stored comparisons; nonestimable and blocked rows retain their original status.'],
['Disease-related cell composition','Annotation level × disease background × denominator','BH','Exact families are recorded in the FDR family field in Supplementary Table 43; only the corresponding stored estimable tests enter each family.'],
['Spatial module coherence','Posterior summary × null method','BH','14 modules × 2 coherence metrics = 28 tests per family; Supplementary Table 40.'],
['Spatial module programme direction','Posterior summary','BH','14 modules × 6 programmes = 84 tests per family; Supplementary Table 56.'],
]
cell('Information','C60','Multiple-testing families and correction boundaries for the retained genetic, expression, composition, CoVarNet and spatial analyses.')
for row,name,desc in [(904,'Analysis','Analysis or test whose original multiplicity family is specified.'),(905,'Multiplicity family','Variables that partition the original tests into separate correction families.'),(906,'Correction','Original correction method and, where relevant, the order of filtering and adjustment.'),(907,'Unit / boundary','Tests included in each family, handling of nonestimable terms, and the corresponding current table.')]:
    cell('Data dictionary',f'C{row}',desc)

versions=[list(row) for row in W['Table 60'].iter_rows(min_row=3,values_only=True)]
evidence=[]
def add(name,value,scope,source):
    versions.append([name,value,scope]);evidence.append({'table_row':len(versions)+2,'name':name,'value':value,'scope':scope,'source':str(source)})
cpu=Path('/Users/taki/Desktop/Codex任务1/Panvascular_reanalysis_takeover_20260908/metadata_provenance_impact_audit_20260913_013529/evidence/06_reanalysis_takeover_20260908/CPU_full_atlas_theta2_20260911_160002_388642/pipeline_parameters.json')
params=json.loads(cpu.read_text())
for name,v in params['packages'].items():add(name,v,'Full atlas CPU theta=2 integration; pipeline_parameters.json',cpu)
add('Random seed',params['seed'],'Full atlas CPU integration, PCA, Harmony, neighbours, Leiden and UMAP',cpu)
mac=R/'analysis_code/native/upstream/L2_CPU/Macrophage_L2_round1_20260912_220350_741672/task.json'
add('Random seed',0,'Macrophage L2 CPU analysis; task.json and executed entry point',mac)
sc=R/'analysis_code/native/production_patches/scdrs_repair_20260914T0910/evidence/NATIVE_VERSIONS_AND_SOURCE.json'
sd=json.loads(sc.read_text())
evidence.append({'record':'scDRS original environment','source':str(sc),'content':sd})
for name,v in [('scDRS','1.0.2'),('Python','3.9.25'),('scanpy','1.9.3'),('anndata','0.8.0'),('numpy','1.23.5'),('scipy','1.10.1'),('pandas','1.5.3')]:add(name,v,'Native scDRS scoring and group enrichment environment',sc)
add('Random seed',0,'scDRS native score_cell default, retained by the CLI; 1,000 matched control sets',R/'analysis_code/native/production_patches/scdrs_repair_20260914T0910/evidence/method.py')
add('MAGMA','1.10','Native gene property analysis; production command log',O/'correction_evidence/magma_spatial_native_evidence.json')
gr=Path('/Users/taki/Desktop/Codex任务1/Panvascular_retained_reproduction_dependencies/f59290f1_Panvas_Nature_V1.9.3_submission_package/Submission_Audit_V1.9.2/genetics_server_reports/software_versions_genetics.tsv')
for n,v in [('R','4.3.1'),('GenomicSEM','0.0.5'),('lavaan','0.6.19')]:add(n,v,'Retained genetics environment snapshot (software_versions_genetics.tsv)',gr)
for n,v in [('R','4.4.3'),('susieR','0.16.6')]:add(n,v,'Adopted mismatch-aware SuSiE-RSS v3; per-region runtime sessionInfo',O/'correction_evidence/genetics_runtime_records.json')
for n,v in [('R','4.5.1'),('coloc','5.2.3')]:add(n,v,'Existing17 and new12 candidate colocalisation; original R_sessionInfo.txt',O/'correction_evidence/server_runtime_records.json')
for n,v in [('edgeR','4.0.16'),('metafor','5.0.1')]:add(n,v,'Manuscript-strategy pseudobulk DE and REML candidate meta-analysis; native preflight record',O/'correction_evidence/server_runtime_records.json')
add('Random seed',0,'Manuscript-strategy pseudobulk native.R',R/'analysis_code/native/downstream_final/pseudobulk_manuscript_strategy_20260916_v1/code/native.R')
for n,v in [('R','4.5.1'),('fgsea','1.36.2'),('msigdbr','26.1.0')]:add(n,v,'Retained MAGMA pathway analysis runtime',Path('/Users/taki/Desktop/Codex任务1/final_revision_work/pathway_enrichment/sessionInfo.txt'))
add('Random seed',20260822,'MAGMA GSEA and retained enrichment implementation',R/'analysis_code/native/downstream_base/external_legacy/pathway_enrichment.R')
add('R','4.3.1','Manuscript-strategy pseudobulk pathway runtime; fgsea 1.28.0 is recorded above',O/'correction_evidence/server_runtime_records.json')
add('Random seed',20260822,'Manuscript-strategy pseudobulk GSEA and ORA native_pathway.R',R/'analysis_code/native/downstream_final/pseudobulk_manuscript_strategy_enrichment_v1/code/native_pathway.R')
for n,v in [('R','4.3.3'),('NMF','0.28')]:add(n,v,'CoVarNet primary nsNMF runtime; stable_NMF/sessionInfo.txt',O/'correction_evidence/server_runtime_records.json')
add('Base random seed',20260806,'CoVarNet primary rank evaluation and nsNMF; deterministic offsets identify individual fits',R/'analysis_code/native/downstream_final/code/run_covarnet.py')
add('Random seed',20260826,'CoVarNet parent-lineage-preserving null',R/'analysis_code/native/downstream_final/code/run_covarnet.py')
add('Base random seed',20260807,'CoVarNet additional rank, bootstrap and depth sensitivities; offsets retained in remaining_refits.R',R/'analysis_code/native/downstream_final/remaining_downstream_20260916_v1/covarnet/code/remaining_refits.R')
gpu=R/'analysis_code/native/downstream_base/evidence/GPU_ENVIRONMENT_SNAPSHOT.json'
for n,v in [('torch','2.10.0+cu128'),('cell2location','0.1.5'),('scvi-tools','1.3.3'),('pyro-ppl','1.9.1'),('anndata','0.11.4')]:add(n,v,'Spatial reference and mapping GPU environment',gpu)
add('Random seed',20260814,'Spatial reference/mapping and target-excluded mapping configs',R/'analysis_code/native/downstream_final/overnight_manuscript_integration_20260917/spatial_target_excluded_v1/configs/config.json')

cell('Table 60','A1','Supplementary Table 60 | Software versions and random seeds')
cell('Table 60','A2','Software, resource or setting')
cell('Table 60','B2','Version or value')
cell('Table 60','C2','Workflow and recorded scope')
cell('Information','B61','Software versions and random seeds')
cell('Information','C61','Workflow-specific software versions, resources and random seeds from retained environment, parameter and execution records.')
cell('Information','D61','One software, resource or random-seed record; environments are distinguished by workflow.')
for row,field,desc,unit in [
    (908,'Software, resource or setting','Software package, reference resource or random-seed setting recorded for a specific workflow.','name / setting'),
    (909,'Version or value','Recorded package/resource version or numeric random seed. Version strings are identifiers, not quantitative measurements.','version identifier or integer seed'),
    (910,'Workflow and recorded scope','Workflow to which the version or seed applies and the type of retained provenance record. Different workflows may use different versions.','workflow / provenance context')]:
    cell('Data dictionary',f'B{row}',field);cell('Data dictionary',f'C{row}',desc);cell('Data dictionary',f'D{row}',unit)

numeric=json.loads((O/'numeric_text_cells.json').read_text())
for x in numeric:cell(x['sheet'],x['cell'],int(x['stored_value']))
spec={'changes':changes,'families':families,'versions':versions,'numeric_cells':numeric,'table24_singleline_columns':['AQ','AS'],'table24_row_height':60}
(O/'correction_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2))
(O/'correction_evidence/software_and_seed_sources.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2))
print(json.dumps({'planned_cells':len(changes),'families':len(families),'version_and_seed_rows':len(versions),'protected_files':len(protected)},ensure_ascii=False))
