from pathlib import Path
from collections import Counter
import json, csv, math, re, hashlib
import pandas as pd
import numpy as np
import openpyxl

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116'); O=P/'图件编辑资料/补充表核查_20261118'; S=P/'Source_Data_Figures'
W=openpyxl.load_workbook(P/'Supplementary_Tables.xlsx',read_only=True,data_only=True)
T={}; D={}
for n in range(1,61):
    vals=list(W[f'Table {n:02d}'].values);hs=vals[1]
    T[n]=pd.DataFrame([r for r in vals[2:] if any(v is not None for v in r)],columns=hs)
checks=[]; differences=[]; contextual={}
def save(name,obj):(O/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str))
def check(name,passed,detail):checks.append({'check':name,'pass':bool(passed),'detail':detail})
def rd(path):return pd.read_csv(S/path,sep='\t')
def norm(v):
    if v is None or (not isinstance(v,(list,dict)) and pd.isna(v)):return None
    if isinstance(v,str) and v in {'TRUE','FALSE'}:return v=='TRUE'
    return v
def equal(a,b):
    a,b=norm(a),norm(b)
    if a is None or b is None:return a is b
    if isinstance(a,(float,int,np.number)) and isinstance(b,(float,int,np.number)) and not isinstance(a,bool) and not isinstance(b,bool):return math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-14 if abs(a)>1e-10 and abs(b)>1e-10 else 0)
    return a==b
def comp(n,source,mp,keys,subset=None,label=None,source_subset=None):
    a=T[n].copy(); b=rd(source)
    if subset is not None:a=subset(a)
    if source_subset is not None:b=source_subset(b)
    b=b[list(mp.values())].rename(columns={v:k for k,v in mp.items()})
    tag=label or f'Table{n}_vs_{source}'
    missing=[];mis=[];compared=0
    if a.duplicated(keys).any() or b.duplicated(keys).any():
        check(tag,False,{'reason':'duplicate comparison key','workbook_duplicates':int(a.duplicated(keys).sum()),'source_duplicates':int(b.duplicated(keys).sum())});return
    aa={tuple(norm(r[k]) for k in keys):r for r in a.to_dict('records')}
    for row in b.to_dict('records'):
        key=tuple(norm(row[k]) for k in keys)
        if key not in aa:missing.append(key);continue
        ar=aa[key]
        for k in mp:
            compared+=1
            if not equal(ar[k],row[k]):mis.append({'key':key,'column':k,'workbook':ar[k],'source':row[k]})
    check(tag,not missing and not mis,{'source_rows':len(b),'workbook_rows_in_scope':len(a),'compared_cells':compared,'missing_keys':missing[:20],'mismatches':mis[:20],'source':str(S/source)})
    differences.extend([dict(table=n,source_file=str(source),**q) for q in mis])

comp(4,'ExtendedData1/benchmark_publication_CPU_only.tsv',{'Method':'method','Biology score':'biology_score','Batch score':'batch_score','Overall score':'overall_score','Rank':'rank'},['Method'])
comp(8,'Figure1/L1_tissue_sample_mean_fraction.tsv',{'Tissue':'tissue',**{h:h for h in T[8].columns if h!='Tissue'}},['Tissue'])
mp15={'Latent factor':'lhs','Parameter':'op','Observed phenotype':'rhs','Unstandardised estimate':'Unstand_Est','Unstandardised standard error':'Unstand_SE','Standardised estimate (genotype scale)':'STD_Genotype','Standard error (genotype scale)':'STD_Genotype_SE','Fully standardised estimate':'STD_All','P':'p_value'}
trait_names={'AA':'aortic dissection','AAA':'abdominal aortic aneurysm','TAAD':'thoracic aortic aneurysm','CAD':'coronary artery disease','Stroke':'stroke','PAD':'peripheral artery disease','CA':'carotid atherosclerosis','CAVS':'calcific aortic valve stenosis'}
comp(15,'Figure2/Figure2_Fig3d_bifactor_loadings.tsv',mp15,['Latent factor','Parameter','Observed phenotype'],source_subset=lambda x:x.assign(rhs=x.rhs.map(trait_names),op=x.op.replace({'=~':'loading'})))
common={'Gene':'gene_symbol','Annotation level':'annotation_level','Cell state':'cell_state','Disease contexts (n)':'n_disease_contexts','GWAS-aligned contexts (n)':'n_GWAS_aligned_contexts','Disease contexts':'disease_contexts','Meta log2FC':'meta_logFC','Meta SE':'meta_SE','Meta Z':'meta_Z','Meta P':'meta_P','Meta 95% CI lower':'meta_CI_low','Meta 95% CI upper':'meta_CI_high','I² (%)':'I2','Direction agreement':'direction_consistency','Consensus direction':'consensus_direction','Significant contexts (n)':'n_context_FDR_lt_0.05','Median context log2FC':'median_context_logFC','Ensembl gene ID':'gene_id','Target-set BH FDR':'target_set_meta_FDR'}
comp(24,'Figure3/Fig3d_e_all29.tsv',common,['Gene','Annotation level','Cell state'])
for lev,filename in [('L1','Figure3/Fig3a_L1.tsv'),('L2','Figure3/Fig3b_L2_selected19.tsv')]:
    mp={'Annotation level':'annotation_level','GWAS scenario':'scenario','Weighting':'weighting','MAGMA variable ID':'VARIABLE','MAGMA β':'BETA','MAGMA standardised β':'BETA_STD','MAGMA SE':'SE','MAGMA P':'P','Group':'group','MAGMA BH FDR':'FDR_BH','BH family':'BH_family','Family tests (n)':'n_tests_in_family','Cells (n)':'n_cell','ScDRS Z':'primary_r0.1__assoc_mcz','ScDRS P':'primary_r0.1__assoc_mcp','ScDRS FDR':'primary_r0.1__assoc_mcp_BH'}
    comp(19,filename,mp,['Annotation level','GWAS scenario','Weighting','Group'])
mp28={'Gene':'gene','GWAS scenario':'scenario','Context':'context','QTL source':'qtl_source','QTL sample size':'qtl_sample_size','Shared variants (n)':'n_overlap','GWAS effective N':'gwas_effective_N','Minimum GWAS P':'min_gwas_p','Minimum eQTL P':'min_eqtl_p','PP.H0':'PP.H0','PP.H1':'PP.H1','PP.H2':'PP.H2','PP.H3':'PP.H3','PP.H4':'PP.H4','PP.H4 / (PP.H3+PP.H4)':'H4_given_H3H4','Evidence class':'classification','Primary prior p12':'primary_prior_p12'}
comp(28,'Figure3/Fig3f_strong5_contexts.tsv',mp28,['Gene','GWAS scenario','Context'])
cols30=rd('Figure4_current/14_final_module_evidence_and_tiers.tsv').columns.tolist()
comp(30,'Figure4_current/14_final_module_evidence_and_tiers.tsv',dict(zip(T[30].columns,cols30[:34])),['Module'])
sp={'Gene':'gene','Cell context':'cell_context','Units (n)':'k','β':'beta','SE':'se','Z':'z','P':'p','τ²':'tau2','I² (%)':'I2','Posterior summary':'posterior','Expression programme':'pathology','Model term':'term','Section set':'sensitivity','Datasets (n)':'n_datasets','95% CI lower':'CI_low','95% CI upper':'CI_high','BH FDR (87 tests)':'FDR_BH_all_29_pairs_x_3_terms'}
sk=['Gene','Cell context','Posterior summary','Expression programme','Model term','Section set']
comp(36,'Figure4_current/all29_original_meta.tsv',sp,sk,source_subset=lambda x:x[(x.posterior=='q05')&(x.sensitivity=='all_sections')])
comp(37,'Figure4_current/all29_original_meta.tsv',sp,sk)
comp(40,'Figure4_current/module_coherence_meta.tsv',dict(zip(T[40].columns,rd('Figure4_current/module_coherence_meta.tsv').columns)),['Posterior summary','Module','Null method','Metric'])
comp(43,'ExtendedData10/effects_with_estimability_and_BH.tsv',{'Annotation level':'level','State':'state','Model id':'model_id','Estimate':'estimate','SE':'SE','Z':'z','P':'P','Model term':'term','FDR family':'FDR_family','Inference status':'inference_status','BH FDR':'FDR_BH'},['Annotation level','State','Model id','Model term'])
comp(44,'ExtendedData10/model_audit_corrected.tsv',{'Annotation level':'level','State':'state','Model id':'model_id','Background':'background','Denominator':'denominator','Samples (n)':'n_samples','Donors (n)':'n_donors','Status':'status','Formula':'formula','Scope kind':'scope_kind'},['Annotation level','State','Model id'])
mp52={'Gene':'gene','Cell context':'cell_context','Posterior summary':'posterior','Expression programme':'pathology','Model term':'term','Study ID':'GSE','Analysis unit ID':'biological_unit_id','Section set':'sensitivity','β':'beta','SE':'se','Sections (n)':'n_sections','Section covariance unmodelled':'serial_section_covariance_unavailable'}
comp(52,'ExtendedData9/target_biological_unit_models.tsv',mp52,['Gene','Cell context','Posterior summary','Expression programme','Model term','Study ID','Analysis unit ID','Section set'])
mp53={'Gene':'gene','Cell context':'cell_context','Study ID':'GSE','Section ID':'section_id','Analysis unit ID':'biological_unit_id','Metadata status':'metadata_status','QC fallback':'QC_fallback','Posterior summary':'posterior','Expression programme':'pathology','Expected direction':'expected_direction','Target removed from program':'target_removed_from_program','N program genes':'n_program_genes','Spots (n)':'n_spots','Positive spots (n)':'n_positive_spots','Gene present in source':'gene_present_in_source','Model status':'model_status','Model term':'term','β':'beta','SE':'se','Z':'z','P':'p','Model':'model'}
comp(53,'ExtendedData9/target_section_models.tsv',mp53,['Gene','Cell context','Study ID','Section ID','Posterior summary','Model term'])
comp(56,'Figure4_current/module_program_direction_meta.tsv',dict(zip(T[56].columns,rd('Figure4_current/module_program_direction_meta.tsv').columns)),['Posterior summary','Module','Program'])

# Internal table consistency: only arithmetic/reconciliation of saved numbers, no refitting or new FDR.
for n,keys in [(2,['Sample ID']),(3,['Sample ID']),(6,['L1 category','L2 state']),(7,['L1 category','L2 state']),(19,['Annotation level','GWAS scenario','Weighting','Group']),(24,['Gene','Annotation level','Cell state']),(26,['Gene','Annotation level','Cell state']),(35,['Section ID']),(36,sk),(37,sk)]:
    check(f'unique_key_table{n}',not T[n].duplicated(keys).any(),{'keys':keys,'duplicates':int(T[n].duplicated(keys).sum())})
for na,nb,keys,fields in [(2,3,['Sample ID'],['Study ID','Cells (n)','Specimen ID','Donor ID','Tissue']), (6,7,['L1 category','L2 state'],[h for h in T[7].columns if h not in ['L1 category','L2 state']]),(24,26,['Gene','Annotation level','Cell state'],list(common)[3:]),(30,32,['Module'],[h for h in T[32].columns if h!='Module'])]:
    a=T[na].set_index(keys);b=T[nb].set_index(keys);diff=[]
    for k,r in a.iterrows():
        if k not in b.index: diff.append((str(k),'missing'));continue
        for h in fields:
            if not equal(r[h],b.loc[k,h]):diff.append((str(k),h,str(r[h]),str(b.loc[k,h])))
    check(f'internal_table{na}_table{nb}',not diff,{'rows':len(a),'fields':len(fields),'differences':diff[:10]})
for n,est,se,lo,hi in [(24,'Meta log2FC','Meta SE','Meta 95% CI lower','Meta 95% CI upper'),(26,'Meta log2FC','Meta SE','Meta 95% CI lower','Meta 95% CI upper'),(36,'β','SE','95% CI lower','95% CI upper'),(37,'β','SE','95% CI lower','95% CI upper')]:
    a=T[n].dropna(subset=[est,se,lo,hi]);bad=a[(a[lo]>a[est])|(a[hi]<a[est])|(a[lo]>a[hi])]
    check(f'CI_order_table{n}',len(bad)==0,{'complete_rows':len(a),'bad_rows':(bad.index+3).tolist()})
    # Report empirical stored CI/SE relationship; does not derive or replace any output.
    ratios=((a[hi]-a[lo])/(2*a[se])).replace([np.inf,-np.inf],np.nan).dropna()
    contextual[f'table{n}_CI_SE_ratio_range']=[float(ratios.min()),float(ratios.max())]
for n,p,q in [(19,'MAGMA P','MAGMA BH FDR'),(19,'ScDRS P','ScDRS FDR'),(21,'P','Adjusted P'),(22,'P','Adjusted P'),(24,'Meta P','Target-set BH FDR'),(26,'Meta P','Target-set BH FDR'),(36,'P','BH FDR (87 tests)'),(37,'P','BH FDR (87 tests)'),(40,'P','BH FDR (28 tests)'),(43,'P','BH FDR'),(56,'P direction positive','BH FDR (84 tests)')]:
    a=T[n].dropna(subset=[p,q]);bad=a[a[q]+1e-12<a[p]]
    check(f'FDR_not_less_than_P_table{n}_{p}',len(bad)==0,{'complete_rows':len(a),'bad_rows':(bad.index+3).tolist()[:15]})

# Reproduce only tabulation of already saved eligibility/significance flags.
contextual['Table19_families']=T[19].groupby(['Annotation level','GWAS scenario','Weighting']).agg(tests=('MAGMA P','size'),magma_sig=('MAGMA BH FDR',lambda s:int((s<.05).sum())),scdrs_sig=('ScDRS FDR',lambda s:int((s<.05).sum()))).reset_index().to_dict('records')
contextual['Table20']=T[20].to_dict('records')
contextual['spatial_counts']=T[36].assign(significant=T[36]['BH FDR (87 tests)']<.05).groupby('Model term')['significant'].sum().to_dict()
ab=T[36][T[36]['Model term']=='cell_abundance']; sig=ab[ab['BH FDR (87 tests)']<.05]
contextual['significant_spatial_abundance']=sig[['Gene','Cell context','β','BH FDR (87 tests)']].to_dict('records')
check('spatial_10_8_1',contextual['spatial_counts']=={'cell_abundance':10,'interaction':1,'pathology':8},contextual['spatial_counts'])
contextual['pooled_FN1_SCARB1']=ab[ab.Gene.isin(['FN1','SCARB1'])].to_dict('records')
u=T[52];u=u[(u['Gene'].isin(['FN1','SCARB1']))&(u['Posterior summary']=='q05')&(u['Section set']=='all_sections')&(u['Model term']=='cell_abundance')]
cnt=u.groupby('Gene').size().to_dict();uniq=u['Analysis unit ID'].nunique()
check('Figure4c_60_estimates_34_units',cnt=={'FN1':34,'SCARB1':26} and uniq==34,{'gene_unit_counts':cnt,'unique_record_ids':uniq})
check('spatial_49_sections_35_units',len(T[35])==49 and T[35]['Analysis unit ID'].nunique()==35,{'sections':len(T[35]),'units':T[35]['Analysis unit ID'].nunique(),'studies':T[35]['Study ID'].nunique()})
contextual['source_overlap_50']=T[50].to_dict('records')
contextual['true_spatial_rows']=len(T[52]);contextual['all_section_model_status']=T[54]['Model status'].value_counts(dropna=False).to_dict()
contextual['row_units_21']=T[21].groupby(['Analysis','Database'],dropna=False).agg(rows=('Pathway','size'),sig=('Adjusted P',lambda s:int((s<.05).sum()))).reset_index().to_dict('records')
contextual['recordtypes_22']=T[22]['Record type'].value_counts(dropna=False).to_dict()

save('source_comparisons.json',checks);save('source_mismatches.json',differences);save('key_result_reconciliation.json',contextual)
OLD=Path('/Users/taki/Desktop/投稿包_最新版_20260919/图件编辑资料/全流程代码审核_20260922/Panvascular_code_review')
comp(8,OLD/'analysis_code/native/downstream_final/descriptive_updated_v2/L1_tissue_donor_mean_fraction.tsv',{'Tissue':'tissue',**{h:h for h in T[8].columns if h!='Tissue'}},['Tissue'],label='Table8_matches_existing_donor_average')
mp5={'L1 category':'level1_annotation','Gene identifier':'gene_id','Mean log1p expression':'mean_log1p','Detected fraction':'detection_fraction','Other-cell mean log1p':'mean_log1p_rest','Other-cell detected fraction':'detection_fraction_rest','Mean expression difference':'mean_difference','Cells (n)':'n_cells','Other cells (n)':'rest_n_cells','Reference scope':'reference_scope','Descriptive rank':'descriptive_rank'}
comp(5,'Current_Display_20260920/marker_audit_20260921/original_all392_ordered.tsv',mp5,['L1 category','Gene identifier'],subset=lambda x:x[x['Record type']=='L1 canonical-union descriptive expression'])
comp(21,'Current_Display_20260920/Figure3c_current/MAGMA_complete_GSEA.tsv',{'Pathway':'pathway','P':'pval','Adjusted P':'padj','NES':'NES','Gene-set size (n)':'size','Database':'database','Ranked genes (n)':'n_ranked_genes'},['Pathway','Database'],subset=lambda x:x[x.Analysis=='preranked_GSEA'])

# Full Table 22 reconciliation against the 93 already generated result folders.
# Detail rows are explicitly BH_FDR_lt_0.05; the ORA index also applies overlap>=3.
# Check those distinct saved views separately. No test or FDR is recomputed.
I=T[22][T[22]['Record type']=='ANALYSIS_INDEX'].set_index('Analysis ID')
R=T[22][T[22]['Record type']=='BH_FDR_lt_0.05']
enrichment_checks=[]; detail_comparisons=0; total_gsea=0;total_ora=0;total_saved=0
for aid, idx in I.iterrows():
    folder=P/'Supplementary_Table22_full_results'/aid
    for typ,filename,mp in [('GSEA','fgsea_all.tsv.gz',{'Pathway':'pathway','P':'pval','Adjusted P':'padj','log2 error estimate':'log2err','ES':'ES','NES':'NES','Gene-set size (n)':'size','Leading-edge genes':'leadingEdge','Database':'database','Ranked genes (n)':'n_ranked_genes'}),('ORA','ora_all_tested_before_display_filter.tsv.gz',{'Pathway':'pathway','P':'pval','Adjusted P':'padj','Pathway genes (n)':'pathway_size','Overlap (n)':'overlap','Database':'database','Direction':'direction','Background genes (n)':'n_background_genes','Selected genes (n)':'n_selected_genes','Display overlap ge3':'display_overlap_ge3'})]:
        f=folder/filename
        if not f.exists():enrichment_checks.append({'analysis_id':aid,'kind':typ,'pass':False,'missing_file':str(f)});continue
        df=pd.read_csv(f,sep='\t'); total_gsea+=len(df) if typ=='GSEA' else 0;total_ora+=len(df) if typ=='ORA' else 0
        recorded=df[df.padj<.05].rename(columns={v:k for k,v in mp.items()})
        index_sig=int(((df.padj<.05)&((df.overlap>=3) if typ=='ORA' else True)).sum())
        wr=R[(R['Analysis ID']==aid)&(R['Analysis kind']==typ)]
        keys=['Pathway','Database']+(['Direction'] if typ=='ORA' else [])
        a={tuple(norm(r[k]) for k in keys):r for r in wr.to_dict('records')}
        b={tuple(norm(r[k]) for k in keys):r for r in recorded.to_dict('records')}
        bad=[]
        for key in a.keys()&b.keys():
            for col in mp:
                detail_comparisons+=1
                if not equal(a[key][col],b[key][col]):bad.append({'key':key,'column':col,'workbook':a[key][col],'source':b[key][col]})
        nfield='GSEA pathways (n)' if typ=='GSEA' else 'ORA pathways (n)'
        sfield='GSEA pathways, FDR<0.05 (n)' if typ=='GSEA' else 'ORA pathways, FDR<0.05 and overlap≥3 (n)'
        same=(len(df)==idx[nfield] and index_sig==idx[sfield] and set(a)==set(b) and len(a)==len(wr) and len(b)==len(recorded) and not bad)
        enrichment_checks.append({'analysis_id':aid,'kind':typ,'pass':same,'native_rows':len(df),'index_native_rows':int(idx[nfield]),'native_BH_significant_rows':len(recorded),'native_index_rule_significant_rows':index_sig,'index_significant_rows':int(idx[sfield]),'workbook_rows':len(wr),'workbook_only_keys':list(set(a)-set(b))[:5],'native_only_keys':list(set(b)-set(a))[:5],'mismatches':bad[:5]})
        total_saved+=len(wr)
    print('Checked pathway files',aid,flush=True)
check('Table22_all93_native_results',all(c['pass'] for c in enrichment_checks),{'analyses':len(I),'native_GSEA':total_gsea,'native_ORA':total_ora,'saved_significant_rows':total_saved,'compared_cells':detail_comparisons,'failed_blocks':[x for x in enrichment_checks if not x['pass']]})
save('enrichment_file_reconciliation.json',enrichment_checks)
save('source_comparisons.json',checks);save('source_mismatches.json',differences);save('key_result_reconciliation.json',contextual)
print(json.dumps({'checks':len(checks),'failed':[c for c in checks if not c['pass']]},ensure_ascii=False,indent=2))
