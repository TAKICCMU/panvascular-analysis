import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[key] = '1'
from pathlib import Path
import csv, json, hashlib, html, collections, sys, re
import numpy as np
import pandas as pd

R = Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
D = Path('/data/anzhen_vhd/scRNAseq/codex/panvascular_downstream_storage_v1')
BASE = D/'covarnet_pseudobulk_launch_review_20260915'
W = Path(__file__).resolve().parents[1]
PRIOR = BASE/'covarnet_group_review_20260915T155047Z'
OLD = BASE/'full_review_v1'
COMM = R/'production_launches/communication_L1_20260914T015040Z'
PHASE = Path(json.loads((R/'pseudobulk/descriptive_final_atlas_v1/config.json').read_text())['cohort_output']).parent.parent

def read(p):
    return pd.read_csv(p, sep='\t', keep_default_na=False)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()
def save(name,df):
    df.to_csv(W/name,sep='\t',index=False)

def main():
    m=read(PRIOR/'covarnet_sample_sidecar.tsv')
    assert m.sample_uid.is_unique and len(m)==247
    m['group_review_required']=m.group_review_required.astype(str).str.lower().eq('true')
    for gsm in ['GSM7048232','GSM7048234']:
        mask=m.sample_uid.eq('GSE224587__'+gsm)
        assert mask.sum()==1
        m.loc[mask,'covarnet_group']='健康组（同供者近端参照）'
        m.loc[mask,'group_evidence']='用户明确确认近端为健康组；原始AAA供者诊断与配对关系保留'
        m.loc[mask,'group_review_required']=False
    m['group_binary']=np.where(m.covarnet_group.str.startswith('健康'),'healthy','disease')
    m['covarnet_include']=m.covarnet_include.astype(str).str.lower().eq('true')
    assert m.covarnet_include.sum()==230
    a=read(OLD/'all247_sample_aggregation_audit_source_updated.tsv')
    meta=pd.read_csv(R/'L2_targeted_followup_metadata_v4_20260913_154850_913989/inputs/metadata_h5ad_patch_v4.csv',keep_default_na=False)
    assert meta.sample_uid.is_unique and set(meta.sample_uid)==set(m.sample_uid)
    m=m.merge(meta[['sample_uid','sample_name','age','sex','assay_type','vascular_bed','cohort_disease_context','donor_disease_context','tissue_disease_status']],on='sample_uid',validate='one_to_one')
    resolved=a.set_index('sample_uid').donor_resolved.astype(str).str.lower().eq('true')
    m['donor_resolved']=m.sample_uid.map(resolved)
    assert not m.donor_resolved.isna().any()
    assert not (m.donor_resolved & m.donor_uid.astype(str).str.lower().isin(['','nr','na','nan','unknown','none'])).any()
    m['donor_block']=np.where(m.donor_resolved,m.donor_uid,'UNRESOLVED_SAMPLE::'+m.sample_uid)
    m['analysis_profile_uid']=m.sample_uid  # user requires separate specimens/captures
    save('sample_sidecar_updated.tsv',m)
    # Review group meaning without changing user-approved group assignments.
    findings=[]
    for x in m.itertuples():
        if not x.covarnet_include: continue
        dx=x.source_diagnosis.lower(); flags=[]
        if x.group_binary=='healthy':
            if x.study=='GSE255895':
                flags.append(('HEALTHY_TISSUE_NOT_HEALTHY_DONOR','无斑块胸廓内动脉，但供者均为冠心病搭桥患者；不能称健康人'))
            elif x.study=='GSE248610':
                flags.append(('DISEASE_SPECIFIC_REFERENCE_NOT_GENERAL_HEALTH','糖尿病大血管病的未钙化段；并非无病组织已获证实'))
            elif 'non-aneurysmal diseased control' in dx:
                flags.append(('DISEASE_SPECIFIC_REFERENCE_NOT_GENERAL_HEALTH','非动脉瘤对照但供者患病；不能自动等同一般健康人'))
            elif any(t in dx for t in ['proximal adjacent','abdominal aortic aneurysm']):
                flags.append(('PAIRED_REGION_REFERENCE','病人自身邻近/近端对照；保留健康分析组，不解释为健康供者'))
            elif 'transplant recipient' in dx or 'heart-transplant control' in dx:
                flags.append(('REFERENCE_DONOR_COMORBIDITY','参照组织来自移植/心脏疾病背景供者；组织状态与供者健康分开记录'))
        if not x.donor_resolved:
            flags.append(('MISSING_DONOR','供者未提供；用独立样本占位，不把多个NR当同一人，不虚报供者人数'))
        if x.sample_uid.startswith('GSE131778__') and 'DONOR_USER_CONFIRMED' in x.donor_uid:
            flags.append(('USER_CONFIRMED_DONOR_OVERRIDE','按用户确认独立供者；原年龄/性别来源仍保留，不称为新完成的供者临床核验'))
        for code,reason in flags:
            findings.append(dict(sample_uid=x.sample_uid,study=x.study,donor_uid=x.donor_uid,group=x.group_binary,source_diagnosis=x.source_diagnosis,flag=code,explanation=reason,automatically_changed=False))
    save('grouping_findings.tsv',pd.DataFrame(findings))
    donor_checks=[]
    for donor,z in m[m.donor_resolved].groupby('donor_uid'):
        for col in ['sex','age']:
            vals=sorted(set(str(v) for v in z[col] if str(v).lower() not in ['','nr','na','nan','unknown']))
            if len(vals)>1: donor_checks.append(dict(donor_uid=donor,field=col,values=';'.join(vals),samples=';'.join(z.sample_uid)))
    save('donor_metadata_conflicts.tsv',pd.DataFrame(donor_checks,columns=['donor_uid','field','values','samples']))
    # Counts are the exact published all-cell sample counts, not old biological aggregates.
    counts=read(PHASE/'counts/sample_state_counts.tsv')
    keys=['level1_annotation','level2_annotation']
    assert not counts.duplicated(['sample_uid']+keys).any()
    assert set(counts.sample_uid)==set(m.sample_uid)
    assert (counts['count']>=0).all() and (counts['count']%1==0).all()
    totals=counts.groupby('sample_uid')['count'].sum()
    assert totals.reindex(m.sample_uid).to_numpy().tolist()==m.n_cells.astype(int).tolist()
    # Exactly reproduce the historical eligibility pathway and completed coverage.
    units=read(COMM/'inputs/reviewed_units.tsv');blocked=read(COMM/'inputs/blocked_units.tsv')
    plan=read(COMM/'inputs/specimen_plan.tsv');selected=read(COMM/'inputs/selected_specimens.tsv')
    good=units[units.aggregate_approved.astype(str).str.lower().eq('true') & units.biological_sample_resolved.astype(str).str.lower().eq('true')]
    assert len(good)==206 and good.biological_sample_uid.nunique()==206
    assert len(blocked)==41 and len(selected)==203
    assert set(selected.biological_sample_uid)==set(plan.loc[plan.eligible.astype(str).str.lower().eq('true'),'biological_sample_uid'])
    summary=read(D/'communication_data_v1/parallel4_r2/completion_review_v1/communication_203_specimens_summary.tsv')
    assert len(summary)==203
    l1=counts.groupby(['sample_uid','level1_annotation'],as_index=False)['count'].sum()
    l1=l1[~l1.level1_annotation.isin(['Unknown','Not_assessed','Pending'])]
    ntypes=l1[l1['count']>=50].groupby('sample_uid').size()
    old=units.set_index('sample_uid');done=set(selected.biological_sample_uid)
    comm=[]
    for x in m.itertuples():
        ou=old.loc[x.sample_uid];finished=ou.biological_sample_uid in done
        n=int(ntypes.get(x.sample_uid,0));eligible=n>=2
        comm.append(dict(sample_uid=x.sample_uid,study=x.study,old_biological_sample_uid=ou.biological_sample_uid,new_biological_sample_uid=x.biological_sample_uid,new_donor_uid=x.donor_uid,old_completed=finished,old_blocked=x.sample_uid in set(blocked.sample_uid),known_L1_at_least_50_cells=n,new_separate_sample_eligible=eligible,CoVarNet_include=x.covarnet_include,communication_scope_not_automatically_changed=True,action='REUSE_EXISTING_MEMBERSHIP' if finished else 'ELIGIBLE_TO_ADD_AFTER_APPROVAL' if eligible else 'INSUFFICIENT_L1_TYPES_AT_50'))
    comm=pd.DataFrame(comm);save('communication_old_vs_new_247.tsv',comm)
    save('communication_old_3_ineligible.tsv',plan[~plan.eligible.astype(str).str.lower().eq('true')])
    save('communication_newly_eligible.tsv',comm[~comm.old_completed&comm.new_separate_sample_eligible])
    save('communication_still_ineligible.tsv',comm[~comm.new_separate_sample_eligible])
    # Simple user-facing table, provenance is available in sidecar.
    tab=pd.DataFrame({'研究':m.study,'样本':m.sample_uid,'供者':m.donor_uid,'分组':m.covarnet_group,'CoVarNet纳入':np.where(m.covarnet_include,'是','否')})
    tab.to_csv(W/'研究_样本_供者_分组.csv',index=False,encoding='utf-8-sig')
    (W/'样本表.html').write_text('<meta charset="utf-8"><style>body{font-family:system-ui}td,th{padding:6px;border:1px solid #ddd}table{border-collapse:collapse}</style>'+tab.to_html(index=False,escape=True))
    conf={'status':'AUDIT_ONLY_NOT_APPROVED_TO_RUN','sample_registry':str(W/'sample_sidecar_updated.tsv'),'sample_registry_sha256':sha(W/'sample_sidecar_updated.tsv'),'profile_unit':'sample_uid','donor_block':'donor_uid; unresolved samples have individual placeholders, not inferred donors','all_three_GSE224587_proximal_samples':'healthy; same donor as matching center','exclude_study':['GSE220774'],'exclude_samples':['GSE224559__GSM7025944','GSE224559__GSM7025945'],'scope':'CoVarNet only; no historical analysis mutation','CoVarNet_run':False}
    (W/'config.json').write_text(json.dumps(conf,indent=2,ensure_ascii=False))
    result={'samples':247,'included_samples':int(m.covarnet_include.sum()),'included_cells':int(m.loc[m.covarnet_include,'n_cells'].astype(int).sum()),'group_counts':m[m.covarnet_include].group_binary.value_counts().to_dict(),'missing_donor_samples':int((m.covarnet_include&~m.donor_resolved).sum()),'missing_donor_cells':int(m.loc[m.covarnet_include&~m.donor_resolved,'n_cells'].sum()),'donor_metadata_conflicts':len(donor_checks),'communication':{'old_total_sample_records':247,'old_blocked_records':41,'old_safe_specimens':206,'old_eligible_and_completed':203,'old_below_type_threshold':3,'new_all_sample_eligible':int(comm.new_separate_sample_eligible.sum()),'new_additions':int((~comm.old_completed&comm.new_separate_sample_eligible).sum()),'new_ineligible':int((~comm.new_separate_sample_eligible).sum()),'covarnet_overlap_eligible':int((comm.new_separate_sample_eligible&comm.CoVarNet_include).sum()),'rerun_started':False},'models_run':False}
    (W/'SCOPE_AUDIT.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    sources=[PRIOR/'covarnet_sample_sidecar.tsv',PHASE/'counts/sample_state_counts.tsv',PHASE/'counts/ontology.tsv',COMM/'code/prepare.py',COMM/'inputs/reviewed_units.tsv',COMM/'inputs/specimen_plan.tsv',COMM/'inputs/selected_specimens.tsv',D/'communication_data_v1/parallel4_r2/completion_review_v1/VALIDATION.json']
    (W/'SOURCE_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in sources},indent=2))
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=='__main__': main()
