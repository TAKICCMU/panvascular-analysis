#!/usr/bin/env python3
"""Read only frozen proposal metadata; build NOT_APPROVED clinical decision tables.
No counts, h5ad, source annotations, approval config, or statistical model is opened.
"""
import argparse, hashlib, json
from pathlib import Path
import pandas as pd

MISSING={'','na','nr','nan','none','null','not_available','not available','unknown','<na>'}
def absent(v):return pd.isna(v) or str(v).strip().lower() in MISSING
def donor_set(g):return set(g.loc[~g.donor_uid.map(absent)&~g.donor_uid.astype(str).str.contains('PROVISIONAL|UNRESOLVED|PLACEHOLDER',case=False,regex=True),'donor_uid'])
def jl(values):return json.dumps(sorted(set(map(str,values))),ensure_ascii=False,separators=(',',':'))
def vals(g,col):return jl(g.loc[~g[col].map(absent),col])
def truth(s):return s.astype(str).str.lower().isin(['true','1'])
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

NOTES={
 'contrast_001':'同一right-coronary标本技术partition须先核查；8条sample不等于8名供者；当前无reference。',
 'contrast_002':'reference为non-aneurysmal diseased control，不是健康人；仅2名reference。',
 'contrast_004':'估计同供者plaque core vs proximal adjacent segment局部差异，不是患病者vs健康者；仅3完整对；年龄/性别主效应不可与donor固定效应并列识别。',
 'contrast_010':'reference是终末期心衰移植受者的非病变主动脉瓣，年龄10–33岁；不能与其他study老年CAVD直接拼成万能健康对照。',
 'contrast_011':'Type-A夹层与heart-transplant donor aorta比较；reference3人，age/sex构成待审；不是所有reference供者均已证明健康。',
 'contrast_012':'Fibrosa/Mixed/Ventricularis等重复组织区段属于同供者，不是技术重复或独立人；4 vs1供者；当前unpaired执行器会拒绝同donor多specimen。',
 'contrast_013':'急性/亚急性/慢性Type-A夹层合并case是异质疾病阶段；年龄当前全缺失，case仅male，reference含female；9 vs5仅state过滤前上界。',
 'contrast_017':'source写healthy donor pulmonary artery，仅按来源含义保留；reference仅3人，年龄/性别构成仍须审。',
 'contrast_018':'肺纤维化相关PH不同于PAH，不合并成同一case；reference仅2人。',
 'contrast_019':'同biological_sample存在2个capture/partition，均未批准合并；无reference。',
 'contrast_021':'study为GSE248610而sample/source是GSE248609；保留真实别名，不重命名。问题为糖尿病患者calcified vs non-calcified段，非diabetes vs healthy；仅2完整对，另2单臂供者不能当已配对。',
 'contrast_024':'internal thoracic artery来自CAD/CABG患者，是reference组织而非健康供者；全snRNA-seq；与其他study病变冠脉scRNA简单合并会混淆study/血管床/assay。',
 'contrast_037':'当前分组均reference，young/adult是潜在另一科学问题；不得本轮重编码成疾病case/reference。',
 'contrast_038':'CAVD vs non-CAVD valve仅4 vs2供者；不能自动加入其他context的儿童reference来凑样本。'
}

def main(root):
 root=Path(root).resolve();inputs=[root/'configs/cohort_design.tsv',root/'configs/contrasts.tsv',root/'configs/pooled_context_contrasts_PROPOSED.tsv',root/'CURRENT_CONTRAST_FEASIBILITY.tsv']
 before={str(p.relative_to(root)):sha(p) for p in inputs}
 m=pd.read_csv(inputs[0],sep='\t',keep_default_na=False);specs=pd.read_csv(inputs[1],sep='\t',keep_default_na=False)
 assert len(m)==247 and m.sample_uid.is_unique and len(specs)==38 and specs.contrast_id.is_unique
 assert not truth(specs.approved).any(), 'Never promote approval'
 rows=[];matches={};members=[]
 for r in specs.itertuples(index=False):
  g=m.copy()
  for col in ['study','vascular_bed','cohort_disease_context','assay_type']:g=g[g[col].astype(str).eq(str(getattr(r,col)))]
  g=g[g.comparison_group.isin([r.case,r.reference])].copy();matches[r.contrast_id]=g
  ca=g[g.comparison_group.eq(r.case)];re=g[g.comparison_group.eq(r.reference)];cd=donor_set(ca);rd=donor_set(re);pairs=cd&rd
  paired=str(r.paired).lower() in ['true','1'];nmin=int(r.min_independent_per_arm);both=bool(cd and rd)
  good=len(cd)>=nmin and len(rd)>=nmin and (not paired or len(pairs)>=nmin)
  grade='PRIMARY_CANDIDATE_PENDING_STATE_AND_REVIEW' if good else ('EXPLORATORY_DESIGN_CANDIDATE_NOT_EXECUTABLE_AT_CURRENT_THRESHOLD' if both else 'NOT_ESTIMABLE_AS_CURRENT_CASE_REFERENCE_CONTRAST')
  notes=[NOTES.get(r.contrast_id,'病例/参照仅按现有cohort原字段保留，最终release须再核对来源和同供者结构。')]
  if r.study.startswith('HCA-'):notes.append('多个动脉床可来自同一供者；当前reference-only，不把9个血管床当9个独立疾病context。')
  if r.study=='GSE198750':notes.append('当前study名称来自metadata追溯，不能据名称推断重新加入了历史重复counts；保留atlas_source_id/sample身份供审核。')
  if not both:notes.append('缺同study/同bed/同assay一臂，不从其他study移入参照消除结构性缺失；仍可用于图谱描述。')
  elif not good:notes.append('探索建议不是执行豁免；原每臂5/配对5对门槛不变。')
  donor_arm_bios=g.groupby(['donor_uid','comparison_group']).biological_sample_uid.nunique()
  repeated=donor_arm_bios[donor_arm_bios.gt(1)]
  if len(repeated):notes.append('同donor/arm多biological_sample，当前simple paired/unpaired模型不可直接当独立单位；需明确真实标本合并或另行实现重复测量模型。')
  bio=g.groupby('biological_sample_uid').agg(n_partitions=('sample_uid','size'),all_aggregate_flags=('aggregate_approved',lambda s:truth(s).all()),n_inference_flags=('inference_eligible',lambda s:int(truth(s).sum())))
  multi=bio[bio.n_partitions.gt(1)]
  covariate_proposal='~ donor_uid + group（仅已审核完整配对；donor固定效应吸收不随配对变化的年龄/性别）' if paired else '~ group；可另审~ sex + group或~ age + sex + group，须完整性/秩/自由度通过，不静默去协变量'
  meta_suitability='NOT_SUITABLE_CURRENT_SINGLE_ARM' if not both else 'WITHIN_CONTEXT_MODEL_CANDIDATE_REQUIRES_REFERENCE_AND_DONOR_REVIEW'
  joint_id='';joint_components=[]
  if r.contrast_id in ['contrast_011','contrast_013']:
   joint_id='pooled_context_000_EXISTING_UNAPPROVED';joint_components=['contrast_011','contrast_013'];meta_suitability='CONDITIONAL_JOINT_ASCENDING_AORTA_AD_NOT_APPROVED'
  if r.contrast_id in ['contrast_015','contrast_023']:
   notes.append('当前vascular_bed=aorta，不静默等同ascending aorta纳入既有联合AD提案。')
  row=dict(contrast_id=r.contrast_id,comparison_question=f'{r.cohort_disease_context}: {r.study}内{r.vascular_bed}的当前case vs reference',study=r.study,vascular_bed=r.vascular_bed,assay_type=r.assay_type,meta_context=r.meta_context,
    case_source_diagnosis=vals(ca,'source_diagnosis'),reference_source_diagnosis=vals(re,'source_diagnosis'),case_donor_disease_context=vals(ca,'donor_disease_context'),reference_donor_disease_context=vals(re,'donor_disease_context'),case_tissue_disease_status=vals(ca,'tissue_disease_status'),reference_tissue_disease_status=vals(re,'tissue_disease_status'),reference_is_healthy_person_assumed=False,
    case_sample_ids=jl(ca.sample_uid),reference_sample_ids=jl(re.sample_uid),case_source_series=vals(ca,'source_series_accession'),reference_source_series=vals(re,'source_series_accession'),atlas_source_ids=vals(g,'atlas_source_id'),case_tissue_details=vals(ca,'tissue_detail'),reference_tissue_details=vals(re,'tissue_detail'),
    n_case_sample_records=len(ca),n_reference_sample_records=len(re),n_case_biological_samples=ca.biological_sample_uid.nunique(),n_reference_biological_samples=re.biological_sample_uid.nunique(),n_case_donors_upper=len(cd),n_reference_donors_upper=len(rd),complete_pair_donors=jl(pairs),n_complete_pairs_upper=len(pairs),case_unpaired_donors=jl(cd-rd),reference_unpaired_donors=jl(rd-cd),paired=paired,
    n_repeated_donor_arm_specimen_groups=len(repeated),n_technical_multi_partition_biological_units=len(multi),n_multi_partition_units_all_aggregate_flags_true=int(multi.all_aggregate_flags.sum()),n_current_aggregate_flag_usable_biological_units=int(bio.all_aggregate_flags.sum()),n_current_inference_eligible_sample_records=int(truth(g.inference_eligible).sum()),technical_flag_is_not_model_approval=True,
    both_arms_nonempty=both,minimum_independent_per_arm=nmin,minimum_complete_pairs=nmin if paired else 'NOT_APPLICABLE',passes_donor_upper_bound=good,proposed_grade=grade,approved=False,de_design_approved=False,meta_se_approved=False,
    source_formula=r.formula,source_covariates=r.covariates,proposed_formula_options=covariate_proposal,case_age_known_samples=int((~ca.age.map(absent)).sum()),reference_age_known_samples=int((~re.age.map(absent)).sum()),case_sex_values=vals(ca,'sex'),reference_sex_values=vals(re,'sex'),
    within_contrast_study_bed_assay_confounding='Both arms share study/bed/assay strata; finer wall/layer,stage,age,sex remain unassessed' if both else 'One arm absent; cross-study/bed/assay borrowed controls not identifiable by this within-study design',
    per_state_estimability_checks='approved complete biological units; >=20 cells/state/unit; >=5 donors/arm or >=5 complete pairs; prespecified complete-case covariates; no duplicate donor/arm; full rank; residual df>=2; gene filtering',
    joint_context_suitability=meta_suitability,joint_context_proposal_id=joint_id,joint_component_contrasts=jl(joint_components),joint_vs_study_meta_rule='If a joint context model is approved, only its single effect enters the chosen meta run; never also its component-study effects; one disease remains one context; shared donor across contexts blocks',
    clinical_interpretation_and_risks=' '.join(notes),input_scope='FROZEN_PROPOSAL_METADATA_ONLY_NOT_FINAL_RELEASE_RECONCILIATION',status='PROPOSED_NOT_APPROVED_NOT_EXECUTED')
  rows.append(row)
  for _,x in g.iterrows():members.append({'contrast_id':r.contrast_id,**{k:x[k] for k in ['sample_uid','source_sample_accession','study','atlas_source_id','biological_sample_uid','donor_uid','comparison_group','vascular_bed','assay_type','source_diagnosis','tissue_detail','age','sex','aggregate_approved','inference_eligible']}})
 out=pd.DataFrame(rows);assert len(out)==38 and out.both_arms_nonempty.sum()==11 and out.passes_donor_upper_bound.sum()==1
 old=pd.read_csv(inputs[3],sep='\t').set_index('contrast_id');new=out.set_index('contrast_id')
 for ours,theirs in [('n_case_donors_upper','n_case_donors_recomputed'),('n_reference_donors_upper','n_reference_donors_recomputed'),('n_complete_pairs_upper','n_complete_pairs_recomputed')]:assert new[ours].equals(old[theirs])
 out.to_csv(root/'CLINICAL_CONTRAST_DECISIONS_PROPOSED.tsv',sep='\t',index=False)
 pd.DataFrame(members).to_csv(root/'CLINICAL_CONTRAST_SAMPLE_MEMBERSHIP_PROPOSED.tsv',sep='\t',index=False)
 combined=pd.concat([matches['contrast_011'],matches['contrast_013']]);jc=donor_set(combined[combined.comparison_group.eq('case')]);jr=donor_set(combined[combined.comparison_group.eq('reference')])
 joint=[dict(proposal_id='pooled_context_000_EXISTING_UNAPPROVED',meta_context='aortic dissection',component_contrasts=jl(['contrast_011','contrast_013']),study_allowlist=jl(['GSE213740','GSE222318']),vascular_bed='ascending aorta',assay_type='scRNA-seq',source_formula='~ study + group',n_case_donors_upper=len(jc),n_reference_donors_upper=len(jr),n_overlap_donor_keys=len(jc&jr),existing_proposal_min_per_arm_per_study=2,current_primary_per_arm_threshold=5,threshold_decision='Existing unapproved pooled proposal has per-study=2; not adopted/changed by this revision. With unchanged per-study5, only GSE222318 remains; no new independent context.',reference_review='Compare actual transplant reference diagnoses, age missingness, sex and AD stage across studies; case-only/control-only pooling prohibited',alternative_estimand='One study-adjusted within-context effect only if clinical comparability and explicit joint support rule separately approved',component_effect_meta_policy='Do not combine pooled effect and contrast_011/013 effects as independent meta inputs',cross_context_gain=0,approved=False,source_proposal_preserved=True,status='NOT_APPROVED_NOT_EXECUTED')]
 pd.DataFrame(joint).to_csv(root/'JOINT_CONTEXT_DECISIONS_PROPOSED.tsv',sep='\t',index=False)
 assert before=={str(p.relative_to(root)):sha(p) for p in inputs}
 result={'status':'READ_ONLY_METADATA_RECOMPUTATION_PASS','n_contrasts':38,'n_both_arms_nonempty':11,'n_meeting_donor_upper_bound':1,'only_donor_upper_bound_candidate':out.loc[out.passes_donor_upper_bound,'contrast_id'].tolist(),'n_approved_contrasts':0,'n_current_inference_eligible_sample_records':int(truth(m.inference_eligible).sum()),'metadata_source_files_unchanged':before,'script_sha256':sha(__file__),'clinical_approval_changed':False,'models_run':False,'real_counts_read':False,'source_boundary':'current proposed cohort table; not independently reauthenticated final-release metadata'}
 (root/'CLINICAL_DECISIONS_RECOMPUTATION.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
 print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);main(p.parse_args().root)
