from common import *
import pandas as pd,numpy as np,shutil
out=W/'communication_disease';out.mkdir(exist_ok=False)
m=pd.read_csv(AUDIT/'sample_sidecar_updated.tsv',sep='\t',keep_default_na=False)
eligible=pd.read_csv(W/'communication/inputs/selected_specimens.tsv',sep='\t').sample_uid
d=m[m.sample_uid.isin(eligible)].rename(columns={'sample_uid':'sample_unit','study':'study_canonical','tissue':'tissue_v2','cohort_disease_context':'disease_context'}).copy()
d['donor_verified']=True;d['condition_class']=d.group_binary.map({'healthy':'control','disease':'disease'})
assert d.condition_class.notna().all() and d.donor_uid.ne('').all()
keys=['study_canonical','tissue_v2','disease_context'];aud=[];kept=[]
for key,z in d.groupby(keys):
 n=z.groupby('condition_class').donor_uid.nunique();ok=n.get('control',0)>=2 and n.get('disease',0)>=2
 aud.append(dict(zip(keys,key),n_disease=int(n.get('disease',0)),n_control=int(n.get('control',0)),eligible=ok,repeated_donor=z.donor_uid.duplicated().any()))
 if ok:kept.extend(z.sample_unit)
d['cohort_eligible']=d.sample_unit.isin(kept);d.to_csv(out/'design_all229.tsv',sep='\t',index=False);d[d.cohort_eligible].to_csv(out/'design_eligible.tsv',sep='\t',index=False);pd.DataFrame(aud).to_csv(out/'cohort_readiness.tsv',sep='\t',index=False)
# Exact native inherited CR2 program; dependency isolated in this version directory.
shutil.copy2(C/'external_legacy/communication/g08_cr2_pattern_reference.R',out/'cr2_native.R')
rng=np.random.default_rng(731);rows=[]
for i in range(8):
 base=rng.normal()
 for g in ['control','disease']:rows.append(dict(pattern_id='paired_donor_fixture',sample_order=len(rows),condition_class=g,resolved_person_key='D'+str(i),validation_y=base+(g=='disease')*.6+rng.normal(0,.15)))
pd.DataFrame(rows).to_csv(out/'cr2_test_input.tsv',sep='\t',index=False)
e=os.environ.copy();e['G08_R_LIBRARY']=str(W/'Rlib')
subprocess.run(['/usr/bin/Rscript',out/'cr2_native.R',out/'cr2_test_input.tsv',out/'cr2_test_output.tsv'],env=e,check=True)
z=pd.read_csv(out/'cr2_test_output.tsv',sep='\t');assert z.cr2_reference_status.eq('passed').all() and z.reference_se.gt(0).all() and z.reference_df_satterthwaite.gt(0).all()
dump(out/'NATIVE_CR2_TEST.json',dict(status='PASS',repeated_people=8,records=16,native_script_sha256=sha(out/'cr2_native.R'),R_library=str(W/'Rlib')))
dump(out/'DESIGN.json',dict(scope='Within each disease context: disease vs its designated reference; no invalid direct comparison across vascular beds',min_independent_donors_per_arm=2,min_cohorts_meta=2,method='Inherited difference of scores; repeated donor CR2/Satterthwaite; context DL meta only across independent cohorts',missing_interactions='missing/unmeasured is NA, never fabricated zero',unpaired_cohort_SE='inherited sample variance',BH_scope='method x disease_context x threshold_scheme',full_L2_communication='NOT_IN_ORIGINAL_PLAN',status='WAITING_FOR_VALIDATED_229_SAMPLE_L1_COMMUNICATION'))
print('Communication disease cohorts',sum(x['eligible'] for x in aud),'of',len(aud),'CR2 native test passed')
