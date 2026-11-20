"""Read-only extraction of Table 02 for a clinical-metadata review graphic.

Run with Python and openpyxl. Deduplicate only by the recorded Donor ID;
do not infer identity, fill missing demographics, or edit the workbook.
"""
from pathlib import Path
import collections,csv,hashlib,json,math
import openpyxl

O=Path(__file__).parent;P=O.parents[2];R=O/'review'
source=P/'Supplementary_Tables.xlsx'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source_hash=sha(source)
w=openpyxl.load_workbook(source,read_only=True,data_only=True)
rows=list(w['Table 02'].values)
records=[dict(zip(rows[2],r)) for r in rows[3:] if r[0] is not None]
w.close()
assert len(records)==247 and len({r['Sample ID'] for r in records})==247
assert sum(r['Cells (n)'] for r in records)==1175115
groups=collections.defaultdict(list)
for row in records:
    key=row['Donor ID']
    assert key is not None and str(key).strip().lower() not in {'','na','nr','none','unknown','nan'}
    assert row['Donor resolved'] in {True,'True','true'}
    assert row['Sex'] in {'male','female',None}
    assert row['Age'] is None or (isinstance(row['Age'],(int,float)) and math.isfinite(row['Age']) and 0<=row['Age']<=110)
    groups[key].append(row)
donors=[]
for donor,rs in sorted(groups.items()):
    age={r['Age'] for r in rs if r['Age'] is not None}
    sex={r['Sex'] for r in rs if r['Sex'] is not None}
    assert len(age)<=1 and len(sex)<=1,(donor,age,sex)
    donors.append({'donor_id':donor,'age_years':next(iter(age),None),'sex':next(iter(sex),None),
        'sample_record_count':len(rs),'sample_ids':' | '.join(r['Sample ID'] for r in rs),
        'study_ids':' | '.join(sorted({r['Study ID'] for r in rs})),
        'metadata_source':'Supplementary_Tables.xlsx / Table 02'})
assert len(donors)==172
def tsv(path,items,columns):
    with path.open('w',newline='') as fh:
        wr=csv.DictWriter(fh,fieldnames=columns,delimiter='\t',extrasaction='ignore')
        wr.writeheader();wr.writerows(items)
tsv(R/'Figure1e_clinical_sample_records.tsv',records,
    ['Sample ID','Study ID','Donor ID','Age','Sex','Cells (n)','Evidence source','Previous donor ID','Decision'])
tsv(R/'Figure1e_clinical_donor_records.tsv',donors,list(donors[0]))
counts=collections.Counter(r['sex'] or 'Not reported' for r in donors)
age_available=sum(r['age_years'] is not None for r in donors)
audit={'source':str(source),'sheet':'Table 02','source_sha256':source_hash,
 'sample_records':247,'recorded_donor_ids':len(donors),'study_ids':len({r['Study ID'] for r in records}),
 'duplicate_sample_ids':0,'invalid_donor_ids':0,'within_donor_age_sex_conflicts':0,
 'donors_with_multiple_sample_records':sum(len(v)>1 for v in groups.values()),
 'age_reported':age_available,'age_not_reported':len(donors)-age_available,'sex_counts':dict(counts),
 'unit':'Recorded Donor ID, one observation per ID; not weighted by cell counts or sample multiplicity.',
 'deduplication':'Take the unique available stored age and sex within each donor ID; fail on conflicts. No cross-donor filling.',
 'limitations':['The figure describes the current Table 02 metadata, not a new independent verification of donor identities or demographics.',
 'Seven GSE131778 sample records have user-confirmed separate donor IDs; their stored age/sex values are inherited in Table 02. This figure does not independently verify those inherited values.',
 'No disease or healthy-donor categories are inferred from analysis-group labels.'],
 'source_unchanged':sha(source)==source_hash}
assert audit['source_unchanged']
(R/'Figure1e_clinical_metadata_check.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False))
print(json.dumps({k:v for k,v in audit.items() if k not in ['source','source_sha256','limitations']},indent=2))
