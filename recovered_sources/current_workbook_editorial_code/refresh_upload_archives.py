from pathlib import Path
import csv,json,hashlib,zipfile,gzip,os
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');R=P/'Public_Release';U=P/'上传文件';U.mkdir(exist_ok=True);A=P/'投稿核查'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
# Explicit archive scopes; never include a full h5ad or a whole internal working directory.
jobs={
'Analysis_Code.zip':[(R/'Analysis_Code','Analysis_Code')],
'Figure_Code.zip':[(R/'Figure_Code','Figure_Code')],
'Atlas_Annotations_Metadata_Coordinates.zip':[(R/'Atlas_Metadata','Atlas_Metadata')],
'Derived_SharedFactor_GWAS.zip':[(R/'Derived_GWAS','Derived_GWAS')],
'Source_Data_Figures.zip':[(P/'Source_Data_Figures','Source_Data_Figures'),(R/'Supporting_Visuals','Supporting_Visuals')],
'Supplementary_Data.zip':[(P/'Supplementary_Tables.xlsx','Supplementary_Tables.xlsx'),(P/'Supplementary_Table20_full_results','Supplementary_Table20_full_results'),(R/'Supplementary_Tables_TSV','Supplementary_Tables_TSV'),(P/'Supplementary_Information.docx','Supplementary_Information.docx')]}
jobs={k:v for k,v in jobs.items() if k in ['Supplementary_Data.zip','Source_Data_Figures.zip']}
rows=[r for r in csv.DictReader((P/'UPLOAD_MANIFEST.tsv').open(),delimiter='\t') if Path(r['file']).name not in jobs]
for name,roots in jobs.items():
 entries=[]
 for root,label in roots:
  files=[root] if root.is_file() else sorted(root.rglob('*'))
  for f in files:
   if not f.is_file() or f.name=='.DS_Store' or '__pycache__' in f.parts:continue
   assert f.suffix.lower() not in ['.h5ad','.h5','.hdf5','.pem','.pth','.pt']
   arc=label if root.is_file() else str(Path(label)/f.relative_to(root));entries.append((f,arc))
 manifest=[{'file':arc,'bytes':f.stat().st_size,'sha256':sha(f)} for f,arc in entries]
 dst=U/name
 with zipfile.ZipFile(dst,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True,strict_timestamps=False) as z:
  for f,arc in entries:z.write(f,arc)
  z.writestr('RELEASE_CONTENT_MANIFEST.json',json.dumps(manifest,indent=2,ensure_ascii=False))
 with zipfile.ZipFile(dst) as z:
  assert z.testzip() is None
 rows.append({'file':str(dst.relative_to(P)),'files':len(entries),'bytes':dst.stat().st_size,'sha256':sha(dst),'status':'LOCAL_ARCHIVE_CRC_VERIFIED_NOT_UPLOADED'})
 print(name,rows[-1]['bytes'],flush=True)
with (P/'UPLOAD_MANIFEST.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(A/'archive_validation.json').write_text(json.dumps({'archives':rows,'full_atlas_included':False,'remote_upload_performed':False},ensure_ascii=False,indent=2))
# Native-file receipt integrity, without scientific re-execution.
manifest=json.loads((R/'Analysis_Code/native/NATIVE_SOURCE_MANIFEST.json').read_text());bad=[]
for item in manifest:
 f=R/'Analysis_Code/native'/item['path']
 if not f.exists() or sha(f)!=item['sha256']:bad.append(item['path'])
assert not bad,bad
(A/'native_code_copy_validation.json').write_text(json.dumps({'native_source_files':len(manifest),'missing_or_changed':bad,'test_scope':'byte-identical source copy; not scientific rerun'},indent=2))
