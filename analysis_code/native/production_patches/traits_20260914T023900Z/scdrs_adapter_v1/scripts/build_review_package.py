"""Package code/review evidence only, <=35MiB archive; no final data/models."""
from pathlib import Path
import hashlib,json,zipfile,datetime
root=Path(__file__).resolve().parents[1]
exclude={'__pycache__','.pytest_cache','runs','work'}
package_name=root.name
files=[]
for p in root.rglob('*'):
 if not p.is_file() or any(x in exclude for x in p.relative_to(root).parts):continue
 if p.name in ['PACKAGE_MANIFEST.json','SHA256SUMS.txt'] or p.suffix in ['.zip','.partial'] or p.name.endswith('.zip.sha256'):continue
 if p.suffix in ['.h5ad','.h5','.loom','.pt','.ckpt']:continue
 if p.stat().st_size>10*1024*1024:raise RuntimeError(f'Review file too large: {p}')
 files.append(p)
records=[{'path':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(files)]
(root/'PACKAGE_MANIFEST.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'PAUSED_FOR_DOWNSTREAM_CODE_REVISION_REVIEW','no_full_scientific_analysis':True,'files':records},indent=2,ensure_ascii=False))
(root/'SHA256SUMS.txt').write_text(''.join(r['sha256']+'  '+r['path']+'\n' for r in records))
z=root/(package_name+'_review.zip')
with zipfile.ZipFile(z.with_suffix('.partial'),mode='w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as f:
 for p in files+[root/'PACKAGE_MANIFEST.json',root/'SHA256SUMS.txt']:f.write(p,package_name+'/'+str(p.relative_to(root)))
z.with_suffix('.partial').replace(z)
if z.stat().st_size>35*1024*1024:raise RuntimeError('Split review modules; archive exceeds35MiB')
h=hashlib.sha256(z.read_bytes()).hexdigest();z.with_suffix(z.suffix+'.sha256').write_text(h+'  '+z.name+'\n')
with zipfile.ZipFile(z) as f:
 for r in records:
  if hashlib.sha256(f.read(package_name+'/'+r['path'])).hexdigest()!=r['sha256']:raise RuntimeError('ZIP readback mismatch')
print(json.dumps({'files':len(records),'zip_bytes':z.stat().st_size,'sha256':h,'status':'ZIP_READBACK_PASS'}))
