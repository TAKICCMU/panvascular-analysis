from pathlib import Path
import subprocess,tarfile,io,json
O=Path('/Users/taki/Desktop/投稿包_最新版_20260919/图件编辑资料/全流程代码审核_20260922');B=O/'Panvascular_code_review'
script=r'''import pathlib,sys,tarfile,hashlib,io,json
root=pathlib.Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908/scib')
rows=[]
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as t:
 for p in sorted(root.rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts:continue
  if p.suffix.lower() not in ['.py','.r','.sh','.json','.yaml','.yml','.md','.tsv','.txt','.log','.diff']:continue
  if p.stat().st_size>3000000:continue
  raw=p.read_bytes();rel=p.relative_to(root).as_posix();info=tarfile.TarInfo(rel);info.size=len(raw);t.addfile(info,io.BytesIO(raw));rows.append({'source_path':str(p),'relative_path':rel,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
 raw=json.dumps(rows,indent=2).encode();info=tarfile.TarInfo('RECOVERY_SOURCE_MANIFEST.json');info.size=len(raw);t.addfile(info,io.BytesIO(raw))
'''
r=subprocess.run(['ssh','-o','BatchMode=yes','codex-shpc-gpu','python3','-'],input=script.encode(),capture_output=True,check=True)
out=B/'recovered_sources/scib';out.mkdir(exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(r.stdout),mode='r:gz') as t:
 for m in t:
  assert m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts
  p=out/m.name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(t.extractfile(m).read())
print(len(json.loads((out/'RECOVERY_SOURCE_MANIFEST.json').read_text())),'files recovered')
