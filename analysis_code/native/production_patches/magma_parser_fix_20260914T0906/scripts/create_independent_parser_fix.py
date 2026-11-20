from pathlib import Path
import shutil,difflib,ast,json,hashlib
from pypdf import PdfReader
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908');OLD=R/'production_launches/magma_category_fix_20260914T0301';W=R/'production_launches/magma_parser_fix_20260914T0906';B=OLD/'magma_adapter_v1';C=W/'magma_adapter_v2'
for d in ['scripts','evidence','configs','approvals','logs','outputs']:(W/d).mkdir(parents=True,exist_ok=True)
for d in ['src','scripts','external_legacy']:shutil.copytree(B/d,C/d,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
shutil.copyfile('/tmp/magma_recovery_20260914.py',C/'src/magma_recovery.py')
f=C/'scripts/external_run_module.py';old=f.read_text();start=old.index("        if module=='magma':");end=old.index("        elif module=='scdrs':",start);new=old[:start]+"        if module=='magma':\n            from magma_recovery import run_recovery\n            run_recovery(a,ac,cfg,ext,out,call)\n"+old[end:];ast.parse(new);f.write_text(new)
(W/'evidence/EXACT_RUNNER_DIFF.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(B/'scripts/external_run_module.py'),tofile=str(f))))
manual=Path('/home/anzhen_vhd/scRNAseq/manual_v1.10.pdf');r=PdfReader(manual);(W/'evidence/MAGMA_MANUAL_v1.10_P21_P24.txt').write_text('\n\n'.join('PAGE '+str(i+1)+'\n'+r.pages[i].extract_text() for i in [20,23]));(W/'evidence/MANUAL_SOURCE.json').write_text(json.dumps({'path':str(manual),'sha256':hashlib.sha256(manual.read_bytes()).hexdigest(),'pages':[21,24],'interpretation':'condition variables have one result per MODEL; MODEL identifies a regression, Average is not a target enrichment'},indent=2));print(str(W))
