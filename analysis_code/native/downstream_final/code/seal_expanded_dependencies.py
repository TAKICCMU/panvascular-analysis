from common import *
import shutil
rec=verify();old=sha(W/'RECEIPT.json');(W/'receipts').mkdir(exist_ok=True);shutil.copy2(W/'RECEIPT.json',W/'receipts'/(old+'.json'))
roots=[W/'code',W/'communication/code',W/'communication/inputs',W/'communication_disease',W/'Rlib',C/'src/external_statistics.py',C/'scripts/external_communication_meta.py',C/'external_legacy/communication/g08_cr2_pattern_reference.R',C/'external_legacy/communication/disease_context_meta.py']
for root in roots:
 if not root.exists():continue
 files=[root] if root.is_file() else list(root.rglob('*'))
 for p in files:
  if p.is_file() and '__pycache__' not in str(p):rec['bound_files'][str(p)]=sha(p)
rec['previous_receipt_sha256']=old;rec['additive_seal_time']=time.time();dump(W/'RECEIPT.json',rec)
print('Extended verified receipt',sha(W/'RECEIPT.json'))
