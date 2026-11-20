from common import *
import shutil
rec=verify();history=W/'receipts';history.mkdir(exist_ok=True);old=sha(W/'RECEIPT.json');shutil.copy2(W/'RECEIPT.json',history/(old+'.json'))
for root in [W/'code',W/'pseudobulk',W/'communication/code',W/'communication/inputs']:
 if not root.exists():continue
 for p in root.rglob('*'):
  if p.is_file() and '__pycache__' not in str(p) and '/results/' not in str(p):rec['bound_files'][str(p)]=sha(p)
for p in [R/'communication/resources/LIANA_consensus.tsv',R/'communication/resources/CellChatDB_human_frozen.rds']:
 rec['bound_files'][str(p)]=sha(p)
rec['previous_receipt_sha256']=old;rec['additive_seal_time']=time.time();rec['note']='Earlier bound contents verified unchanged; added authorized code and exact new inputs; no old checksum replaced'
dump(W/'RECEIPT.json',rec);print('PASS additive receipt',sha(W/'RECEIPT.json'))
