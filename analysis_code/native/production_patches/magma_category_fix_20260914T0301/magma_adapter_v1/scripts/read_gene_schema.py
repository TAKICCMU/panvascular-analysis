"""Read only var metadata from current working atlas for exact-ID code preparation."""
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.core import read_frame,atomic_json
root=Path(__file__).resolve().parents[1];probe=json.load(open(root/'evidence/READ_ONLY_INPUT_CHECK.json'));p=probe['source']
v=read_frame(p,'var');v.to_csv(root/'evidence/working_gene_metadata.tsv.gz',sep='\t')
atomic_json(root/'evidence/GENE_SCHEMA.json',{'path':p,'role':'working var schema only','n_genes':len(v),'columns':list(v),'first_ids':list(v.index[:10]),'index_unique':v.index.is_unique})
print(json.dumps({'columns':list(v),'first_ids':list(v.index[:10])}))
