"""Resolve only existing var fields; keep ambiguous/missing mappings explicit."""
from pathlib import Path
import pandas as pd
p=Path(__file__).resolve().parents[1];v=pd.read_csv(p/'evidence/working_gene_metadata.tsv.gz',sep='\t',index_col=0)
m=pd.DataFrame({'feature_id':v.index,'gene_symbol':v.gene_symbol.astype('string').to_numpy(),'ensembl_id':v.stable_ensembl_id.astype('string').to_numpy()})
m['ensembl_id']=m.ensembl_id.str.replace(r'\.\d+$','',regex=True)
m['valid_ensembl']=m.ensembl_id.str.match(r'^ENSG\d+$',na=False)
m['duplicate_ensembl']=m.ensembl_id.duplicated(keep=False)&m.valid_ensembl
m['duplicate_symbol']=m.gene_symbol.duplicated(keep=False)&m.gene_symbol.notna()
m['status']='PROPOSED_FROM_WORKING_VAR_NOT_FINAL_RELEASE'
m.to_csv(p/'configs/gene_map_PROPOSED.tsv',sep='\t',index=False)
print({'n_features':len(m),'valid_ensembl':int(m.valid_ensembl.sum()),'duplicate_ensembl_rows':int(m.duplicate_ensembl.sum()),'duplicate_symbol_rows':int(m.duplicate_symbol.sum())})
