"""Exact identity closure over supplied Ensembl IDs, symbols and verified gene maps."""
import re
IDENTITY_FIELDS=('gene_ids','gene_symbol','gene_id_clean','gene_symbol_std','merge_gene','stable_ensembl_id','display_gene_symbol','stable_feature_id','merge_key','source_id','source_symbol','atlas_gene')
def token(v):
 s=str(v).strip().upper()
 if s in ('','NAN','NONE','NA'):return ''
 return s.split('.')[0] if re.fullmatch(r'ENSG\d+(?:\.\d+)?',s) else s
def row_tokens(index,row):
 return {t for t in [token(index)]+[token(row[k]) for k in IDENTITY_FIELDS if k in row] if t}
def excluded_mask(var,excluded):
 import numpy as np
 return np.array([bool(row_tokens(i,r)&excluded) for i,r in var.iterrows()],dtype=bool)
def assert_excluded(var,excluded):
 assert not excluded_mask(var,excluded).any(), 'Target identity or alias remains in model feature matrix'
