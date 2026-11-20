"""Resolve approved final-L2 anchors against exact atlas gene IDs/symbols."""
import json,re
import pandas as pd
from .core import Blocked

def final_anchor_panel(ontology,genes,symbol_column='gene_symbol',source='final_ontology'):
    required={'level1_annotation','level2_annotation','marker_anchor','approved'}
    if not required.issubset(ontology):raise Blocked('Final ontology lacks approved L2 marker_anchor contract')
    named=ontology.loc[~ontology.level2_annotation.isin(['Unknown','Not_assessed'])].copy()
    if not named.approved.astype(str).str.lower().isin(['true','1']).all():raise Blocked('Final L2 anchors not approved')
    ids=set(genes.index.astype(str));symbols={}
    if symbol_column in genes:
        for gene,value in genes[symbol_column].items():
            if pd.notna(value) and str(value).strip():symbols.setdefault(str(value),set()).add(str(gene))
    rows=[]
    for _,r in named.iterrows():
        val=r.marker_anchor
        if pd.isna(val) or not str(val).strip():raise Blocked('Missing marker anchors for '+str(r.level2_annotation))
        val=str(val).strip()
        tokens=json.loads(val) if val.startswith('[') else re.split(r'[,;/|]',val)
        if not isinstance(tokens,list) or not tokens:raise Blocked('Anchor must be gene list')
        for token in tokens:
            token=str(token).strip();matches=({token} if token in ids else set())|symbols.get(token,set())
            if len(matches)!=1:raise Blocked(f'Anchor must resolve exactly once: {token}, matches={sorted(matches)}')
            rows.append({'level1_annotation':r.level1_annotation,'level2_annotation':r.level2_annotation,'gene_id':next(iter(matches)),'anchor_original':token,'panel':r.level2_annotation,'approved':True,'evidence_source':source})
    return pd.DataFrame(rows).drop_duplicates(['level1_annotation','level2_annotation','gene_id'])
