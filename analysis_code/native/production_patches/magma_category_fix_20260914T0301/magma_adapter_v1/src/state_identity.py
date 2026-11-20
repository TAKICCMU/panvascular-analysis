"""Exact, reversible final-ontology identities; never infer states from spelling."""
import hashlib,json
import pandas as pd
MISSING={'','Unknown','Not_assessed','NR','nan','None','not_available'}

def build_state_identity_map(ontology):
    cols=['level1_annotation','level2_annotation']
    if not set(cols)<=set(ontology):raise ValueError('Final ontology columns missing')
    if ontology[cols].isna().any().any():raise ValueError('Null ontology identity')
    # Reserved denominator rows are not CoVarNet states and need no assigned state ID.
    ontology=ontology.loc[~ontology[cols[0]].astype(str).isin(MISSING)&~ontology[cols[1]].astype(str).isin(MISSING)].copy()
    if 'level2_state_id' in ontology:
        if 'state_id' in ontology and not ontology.state_id.astype(str).equals(ontology.level2_state_id.astype(str)):raise ValueError('Conflicting state_id and level2_state_id')
        ontology['state_id']=ontology.level2_state_id
    x=ontology[cols+(['state_id'] if 'state_id' in ontology else [])].copy()
    if x.isna().any().any():raise ValueError('Null named-state ontology identity')
    x=x.astype(str).drop_duplicates()
    if 'state_id' not in x:
        if x.level2_annotation.duplicated().any():raise ValueError('Final L2 label has multiple parents without explicit state_id')
        x['state_id']=x.level2_annotation
    if x.state_id.isin(MISSING).any():raise ValueError('Missing stable state_id')
    if x.duplicated(cols).any():raise ValueError('Parent/L2 pair has multiple state IDs')
    x['parent_qualified']=x.level1_annotation+' | '+x.level2_annotation
    x['safe_label']=['S_'+hashlib.sha256(json.dumps([p,s],ensure_ascii=False).encode()).hexdigest()[:24] for p,s in x[cols].itertuples(index=False,name=None)]
    for col in ['state_id','parent_qualified','safe_label']:
        if x[col].duplicated().any():raise ValueError('Non-bijective state identity: '+col)
    x['spatial_label']=x.parent_qualified
    x['spatial_reference_state_original']=x.parent_qualified
    x['spatial_reference_state_safe']=x.safe_label
    return x.reset_index(drop=True)

def map_state_axis(values,state_map,source,target='state_id'):
    if source not in state_map or target not in state_map:raise ValueError('Unknown state identity namespace')
    if state_map[source].duplicated().any() or state_map[target].duplicated().any():raise ValueError('State map is not bijective')
    lookup=dict(zip(state_map[source].astype(str),state_map[target].astype(str)))
    vals=list(map(str,values));missing=sorted(set(vals)-set(lookup))
    if missing:raise ValueError('States missing from final ontology: '+repr(missing))
    return pd.Index([lookup[v] for v in vals])
