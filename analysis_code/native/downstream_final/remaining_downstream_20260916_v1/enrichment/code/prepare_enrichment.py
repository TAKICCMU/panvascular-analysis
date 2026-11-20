import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1]
W=P.parents[1]
R=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
out=P/'inputs';out.mkdir(exist_ok=False)
gs=R/'resources/pathway_enrichment/MSigDB_2026.1_Hs_gene_sets_used.tsv.gz'
g=pd.read_csv(gs,sep='\t');assert g.db_version.nunique()==1 and str(g.db_version.iloc[0])=='2026.1.Hs'
sources={str(gs):sha(gs)};rows=[];audits=[]
def save(d,aid,source,kind,metadata):
    d=d.copy();d['gene']=d['gene'].fillna('').astype(str)
    good=d.gene.ne('') & np.isfinite(d.statistic)
    ambiguous=d.duplicated('gene',keep=False)
    for _,z in d[~good|ambiguous].iterrows():audits.append({'analysis_id':aid,'gene':z['gene'],'reason':'ambiguous_symbol' if z['gene'] in set(d[ambiguous].gene) else 'missing_symbol_or_statistic'})
    d=d[good & ~ambiguous];assert not d.gene.duplicated().any() and len(d)>100
    d['analysis_id']=aid
    p=out/(aid+'.tsv');d[['analysis_id','gene','statistic','selected_up','selected_down']].to_csv(p,sep='\t',index=False)
    sources[str(source)]=sha(source);sources[str(p)]=sha(p)
    rows.append(dict(id=aid,kind=kind,ranks=str(p),source=str(source),n_genes=len(d),metadata=metadata))
for scenario,folder in [('primary_r0.1','r0.1'),('primary_r0.6','r0.6'),('leave_one_CAD_r0.1','LOTr0.1')]:
    p=Path('/home/anzhen_vhd/更新版MAGMA_scDRS')/folder/'magma.genes.out'
    d=pd.read_csv(p,sep='\t');assert d.GENE.is_unique and d.P.between(0,1).all()
    # Genetic Z direction does not represent expression direction. The two ORA sets are explicitly renamed in reporting.
    from scipy.stats import false_discovery_control
    d['gene']=d.SYMBOL;d['statistic']=d.ZSTAT
    d['selected_up']=d.P < .05/len(d);d['selected_down']=false_discovery_control(d.P.to_numpy()) < .05
    save(d,'MAGMA_'+scenario,p,'genetic',{'ORA_up':'Bonferroni_associated','ORA_down':'FDR_associated','direction_not_expression':True})
pb=W/'pseudobulk_3backgrounds_20260916_v1'
roster=pd.read_csv(pb/'tables/model_roster.tsv',sep='\t').set_index('model_id')
for p in sorted((pb/'results').glob('*/DE.tsv')):
    status=json.loads((p.parent/'MODEL.json').read_text());assert status['status']=='COMPLETE'
    d=pd.read_csv(p,sep='\t');assert d.gene_id.is_unique and d.PValue.between(0,1).all()
    negative=d.F.lt(0)
    assert np.allclose(d.loc[negative,'PValue'],1), 'Negative F with P<1 requires review'
    for _,z in d[negative].iterrows():audits.append({'analysis_id':'PB_'+p.parent.name,'gene':z.gene_symbol,'reason':'native_negative_F_with_P1_ranked_zero_only; original_F='+str(z.F)})
    d['gene']=d.gene_symbol;d['statistic']=np.sign(d.logFC)*np.sqrt(d.F.clip(lower=0))
    d['selected_up']=(d.FDR<.05)&(d.logFC>0);d['selected_down']=(d.FDR<.05)&(d.logFC<0)
    aid='PB_'+p.parent.name
    save(d,aid,p,'disease_context',roster.loc[p.parent.name].fillna('').to_dict())
    sources[str(p.parent/'MODEL.json')]=sha(p.parent/'MODEL.json')
assert sum(x['kind']=='disease_context' for x in rows)==197
pd.DataFrame(audits).to_csv(P/'gene_symbol_exclusion_audit.tsv',sep='\t',index=False)
side=W/'disease_background_3groups_review_v2/sample_background_sidecar_v2.tsv'
assert sha(side)=='852c08353c53a1b4d32b0b5928c9d1ac6a0278e365748912adfaf2ff30cf2ce4'
sources[str(side)]=sha(side)
(P/'ANALYSES.json').write_text(json.dumps(rows,indent=2))
(P/'PREPARATION.json').write_text(json.dumps(dict(status='RANKS_PREPARED',analyses=len(rows),gene_sets=str(gs),bound_files=sources,method='197 approved three-background single-context DE signed sqrt(QL F) ranks, not SE; three genetic native ZSTAT ranks; ambiguous symbols excluded wholesale, no minimum-P selection; ORA BH over all nonempty testable pathways before display filter; no unapproved approximate cross-context meta'),indent=2))
print('PREPARED',len(rows),flush=True)
