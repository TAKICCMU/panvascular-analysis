"""Opt-in CPU-only nsNMF integration test on synthetic counts; never production data."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from src.covarnet import prepare, write_prepared, run_nmf
from src.counts import state_counts


def main(out):
    out=Path(out)
    if out.exists():raise FileExistsError('Fresh synthetic output required')
    rng=np.random.default_rng(31)
    samples=[f'synthetic_s{i}' for i in range(36)]
    md=pd.DataFrame(dict(sample_uid=samples,biological_sample_uid=samples,
        donor_uid=[f'synthetic_person_{i%18}' for i in range(36)],donor_resolved=True,
        study=[f'synthetic_study_{i//12}' for i in range(36)],
        tissue=['artery' if i<18 else 'valve' for i in range(36)],
        condition='synthetic_group',aggregate_approved=False))
    records=[]
    for i,s in enumerate(samples):
        for j in range(6):
            n=30+int(rng.integers(0,30))
            records.extend([(s,'Fibroblast' if j<3 else 'Macrophage',f'synthetic_state_{j}')]*n)
        records.extend([(s,'Fibroblast','Unknown')]*20)
        records.extend([(s,'Unknown','Not_assessed')]*20)
    obs=pd.DataFrame(records,columns=['sample_uid','level1_annotation','level2_annotation'])
    obs.index=[f'synthetic_cell_{i}' for i in range(len(obs))]
    counts=state_counts(obs)
    ontology=obs[['level1_annotation','level2_annotation']].drop_duplicates().copy()
    ontology['covarnet_eligible']=~ontology.level2_annotation.isin(['Unknown','Not_assessed'])
    options=dict(rank_min=2,rank_max=4,rank_nrun=2,final_nrun=2,stability_nrun=2,
                 n_bootstrap=2,parent_null_permutations=9,workers=1,threads_per_worker=1)
    prepared=prepare(counts,md,ontology,options)
    write_prepared(prepared,out,{'fixture_only':True,'synthetic':True,'production_parameters_modified':False})
    result=run_nmf(out,options)
    assert result['status']=='COMPLETE'
    verify(out)

def verify(out):
    out=Path(out)
    manifest=json.loads((out/'INPUT_MANIFEST.json').read_text())
    assert manifest['source_identity']['fixture_only'] is True
    result=json.loads((out/'COVARNET_COMPLETE.json').read_text())
    assert result['status']=='COMPLETE'
    options=result['parameters']
    stable=out/'stable_NMF'
    network=pd.read_csv(stable/'stable_node_pairwise_spearman_network.tsv.gz',sep='\t')
    assert not {'P','FDR','pvalue'}&set(network)
    assert network.inference.eq('DESCRIPTIVE_NO_P_OR_FDR').all()
    assert network.n_profiles.eq(36).all()
    assert network.n_resampling_identity_blocks.eq(18).all()
    raw=network[network['transform'].eq('raw_proportion')].sort_values(['node1','node2']).rho.to_numpy()
    logit=network[network['transform'].eq('logit_proportion')].sort_values(['node1','node2']).rho.to_numpy()
    assert np.allclose(raw,logit,equal_nan=True)
    hfiles=list(stable.glob('sample_module_relative_coefficient_K*.tsv'));assert len(hfiles)==1
    h=pd.read_csv(hfiles[0],sep='\t',index_col=0)
    assert np.allclose(h.sum(axis=1),1)
    assert not list(stable.glob('sample_module_abundance*'))
    audit=pd.read_csv(next(stable.glob('module_matching_collapse_audit_K*.tsv')),sep='\t')
    assert set(audit.matching_mode)=={'inherited_many_to_one_primary','one_to_one_hungarian_sensitivity'}
    assert (audit.loc[audit.matching_mode.eq('one_to_one_hungarian_sensitivity'),'n_reference_rows_in_duplicate_matches']==0).all()
    checks={'status':'PASS','synthetic_only':True,'n_cells':manifest['n_cells'],'n_profiles':36,'n_nodes':6,
            'native_nmf':True,'native_parent_null':True,'production_run':False,
            'assertions':['real cell ledger->counts->prepare->nsNMF->parent null',
                          'descriptive rho no naive P/FDR','raw/logit are same rank relation',
                          'H relative coefficient column closure','two matching modes separate',
                          'one-to-one no duplicate matches'],'synthetic_parameters':options}
    (out/'REVISION_ASSERTIONS.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':
    if len(sys.argv)>2 and sys.argv[2]=='--verify-only':verify(sys.argv[1])
    else:main(sys.argv[1])
