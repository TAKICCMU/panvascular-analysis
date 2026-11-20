"""One streamed raw pass for raw-sum pseudobulk and marker detection statistics."""
from pathlib import Path
import json,gzip
import numpy as np
import pandas as pd
from scipy import sparse,io
from .core import Blocked,expression_blocks,read_frame,atomic_json

def factor_pairs(frame,cols):
    # MultiIndex avoids parsing pipes/spaces in biological labels.
    idx=pd.MultiIndex.from_frame(frame[cols].astype(str)); codes,uniques=pd.factorize(idx,sort=True)
    return codes,uniques.to_frame(index=False).set_axis(cols,axis=1)

def aggregate_sparse(path,obs,cols,slot='X',block_rows=4096,detection=False,integer=False):
    expected=read_frame(path,columns=[]).index
    if not obs.index.equals(expected):raise Blocked('Expression/metadata cell order differs')
    codes,units=factor_pairs(obs,cols);genes=read_frame(path,'var');ng=len(genes)
    out=sparse.csr_matrix((len(units),ng),dtype='int64' if integer or detection else 'float64')
    # obs must have complete atlas rows; excluded records encoded with a keep column.
    keep=obs['_keep'].to_numpy(bool) if '_keep' in obs else np.ones(len(obs),bool)
    for start,stop,x in expression_blocks(path,slot,block_rows):
        if stop>len(obs):raise Blocked('Expression/metadata row count mismatch')
        if not np.isfinite(x.data).all() or (x.data<0).any():raise Blocked('Invalid expression values')
        if (integer or detection) and not np.equal(x.data,np.floor(x.data)).all():raise Blocked('Raw matrix contains fractional counts')
        x=x.astype('int64' if integer or detection else 'float64')
        if detection:x.data=(x.data>0).astype('int64');x.eliminate_zeros()
        loc=np.flatnonzero(keep[start:stop]);weights=sparse.csr_matrix((np.ones(len(loc),dtype='int64'),(codes[start:stop][loc],loc)),shape=(len(units),stop-start))
        out=out+weights@x
    units['n_cells']=np.bincount(codes[keep],minlength=len(units));valid=units.n_cells.gt(0).to_numpy()
    return out[valid].tocsr(),units.loc[valid].reset_index(drop=True),genes

def export_pseudobulk(config,obs,cohort,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    c=cohort.set_index('sample_uid');aligned=c.reindex(obs.sample_uid)
    if aligned.index.isna().any() or aligned.biological_sample_uid.isna().any():raise Blocked('Missing biological metadata')
    xobs=obs.copy();xobs['biological_sample_uid']=aligned.biological_sample_uid.to_numpy()
    xobs['_keep']=aligned.inference_eligible.to_numpy(bool)
    records={}
    for name,key in [('L1',config['atlas']['level1_key']),('L2',config['atlas']['level2_key'])]:
        xobs['cell_state']=obs[key].astype(str)
        xobs['_keep']=aligned.inference_eligible.to_numpy(bool) & ~xobs.cell_state.isin(['Unknown','Not_assessed']).to_numpy()
        # L2 states are keyed by parent too; no ontology string parsing.
        cols=['biological_sample_uid','cell_state']
        count,units,genes=aggregate_sparse(config['atlas']['path'],xobs,cols,integer=True,block_rows=config['execution']['block_rows'])
        units['pseudobulk_id']=[f'{name}_PB{i:06d}' for i in range(len(units))]
        units['eligible_min_cells']=units.n_cells.ge(config['pseudobulk']['min_cells'])
        one=cohort.loc[cohort.inference_eligible].drop_duplicates('biological_sample_uid')
        for bio,g in cohort.loc[cohort.inference_eligible].groupby('biological_sample_uid'):
            for field in ['donor_uid','study','comparison_group','vascular_bed','tissue','assay_type','cohort_disease_context','age','sex']:
                if g[field].fillna('<MISSING>').astype(str).nunique()>1:raise Blocked(f'Cannot aggregate inconsistent {field} for {bio}')
        units=units.merge(one.drop(columns=['n_cells','sample_uid'],errors='ignore'),on='biological_sample_uid',how='left',validate='many_to_one')
        folder=out/name;folder.mkdir()
        sparse.save_npz(folder/'counts_units_by_genes.npz',count)
        with gzip.open(folder/'counts_genes_by_units.mtx.gz','wb') as f:io.mmwrite(f,count.T,field='integer')
        units.to_csv(folder/'units.tsv',sep='\t',index=False)
        genes.reset_index().to_csv(folder/'genes.tsv',sep='\t',index=False)
        records[name]={'shape':list(count.shape),'raw_count_sum':int(count.sum()),'n_units':len(units),'min_cells':config['pseudobulk']['min_cells'],'low_coverage_units_retained_but_ineligible':int((~units.eligible_min_cells).sum())}
    atomic_json(out/'EXPORT.json',records)
    return {'status':'COMPLETE','outputs':records}

def marker_summary(config,obs,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    records={}
    for name,key in [('L1',config['atlas']['level1_key']),('L2',config['atlas']['level2_key'])]:
        cols=[config['atlas']['level1_key']] if name=='L1' else [config['atlas']['level1_key'],key]
        sums,units,genes=aggregate_sparse(config['atlas']['path'],obs,cols,slot='layers/log1p',block_rows=config['execution']['block_rows'])
        detections,du,_=aggregate_sparse(config['atlas']['path'],obs,cols,slot='X',detection=True,block_rows=config['execution']['block_rows'])
        if not units.equals(du):raise Blocked('Raw/log group alignment mismatch')
        total_n=units.n_cells.sum();total=sums.sum(axis=0).A1;total_d=detections.sum(axis=0).A1
        scopes=['global_rest'] if name=='L1' else ['global_rest','within_parent_rest']
        for scope in scopes:
            filename=f'{name}_all_gene_marker_evidence.tsv.gz' if scope=='global_rest' else f'{name}_within_parent_marker_evidence.tsv.gz'
            with gzip.open(out/filename,'wt') as f:
                for i,row in units.iterrows():
                    n=int(row.n_cells);v=sums.getrow(i).toarray().ravel();vd=detections.getrow(i).toarray().ravel()
                    mean=v/n;det=vd/n
                    if scope=='within_parent_rest':
                        parent_mask=units[cols[0]].eq(row[cols[0]]).to_numpy()
                        scope_n=int(units.loc[parent_mask,'n_cells'].sum())
                        scope_sum=sums[parent_mask].sum(axis=0).A1;scope_det=detections[parent_mask].sum(axis=0).A1
                    else:scope_n,scope_sum,scope_det=total_n,total,total_d
                    rest_n=scope_n-n
                    rest_mean=(scope_sum-v)/rest_n if rest_n else np.full(len(genes),np.nan)
                    rest_det=(scope_det-vd)/rest_n if rest_n else np.full(len(genes),np.nan)
                    t=pd.DataFrame({'gene_id':genes.index,'mean_log1p':mean,'detection_fraction':det,'mean_log1p_rest':rest_mean,'detection_fraction_rest':rest_det,'mean_difference':mean-rest_mean,'n_cells':n,'rest_n_cells':rest_n,'reference_scope':scope})
                    for col in cols:t[col]=row[col]
                    t['descriptive_rank']=t.mean_difference.rank(method='min',ascending=False)
                    t.to_csv(f,sep='\t',index=False,header=(i==0))
        records[name]={'states':len(units),'genes':len(genes),'reference_scopes':scopes,'method':'descriptive mean/detection, no cell-level disease P values'}
    atomic_json(out/'MARKER_METHOD.json',records);return {'status':'COMPLETE','outputs':records}
