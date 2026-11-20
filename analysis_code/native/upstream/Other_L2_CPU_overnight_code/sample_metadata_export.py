# -*- coding: utf-8 -*-
"""Read ONLY obs from the frozen atlas, export sample/clinical metadata for manual review.
No counts/log1p/coordinates read, no h5ad write, no missing-value inference.
"""
from __future__ import annotations
from pathlib import Path
from datetime import datetime
import hashlib, json, re

ROOT=Path('/home/anzhen_vhd/scRNAseq/codex/panvascular_atlas_rebuild_v2_20260907/06_reanalysis_takeover_20260908')
DEFAULT_SOURCE=ROOT/'panvascular_CPU_theta2_level1_frozen_v1.h5ad'
CORE_FIELDS=['study','study_uid','sample_id','sample_name','donor_uid_v2','donor_uid','donor_id',
             'age','sex','tissue','condition','disease','disease_context']
MISSING_TOKENS={'','nr','na','n/a','nan','none','null','unknown','not reported','missing','<na>'}
PREFIX=re.compile(r'^(?:study|sample|donor|patient|subject|age|sex|gender|tissue|organ|disease|condition|'
                  r'clinical|phenotype|diagnosis|assay|platform|technology|chemistry|library|batch|'
                  r'geo|gse|gsm|accession|source|provenance|orig\.ident)(?:$|[_ .-])',re.I)
EXACT={'sample','study','donor','age','sex','gender','tissue','organ','disease','condition','disease_context',
       'assay','platform','technology','chemistry','batch','orig.ident'}
EXCLUDE=re.compile(r'cell[ _-]?(?:id|type|state)|barcode|annotation|leiden|louvain|cluster|^level[123]',re.I)


def _require(ok,msg):
    if not bool(ok):raise RuntimeError(msg)


def _scalar(x):
    import numpy as np
    import pandas as pd
    if x is None or x is pd.NA or x is pd.NaT:return None
    if isinstance(x,np.generic):x=x.item()
    if isinstance(x,float) and pd.isna(x):return None
    if isinstance(x,bytes):return x.decode('utf-8',errors='replace')
    if isinstance(x,(datetime,pd.Timestamp)):return x.isoformat()
    if isinstance(x,(str,int,float,bool)):return x
    return str(x)


def _json(x):return json.dumps(x,ensure_ascii=False,separators=(',',':'),allow_nan=False)


def _missing(x):return x is None or (isinstance(x,str) and x.strip().lower() in MISSING_TOKENS)


def _value_records(series):
    # Conserve original values. Explicitly discard zero-count categorical levels.
    records=[]
    for v,n in series.value_counts(dropna=False,sort=False).items():
        if int(n)==0:continue
        value=_scalar(v)
        records.append({'value':value,'value_type':type(value).__name__,'n_cells':int(n),
                        'missing_like':_missing(value),
                        'provisional':isinstance(value,str) and 'provisional' in value.lower()})
    # Coalesce equivalent actual-null representations, but never sex/age/study naming conventions.
    merged={}
    for r in records:
        k=(r['value_type'],_json(r['value']))
        if k in merged:merged[k]['n_cells']+=r['n_cells']
        else:merged[k]=r
    return list(merged.values())


def summarize_obs(obs,extra_fields=None):
    """Pure-dataframe implementation, also usable in tests; no mutations of obs."""
    import numpy as np
    import pandas as pd
    _require(obs.index.is_unique,'obs_names重复。')
    _require('sample_uid' in obs and obs.sample_uid.notna().all(),'缺少有效sample_uid；不猜测补填。')
    _require(not obs.sample_uid.astype(str).str.strip().str.lower().isin(MISSING_TOKENS).any(),'sample_uid中有空/占位值。')
    fields=[f for f in CORE_FIELDS if f in obs]
    for f in obs.columns:
        if f=='sample_uid' or f in fields:continue
        if (PREFIX.search(str(f)) or str(f).lower() in EXACT) and not EXCLUDE.search(str(f)):
            fields.append(f)
    for f in (extra_fields or []):
        _require(f in obs,f'指定额外字段不存在：{f}')
        if f not in fields and f!='sample_uid':fields.append(f)
    missing_fields=[f for f in CORE_FIELDS if f not in obs]
    n=len(obs);sample_strings=obs.sample_uid.astype(str)
    # Do not silently collapse numeric versus string sample_uid values that stringify alike.
    _require(obs.sample_uid.nunique(dropna=False)==sample_strings.nunique(),'sample_uid字符串转换造成重复。')
    groups=sample_strings.groupby(sample_strings,sort=True,observed=True).indices
    main=[];long=[];conflicts=[]
    for sid,ix in groups.items():
        block=obs.iloc[ix];row={'sample_uid':sid,'n_cells':len(ix)}
        if 'level1_annotation' in block:
            row['n_L1_Unknown']=int(block.level1_annotation.astype(str).eq('Unknown').sum())
            row['n_L1_known']=len(ix)-row['n_L1_Unknown']
        conflict_fields=[];partial_fields=[];allmissing_fields=[]
        for f in fields:
            vals=_value_records(block[f]);valid=[r for r in vals if not r['missing_like']]
            n_missing=sum(r['n_cells'] for r in vals if r['missing_like'])
            actual_null=sum(r['n_cells'] for r in vals if r['value'] is None)
            provisional=sum(r['n_cells'] for r in vals if r['provisional'])
            if len(valid)>1:
                value='[CONFLICT]';status='MULTIPLE_REPORTED_VALUES';conflict_fields.append(f)
                conflicts.append({'sample_uid':sid,'field':f,'n_cells':len(ix),
                                  'reported_values_counts_json':_json(vals),
                                  'action':'REVIEW_ONLY_NO_VALUE_CHOSEN'})
            elif len(valid)==1:
                value=valid[0]['value'];status='PARTIAL_MISSING' if n_missing else 'ONE_REPORTED_VALUE'
                if n_missing:partial_fields.append(f)
            else:
                # Preserve NR/other placeholder if there is one; true missing remains blank.
                literals=[r['value'] for r in vals if r['value'] is not None]
                value=literals[0] if len(literals)==1 else ('' if not literals else '[MISSING_TOKENS]')
                status='ALL_MISSING_OR_PLACEHOLDER';allmissing_fields.append(f)
            row[f]=value;row[f+'__status']=status
            row[f+'__n_reported_values']=len(valid)
            row[f+'__n_missing_like_cells']=n_missing
            row[f+'__n_provisional_cells']=provisional
            row[f+'__values_counts_json']=_json(vals)
            long.append({'sample_uid':sid,'field':f,'n_cells':len(ix),'n_reported_values':len(valid),
                         'n_actual_null_cells':actual_null,'n_missing_like_cells':n_missing,
                         'n_provisional_cells':provisional,'status':status,'values_counts_json':_json(vals)})
        for f in missing_fields:row[f]='';row[f+'__status']='SOURCE_FIELD_ABSENT'
        row['has_metadata_conflict']=bool(conflict_fields)
        row['conflicting_fields']=';'.join(conflict_fields)
        row['partially_missing_fields']=';'.join(partial_fields)
        row['entirely_missing_fields']=';'.join(allmissing_fields+missing_fields)
        # Empty correction fields: never auto-fill them with current values.
        for f in ('study','donor_uid','age','sex','tissue','disease_context'):
            row['verified_'+f]=''
        row['verification_reference']='';row['verification_note']=''
        main.append(row)
    main=pd.DataFrame(main)
    front=['sample_uid','study','n_cells','n_L1_known','n_L1_Unknown','donor_uid_v2','age','sex','tissue',
           'condition','disease','disease_context','has_metadata_conflict','conflicting_fields']
    front=[c for c in front if c in main]
    main=main[front+[c for c in main if c not in front]]
    long=pd.DataFrame(long)
    conflict=pd.DataFrame(conflicts,columns=['sample_uid','field','n_cells','reported_values_counts_json','action'])
    # Distinct joint clinical records preserve covariation (e.g. two donors with two ages).
    joint_fields=[f for f in ('study','study_uid','donor_uid_v2','donor_uid','age','sex','tissue',
                              'condition','disease','disease_context') if f in obs]
    joint=obs[['sample_uid']+joint_fields].copy()
    for c in joint:joint[c]=joint[c].astype(object).where(joint[c].notna(),None)
    joint=(joint.groupby(list(joint.columns),dropna=False,observed=True,sort=False).size()
           .rename('n_cells').reset_index())
    if 'level1_annotation' in obs:
        l1counts=pd.crosstab(sample_strings,obs.level1_annotation.astype(str),dropna=False)
        l1counts.index.name='sample_uid';l1counts=l1counts.reset_index()
    else:l1counts=pd.DataFrame({'sample_uid':list(groups)})
    inventory=pd.DataFrame({'source_field':list(obs.columns),'source_dtype':[str(obs[c].dtype) for c in obs],
                           'included_in_metadata_export':[c in fields or c=='sample_uid' for c in obs]})
    _require(int(main.n_cells.sum())==n and main.sample_uid.is_unique,'导出样本计数不守恒。')
    _require(int(joint.n_cells.sum())==n,'联合临床信息组合计数不守恒。')
    # Main CSV stays readable. Expanded fields/frequencies are preserved in companion tables.
    primary=['sample_uid','study','study_uid','sample_id','sample_name','n_cells','n_L1_known','n_L1_Unknown',
             'donor_uid_v2','donor_uid','donor_id','age','sex','tissue','condition','disease','disease_context',
             'has_metadata_conflict','conflicting_fields','partially_missing_fields','entirely_missing_fields',
             'donor_uid_v2__n_provisional_cells','donor_uid__n_provisional_cells',
             'verified_study','verified_donor_uid','verified_age','verified_sex','verified_tissue',
             'verified_disease_context','verification_reference','verification_note']
    compact=main[[c for c in primary if c in main]].copy()
    return {'sample_study_clinical_metadata':compact,'sample_metadata_wide_detail':main,'sample_metadata_values_long':long,
            'sample_metadata_conflicts':conflict,'sample_clinical_distinct_combinations':joint,
            'sample_level1_cell_counts':l1counts,'metadata_field_inventory':inventory}


def export_sample_metadata(source=DEFAULT_SOURCE,output_parent=None,extra_fields=None):
    """CSV only; requires anndata/h5py from the existing omicverse environment."""
    import h5py
    try:from anndata.io import read_elem
    except ImportError:from anndata._io.specs import read_elem
    source=Path(source).expanduser().resolve();_require(source.is_file(),f'找不到冻结文件：{source}')
    initial=(source.stat().st_size,source.stat().st_mtime_ns)
    # Direct HDF5 element read avoids even the heavy layers that backed AnnData may load.
    with h5py.File(source,'r') as f:
        obs=read_elem(f['obs'])
        _require('uns/level1_frozen' in f and bool(read_elem(f['uns/level1_frozen'])),'输入不是冻结L1。')
    tables=summarize_obs(obs,extra_fields)
    parent=Path(output_parent) if output_parent else source.parent
    out=parent/('sample_metadata_review_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));out.mkdir(parents=True,exist_ok=False)
    for name,frame in tables.items():frame.to_csv(out/(name+'.csv'),index=False,encoding='utf-8-sig',na_rep='')
    report={'source_h5ad':str(source),'source_stat':{'size':initial[0],'mtime_ns':initial[1]},
            'n_cells':len(obs),'n_samples':int(tables['sample_study_clinical_metadata'].shape[0]),
            'n_conflicting_sample_fields':len(tables['sample_metadata_conflicts']),
            'source_h5ad_modified':False,'expression_read':False,'missing_values_inferred':False,
            'meaning':'condition/disease are coarse as recorded; disease_context preserves detailed context; no external lookup.'}
    _require(initial==(source.stat().st_size,source.stat().st_mtime_ns),'导出期间源文件发生变化。')
    (out/'EXPORT_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'README_CN.txt').write_text(
        '主文件：sample_study_clinical_metadata.csv，每个sample_uid一行。\n'
        '不读取raw count/log1p矩阵，不修改h5ad。空值、NR和PROVISIONAL均按来源保留，不推断。\n'
        '同一样本一个字段有多个不同的非缺失值，主值显示[CONFLICT]，所有原值和对应细胞数在values_counts_json和conflicts表。\n'
        '一个字段仅部分细胞缺失时，ONE值只是样本级汇总，不代表已给缺失细胞回填。\n'
        '联合组合表保留study/donor/age/sex/疾病背景之间的原始对应组合。\n'
        'source的原字段及其单位优先，不将数字sex代码转换为男女，不平均年龄，不把sample推断成donor。\n'
        '补充信息请写verified_*、verification_reference和verification_note列。人工改CSV不会自动影响本轮L2结果或h5ad。\n'
        '新字段合并必须另行审核、按sample_uid/donor核对并新建metadata版本。禁止在运行期间覆盖冻结h5ad。\n',encoding='utf-8')
    print('临床元数据已导出：',out/'sample_study_clinical_metadata.csv')
    print(f"样本{report['n_samples']}；细胞{report['n_cells']:,}；样本-字段冲突{report['n_conflicting_sample_fields']}")
    return out


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',default=str(DEFAULT_SOURCE));args=p.parse_args()
    export_sample_metadata(args.source)
