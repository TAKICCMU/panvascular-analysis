import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[key]='1'
from pathlib import Path
import json,re,hashlib
import numpy as np
import pandas as pd
from scipy import sparse
from audit_scope import R,D,BASE,OLD,W,PHASE,read,save,sha

def good_gene(g):
    return not re.match(r'^(MT-|MTRNR|RPL|RPS|MRPL|MRPS)',g)

def main():
    m=read(W/'sample_sidecar_updated.tsv')
    m=m[m.covarnet_include.astype(str).str.lower().eq('true')].copy()
    # Same established minimum total cells, evaluated on distinct user-defined samples.
    small=m[m.n_cells.astype(int)<100].copy();save('profiles_below_100_cells.tsv',small)
    m=m[m.n_cells.astype(int)>=100].copy()
    assert m.sample_uid.is_unique
    counts=read(PHASE/'counts/sample_state_counts.tsv')
    counts=counts[counts.sample_uid.isin(m.sample_uid)].copy()
    ont=read(PHASE/'counts/ontology.tsv')
    keys=['level1_annotation','level2_annotation']
    named=ont[~ont.level1_annotation.isin(['Unknown','Not_assessed','Pending'])&~ont.level2_annotation.isin(['Unknown','Not_assessed','Pending'])].copy()
    assert len(named)==116 and not named.duplicated(keys).any()
    oldnodes=read(OLD/'all116_CoVarNet_node_evidence_review.tsv').set_index(keys)
    top=read(OLD/'L2_top20_within_parent_descriptive_markers.tsv')
    anchors=read(OLD/'L2_anchor_within_parent_evidence.tsv')
    pb=R/'pseudobulk/descriptive_final_atlas_v1/outputs'
    genes=read(pb/'genes.tsv');assert genes.gene_id.is_unique
    gi={str(g):i for i,g in enumerate(genes.gene_id)}
    # Panels are inherited descriptive evidence, not an independent validation set.
    panels={};panelrows=[]
    for n in named.itertuples():
        ident=(n.level1_annotation,n.level2_annotation)
        t=top[(top.level1_annotation==ident[0])&(top.level2_annotation==ident[1])].sort_values('descriptive_rank')
        gs=[str(g) for g in t.loc[(t.mean_difference>0)&(t.detection_fraction>t.detection_fraction_rest),'gene_id'] if good_gene(str(g)) and str(g) in gi][:10]
        anchor=str(n.marker_anchor)
        if anchor in gi and anchor not in gs:gs.append(anchor)
        panels[ident]=gs
        for g in gs:panelrows.append(dict(level1_annotation=ident[0],level2_annotation=ident[1],gene=g,is_anchor=g==anchor,source='inherited top20 within-parent expression evidence plus annotation anchor; no selection using disease associations'))
    save('node_evidence_gene_panels.tsv',pd.DataFrame(panelrows))
    allgenes=sorted({g for p in panels.values() for g in p});sel=[gi[g] for g in allgenes];gix={g:i for i,g in enumerate(allgenes)}
    u1=read(pb/'sample_L1/units.tsv.gz');u2=read(pb/'sample_L2/units.tsv.gz')
    assert not u1.duplicated(['sample_uid','level1_annotation']).any()
    assert not u2.duplicated(['sample_uid']+keys).any()
    expected=counts.set_index(['sample_uid']+keys)['count'].sort_index()
    actual=u2[u2.sample_uid.isin(m.sample_uid)].set_index(['sample_uid']+keys).n_cells
    # Unit files can omit zero-cell states; positive counts must match exactly.
    actual=actual[actual>0].sort_index();expected=expected[expected>0]
    assert actual.index.equals(expected.index) and np.array_equal(actual.to_numpy(),expected.to_numpy())
    print('Validated pseudobulk unit-to-final-state counts; loading inherited sparse matrices',flush=True)
    x1=sparse.load_npz(pb/'sample_L1/counts_units_by_genes.npz')
    assert x1.shape==(len(u1),len(genes))
    lib1=np.asarray(x1.sum(axis=1)).ravel();v1=x1[:,sel].toarray();del x1
    x2=sparse.load_npz(pb/'sample_L2/counts_units_by_genes.npz')
    assert x2.shape==(len(u2),len(genes))
    lib2=np.asarray(x2.sum(axis=1)).ravel();v2=x2[:,sel].toarray();del x2
    parent={(r.sample_uid,r.level1_annotation):i for i,r in enumerate(u1.itertuples())}
    mm=m.set_index('sample_uid');cx=counts[counts['count']>0].merge(m[['sample_uid','study','donor_uid','donor_resolved','donor_block']],on='sample_uid',validate='many_to_one')
    sampledetails=[];genestats=[];reviews=[];studyrows=[]
    for n in named.itertuples():
        ident=(n.level1_annotation,n.level2_annotation);panel=panels[ident]
        z=cx[(cx.level1_annotation==ident[0])&(cx.level2_annotation==ident[1])]
        bs=z.groupby('sample_uid')['count'].sum();st=z.groupby('study')['count'].sum()
        cells=int(z['count'].sum());ns=len(bs);nst=len(st)
        frequency=cells>=500 and ns>=10 and nst>=3
        ix=np.flatnonzero(u2.sample_uid.isin(m.sample_uid)&u2.level1_annotation.eq(ident[0])&u2.level2_annotation.eq(ident[1])&(u2.n_cells>=20))
        local=[]
        for i in ix:
            unit=u2.iloc[i];j=parent[(unit.sample_uid,unit.level1_annotation)]
            restcells=int(u1.iloc[j].n_cells-unit.n_cells);restlib=lib1[j]-lib2[i]
            if restcells<20 or restlib<=0 or lib2[i]<=0:continue
            for g in panel:
                k=gix[g];target=float(v2[i,k]);rest=float(v1[j,k])-target
                assert rest>=0
                tcpm=target/lib2[i]*1e6;rcpm=rest/restlib*1e6
                local.append(dict(level1_annotation=ident[0],level2_annotation=ident[1],sample_uid=unit.sample_uid,study=mm.loc[unit.sample_uid,'study'],donor_block=mm.loc[unit.sample_uid,'donor_block'],gene=g,state_cells=int(unit.n_cells),rest_parent_cells=restcells,state_CPM=tcpm,rest_CPM=rcpm,difference_CPM=tcpm-rcpm))
        sampledetails.extend(local)
        supported=[];studysupported={};validsamples=0;validdonors=0
        if local:
            ld=pd.DataFrame(local);validsamples=ld.sample_uid.nunique();validdonors=ld.donor_block.nunique()
            # No repeated capture is counted as an independent donor replicate.
            donor=ld.groupby(['study','donor_block','gene'],as_index=False).difference_CPM.median()
            study=donor.groupby(['study','gene'],as_index=False).difference_CPM.median()
            for gene,gg in study.groupby('gene'):
                npos=int(gg.difference_CPM.gt(0).sum())
                genestats.append(dict(level1_annotation=ident[0],level2_annotation=ident[1],gene=gene,studies_evaluable=len(gg),studies_positive=npos,donor_median_then_study_median=True))
                if npos>=3:supported.append(gene)
            for stu,ss in study.groupby('study'):
                positive=ss[ss.difference_CPM>0].gene.tolist();studysupported[stu]=positive
                studyrows.append(dict(level1_annotation=ident[0],level2_annotation=ident[1],study=stu,positive_genes=';'.join(positive),n_positive_genes=len(positive)))
        coherent_studies=sum(len(v)>=2 for v in studysupported.values())
        program_support=len(supported)>=2 and coherent_studies>=3
        broad=str(n.label_scope)=='marker_parent_group'
        if not frequency:decision='EXCLUDE_LOW_REPRESENTATION';reason='未达到既定500细胞、10样本、3研究门槛'
        elif broad:decision='SENSITIVITY_ONLY_BROAD_PARENT';reason='标签是宽泛L1父类剩余群；不当作已验证独立L2状态，建议仅补充敏感性分析'
        elif not program_support:decision='HOLD_LIMITED_PROGRAM_REPLICATION';reason='当前可评估样本不足或多基因跨研究证据不够；保留L2标签，暂不列核心节点'
        else:decision='RECOMMEND_CORE';reason='通过既定频率门槛；现有已批准状态有多基因跨研究表达支持，建议纳入'
        dominant=st.idxmax() if len(st) else '';domfrac=float(st.max()/cells) if cells else 0
        remaining=z[z.study.ne(dominant)]
        old=oldnodes.loc[ident]
        reviews.append(dict(level2_state_id=n.level2_state_id,level1_annotation=ident[0],level2_annotation=ident[1],label_scope=n.label_scope,marker_anchor=n.marker_anchor,n_cells=cells,n_samples=ns,n_studies=nst,n_resolved_donors=z.loc[z.donor_resolved.astype(str).str.lower().eq('true'),'donor_uid'].nunique(),n_unresolved_sample_blocks=z.loc[~z.donor_resolved.astype(str).str.lower().eq('true'),'donor_block'].nunique(),frequency_pass=frequency,old_n_cells=int(old.n_cells),old_frequency_pass=old.frequency_pass_full_provisional,dominant_study=dominant,dominant_study_fraction=domfrac,dominant_sample_fraction=float(bs.max()/cells) if cells else 0,n_cells_without_largest_study=int(remaining['count'].sum()),n_samples_without_largest_study=remaining.sample_uid.nunique(),n_studies_without_largest_study=remaining.study.nunique(),expression_evaluable_samples=validsamples,expression_evaluable_donor_blocks=validdonors,replicated_positive_genes=';'.join(supported),n_replicated_positive_genes=len(supported),n_studies_multigene_positive=coherent_studies,evidence_panel=';'.join(panel),recommendation=decision,reason=reason,core_candidate=decision=='RECOMMEND_CORE',user_approved=False))
    review=pd.DataFrame(reviews);save('L2_node_eligibility_review.tsv',review)
    save('node_marker_sample_evidence.tsv',pd.DataFrame(sampledetails))
    save('node_marker_study_replication.tsv',pd.DataFrame(genestats))
    save('node_program_by_study.tsv',pd.DataFrame(studyrows))
    save('core_node_candidates.tsv',review[review.core_candidate])
    save('nodes_not_recommended_for_core.tsv',review[~review.core_candidate])
    save('analysis_profile_ledger.tsv',m)
    # Exact denominators preserved, even for Unknown and non-retained states.
    totals=counts.groupby('sample_uid')['count'].sum()
    mapping=named[keys+['level2_state_id']]
    c=counts.merge(mapping,on=keys,how='left',validate='many_to_one')
    core=set(review.loc[review.core_candidate,'level2_state_id'])
    matrix=c[c.level2_state_id.isin(core)].pivot(index='level2_state_id',columns='sample_uid',values='count').fillna(0).reindex(index=sorted(core),columns=m.sample_uid,fill_value=0)
    matrix.to_csv(W/'candidate_node_counts.tsv',sep='\t')
    matrix.divide(totals.reindex(m.sample_uid),axis=1).to_csv(W/'candidate_node_fraction_all_cells.tsv',sep='\t')
    zero=matrix.sum(axis=0).eq(0)
    save('zero_retained_node_profiles.tsv',m[m.sample_uid.isin(zero[zero].index)])
    rules={'frequency':'Existing thresholds unchanged: >=100 all-atlas cells per profile; >=500 state cells, >=10 nonzero profiles and >=3 studies per node','profile':'One separate input sample per profile, as user requested; not old automatic donor+tissue+condition pooling','expression_review':'Descriptive project-level triage, not official CoVarNet criteria or independent biological validation: inherited top20 non-ribosomal/non-mitochondrial evidence up to10 genes plus anchor; state and rest-parent >=20cells; full-library CPM; median sample differences within donor, then median donors per study. Recommend core when >=2 genes positive in >=3 studies and >=3 studies support >=2 genes. No P-values used. Frequency alone does not prove stability.','parent_groups':'Broad marker_parent_group labels are sensitivity-only suggestions; never relabelled or removed from atlas/denominators','concentration':'Report all source concentrations continuously; no automatic dominant-study cutoff','sensitivity':'Proposed donor-block bootstrap, leave-one-study-out and existing denominator/within-parent/parent-null checks. No sensitivity model executed here.','circularity':'State panels and prior annotations derive from the same atlas. This is descriptive reproducibility evidence, not independent validation. NMF stability remains untested.','authorization':'All node decisions are assistant recommendations awaiting user review.'}
    (W/'NODE_REVIEW_RULES.json').write_text(json.dumps(rules,indent=2,ensure_ascii=False))
    summary={'status':'ALL_116_NODES_AUDITED_RECOMMENDATIONS_PENDING_USER_REVIEW','profiles':len(m),'small_profiles_excluded':len(small),'cells':int(counts['count'].sum()),'nodes_audited':len(review),'frequency_pass':int(review.frequency_pass.sum()),'recommendations':review.recommendation.value_counts().to_dict(),'evidence_rows':len(sampledetails),'zero_core_profiles':int(zero.sum()),'full_atlas_expression_loaded':False,'NMF_run':False,'sensitivity_run':False,'annotation_changed':False,'pandas':pd.__version__,'numpy':np.__version__}
    (W/'NODE_AUDIT.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
    print(json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=='__main__': main()
