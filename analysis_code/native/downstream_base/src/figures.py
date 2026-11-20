"""Data-derived figures and per-panel provenance. No historical numerical targets."""
from pathlib import Path
import json
import h5py,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .core import sha,atomic_json,Blocked

def export_figures(config,obs,run_dir):
    run=Path(run_dir);out=run/'figures';out.mkdir(exist_ok=False)
    pf=json.load(open(run/'preflight.json'));key=pf['global_key_map']['umap']
    with h5py.File(config['atlas']['path'],'r') as f:xy=f['obsm'][key][:]
    if xy.shape!=(len(obs),2):raise Blocked('Registered UMAP shape invalid')
    panels=[]
    for col in [config['atlas']['level1_key'],config['atlas']['level2_key'],'study','tissue','comparison_group']:
        labels=obs[col].astype(str);cats=sorted(labels.unique());colors=plt.cm.tab20(np.linspace(0,1,len(cats))) if len(cats)<=20 else plt.cm.gist_ncar(np.linspace(0,.95,len(cats)))
        rng=np.random.default_rng(config['execution']['seed']);order=rng.permutation(len(obs))
        fig,ax=plt.subplots(figsize=(12,8));lut={c:colors[i] for i,c in enumerate(cats)}
        ax.scatter(xy[order,0],xy[order,1],c=np.array([lut[v] for v in labels.iloc[order]]),s=.12,linewidths=0,rasterized=True)
        for c in cats:ax.scatter([],[],color=lut[c],s=12,label=c)
        ax.legend(bbox_to_anchor=(1,1),loc='upper left',fontsize=5,ncol=1 if len(cats)<30 else 2,frameon=False)
        ax.set(xlabel='Frozen global UMAP 1',ylabel='Frozen global UMAP 2',title=col);fig.tight_layout()
        for ext in ['png','pdf']:fig.savefig(out/f'UMAP_{col}.{ext}',dpi=300,bbox_inches='tight')
        plt.close(fig)
        pd.DataFrame({'category':cats,'RGBA':[json.dumps(lut[c].tolist()) for c in cats]}).to_csv(out/f'palette_{col}.tsv',sep='\t',index=False)
        panels.append({'panel':f'UMAP_{col}','input':config['atlas']['path'],'input_sha256':pf['sha256'],'embedding_key':key,'recomputed':False,'all_cells':len(obs),'shuffle_seed':config['execution']['seed'],'claim':'descriptive location only'})
    s=obs.groupby([config['atlas']['level1_key'],config['atlas']['level2_key']],observed=True).size().rename('n_cells').reset_index();s['fraction_all_QC']=s.n_cells/len(obs);s.to_csv(out/'composition_source.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(12,max(4,len(s)*.18)));ax.barh(np.arange(len(s)),s.n_cells);ax.set_yticks(np.arange(len(s)),s.iloc[:,1],fontsize=5);ax.set_xlabel('Cells');fig.tight_layout();fig.savefig(out/'L2_composition.png',dpi=250);plt.close(fig)
    pd.DataFrame(panels).to_csv(out/'panel_manifest.tsv',sep='\t',index=False)
    # Source linkage is generated only from actual available result files.
    files=[{'path':str(p.relative_to(run)),'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()]
    atomic_json(out/'FIGURE_MANIFEST.json',files)
    return {'status':'COMPLETE','n_panels':len(panels)+1}

def dotplot(config,run_dir,level='L2'):
    run=Path(run_dir)
    if level=='L2':
        panel_path=config['figures'].get('l2_anchor_panel')
        if panel_path:
            if not Path(panel_path).is_file():raise Blocked('Final L2 anchor panel not released')
            ont=pd.read_csv(panel_path,sep='\t')
        else:
            from .anchors import final_anchor_panel
            from .core import read_frame
            op=Path(config['outputs']['release_support'])/'level2_ontology.tsv'
            ont=final_anchor_panel(pd.read_csv(op,sep='\t'),read_frame(config['atlas']['path'],'var'),config.get('pseudobulk',{}).get('gene_symbol_column','gene_symbol'),str(op)+'#sha256='+sha(op))
            anchor_out=run/'figures';anchor_out.mkdir(exist_ok=True);ont.to_csv(anchor_out/'resolved_final_L2_anchors.tsv',sep='\t',index=False)
        required={'gene_id','panel','level1_annotation','level2_annotation','approved','evidence_source'}
        if not required.issubset(ont):raise Blocked('L2 anchor panel missing provenance/parent/state/approval')
        if not ont.approved.astype(str).str.lower().eq('true').all():raise Blocked('L2 anchors await approval')
        if ont.evidence_source.isna().any() or ont.evidence_source.astype(str).str.strip().eq('').any():raise Blocked('Missing L2 anchor evidence')
        ontology=pd.read_csv(Path(config['outputs']['release_support'])/'level2_ontology.tsv',sep='\t')
        keys=['level1_annotation','level2_annotation']
        expected=set(map(tuple,ontology.loc[~ontology.level2_annotation.isin(['Unknown','Not_assessed']),keys].to_numpy()))
        observed=set(map(tuple,ont[keys].to_numpy()))
        if expected!=observed:raise Blocked('L2 anchor state coverage differs from final ontology')
    else:ont=pd.read_csv(config['figures']['marker_panel'],sep='\t')
    if not {'gene_id','panel'}.issubset(ont):raise Blocked('Marker panel needs gene_id,panel (resolve gene symbols explicitly)')
    from .core import read_frame
    atlas_gene_ids=set(read_frame(config['atlas']['path'],'var').index.astype(str))
    missing_genes=set(ont.gene_id.astype(str))-atlas_gene_ids
    if missing_genes:raise Blocked('Marker panel has missing exact gene IDs: '+repr(sorted(missing_genes)))
    parts=[]
    for block in pd.read_csv(run/'markers'/f'{level}_all_gene_marker_evidence.tsv.gz',sep='\t',chunksize=100000):parts.append(block[block.gene_id.isin(ont.gene_id)])
    d=pd.concat(parts,ignore_index=True)
    if d.empty:raise Blocked('No exact marker IDs matched; no silent symbol substitution')
    key=config['atlas']['level2_key'] if level=='L2' else config['atlas']['level1_key'];states=sorted(d[key].unique());genes=list(dict.fromkeys(ont.gene_id))
    fig,ax=plt.subplots(figsize=(max(9,len(genes)*.25),max(5,len(states)*.18)))
    sc=ax.scatter(d.gene_id.map({g:i for i,g in enumerate(genes)}),d[key].map({v:i for i,v in enumerate(states)}),s=d.detection_fraction*70,c=d.mean_log1p,cmap='Reds',linewidths=.1,edgecolors='grey')
    ax.set_xticks(range(len(genes)),genes,rotation=90,fontsize=6);ax.set_yticks(range(len(states)),states,fontsize=6);fig.colorbar(sc,ax=ax,label='Mean CP10k log1p');ax.set_title('Dot area: raw-count detection fraction');fig.tight_layout()
    out=run/'figures';out.mkdir(exist_ok=True);fig.savefig(out/f'{level}_DotPlot.png',dpi=300);fig.savefig(out/f'{level}_DotPlot.pdf');plt.close(fig)
    d.to_csv(out/f'{level}_DotPlot_source.tsv',sep='\t',index=False)
    return {'status':'COMPLETE','n_rows':len(d)}
