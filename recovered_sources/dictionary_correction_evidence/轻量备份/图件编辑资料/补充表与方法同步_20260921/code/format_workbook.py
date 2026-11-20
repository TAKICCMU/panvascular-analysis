from pathlib import Path
from openpyxl import load_workbook
from openpyxl.styles import Font,PatternFill,Alignment,Border,Side
from openpyxl.utils import get_column_letter as L
from openpyxl.worksheet.page import PageMargins
import csv,json,re,shutil,hashlib,math,sys
from table_metadata import META
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/补充表与方法同步_20260921'; B=O/'轻量备份';B.mkdir(exist_ok=True); S=O/'待替换';S.mkdir(exist_ok=True)
src=P/'Supplementary_Tables.xlsx'
if not (B/src.name).exists():shutil.copy2(src,B/src.name)
w=load_workbook(B/src.name)
# Explicit display names. Source fields remain reversible through column_name_map.tsv.
D={
'gene symbol':'Gene','Gene symbol':'Gene','gene':'Gene','Gene':'Gene','gene id':'Ensembl gene ID','level1 annotation':'L1 category','level2 annotation':'L2 state','state':'State','state id':'State ID','cell state':'Cell state','cell context':'Cell context','annotation level':'Annotation level','level':'Annotation level','study':'Study ID','study canonical':'Study ID','GSE':'Study ID','sample uid':'Sample ID','analysis sample uid':'Analysis sample ID','biological sample uid':'Specimen ID','donor uid':'Donor ID','biological unit id':'Analysis unit ID','sample id':'Sample ID','section id':'Section ID','analysis profile uid':'Profile ID','old biological sample uid':'Previous specimen ID','old donor uid':'Previous donor ID','n cells':'Cells (n)','n cell':'Cells (n)','samples':'Samples (n)','studies':'Studies (n)','donors':'Donors (n)','tissues':'Tissues (n)','n samples':'Samples (n)','n donors':'Donors (n)','n studies':'Studies (n)','n tissues':'Tissues (n)','n sections':'Sections (n)','n spots':'Spots (n)','n positive spots':'Positive spots (n)','posterior':'Posterior summary','pathology':'Expression programme','term':'Model term','sensitivity':'Section set','k':'Units (n)','n datasets':'Datasets (n)','datasets':'Datasets (n)','biological units':'Analysis units (n)','biological units, n':'Analysis units (n)','Biological units, n':'Analysis units (n)',
'beta':'β','se':'SE','z':'Z','p':'P','pval':'P','P value':'P','Adjusted P value':'Adjusted P','padj':'Adjusted P','FDR BH':'BH FDR','tau2':'τ²','I2':'I² (%)','CI low':'95% CI lower','CI high':'95% CI upper','ci low':'95% CI lower','ci high':'95% CI upper','direction consistency':'Direction agreement','median I2 percent':'Median I² (%)',
'meta logFC':'Meta log2FC','meta SE':'Meta SE','meta Z':'Meta Z','meta P':'Meta P','meta CI low':'Meta 95% CI lower','meta CI high':'Meta 95% CI upper','target set meta FDR':'Target-set BH FDR','median context logFC':'Median context log2FC','n disease contexts':'Disease contexts (n)','n GWAS aligned contexts':'GWAS-aligned contexts (n)','n context FDR lt 0.05':'Significant contexts (n)',
'Pathway or gene set':'Pathway','pathway':'Pathway','Normalised enrichment score':'NES','Enrichment score':'ES','log2 error estimate':'log2 error estimate','log2err':'log2 error estimate','Gene-set size':'Gene-set size (n)','size':'Gene-set size (n)','leadingEdge':'Leading-edge genes','pathway size':'Pathway genes (n)','n background genes':'Background genes (n)','n selected genes':'Selected genes (n)','n ranked genes':'Ranked genes (n)','overlap':'Overlap (n)','Odds ratio':'OR','Analysis identifier':'Analysis ID','analysis id':'Analysis ID','record type':'Record type','analysis kind':'Analysis kind',
'mean log1p':'Mean log1p expression','detection fraction':'Detected fraction','mean log1p rest':'Other-cell mean log1p','detection fraction rest':'Other-cell detected fraction','mean difference':'Mean expression difference','rest n cells':'Other cells (n)','mean gene zscore':'Mean gene Z score','canonical program':'Canonical programme','reference scope':'Reference scope',
'pooled effect disease minus control':'Effect (disease − control)','pooled se':'Pooled SE','pooled p':'Pooled P','n cohorts':'Cohorts (n)','cohort effects':'Cohort contributions','source':'Sender state','target':'Receiver state','top cohort abs effect share':'Largest cohort effect share','context meta supported':'Context support','interaction key':'Interaction ID','threshold scheme':'Eligibility scheme','signaling category':'Signalling category',
'QTL status':'QTL status','qtl source':'QTL source','qtl sample size':'QTL sample size','n overlap':'Shared variants (n)','gwas effective N':'GWAS effective N','min gwas p':'Minimum GWAS P','min eqtl p':'Minimum eQTL P','H4 given H3H4':'PP.H4 / (PP.H3+PP.H4)','top shared snp':'Top shared variant','top shared SNP PP.H4':'Top-variant SNP.PP.H4','classification':'Evidence class','n contexts requested':'Requested contexts (n)','n contexts with qtl':'Available QTL contexts (n)','n contexts no qtl':'Missing QTL contexts (n)','n remote failures':'Retrieval failures (n)','n contexts tested':'Tested contexts (n)','n strong':'Strong contexts (n)','n suggestive':'Suggestive contexts (n)','n distinct H3':'Distinct-signal contexts (n)','maximum PP.H4':'Maximum PP.H4',
'module':'Module','module tier':'Module tier','分级':'Module tier (Chinese)','n total atlas cells':'Total atlas cells (n)','n technical samples':'Technical samples (n)','covarnet sample uid':'CoVarNet profile ID','covarnet donor key':'CoVarNet donor key','tissue v2':'Tissue','condition v2':'Condition','donor':'Association donor ID','n tests':'Tests (n)','FDR05':'FDR <0.05 (n)','n tests in family':'Family tests (n)','FDR BH all 29 pairs x 3 terms':'BH FDR (87 tests)','FDR BH 28':'BH FDR (28 tests)','FDR BH 84':'BH FDR (84 tests)',
'bras':'BRAS','clisi knn':'cLISI','graph connectivity':'Graph connectivity','ilisi knn':'iLISI','isolated labels':'Isolated-label score','kbet per label':'kBET','nmi ari cluster labels kmeans ari':'KMeans ARI','nmi ari cluster labels kmeans nmi':'KMeans NMI','pcr comparison':'PCR comparison','silhouette label':'Cell-type ASW','Bio conservation':'Biology score','Batch correction':'Batch score','Total':'Overall score','rank':'Rank',
'parent qualified':'Parent-qualified state','safe label':'Safe state label','spatial label':'Spatial factor label','spatial reference state original':'Reference state name','spatial reference state safe':'Reference safe label','selected genes':'Marker genes','analysis approved':'Analysis approval record','inferential design approved':'Design approval record','same donor samples are independent donors':'Treat same-donor samples as independent','automatic cross specimen pooling':'Automatic specimen pooling','group review required':'Grouping-review flag',
'primary recomputed n contexts tested':'Primary tested contexts (n)','primary recomputed n strong':'Primary strong contexts (n)','primary recomputed maximum PP.H4':'Primary maximum PP.H4','primary recomputed best context':'Primary best context','primary recomputed best classification':'Primary best evidence class','leave one CAD recomputed n contexts tested':'Leave-CAD tested contexts (n)','leave one CAD recomputed n strong':'Leave-CAD strong contexts (n)','leave one CAD recomputed maximum PP.H4':'Leave-CAD maximum PP.H4','leave one CAD recomputed best context':'Leave-CAD best context','leave one CAD recomputed best classification':'Leave-CAD best evidence class','primary FUMA eqtl mapping rows':'Primary FUMA eQTL rows (n)','leave CAD FUMA eqtl mapping rows':'Leave-CAD FUMA eQTL rows (n)','primary FUMA eqtl tissues':'Primary FUMA eQTL tissues','leave CAD FUMA eqtl tissues':'Leave-CAD FUMA eQTL tissues',
'Cells expressing in cluster (%)':'Within-cluster detected fraction','Cells expressing outside cluster (%)':'Outside-cluster detected fraction','n cells x':'Unknown cells (n)','n cells y':'Total sample cells (n)',
'VARIABLE':'MAGMA variable ID','variable':'Mapped variable ID','NGENES':'MAGMA genes (n)','BETA':'MAGMA β','BETA STD':'MAGMA standardised β','SE':'SE','P':'P',
'posterior':'Posterior summary','serial section covariance unavailable':'Section covariance unmodelled','mean outside q05 q95 entries':'Mean outside q05–q95 (n)','original source':'Source path','quantile order':'Quantile-order check','finite nonnegative':'Finite/nonnegative check','exact original spots':'Spot-identity check',
'legacy fallback  within  p':'Fallback within-spot P','legacy fallback  within  nullmean':'Fallback within-spot null mean','legacy fallback  lag  p':'Fallback lag P','legacy fallback  lag  nullmean':'Fallback lag null mean','strict parent  within  p':'Parent-matched within-spot P','strict parent  within  nullmean':'Parent-matched within-spot null mean','strict parent  lag  p':'Parent-matched lag P','strict parent  lag  nullmean':'Parent-matched lag null mean'
}
# Field semantics documented only where supported; no numeric rescaling.
def display(h,n):
 if n==26 and h=='gene id':return 'Gene identifier'
 if n==18 and h=='QC fallback':return 'Median QC-fallback indicator'
 if n==18 and h=='n spots':return 'Median spots per unit'
 if n==18 and h.endswith('  p'):return 'Median '+D[h]+' (descriptive)'
 if n==53 and h.startswith('FDR05 '):return ('Original-reference ' if h.endswith('old') else 'Target-excluded ')+'FDR <0.05'
 if n==19 and h=='Overlapping genes':return 'Overlap (n)'
 if n==5 and h in ['SE','P','FDR BH']:return {'SE':'MAGMA SE','P':'MAGMA P','FDR BH':'MAGMA BH FDR'}[h]
 if h.lower()=='scenario':return 'GWAS scenario'
 if n==53 and h.endswith((' old',' excluded')):
  core,suffix=h.rsplit(' ',1);return ('Original-reference ' if suffix=='old' else 'Target-excluded ')+D.get(core,core)
 if h.endswith((' LIANA',' CellChat')):
  core,suffix=h.rsplit(' ',1);return suffix+' '+D.get(core,core)
 if h.startswith('rho  '):return 'ρ: '+h[5:].replace('fibrous cap','fibrous-cap')
 if h in D:return D[h]
 h=re.sub(r', n$', ' (n)',h)
 h=h.replace('systemic','atherosclerotic').replace('Systemic','Atherosclerotic') if n in [7,8,9] else h
 h=re.sub(' +',' ',h)
 return h[0].upper()+h[1:]
SEM={
'Gene':('Reported gene symbol; original gene identities are unchanged.','gene symbol'),
'Ensembl gene ID':('Stable source gene identifier; not a row number or expression value.','Ensembl identifier'),
'L1 category':('Final Level 1 annotation. Unknown is a valid annotation, not missing data.','category'),
'L2 state':('Final Level 2 annotation within the stated parent lineage.','category'),
'Analysis unit ID':('Recorded grouping used to combine section estimates before between-unit meta-analysis; the grouping is retained as approved.','identifier'),
'Specimen ID':('Biological specimen identifier, distinct from a technical library and from donor identity.','identifier'),
'Sample ID':('Original technical sample identifier.','identifier'),
'Donor ID':('Recorded donor identifier; NR/blank is not a new independent donor.','identifier'),
'Study ID':('Source-study accession or its canonical identifier.','identifier'),
'Profile ID':('Retained analysis profile identifier; not interchangeable with donor count.','identifier'),
'Direction agreement':('Fraction of contributing effects with the recorded consensus direction; the qualifying units are defined by this table.','fraction (0–1)'),
'Within-cluster detected fraction':('Source fraction of cells with detectable expression within the source cluster; includes zero-expression cells in the denominator.','fraction (0–1)'),
'Outside-cluster detected fraction':('Source fraction of cells with detectable expression outside the source cluster.','fraction (0–1)'),
'Mean log1p expression':('Mean full-library log1p(CP10k) expression in the indicated group, including zero-expression cells.','log1p(CP10k)'),
'Detected fraction':('Within-group fraction of cells with raw counts >0. Missing measurement is not converted to zero.','fraction (0–1)'),
'Mean gene Z score':('Canonical programme summary on the recorded gene-standardised scale; not raw mean expression or Figure 1g min–max colour.','gene Z score'),
'Model term':('cell_abundance: conditional abundance coefficient; pathology: expression-programme coefficient; interaction: abundance × programme coefficient.','model term'),
'Expression programme':('Molecular expression signature used in the spatial model; not an independently assigned histopathological region.','programme name'),
'Posterior summary':('cell2location posterior summary used for abundance: q05 or posterior mean.','summary label'),
'Section set':('all_sections or the recorded exclusion of QC-fallback sections; not an expression-based exclusion.','setting'),
'Section covariance unmodelled':('Recorded flag that within-unit section residual covariance was not explicitly estimated. It is not a grouping-error flag.','boolean'),
'BH FDR (87 tests)':('BH-adjusted P across 29 candidate pairs × 3 model terms within the stated posterior/QC setting. Non-estimable terms retain their recorded P=1 accounting.','probability (0–1)'),
'Target-set BH FDR':('BH-adjusted meta-analysis P within cell state over the prespecified priority-gene set; distinct from all-gene FDR.','probability (0–1)'),
'PP.H4 / (PP.H3+PP.H4)':('Posterior probability of a shared signal conditional on either a distinct or shared two-trait signal.','probability (0–1)'),
'Top-variant SNP.PP.H4':('Variant-level posterior conditional on H4; not the region-level PP.H4.','probability (0–1)'),
'NES':('Normalised enrichment score from the indicated preranked GSEA; positive and negative values preserve direction.','signed enrichment score'),
'OR':('MAGMA ORA odds ratio with 0.5 continuity correction.','odds ratio'),
'Leading-edge genes':('GSEA leading-edge subset supplied by fgsea, not the full pathway membership.','gene list'),
'Overlap (n)':('Number of selected genes overlapping the tested pathway.','count'),
'Gene-set size (n)':('Number of pathway genes represented in the ranked GSEA input.','count'),
'Adjusted P':('Original adjusted P for the stated enrichment analysis; adjustment families follow the Methods.','probability (0–1)'),
'GSEA FDR05':('Number of GSEA pathways meeting native adjusted P <0.05 for the indexed analysis.','count'),
'ORA FDR05 overlap3':('Number of ORA pathways meeting native adjusted P <0.05 and recorded overlap ≥3 for the indexed analysis.','count'),
'Native NA':('Count of unavailable native enrichment outputs for this analysis index; not a zero P value.','count'),
'BRAS':('Batch-removal adjusted silhouette metric on its original scIB scale.','score'),
'Overall score':('0.6 × biology score + 0.4 × batch score on the implemented metric scales.','score'),
'Meta log2FC':('REML random-effects meta-analysis log2 disease-versus-reference fold change.','log2 fold change'),
'Meta SE':('Standard error of the REML meta-analysis log2 fold change.','log2 fold change'),
'Median direction consistency':('Median fraction of contributing studies with positive within-study correlation across the top-five-node pairs.','fraction (0–1)'),
'Positive pair fraction':('Fraction of evaluated top-five-node pairs with positive pooled correlation.','fraction (0–1)'),
'Scaled mean expression':('Recorded scaled expression, distinct from unscaled means; scaling definition follows the corresponding display.','scaled value'),
'Scope':('Analysis-specific scope of the associated count or result.','text')
}
def semantics(n,h,new):
 if n==26 and h=='gene id':return ('Gene symbol or source gene identifier present in the marker export; not uniformly an Ensembl ID.','source identifier')
 if n==18 and new.startswith('Median '):return ('Median of the recorded section-to-unit-to-dataset summary; P medians are descriptive, not a combined significance test.','original scale, median summary')
 if n==53 and h.startswith('FDR05 '):return ('Boolean indicator of BH FDR <0.05 for the stated reference implementation, not a count.','boolean')
 if n==35 and h=='term':return ('Disease-versus-reference coefficient in the stated beta-binomial proportion model.','model term')
 if new=='Metadata status' and n in [21,58,59]:return ('Source metadata status retained from model execution. It does not supersede author confirmation of the current analysis-unit grouping.','source status')
 if new in SEM:return SEM[new]
 if new.startswith('PP.H'):return ('Colocalisation hypothesis posterior; H0 neither, H1 GWAS only, H2 QTL only, H3 distinct signals, H4 shared signal.','probability (0–1)')
 if 'I²' in new:return ('Between-unit heterogeneity statistic on its original percentage scale.','percent (0–100)')
 if 'τ²' in new:return ('Estimated between-unit effect variance.','squared effect units')
 if new in ['β','SE','95% CI lower','95% CI upper'] and n in [32,48,53,57,58]:return ({'β':'Conditional spatial log-link coefficient for the stated model term.','SE':'Standard error of the conditional spatial coefficient; not a confidence interval.','95% CI lower':'Lower endpoint of the recorded 95% confidence interval.','95% CI upper':'Upper endpoint of the recorded 95% confidence interval.'}[new],'natural-log link coefficient')
 if new=='Units (n)':return ('Number of contributing analysis units in the spatial meta-analysis, not number of spots or gene–unit estimates.','count')
 if new=='Direction agreement' and n in [13,43]:return ('Fraction of disease contexts agreeing with the consensus log2-fold-change direction.','fraction (0–1)')
 if '(n)' in new or new.startswith(('N GSEA','N ORA')):return (new.replace(' (n)','')+' represented by this row.','count')
 if new in ['P','MAGMA P','scDRS P','Meta P','Pooled P'] or new.endswith(' P'):return ('Original P for '+new.replace(' P','')+' in the stated model; source numeric underflow zeros and eligibility flags are retained.','probability (0–1)')
 if 'FDR' in new and not any(s in new for s in ['family','Family','05','<0.05','fraction']):return ('Original adjusted P for the stated testing family; no correction has been recomputed.','probability (0–1)')
 if 'fraction' in new.lower() or 'share' in new.lower():return ('Recorded '+new.lower()+', with its analysis-specific denominator; no conversion to percent.','fraction (0–1)')
 if 'ID' in new or 'accession' in new.lower():return ('Recorded '+new+'; preserved as an identifier, not recoded into another unit.','identifier')
 if new.startswith(('Original-reference ','Target-excluded ')):
  root=new.split(' ',1)[1]; d,u=semantics(n,h.rsplit(' ',1)[0],root);return (('Original-reference model. ' if new.startswith('Original') else 'Target-excluded model. ')+d,u)
 if any(x in new for x in ['LIANA ','CellChat ']):
  method,root=new.split(' ',1);d,u=semantics(n,root,root);return (method+': '+d,u)
 if 'SE' in new or 'standard error' in new.lower():return ('Standard error on the corresponding effect/estimate scale; not a confidence interval.','effect units')
 if 'CI' in new:return ('Recorded confidence bound on the corresponding effect scale; no bounds were reconstructed.','effect units')
 return ('Recorded '+new[0].lower()+new[1:]+' for '+META[n]['row_unit'][0].lower()+META[n]['row_unit'][1:], 'source scale / category')
font=Font(name='Arial',size=10,color='000000');bold=Font(name='Arial',size=10,bold=True,color='000000');fill=PatternFill('solid',fgColor='E9ECEF');white=PatternFill(fill_type=None);line=Side(style='thin',color='AEB4BA')
maps=[];dictionary=[];conversions=[];text_changes=[];structures=[]
for n in range(1,65):
 s=w[f'Table {n:02}'];cols=s.max_column;rows=s.max_row;m=META[n];old_headers=[s.cell(3,j).value for j in range(1,cols+1)];new_headers=[display(str(h),n) for h in old_headers]
 assert len(new_headers)==len(set(new_headers)),(n,new_headers)
 last=rows
 while last>3 and (sum(c.value is not None for c in s[last])==0 or (sum(c.value is not None for c in s[last])==1 and isinstance(s.cell(last,1).value,str))):last-=1
 s.cell(1,1,f'Supplementary Table {n} | {m["title"]}')
 s.cell(2,1,m['contents']+' '+m['row_unit']+' Blank/NA values retain missing or inapplicable status, never numeric zero.')
 for j,(h,new) in enumerate(zip(old_headers,new_headers),1):
  d,u=semantics(n,h,new)
  if n==9 and h=='Latent factor': d+=' Fixed model symbols are retained; AORTIC denotes Aortic and ATHERO denotes Atherosclerotic.'
  s.cell(3,j,new);maps.append([s.title,'Main data',h,new,d,'Semantic clarification' if h!=new else 'Retained accurate name'])
  dictionary.append([n,new,d,u,'Blank/NA: not supplied or not applicable to this record; consult Record type / model status / eligibility. Explicit zero and FALSE are retained.'])
  vals=[s.cell(i,j).value for i in range(4,min(last,30)+1)]
  isid= ('ID' in new or new in ['PMID','Year','Chromosome','Position (GRCh37)','Start','End','Seed','Rank','Matched module','No. current catalog GWS matches'])
  isnumeric=any(isinstance(v,(int,float)) and not isinstance(v,bool) for v in vals) or all(v is None or (isinstance(v,str) and re.fullmatch(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',v)) for v in vals)
  width=14 if isnumeric and not isid else 24
  if any(x in new.lower() for x in ['pathway','description','interpretation','evidence','publication','context','gene list','marker genes','contributions','tissue','cell state']):width=38
  if new in ['Gene','PMID','Year','Rank','Chromosome']:width=15
  if new in ['Cell state','Reference state name','Reference safe label','Parent-qualified state']:width=42
  width=max(width,min(36,len(new)*.64))
  s.column_dimensions[L(j)].width=width
  s.column_dimensions[L(j)].hidden=False
  for i in range(4,last+1):
   c=s.cell(i,j);v=c.value
   # Convert only entirely numeric measures, never IDs or mixed-status strings.
   if not isid and isnumeric and isinstance(v,str) and re.fullmatch(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',v):
    f=float(v)
    if math.isfinite(f) and not (f==0 and any(ch in '123456789' for ch in v.split('e')[0].split('E')[0])):
     c.value=f;conversions.append([s.title,c.coordinate,v,repr(f)])
   # The only data-text changes: genetic-domain labels, never IDs/paths/publication titles.
   if n in [7,8,9] and isinstance(c.value,str) and 'systemic' in c.value.lower() and not any(x in str(h).lower() for x in ['id','model','path','publication']):
    v=c.value;c.value=re.sub('systemic',lambda ma:'Atherosclerotic' if ma.group()[0].isupper() else 'atherosclerotic',v,flags=re.I);text_changes.append([s.title,c.coordinate,v,c.value])
   if n==47 and h=='Scope' and '35 recorded biological units' in str(c.value):
    old=c.value;c.value='35 analysis units; section estimates pooled within units before between-unit meta-analysis';text_changes.append([s.title,c.coordinate,old,c.value])
   if n==51 and h=='Analysis' and c.value=='Pathway enrichment':
    old=c.value;c.value='Pseudobulk pathway enrichment';text_changes.append([s.title,c.coordinate,old,c.value])
   v=c.value;c.font=font;c.fill=white;c.border=Border();c.alignment=Alignment(horizontal='right' if isinstance(v,(int,float)) and not isinstance(v,bool) else 'left',vertical='top',wrap_text=True,shrink_to_fit=False)
   if isinstance(v,(int,float)) and not isinstance(v,bool):
    if new in ['PMID','Year']:c.number_format='0'
    elif v==int(v) and not ('P'==new or 'FDR' in new or 'PP.' in new):c.number_format='#,##0'
    elif v!=0 and (abs(v)<.001 or abs(v)>=1e6):c.number_format='0.000E+00'
    else:c.number_format='0.0000'
  for i in [1,2,3]:s.cell(i,j).alignment=Alignment(horizontal='left',vertical='center',wrap_text=True,shrink_to_fit=False)
  s.cell(3,j).font=bold;s.cell(3,j).fill=fill;s.cell(3,j).border=Border(bottom=line)
 s.cell(1,1).font=Font(name='Arial',size=12,bold=True);s.cell(2,1).font=font
 s.row_dimensions[1].height=30;s.row_dimensions[2].height=58;s.row_dimensions[3].height=max(42,min(76,max(len(h) for h in new_headers)*1.1))
 for i in range(4,rows+1):
  s.row_dimensions[i].hidden=False
  # Permit long prose to remain fully visible; no shrink-to-fit or data truncation.
  lines=max([max(1,math.ceil(len(str(s.cell(i,j).value or ''))/(s.column_dimensions[L(j)].width*.9))) for j in range(1,cols+1)])
  s.row_dimensions[i].height=min(409,max(15,lines*12))
 for i in range(last+1,rows+1):
  c=s.cell(i,1);c.font=font;c.alignment=Alignment(wrap_text=True,vertical='top',shrink_to_fit=False)
  if isinstance(c.value,str) and n in [7,8,9]:
   old=c.value;c.value=re.sub('systemic','atherosclerotic',old,flags=re.I)
   if old!=c.value:text_changes.append([s.title,c.coordinate,old,c.value])
 s.freeze_panes='C4' if cols>12 else 'B4';s.auto_filter.ref=f'A3:{L(cols)}{last}';s.sheet_view.showGridLines=False
 s.print_title_rows='1:3';s.print_options.horizontalCentered=False;s.page_setup.orientation='landscape';s.page_setup.paperSize=s.PAPERSIZE_A4;s.page_setup.scale=100;s.sheet_properties.pageSetUpPr.fitToPage=False;s.page_setup.fitToWidth=None;s.page_setup.fitToHeight=None;s.page_margins=PageMargins(left=.25,right=.25,top=.35,bottom=.35,header=.15,footer=.15)
 s.print_area=f'A1:{L(cols)}{s.max_row}'
 structures.append({'sheet':s.title,'data_start':4,'data_end':last,'rows':last-3,'columns':cols,'title':m['title']})
 print('formatted',n,flush=True)
info=w['Information'];info.delete_rows(1,info.max_row);info.append(['Table','Title','Contents','Row unit','Worksheet'])
for n,m in META.items():
 info.append([n,m['title'],m['contents'],m['row_unit'],f'Table {n:02}']);info.cell(n+1,5).hyperlink=f"#'Table {n:02}'!A1"
if 'Data dictionary' in w:del w['Data dictionary']
s=w.create_sheet('Data dictionary');s.append(['Table','Display column','Definition','Unit/scale','Missing-value meaning'])
for r in dictionary:s.append(r)
for s,widths in [(info,[9,42,80,65,16]),(w['Data dictionary'],[9,42,90,25,72]),(w['Abbreviations'],[24,80])]:
 for row in s:
  for c in row:c.font=bold if c.row==1 else font;c.fill=fill if c.row==1 else white;c.alignment=Alignment(horizontal='left',vertical='top',wrap_text=True,shrink_to_fit=False)
 for j,wd in enumerate(widths,1):s.column_dimensions[L(j)].width=wd
 s.freeze_panes='B2';s.auto_filter.ref=s.dimensions;s.sheet_view.showGridLines=False
 for i in range(1,s.max_row+1):s.row_dimensions[i].height=48 if s==info else 60
 s.row_dimensions[1].height=30
 s.print_title_rows='1:1';s.page_setup.scale=100;s.sheet_properties.pageSetUpPr.fitToPage=False
w.save(S/src.name)
for name,header,rows in [
 ('column_name_map.tsv',['Worksheet','Block','Original column','Display column','Definition','Change reason'],maps),
 ('numeric_text_conversions.tsv',['Worksheet','Cell','Original text','Numeric value'],conversions),
 ('label_text_changes.tsv',['Worksheet','Cell','Original','Updated'],text_changes)]:
 with (O/name).open('w') as f:ww=csv.writer(f,delimiter='\t');ww.writerow(header);ww.writerows(rows)
(O/'table_metadata.json').write_text(json.dumps(META,ensure_ascii=False,indent=2));(O/'formatted_structure.json').write_text(json.dumps(structures,indent=2))
# Read-back comparison against every original scientific cell. No skipped failures or filtered rows.
a=load_workbook(B/src.name,read_only=True);b=load_workbook(S/src.name,read_only=True)
conv={(s,c):v for s,c,t,v in conversions};tc={(s,c):v for s,c,t,v in text_changes};counts=0
for n in range(1,65):
 x=a[f'Table {n:02}'];y=b[x.title];assert(x.max_row,x.max_column)==(y.max_row,y.max_column)
 for ri,(rx,ry) in enumerate(zip(x.iter_rows(min_row=4),y.iter_rows(min_row=4)),4):
  for cj,(cx,cy) in enumerate(zip(rx,ry),1):
   key=(x.title,f"{L(cj)}{ri}");old=cx.value;new=cy.value
   if key in conv:
    assert float(old)==new or math.isclose(float(old),new,rel_tol=2e-15,abs_tol=0),(key,old,new)
   elif key in tc:assert new==tc[key]
   elif isinstance(old,float):assert old==new or math.isclose(old,new,rel_tol=2e-15,abs_tol=0),(key,old,new)
   else:assert old==new or (old=="" and new is None),(key,old,new)
   counts+=1
assert b.sheetnames[:66]==a.sheetnames
(O/'workbook_validation.json').write_text(json.dumps({'numbered_tables':64,'scientific_cells_compared':counts,'data_rows_and_columns_unchanged':True,'row_order_unchanged':True,'numeric_text_conversions':len(conv),'documented_domain_label_updates':len(tc),'raw_numeric_precision':'openpyxl round-trip; equality or relative tolerance 2e-15, no rounding or rescaling','formulas_source':0,'missing_and_boolean_values_unchanged':True,'blank_text_serialization':'Original empty strings may reload as None; no non-empty missing/status token changed','all_data_visible':True,'all_headers_nonempty_unique':True},indent=2))
print('VALIDATED',counts,'cells; numeric texts:',len(conv))
