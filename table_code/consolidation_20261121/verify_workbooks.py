"""Read-only verification of editorially consolidated workbooks."""
from pathlib import Path
import json,math,hashlib,collections,sys,gzip
from openpyxl import load_workbook
O=Path(__file__).parent
plan=json.loads((O/'CONSOLIDATION_PLAN.json').read_text() if (O/'CONSOLIDATION_PLAN.json').exists() else gzip.decompress((O/'CONSOLIDATION_PLAN.json.gz').read_bytes()).decode())
all_report=json.loads((O/'WORKBOOK_VALUE_VERIFICATION.json').read_text()) if (O/'WORKBOOK_VALUE_VERIFICATION.json').exists() else {};problems=[]
for mode,filename in [('core','Supplementary_Tables.xlsx'),('data','Supplementary_Data.xlsx')]:
    if len(sys.argv)>1 and sys.argv[1]!=mode:continue
    path=O/filename
    if not path.exists():continue
    w=load_workbook(path,read_only=True,data_only=False)
    lay={x['sheet']:x for x in json.loads((O/f'LAYOUT_{mode}.json').read_text())}
    report=[]
    expected_sheets=['Information']+[s['sheet'] for s in plan[mode]]+['Data dictionary','Abbreviations']
    assert w.sheetnames==expected_sheets,(w.sheetnames,expected_sheets)
    for spec in plan[mode]:
        sheet=w[spec['sheet']];tracks=lay[spec['sheet']]['blocks']
        bounds={b['header_row']:('header',b,spec['blocks'][i]) for i,b in enumerate(tracks)}
        count=collections.Counter();rounding=[];checks=0
        for i,track in enumerate(tracks):
            b=spec['blocks'][i]
            rows=sheet.iter_rows(min_row=track['header_row'],max_row=track['end_row'],max_col=len(b['h']))
            actual=next(rows)
            assert [c.value for c in actual]==b['h'],(spec['sheet'],'headers')
            for ri,(cells,values,fmts) in enumerate(zip(rows,b['r'],b['cell_formats']),track['start_row']):
                for ci,(c,v,nf) in enumerate(zip(cells,values,fmts),1):
                    checks+=1;count[c.data_type]+=1
                    if c.value!=v or type(c.value)==bool and type(v)!=bool:
                        numeric=all(isinstance(x,(int,float)) and not isinstance(x,bool) for x in (c.value,v))
                        if numeric and math.isclose(c.value,v,rel_tol=5e-15,abs_tol=0):rounding.append([ri,ci,v,c.value])
                        else:problems.append([mode,spec['sheet'],ri,ci,v,c.value])
                    if isinstance(v,(int,float)) and not isinstance(v,bool) and c.data_type!='n':problems.append([spec['sheet'],c.coordinate,'numeric type lost',c.data_type])
                    if isinstance(v,bool) and c.data_type!='b':problems.append([spec['sheet'],c.coordinate,'boolean type lost',c.data_type])
                    if isinstance(v,str) and c.data_type not in ('s','inlineStr'):problems.append([spec['sheet'],c.coordinate,'text type lost',c.data_type])
                    if c.number_format!=nf:problems.append([spec['sheet'],c.coordinate,'number format',nf,c.number_format])
        report.append({'sheet':spec['sheet'],'source_tables':[b['old'] for b in tracks],'data_rows':sum(b['rows'] for b in tracks),'data_cells_checked':checks,'cell_types':dict(count),'floating_point_serialization_differences':rounding})
    # Native panes, filter tables and merges are inspected in the saved XML by openpyxl.
    w.close();native=load_workbook(path,read_only=False,data_only=False)
    for spec in plan[mode]:
        s=native[spec['sheet']];tracks=lay[s.title]['blocks']
        assert len(s.tables)==len(tracks),(s.title,'filters')
        assert s.freeze_panes is not None,(s.title,'frozen headers')
        for mr in s.merged_cells.ranges:
            assert not any(mr.min_row<=b['end_row'] and mr.max_row>=b['start_row'] for b in tracks),(s.title,'merged data cell',str(mr))
        assert not any(c.data_type=='f' for row in s for c in row),'No formulas should be introduced'
    all_report[mode]={'file':filename,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'numbered_tables':len(plan[mode]) if mode=='core' else 0,'detail_sheets':len(plan[mode]) if mode=='data' else 0,'sheets':report,'native_filters_and_frozen_headers':True,'no_merged_data_cells':True}
    del native
all_report['problems']=problems
(O/'WORKBOOK_VALUE_VERIFICATION.json').write_text(json.dumps(all_report,ensure_ascii=False,indent=2))
assert not problems,problems[:20]
print('VERIFIED',[(mode,sum(x['data_cells_checked'] for x in r['sheets'])) for mode,r in all_report.items() if mode!='problems'])
