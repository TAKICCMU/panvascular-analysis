import fs from 'node:fs/promises';
import {gunzipSync} from 'node:zlib';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const O=decodeURIComponent(new URL('.',import.meta.url).pathname);
let planText;try{planText=await fs.readFile(O+'CONSOLIDATION_PLAN.json','utf8');}catch(e){if(e.code!=='ENOENT')throw e;planText=gunzipSync(await fs.readFile(O+'CONSOLIDATION_PLAN.json.gz')).toString('utf8');}
const plan=JSON.parse(planText);
const mode=process.argv[2]||'core';const specs=plan[mode];
const wb=Workbook.create();const layout=[];
function col(n){let s='';for(n++;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;}
function put(sheet,row,data){const m=data.length,n=data[0].length;sheet.getRange(`A${row}:${col(n-1)}${row+m-1}`).values=data;}
function base(sheet,last,n){let r=sheet.getRange(`A1:${col(n-1)}${last}`);r.format.font={name:'Arial',size:10,color:'#111111'};r.format.verticalAlignment='top';sheet.showGridLines=false;}
function header(sheet,row,n){let r=sheet.getRange(`A${row}:${col(n-1)}${row}`);r.format={fill:'#E9ECEF',font:{name:'Arial',size:10,bold:true,color:'#111111'},wrapText:true};r.format.rowHeight=64;}
function smallSheet(name,data,widths){const s=wb.worksheets.add(name);put(s,1,data);base(s,data.length,data[0].length);header(s,1,data[0].length);s.getRange(`A2:${col(data[0].length-1)}${data.length}`).format.wrapText=true;for(let r=1;r<data.length;r++){const lines=Math.max(...data[r].map((v,c)=>String(v??'').split('\n').reduce((a,x)=>a+Math.max(1,Math.ceil(x.length/(widths[c]*0.95))),0)));s.getRange(`A${r+1}:${col(data[0].length-1)}${r+1}`).format.rowHeight=Math.max(28,15*lines+8);}for(let i=0;i<widths.length;i++)s.getRange(`${col(i)}1:${col(i)}${data.length}`).format.columnWidth=widths[i];s.tables.add(`A1:${col(data[0].length-1)}${data.length}`,true,name.replace(/\W/g,'')+'Index').style='TableStyleLight1';s.freezePanes.freezeRows(1);return s;}
const index=[['Table / data sheet','Title','Contents / sections','Row units','Worksheet']];
for(const s of specs)index.push([s.number??s.sheet,s.title,s.blocks.map(b=>(b.block?b.block+'. ':'')+b.title).join('\n'),s.blocks.map(b=>(b.block?b.block+': ':'')+b.row_unit).join('\n'),s.sheet]);
smallSheet('Information',index,[19,54,72,62,36]);
for(const s of specs){
 const sheet=wb.worksheets.add(s.sheet);const n=Math.max(...s.blocks.map(b=>b.h.length));let row=1;const track={sheet:s.sheet,number:s.number,blocks:[]};
 sheet.getRange('A1').values=[[(s.number?`Supplementary Table ${s.number} | `:'')+s.title]];
 const titleEnd=Math.min(n,8);sheet.getRange(`A1:${col(titleEnd-1)}1`).merge();sheet.getRange(`A1:${col(titleEnd-1)}1`).format={font:{name:'Arial',size:12,bold:true,color:'#111111'},wrapText:true};sheet.getRange(`A1:${col(n-1)}1`).format.rowHeight=44;
 row=2;
 for(const b of s.blocks){
  console.log('WRITING',s.sheet,b.old,b.r.length);
  if(s.blocks.length>1){sheet.getRange(`A${row}`).values=[[b.block+'. '+b.title]];sheet.getRange(`A${row}:${col(Math.min(n,8)-1)}${row}`).merge();sheet.getRange(`A${row}:${col(n-1)}${row}`).format={fill:'#F3F5F7',font:{name:'Arial',size:11,bold:true,color:'#111111'},wrapText:true};sheet.getRange(`A${row}:${col(n-1)}${row}`).format.rowHeight=34;row++;}
  const hr=row;put(sheet,row,[b.h]);header(sheet,row,b.h.length);row++;
  for(let i=0;i<b.r.length;i+=3000){const chunk=b.r.slice(i,i+3000);sheet.getRange(`A${row+i}:${col(b.h.length-1)}${row+i+chunk.length-1}`).values=chunk;}
  const dr=sheet.getRange(`A${row}:${col(b.h.length-1)}${row+b.r.length-1}`);
  dr.format.font={name:'Arial',size:10,color:'#111111'};dr.format.verticalAlignment='top';dr.format.wrapText=true;dr.format.rowHeight=Math.min(84,Math.max(22,b.row_height||24));
  dr.format.numberFormat=b.cell_formats;
  sheet.tables.add(`A${hr}:${col(b.h.length-1)}${row+b.r.length-1}`,true,`Data_${mode}_${b.old}`).style='TableStyleLight1';
  track.blocks.push({old:b.old,block:b.block,header_row:hr,start_row:row,end_row:row+b.r.length-1,rows:b.r.length,cols:b.h.length});
  row+=b.r.length+2;
 }
 const last=row-3;
 for(let i=0;i<n;i++){const widths=s.blocks.filter(b=>i<b.h.length).map(b=>b.widths[i]||18);sheet.getRange(`${col(i)}1:${col(i)}${last}`).format.columnWidth=Math.min(65,Math.max(13,...widths));}
 sheet.showGridLines=false;sheet.freezePanes.freezeRows(track.blocks[0].header_row);sheet.freezePanes.freezeColumns(s.number===14?2:1);
 layout.push(track);console.log('BUILT',s.sheet,track.blocks.map(b=>`${b.old}:${b.rows}`).join(','));
}
const dd=[['Worksheet','Section','Display column','Definition','Unit/scale','Missing-value meaning'],...specs.flatMap(s=>s.dictionary)];smallSheet('Data dictionary',dd,[35,10,44,88,28,75]);
const abr=plan.abbreviations.slice(1);smallSheet('Abbreviations',abr,[24,85]);
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:`${specs[0].sheet}!A1:C5`,include:'values,formulas',tableMaxRows:5,tableMaxCols:3})).ndjson);
await fs.writeFile(O+`LAYOUT_${mode}.json`,JSON.stringify(layout,null,2));
console.log('EXPORTING');
const output=await SpreadsheetFile.exportXlsx(wb);await output.save(O+(mode==='core'?'Supplementary_Tables.xlsx':'Supplementary_Data.xlsx'));
console.log('EXPORTED');
for(const sn of ['Information','Data dictionary','Abbreviations']){const p=await wb.render({sheetName:sn,range:sn==='Data dictionary'?'A1:F7':sn==='Information'?'A1:E6':'A1:B9',scale:1,format:'png'});await fs.writeFile(O+`previews/${mode}_${sn.replace(/ /g,'_')}.png`,new Uint8Array(await p.arrayBuffer()));}
// Uniform bounded views of every sheet's first data block; second blocks are also rendered for the four most complex tables.
for(const t of layout){
 const ranges=[{label:'top',start:1,end:Math.min(t.blocks[0].end_row,t.blocks[0].start_row+4)}];
 if([6,9,14,20,24].includes(t.number))for(const b of t.blocks.slice(1))ranges.push({label:b.block,start:b.header_row-1,end:Math.min(b.end_row,b.start_row+3)});
 for(const r of ranges){const endCol=Math.min(5,t.blocks[0].cols);const p=await wb.render({sheetName:t.sheet,range:`A${r.start}:${col(endCol-1)}${r.end}`,scale:1,format:'png'});await fs.writeFile(O+`previews/${mode}_${t.sheet.replace(/ /g,'_')}_${r.label}.png`,new Uint8Array(await p.arrayBuffer()));}
}
console.log('RENDERED',mode);
