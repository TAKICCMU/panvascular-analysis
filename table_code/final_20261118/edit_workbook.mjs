import fs from 'node:fs/promises';
import {FileBlob, SpreadsheetFile} from '@oai/artifact-tool';
const base='/Users/taki/Desktop/泛血管稿件_批注修订_20261116';
const out=base+'/图件编辑资料/补充表核查_20261118';
const spec=JSON.parse(await fs.readFile(out+'/correction_spec.json','utf8'));
console.log('Importing existing workbook');
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(base+'/Supplementary_Tables.xlsx'));
for(const ch of spec.changes)wb.worksheets.getItem(ch.sheet).getRange(ch.cell).values=[[ch.after]];
for(const [name,rows,last] of [['Table 59',spec.families,'D'],['Table 60',spec.versions,'C']]){
  const s=wb.worksheets.getItem(name);
  for(let i=0;i<rows.length;i++){
    if(i+3>(name==='Table 59'?7:15))s.getRange(`A${i+3}:${last}${i+3}`).copyFrom(s.getRange(`A${i%2+3}:${last}${i%2+3}`),'all');
  }
  s.getRange(`A3:${last}${rows.length+2}`).values=rows;
  s.getRange(`A3:${last}${rows.length+2}`).format.wrapText=true;
  s.getRange(`A3:${last}${rows.length+2}`).format.verticalAlignment='center';
  s.getRange(`A3:${last}${rows.length+2}`).format.rowHeight=name==='Table 59'?66:42;
  s.getRange(`A2:${last}2`).format.rowHeight=34;
  if(name==='Table 59'){
    for(const [col,width] of [['A',32],['B',49],['C',33],['D',88]])s.getRange(`${col}2:${col}${rows.length+2}`).format.columnWidth=width;
  }else{
    for(const [col,width] of [['A',35],['B',22],['C',100]])s.getRange(`${col}2:${col}${rows.length+2}`).format.columnWidth=width;
    s.getRange(`B3:B${rows.length+2}`).setNumberFormat('General');
  }
}
for(const n of spec.numeric_cells)wb.worksheets.getItem(n.sheet).getRange(n.cell).setNumberFormat('#,##0');
const t=wb.worksheets.getItem('Table 24');
t.getRange('A3:BC31').format.rowHeight=spec.table24_row_height;
for(const c of spec.table24_singleline_columns){
  t.getRange(`${c}3:${c}31`).format.wrapText=false;
  t.getRange(`${c}3:${c}31`).format.verticalAlignment='center';
}
wb.worksheets.getItem('Information').getRange('A9:E9').format.rowHeight=72;
wb.worksheets.getItem('Information').getRange('A60:E61').format.rowHeight=54;
wb.worksheets.getItem('Data dictionary').getRange('A114:E129').format.rowHeight=86;
wb.worksheets.getItem('Data dictionary').getRange('A904:E910').format.rowHeight=60;
wb.recalculate();
await fs.mkdir(out+'/staging/previews',{recursive:true});
console.log((await wb.inspect({kind:'table',range:"'Table 59'!A2:D5",tableMaxRows:4,tableMaxCols:4,maxChars:2400})).ndjson);
console.log((await wb.inspect({kind:'table',range:"'Table 60'!A53:C65",tableMaxRows:13,tableMaxCols:3,maxChars:2600})).ndjson);
for(const [name,sheetName,range] of [
 ['table59','Table 59','A1:D7'],['table59_families','Table 59','A14:D20'],
 ['table60','Table 60','A1:C8'],['table60_seeds','Table 60','A53:C65'],
 ['table24','Table 24','A1:L7'],['table24_long','Table 24','AO2:AS5'],
 ['table17','Table 17','A1:F6']]){
 const b=await wb.render({sheetName,range,scale:1.25});
 await fs.writeFile(out+`/staging/previews/${name}.png`,new Uint8Array(await b.arrayBuffer()));
}
console.log('Exporting workbook');
const x=await SpreadsheetFile.exportXlsx(wb);await x.save(out+'/staging/artifact_edited.xlsx');
console.log('Export completed');
