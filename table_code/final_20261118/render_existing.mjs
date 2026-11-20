import fs from 'node:fs/promises';
import { FileBlob, SpreadsheetFile } from '@oai/artifact-tool';
const root='/Users/taki/Desktop/泛血管稿件_批注修订_20261116';
const out=root+'/图件编辑资料/补充表核查_20261118/previews';
await fs.mkdir(out,{recursive:true});
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(root+'/Supplementary_Tables.xlsx'));
console.log('Imported current workbook; read-only render, no export.');
for (const [sheetName,range,name] of [
 ['Information','A1:E8','index'], ['Table 05','A1:H8','markers'],
 ['Table 17','A1:I7','loci'], ['Table 24','A1:S8','candidates'],
 ['Table 36','A1:Q10','spatial'], ['Table 60','A1:C15','software']
]) {
 try { const b=await wb.render({sheetName,range,scale:1.4,format:'png'});
 await fs.writeFile(out+'/'+name+'.png',new Uint8Array(await b.arrayBuffer()));
 console.log('Rendered',sheetName,range);
 } catch(e) {console.log('Render failed',sheetName,String(e));}
}
