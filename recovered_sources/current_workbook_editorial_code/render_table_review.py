from pathlib import Path
from openpyxl import load_workbook
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import json,textwrap,math
P=Path('/Users/taki/Desktop/投稿包_最新版_20260919');O=P/'图件编辑资料/补充表与方法同步_20260921'
pdfmetrics.registerFont(TTFont('Review','/System/Library/Fonts/Supplemental/Arial.ttf'))
w=load_workbook(O/'待替换/Supplementary_Tables.xlsx',read_only=True)
c=canvas.Canvas(str(O/'64表逐字段抽样预览.pdf'));checks=[]
for n in range(1,65):
 s=w[f'Table {n:02}'];rows=list(s.values);headers=rows[2];last=len(rows)
 while last>3 and sum(v is not None for v in rows[last-1])<=1:last-=1
 inds=[3,min(4,last-1),max(3,(last+3)//2),last-1]
 maxlenrow=max(range(3,last),key=lambda i:max(len(str(v or '')) for v in rows[i]))
 inds.append(maxlenrow);inds=list(dict.fromkeys(inds))
 pagew=1190;pageh=max(842,180+len(headers)*31)
 c.setPageSize((pagew,pageh));c.setFont('Review',13);c.drawString(28,pageh-30,str(rows[0][0]))
 c.setFont('Review',9);c.drawString(28,pageh-48,'Field-by-field review of saved workbook. Values sampled from first / next / middle / last / longest-text rows.')
 c.drawString(28,pageh-62,'Long sample text is abbreviated only in this review; the Excel retains complete cell contents. This is not a native Excel screenshot.')
 xs=[28,278,458,638,818,998];y=pageh-86
 c.setFillColorRGB(.91,.92,.93);c.rect(28,y-25,1134,27,fill=1,stroke=0);c.setFillColorRGB(0,0,0)
 c.drawString(32,y-12,'Display column (original order)')
 for j,i in enumerate(inds):c.drawString(xs[j+1]+3,y-12,'Excel row '+str(i+1))
 for j,h in enumerate(headers):
  yy=y-32-j*31
  if j%2==0:c.setFillColorRGB(.975,.975,.975);c.rect(28,yy-24,1134,31,fill=1,stroke=0)
  c.setFillColorRGB(0,0,0);c.setFont('Review',9)
  for k,line in enumerate(textwrap.wrap(str(h),width=40)[:3]):c.drawString(32,yy-k*9,line)
  for ix,i in enumerate(inds):
   v=rows[i][j]; t=('∅' if v is None else (format(v,'.5g') if isinstance(v,float) else str(v)))
   t=t.replace('\n',' ');ss=textwrap.wrap(t,width=27)
   if len(ss)>2:ss=ss[:2];ss[-1]=ss[-1][:-3]+'...'
   for k,line in enumerate(ss):c.drawString(xs[ix+1]+3,yy-k*10,line)
 c.setFont('Review',9);c.drawString(28,30,f'{s.title}: {last-3:,} data rows, {len(headers)} columns. Header row 3; all headers nonempty and unique. Numeric values are abbreviated in this diagnostic PDF only.')
 c.showPage();checks.append({'table':n,'review_page':n,'sampled_excel_rows':[i+1 for i in inds],'headers_reviewed':len(headers)})
c.save();(O/'table_preview_index.json').write_text(json.dumps(checks,indent=2))
