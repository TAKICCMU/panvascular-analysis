import re
Q=lambda x:'{http://www.w3.org/2000/svg}'+x
find=lambda r,id:next(e for e in r.iter() if e.get('id')==id)
def apply(root,name,config):
 changes=[]
 if name=='Figure5':
  for id in ['FigureYa6_native_forest','legend_2']:
   for el in find(root,id).iter():
    tag=el.tag.split('}')[-1];old=el.get('style','')
    if tag not in ['circle','polyline','path','use'] or '#d9dde1' not in old.lower():continue
    role='NO_SUPPORT_CI' if tag=='polyline' else 'NO_SUPPORT_POINT';new=re.sub('#d9dde1',config[role],old,flags=re.I)
    assert 'opacity' not in old
    el.set('style',new);changes.append((el,old,new,id,role))
 if name in [f'ExtendedData_Figure{i}' for i in [2,3,4,5]]:
  n=int(name.split('Figure')[-1])
  for ai,li in ([(1,2),(4,5)] if n==4 else [(1,2)]):
   ax=find(root,f'newED{n}_axes_{ai}');main=next(e for e in ax if 'PathCollection' in e.get('id',''))
   keys=[e for e in find(root,f'newED{n}_axes_{li}') if 'PathCollection' in e.get('id','')]
   for group in [main]+keys:
    for e in group.iter():
     old=e.get('style','')
     if re.search(r'stroke:\s*#555555',old,re.I):
      new=re.sub(r'stroke:\s*#555555','stroke: none',old,flags=re.I);new=re.sub(r'stroke-width:\s*[\d.]+','stroke-width: 0',new);e.set('style',new);changes.append((e,old,new,group.get('id'),'marker outline only'))
 return changes
