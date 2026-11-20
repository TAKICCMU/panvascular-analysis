exec(open(str(__import__('pathlib').Path(__file__).with_name('inventory.py'))).read().split('for d in')[0])
import re,collections
q=lambda s:'{http://www.w3.org/2000/svg}'+s
for name in ['Figure3','Figure4','Figure5','ExtendedData_Figure2','ExtendedData_Figure6','ExtendedData_Figure8','ExtendedData_Figure9','ExtendedData_Figure10']:
 r=E.parse(W/'backup'/f'{name}.svg').getroot();print('\n',name)
 for g in r.iter(q('g')):
  id=g.get('id','')
  if 'axes_' in id or any(x in id for x in ['compact_forest','content','key','matrix','legend']):
   c=collections.Counter(x.lower() for e in g.iter() if e.tag!=q('text') for k,v in e.attrib.items() if k in ['style','fill','stroke'] for x in re.findall(r'#[0-9a-fA-F]{6}',v));print(id,dict(c))
