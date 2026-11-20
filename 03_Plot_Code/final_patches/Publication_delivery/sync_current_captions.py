"""Synchronise only adopted caption text in the three existing DOCX files.

Only authorised w:t text values in word/document.xml are changed. The rest of
each OOXML part is preserved byte for byte, including formatting and fields.
"""
from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
from xml.parsers import expat
from xml.sax.saxutils import escape
import copy, hashlib, json, os, re, shutil, stat

P=Path('/Users/taki/Desktop/泛血管稿件_批注修订_20261116')
R=P/'图件编辑资料/全图核对与交付_20261118'
B=R/'backup_before_caption_sync'
B.mkdir(exist_ok=True)
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
WT='{'+NS['w']+'}t'
sha=lambda x:hashlib.sha256(x).hexdigest()
files=['Manuscript.docx','Figure_Legends.docx','Supplementary_Information.docx']
source=P/'Figure/Figure_Legends.docx'
with ZipFile(source) as z:
    masterxml=z.read('word/document.xml')
masterroot=E.fromstring(masterxml)
def para_text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def caption_key(t):
    m=re.match(r'^(Fig\. \d+|Extended Data Fig\. \d+)\s',t)
    return m.group(1) if m else None
master={caption_key(para_text(p)):para_text(p) for p in masterroot.xpath('//w:p',namespaces=NS) if caption_key(para_text(p))}
assert len(master)==14
edits=json.loads((R/'caption_changes.json').read_text())

def xml_text_spans(raw):
    parser=expat.ParserCreate();spans=[];stack=[]
    def start(name,attrs):
        if name=='w:t':
            a=parser.CurrentByteIndex;b=raw.index(b'>',a)+1
            stack.append((a,b,raw[b-2:b]==b'/>'))
    def end(name):
        if name=='w:t':
            a,b,empty=stack.pop()
            if empty:spans.append((a,b,b,b,True))
            else:
                c=parser.CurrentByteIndex;e=raw.index(b'>',c)+1
                spans.append((a,b,c,e,False))
    parser.StartElementHandler=start;parser.EndElementHandler=end
    parser.Parse(raw,True)
    return sorted(spans)

def patch_phrase(p,old,new):
    s=para_text(p);a=s.index(old);b=a+len(old)
    pos=0;nodes=[]
    for t in p.xpath('.//w:t',namespaces=NS):
        text=t.text or '';nodes.append((t,pos,pos+len(text)));pos+=len(text)
    hits=[x for x in nodes if x[1]<b and x[2]>a]
    assert hits
    for i,(node,start,end) in enumerate(hits):
        t=node.text or ''
        prefix=t[:max(0,a-start)] if i==0 else ''
        suffix=t[max(0,b-start):] if i==len(hits)-1 else ''
        node.text=prefix+(new if i==0 else '')+suffix
    assert para_text(p)==s[:a]+new+s[b:]

protected={str(f.relative_to(P)):sha(f.read_bytes()) for base in [P/'Figure',P/'Figures',P/'Extended_Data',P/'Source_Data_Figures'] for f in base.rglob('*') if f.is_file()}
protected['Supplementary_Tables.xlsx']=sha((P/'Supplementary_Tables.xlsx').read_bytes())
staged=[];report=[]
for name in files:
    file=P/name;original=file.read_bytes()
    with ZipFile(file) as z:
        entries=z.infolist();parts={i.filename:z.read(i.filename) for i in entries};comment=z.comment
    raw=parts['word/document.xml'];root=E.fromstring(raw);oldroot=copy.deepcopy(root)
    paragraphs=root.xpath('//w:p',namespaces=NS);before_paras=oldroot.xpath('//w:p',namespaces=NS)
    changed=[]
    for i,p in enumerate(paragraphs):
        t=para_text(p);key=caption_key(t)
        if not key:continue
        assert key in master
        if t==master[key]:continue
        expected=t;applied=[]
        for edit in edits:
            old,new=edit['original'],edit['revised']
            if old in expected and new not in expected:
                assert expected.count(old)==1
                expected=expected.replace(old,new,1)
                patch_phrase(p,old,new)
                applied.append({'original':old,'adopted':new})
        assert expected==master[key],(name,key,'Unapproved caption differences')
        assert para_text(p)==master[key]
        changed.append({'paragraph_index':i,'caption':key,'phrases':applied})
    expected_count=1 if name=='Supplementary_Information.docx' else 5
    assert len(changed) in [0,expected_count],(name,len(changed))

    oldnodes=oldroot.xpath('//w:t',namespaces=NS);newnodes=root.xpath('//w:t',namespaces=NS)
    assert len(oldnodes)==len(newnodes)
    spans=xml_text_spans(raw);assert len(spans)==len(oldnodes)
    replacements=[]
    for i,(a,b) in enumerate(zip(oldnodes,newnodes)):
        assert a.attrib==b.attrib
        if a.text==b.text:continue
        start,begin,end,close,empty=spans[i]
        assert not empty,'No new text element is required for this bounded change'
        value=escape(b.text or '').encode('utf-8')
        replacements.append((begin,end,value))
    newraw=raw
    for a,b,value in reversed(replacements):newraw=newraw[:a]+value+newraw[b:]
    parsed=E.fromstring(newraw)
    assert E.tostring(parsed)==E.tostring(root)
    # All XML structure, properties, fields and non-text contents are identical.
    aa=copy.deepcopy(oldroot);bb=copy.deepcopy(parsed)
    for x in aa.xpath('//w:t',namespaces=NS):x.text='__TEXT__'
    for x in bb.xpath('//w:t',namespaces=NS):x.text='__TEXT__'
    assert E.tostring(aa)==E.tostring(bb)
    modified_indices={x['paragraph_index'] for x in changed}
    for i,(a,b) in enumerate(zip(before_paras,parsed.xpath('//w:p',namespaces=NS))):
        if i not in modified_indices:assert E.tostring(a)==E.tostring(b),(name,i)
    # Save to staging first; preserve each existing ZIP member's metadata.
    tmp=R/(name+'.staging')
    with ZipFile(tmp,'w') as z:
        z.comment=comment
        for entry in entries:
            z.writestr(copy.copy(entry),newraw if entry.filename=='word/document.xml' else parts[entry.filename])
    with ZipFile(tmp) as z:
        assert z.testzip() is None
        assert z.namelist()==[x.filename for x in entries]
        for part,content in parts.items():
            if part!='word/document.xml':assert z.read(part)==content,(name,part)
        captions={caption_key(para_text(p)):para_text(p) for p in parsed.xpath('//w:p',namespaces=NS) if caption_key(para_text(p))}
        assert all(value==master[key] for key,value in captions.items())
    staged.append((file,tmp,sha(original)))
    report.append({'file':name,'before_sha256':sha(original),'after_sha256':sha(tmp.read_bytes()),'changed_captions':changed,'caption_count_checked':len(captions),'changed_ooxml_parts':['word/document.xml'] if changed else [],'unchanged_other_package_parts':len(parts)-1,'formatting_and_field_structure_unchanged':True,'all_other_paragraphs_identical':True})

# Do not overwrite a document if it changed while preparing this update.
for file,tmp,h in staged:assert sha(file.read_bytes())==h,('Concurrent edit',str(file))
for file,tmp,h in staged:
    backup=B/file.name
    if not backup.exists():shutil.copy2(file,backup)
    os.chmod(tmp,stat.S_IMODE(file.stat().st_mode))
    os.replace(tmp,file)
for x in report:assert sha((P/x['file']).read_bytes())==x['after_sha256']
for rel,h in protected.items():assert sha((P/rel).read_bytes())==h,rel
result={'source':str(source),'source_sha256':sha(source.read_bytes()),'status':'SYNCED_IN_PLACE','documents':report,'protected_files_unchanged':len(protected),'scientific_or_figure_edits':False,'new_desktop_version_created':False}
(R/'caption_sync_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps({'status':result['status'],'changes':{x['file']:[a['caption'] for a in x['changed_captions']] for x in report},'all_other_document_parts_unchanged':True,'protected_files_unchanged':len(protected)},ensure_ascii=False))
