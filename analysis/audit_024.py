"""Inventory literal 0.24 variants, separating stored values from rendered SVG paths."""
import csv
import json
from pathlib import Path
import re
import subprocess
from zipfile import ZipFile
from xml.etree import ElementTree as E
from openpyxl import load_workbook

ROOT=Path(__file__).resolve().parents[1]
PATTERN=re.compile(r'0?\.24\d*|\b24\s*%')


def snappy(block):
    i=0;shift=0;size=0
    while True:
        b=block[i];i+=1;size|=(b&127)<<shift;shift+=7
        if b<128:break
    out=bytearray()
    while i<len(block):
        tag=block[i];i+=1;kind=tag&3
        if kind==0:
            n=tag>>2
            if n<60:n+=1
            else:
                k=n-59;n=int.from_bytes(block[i:i+k],'little')+1;i+=k
            out.extend(block[i:i+n]);i+=n
        else:
            if kind==1:
                n=4+((tag>>2)&7);offset=((tag&224)<<3)|block[i];i+=1
            elif kind==2:n=1+(tag>>2);offset=int.from_bytes(block[i:i+2],'little');i+=2
            else:n=1+(tag>>2);offset=int.from_bytes(block[i:i+4],'little');i+=4
            for _ in range(n):out.append(out[-offset])
    assert len(out)==size
    return out


def texts(p):
    suffix=p.suffix.lower()
    if suffix=='.pdf':
        text=subprocess.check_output(['pdftotext','-layout',str(p),'-'],stderr=subprocess.DEVNULL).decode()
        for i,page in enumerate(text.split('\f'),1):yield f'PDF page {i}',page,'text'
    elif suffix=='.docx':
        with ZipFile(p) as z:
            for f in z.namelist():
                if f.startswith('word/') and f.endswith('.xml'):
                    root=E.fromstring(z.read(f));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
                    for i,para in enumerate(root.findall('.//w:p',ns)):
                        yield f'{f} paragraph {i}',''.join(t.text or '' for t in para.findall('.//w:t',ns)),'text'
    elif suffix=='.xlsx':
        w=load_workbook(p,data_only=True)
        for sheet in w:
            for row in sheet:
                for c in row:
                    if c.value is not None:yield f'{sheet.title}!{c.coordinate}',str(c.value),'cell'
    elif suffix=='.pages':
        with ZipFile(p) as z:
            for f in z.namelist():
                if f.endswith('.iwa'):
                    raw=z.read(f);parts=[];i=0
                    while i<len(raw):
                        kind=raw[i];n=int.from_bytes(raw[i+1:i+4],'little');i+=4
                        chunk=raw[i:i+n];i+=n
                        if kind==0:parts.append(snappy(chunk))
                    yield f,bytes().join(parts).decode('utf8',errors='replace'),'iwa_text'
    elif suffix=='.ipynb':
        obj=json.loads(p.read_text())
        for i,c in enumerate(obj['cells']):
            yield f'cell {i} source',''.join(c.get('source',[])),'text'
            for j,o in enumerate(c.get('outputs',[])):
                yield f'cell {i} output {j}',''.join(o.get('text',[]))+' '.join(''.join(v) for k,v in o.get('data',{}).items() if k.startswith('text/')),'text'
    elif suffix in ('.md','.txt','.tex','.py','.r','.json','.csv','.svg','.chat'):
        for i,line in enumerate(p.read_text(errors='replace').splitlines(),1):
            yield f'line {i}',line,'svg_geometry' if suffix=='.svg' and (' d=' in line or re.match(r'\s*[MLCQZmlcqz] ',line)) else 'text'


def main():
    records=[];files=[];failures=[]
    # Scan baseline file list, not newly generated audit files (avoids self-matching).
    manifest=json.loads((ROOT/'qc/baseline_manifest.json').read_text())
    for entry in manifest['files']:
        p=ROOT/entry['path']
        if not p.exists():continue
        try:
            for location,text,kind in texts(p):
                for match in PATTERN.finditer(text):
                    context=text[max(0,match.start()-110):match.end()+110].replace('\n',' ')
                    records.append({'file':entry['path'],'location':location,'match':match.group(),
                        'context':context,'content_kind':kind,'quantity':'PENDING REVIEW','wrong_cutoff_lr_minus':False,'action':'unchanged'})
            files.append(entry['path'])
        except Exception as e:failures.append({'file':entry['path'],'error':str(e)})
    target=ROOT/'qc/private';target.mkdir(parents=True,exist_ok=True)
    (target/'audit_024_full.json').write_text(json.dumps({'occurrences':records,'files_scanned':files,'extraction_errors':failures},indent=2)+'\n')
    # Private workbook/notebook cell context is never included in committed summaries.
    for r in records:
        if r['content_kind']=='svg_geometry':r['quantity']='SVG geometry; not a likelihood ratio'
        elif r['content_kind']=='cell':r['quantity']='Workbook measurement or cell text; inspect locally'
    print('occurrences',len(records),'extraction errors',failures)
    for i,r in enumerate(records):
        if r['content_kind'] not in ('svg_geometry','cell'):print(i,r['file'],r['location'],r['match'],r['context'][:150])

if __name__=='__main__':main()
