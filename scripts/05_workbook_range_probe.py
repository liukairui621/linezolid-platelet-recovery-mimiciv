"""Inspect XLSX structure through standard public HTTP byte ranges; no DE analysis."""
from pathlib import Path
import io
import json
import hashlib
import datetime as dt
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'public_data'
OUT.mkdir(exist_ok=True)
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252275/suppl/GSE252275_raw_count_data_GeneLevelExpression_.xlsx'

class HTTPRangeFile(io.RawIOBase):
    def __init__(self,url,size,blocksize=65536):
        self.url=url; self.size=size; self.blocksize=blocksize; self.pos=0; self.blocks={}; self.records=[]
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self,offset,whence=0):
        self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
        return self.pos
    def read(self,n=-1):
        if n<0:n=self.size-self.pos
        n=min(n,self.size-self.pos)
        chunks=[]
        while n>0:
            b=self.pos//self.blocksize
            if b not in self.blocks:
                start=b*self.blocksize; end=min(self.size-1,start+self.blocksize-1)
                req=urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{end}'})
                with urllib.request.urlopen(req,timeout=45) as r:
                    if r.status!=206: raise RuntimeError(f'Byte-range request returned {r.status}')
                    data=r.read()
                    if len(data)!=end-start+1: raise RuntimeError('Truncated range')
                self.blocks[b]=data
                p=OUT/f'GSE252275_range_{start}_{end}.bin';p.write_bytes(data)
                self.records.append({'range':[start,end],'path':str(p),'sha256':hashlib.sha256(data).hexdigest()})
            off=self.pos%self.blocksize;take=min(n,len(self.blocks[b])-off)
            chunks.append(self.blocks[b][off:off+take]);self.pos+=take;n-=take
        return b''.join(chunks)

head=urllib.request.Request(URL,method='HEAD')
with urllib.request.urlopen(head,timeout=30) as r:
    size=int(r.headers['Content-Length']);modified=r.headers.get('Last-Modified')
f=HTTPRangeFile(URL,size)
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
result={'url':URL,'content_length':size,'last_modified':modified,'accessed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'sheets':[],'previews':{}}
with zipfile.ZipFile(f) as z:
    w=ET.fromstring(z.read('xl/workbook.xml'))
    result['sheets']=[dict(s.attrib) for s in w.findall('m:sheets/m:sheet',ns)]
    needed=set()
    for name in [n for n in z.namelist() if n.startswith('xl/worksheets/sheet') and n.endswith('.xml')]:
        rows=[]
        with z.open(name) as stream:
            for _,element in ET.iterparse(stream,events=['end']):
                if element.tag.endswith('}row'):
                    cells=[]
                    for cell in element:
                        v=cell.find('m:v',ns);value=v.text if v is not None else None
                        if cell.attrib.get('t')=='s' and value is not None:needed.add(int(value))
                        if cell.attrib.get('t')=='inlineStr':value=''.join(x.text or '' for x in cell.iter() if x.tag.endswith('}t'))
                        cells.append({'cell':cell.attrib.get('r'),'type':cell.attrib.get('t'),'value':value})
                    rows.append(cells);element.clear()
                    if len(rows)>=3:break
        result['previews'][name]=rows
    strings={}
    if needed and 'xl/sharedStrings.xml' in z.namelist():
        idx=0
        with z.open('xl/sharedStrings.xml') as stream:
            for _,element in ET.iterparse(stream,events=['end']):
                if element.tag.endswith('}si'):
                    if idx in needed:strings[idx]=''.join(x.text or '' for x in element.iter() if x.tag.endswith('}t'))
                    element.clear();idx+=1
                    if len(strings)==len(needed):break
        for rows in result['previews'].values():
            for cells in rows:
                for cell in cells:
                    if cell['type']=='s' and cell['value'] is not None:cell['value']=strings.get(int(cell['value']),'<unresolved>')
result['downloaded_bytes']=sum(len(x) for x in f.blocks.values())
result['range_records']=f.records
(OUT/'GSE252275_workbook_range_probe.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({'sheet_names':[s['name'] for s in result['sheets']],
                  'downloaded_bytes':result['downloaded_bytes'],
                  'inspection_file':str(OUT/'GSE252275_workbook_range_probe.json')},indent=2))
