"""Resumable bounded parallel byte-range download from NCBI; validates all ZIP CRCs."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request, hashlib, json, time, zipfile, datetime

R=Path(__file__).resolve().parents[1]
O=R/'public_data'; C=O/'GSE252275_chunks'; C.mkdir(exist_ok=True)
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252275/suppl/GSE252275_raw_count_data_GeneLevelExpression_.xlsx'
with urllib.request.urlopen(urllib.request.Request(URL,method='HEAD'),timeout=60) as r:
    SIZE=int(r.headers['Content-Length']); MOD=r.headers.get('Last-Modified')
BLOCK=2*1024*1024
def fetch(start):
    end=min(SIZE-1,start+BLOCK-1); p=C/f'{start}_{end}.bin'
    if p.exists() and p.stat().st_size==end-start+1:return p
    for trial in range(5):
        try:
            req=urllib.request.Request(URL,headers={'Range':f'bytes={start}-{end}','If-Unmodified-Since':MOD,'Accept-Encoding':'identity'})
            with urllib.request.urlopen(req,timeout=90) as r:
                assert r.status==206, r.status
                assert r.headers.get('Content-Range')==f'bytes {start}-{end}/{SIZE}',r.headers
                data=r.read()
                assert len(data)==end-start+1, len(data)
            partial=p.with_suffix('.partial');partial.write_bytes(data);partial.replace(p)
            return p
        except Exception:
            if trial==4:raise
            time.sleep(2*(trial+1))
starts=list(range(0,SIZE,BLOCK)); done=0; t=time.monotonic()
with ThreadPoolExecutor(max_workers=12) as ex:
    for f in as_completed([ex.submit(fetch,s) for s in starts]):
        f.result();done+=1
        if done%10==0:print(json.dumps({'chunks_done':done,'chunks_total':len(starts),'seconds':round(time.monotonic()-t)}),flush=True)
out=O/'GSE252275_raw_count_data_GeneLevelExpression_.xlsx';partial=out.with_suffix('.xlsx.assembling')
h=hashlib.sha256()
with partial.open('wb') as w:
    for s in starts:
        data=(C/f'{s}_{min(SIZE-1,s+BLOCK-1)}.bin').read_bytes();w.write(data);h.update(data)
assert partial.stat().st_size==SIZE
with zipfile.ZipFile(partial) as z:
    bad=z.testzip();assert bad is None,bad
partial.replace(out)
record={'url':URL,'accessed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'last_modified':MOD,'size':SIZE,'sha256':h.hexdigest(),'zip_crc_all_members_pass':True,'chunks':len(starts),'workers':12,'seconds':round(time.monotonic()-t)}
(O/'GSE252275_workbook_download.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record),flush=True)
