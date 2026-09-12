"""Download public GEO metadata and inspect the platelet count workbook without DE testing."""
from pathlib import Path
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import time
import xml.etree.ElementTree as ET
import zipfile
import requests
import argparse

os.umask(0o077)
ROOT=Path('/root/projects/linezolid_platelet_recovery')
DEST=ROOT/'public_data'
records=[]
parser=argparse.ArgumentParser()
parser.add_argument('--download-workbook',action='store_true')
args=parser.parse_args()

def download(url,path):
    start=time.time()
    if not path.exists():
        with requests.get(url,stream=True,timeout=(20,60)) as response:
            response.raise_for_status()
            temp=Path(str(path)+'.partial')
            with temp.open('wb') as out:
                for chunk in response.iter_content(1024*1024):
                    if chunk: out.write(chunk)
            temp.replace(path)
    digest=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
    records.append({'url':url,'path':str(path),'bytes':path.stat().st_size,'sha256':digest,
                    'accessed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'access':'public anonymous FTP-over-HTTPS'})
    print('DOWNLOADED',path.name,path.stat().st_size,round(time.time()-start,1),flush=True)

summaries={}
for accession in ['GSE252275','GSE310202']:
    family=accession[:-3]+'nnn'
    metadata=DEST/f'{accession}_family.soft.gz'
    download(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{family}/{accession}/soft/{accession}_family.soft.gz',metadata)
    series={}
    samples=[]
    current=None
    with gzip.open(metadata,'rt') as stream:
        for line in stream:
            line=line.rstrip('\n')
            if line.startswith('^SAMPLE = '):
                current={'accession':line.split(' = ',1)[1]}
                samples.append(current)
            elif line.startswith('!Series_') and ' = ' in line:
                key,value=line.split(' = ',1)
                if key in ['!Series_title','!Series_summary','!Series_overall_design','!Series_pubmed_id','!Series_supplementary_file','!Series_relation','!Series_last_update_date']:
                    series.setdefault(key,[]).append(value)
            elif current is not None and line.startswith('!Sample_') and ' = ' in line:
                key,value=line.split(' = ',1)
                if key in ['!Sample_title','!Sample_source_name_ch1','!Sample_organism_ch1','!Sample_characteristics_ch1','!Sample_relation','!Sample_data_processing','!Sample_supplementary_file']:
                    current.setdefault(key,[]).append(value)
    summaries[accession]={'series':series,'sample_records':len(samples),'samples':samples}
    (DEST/f'{accession}_metadata.json').write_text(json.dumps(summaries[accession],indent=2))

(ROOT/'manifests/public_assets.json').write_text(json.dumps(records,indent=2))
if not args.download_workbook:
    print('PUBLIC_METADATA_COMPLETE; full workbook deferred, HTTP range inspection is separate',flush=True)
    raise SystemExit(0)

files=summaries['GSE252275']['series'].get('!Series_supplementary_file',[])
workbooks=[url for url in files if url.endswith('.xlsx')]
if len(workbooks)!=1:
    raise ValueError(f'Expected one original GEO workbook, found {len(workbooks)}')
url=workbooks[0].replace('ftp://','https://')
path=DEST/url.rsplit('/',1)[1]
download(url,path)

ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
inspection={'file':str(path),'sheets':[],'worksheet_previews':{}}
with zipfile.ZipFile(path) as archive:
    workbook=ET.fromstring(archive.read('xl/workbook.xml'))
    inspection['sheets']=[dict(s.attrib) for s in workbook.findall('m:sheets/m:sheet',ns)]
    # Header and a few rows only. No gene selection, group comparison, or inference.
    shared=[]
    if 'xl/sharedStrings.xml' in archive.namelist():
        for event,element in ET.iterparse(archive.open('xl/sharedStrings.xml'),events=['end']):
            if element.tag.endswith('}si'):
                shared.append(''.join(t.text or '' for t in element.iter() if t.tag.endswith('}t')))
                element.clear()
    for sheet in [n for n in archive.namelist() if n.startswith('xl/worksheets/sheet') and n.endswith('.xml')]:
        rows=[]
        for event,element in ET.iterparse(archive.open(sheet),events=['end']):
            if element.tag.endswith('}row'):
                values=[]
                for cell in element:
                    raw=cell.find('m:v',ns)
                    value=raw.text if raw is not None else None
                    if cell.attrib.get('t')=='s' and value is not None: value=shared[int(value)]
                    if cell.attrib.get('t')=='inlineStr': value=''.join(t.text or '' for t in cell.iter() if t.tag.endswith('}t'))
                    values.append({'cell':cell.attrib.get('r'),'value':value})
                rows.append(values)
                element.clear()
                if len(rows)>=4: break
        inspection['worksheet_previews'][sheet]=rows
(DEST/'GSE252275_workbook_inspection.json').write_text(json.dumps(inspection,indent=2))
(ROOT/'manifests/public_assets.json').write_text(json.dumps(records,indent=2))
print('PUBLIC_ASSET_AUDIT_COMPLETE',flush=True)
