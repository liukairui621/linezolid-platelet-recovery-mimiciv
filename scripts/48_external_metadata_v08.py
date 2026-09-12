"""Download public annotations and external cohort metadata; no outcome contrasts."""
from pathlib import Path
import urllib.request,json,hashlib,gzip,re,datetime
R=Path(__file__).resolve().parents[1];O=R/'public_data/reference_sources_v0.8'
sources={
 'hgnc_complete_set.txt':'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt',
 'GSE210797_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210797/soft/GSE210797_family.soft.gz',
 'GSE273700_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273700/soft/GSE273700_family.soft.gz',
 'ventilation.sql':'https://raw.githubusercontent.com/MIT-LCP/mimic-code/main/mimic-iv/concepts/treatment/ventilation.sql',
 'ventilator_setting.sql':'https://raw.githubusercontent.com/MIT-LCP/mimic-code/main/mimic-iv/concepts/measurement/ventilator_setting.sql',
 'oxygen_delivery.sql':'https://raw.githubusercontent.com/MIT-LCP/mimic-code/main/mimic-iv/concepts/measurement/oxygen_delivery.sql',
 'rrt.sql':'https://raw.githubusercontent.com/MIT-LCP/mimic-code/main/mimic-iv/concepts/treatment/rrt.sql'}
log=[]
for name,url in sources.items():
 p=O/name
 try:
  if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=90).read())
  log.append(dict(file=name,url=url,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
  print(name,p.stat().st_size,flush=True)
  if name.endswith('soft.gz'):
   txt=gzip.open(p,'rt').read();samples=[];cur=None
   for line in txt.splitlines():
    if line.startswith('^SAMPLE = '):cur={'GSM':line.split(' = ')[1]};samples.append(cur)
    elif line.startswith('!Sample_') and cur is not None:
     k,v=line.split(' = ',1);cur.setdefault(k,[]).append(v)
   (O/(name+'.json')).write_text(json.dumps(samples,indent=2),encoding='utf-8')
   print('N samples',len(samples),'first samples',[{k:v for k,v in x.items() if k in ['GSM','!Sample_title','!Sample_characteristics_ch1']} for x in samples[:4]])
   print('SUPPLEMENT',[x for x in txt.splitlines() if x.startswith('!Series_supplementary_file')])
 except Exception as e:log.append(dict(file=name,url=url,error=str(e)));print(name,str(e))
(O/'external_download_manifest.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=log),indent=2))
