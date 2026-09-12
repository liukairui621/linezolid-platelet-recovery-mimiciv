"""Acquire public phenotype metadata and matrices, without new contrasts."""
from pathlib import Path
import urllib.request,gzip,json,hashlib,datetime
R=Path(__file__).resolve().parents[1];O=R/'public_data/phenotype_evidence_v0.9';O.mkdir(parents=True,exist_ok=True)
sources={
 'GSE65682_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/soft/GSE65682_family.soft.gz',
 'GSE273700_All_read_count_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273700/suppl/GSE273700_All_read_count_matrix.txt.gz'}
log=[]
for name,url in sources.items():
 p=O/name
 try:
  if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=120).read())
  log.append(dict(file=name,url=url,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
  print(name,p.stat().st_size,flush=True)
  if name.endswith('soft.gz'):
   txt=gzip.open(p,'rt').read();samples=[];s=None
   for line in txt.splitlines():
    if line.startswith('^SAMPLE = '):s={'GSM':line.split(' = ')[1]};samples.append(s)
    elif line.startswith('!Sample_') and s is not None:
     k,v=line.split(' = ',1);s.setdefault(k,[]).append(v)
   (O/'GSE65682_samples.json').write_text(json.dumps(samples,indent=2))
   print('N',len(samples),'first',[{k:v for k,v in x.items() if k in ['GSM','!Sample_title','!Sample_characteristics_ch1']} for x in samples[:2]])
   print('Sample characteristic keys',sorted(set(v.split(':',1)[0] for x in samples for v in x.get('!Sample_characteristics_ch1',[]))))
   print('Series', [x for x in txt.splitlines() if x.startswith(('!Series_summary','!Series_overall_design','!Series_pubmed_id','!Series_supplementary_file'))][:12])
  else:
   with gzip.open(p,'rt') as f:print('Header',next(f)[:1300]);print('First row',next(f)[:350])
 except Exception as e:log.append(dict(file=name,url=url,error=str(e)));print(type(e).__name__,str(e),flush=True)
(O/'download_manifest.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=log),indent=2))
