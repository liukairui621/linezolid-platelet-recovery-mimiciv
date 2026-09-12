"""Public reference acquisition only; no expression or clinical outcome tests."""
from pathlib import Path
import urllib.request, hashlib, json, datetime, zipfile, io, re
import xml.etree.ElementTree as ET
R=Path(__file__).resolve().parents[1]
O=R/'public_data/reference_sources_v0.8';O.mkdir(parents=True,exist_ok=True)
sources={
 'ReactomePathways.gmt.zip':'https://reactome.org/download/current/ReactomePathways.gmt.zip',
 'PMC6753624.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6753624/fullTextXML',
 'PMC13054457.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13054457/fullTextXML',
 'PMC12820328.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12820328/fullTextXML',
 'sepsis3.sql':'https://raw.githubusercontent.com/MIT-LCP/mimic-code/main/mimic-iv/concepts/sepsis/sepsis3.sql'}
manifest=[]
for name,url in sources.items():
 p=O/name
 try:
  if not p.exists():
   req=urllib.request.Request(url,headers={'User-Agent':'research-reference-audit/0.8'})
   data=urllib.request.urlopen(req,timeout=90).read();p.write_bytes(data)
  manifest.append(dict(file=name,url=url,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
  if name.endswith('.xml'):
   root=ET.fromstring(p.read_bytes());txt=' '.join(root.itertext())
   (O/(name+'.txt')).write_text(txt,encoding='utf-8')
   print(name,'GEO',sorted(set(re.findall(r'GSE\d+',txt))))
   for node in root.iter('sec'):
    title=node.find('title')
    if title is not None and re.search('availability|accession', ''.join(title.itertext()),re.I):
     print('DATA:', ' '.join(node.itertext())[:1800])
  print(name, 'OK', p.stat().st_size)
 except Exception as e:manifest.append(dict(file=name,url=url,error=str(e)));print(name,type(e).__name__,str(e))
(O/'download_manifest.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sources':manifest},indent=2),encoding='utf-8')
p=O/'ReactomePathways.gmt.zip'
if p.exists():
 z=zipfile.ZipFile(p)
 print('GMT MEMBERS',z.namelist())
 for member in z.namelist():
  if member.endswith('.gmt'):
   txt=z.read(member).decode();(O/'ReactomePathways.gmt').write_text(txt,encoding='utf-8')
   for line in txt.splitlines():
    a=line.split('\t')
    if re.search('mitochondrial translation$|respiratory electron transport, ATP|platelet degranulation$|megakaryocyte|erythrocyte differentiation|eukaryotic translation initiation$|heme biosynthesis$|platelet activation, signaling',a[0],re.I):print(a[:2],len(a)-2)
