"""Verify archived official references; recover a missing archive only if its hash matches."""
from pathlib import Path
import hashlib,json,urllib.request
R=Path(__file__).resolve().parents[1];D=R/'public_data/reference_sources_v0.10'
m=json.loads((D/'source_manifest.json').read_text())
for s in m['sources']:
 p=D/s['file']
 if p.exists():b=p.read_bytes()
 else:
  req=urllib.request.Request(s['url'],headers={'User-Agent':'Mozilla/5.0 research source archive'})
  with urllib.request.urlopen(req,timeout=60) as response:b=response.read()
  assert hashlib.sha256(b).hexdigest()==s['sha256'],'Remote reference changed; preserve original manifest and fetch a separately versioned reference'
  p.write_bytes(b)
 assert len(b)==s['bytes'] and hashlib.sha256(b).hexdigest()==s['sha256']
print('Three official source archives match recorded hashes.')
