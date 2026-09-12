from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parent
m=json.loads((r/'SHA256SUMS.json').read_text())
for p,h in m.items():
 assert hashlib.sha256((r/p).read_bytes()).hexdigest()==h,p
print(f'Verified {len(m)} release files')
