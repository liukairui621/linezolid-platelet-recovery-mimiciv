"""Resolve human drug-matrix columns against original FASTQ names in official SRA XML."""
from pathlib import Path
import json,urllib.request,urllib.parse,hashlib
from lxml import etree as E
import pandas as pd
R=Path(__file__).resolve().parents[1];O=R/'public_data'
m=json.loads((O/'GSE310202_metadata.json').read_text())
samples=[s for s in m['samples'] if s['!Sample_organism_ch1']==['Homo sapiens']]
ids=[next(t.split('term=')[1] for t in s['!Sample_relation'] if 'term=SRX' in t) for s in samples]
p=O/'GSE310202_human_SRA.xml'
if not p.exists():
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode({'db':'sra','id':','.join(ids)})
    with urllib.request.urlopen(url,timeout=90) as f:data=f.read()
    E.fromstring(data);p.write_bytes(data)
root=E.parse(str(p));records=[]
for pkg in root.findall('EXPERIMENT_PACKAGE'):
    ex=pkg.find('EXPERIMENT');gsm=ex.findtext('DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME')
    filenames=sorted(set(f.get('filename') for f in pkg.findall('.//SRAFile') if f.get('supertype')=='Original'))
    records.append({'GSM':gsm,'SRX':ex.get('accession'),'title':ex.findtext('TITLE'),'original_fastq_files':filenames})
assert len(records)==len(samples),'Missing or duplicate SRA packages'
mapping=[]
for label in ['fibroblast','macrophage','pcss']:
    p2=next(O.glob(f'GSE310202_{label}_raw_count_data.*.gz'))
    x=pd.read_csv(p2,sep='\t',index_col=0)
    for col in x.columns:
        found=[r for r in records if any(n.startswith(col+'_') or n.startswith(col+'.') for n in r['original_fastq_files'])]
        mapping.append({'matrix':label,'column':col,'GSM':found[0]['GSM'] if len(found)==1 else '',
            'SRX':found[0]['SRX'] if len(found)==1 else '',
            'GEO_title':found[0]['title'] if len(found)==1 else '',
            'exact_filename_prefix_matches':len(found)})
out={'sra_xml_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'experiment_packages':len(records),
     'matrix_mapping':mapping,'source_records':records,
     'limits':'Original filename mapping establishes deposited sample identity, not independent donor identity or absence of mislabeled biological treatment. Replicate/donor blocking needs source interpretation.'}
(R/'reports/GSE310202_SRA_SAMPLE_MAPPING_2026-09-11.json').write_text(json.dumps(out,indent=2))
pd.DataFrame(mapping).to_csv(O/'GSE310202_processed/matrix_to_GSM_mapping.tsv',sep='\t',index=False)
design_path=O/'GSE310202_processed/macrophage_design_provisional.tsv'
design=pd.read_csv(design_path,sep='\t').drop(columns=['GSM_mapping'])
design=design.merge(pd.DataFrame(mapping).query("matrix == 'macrophage'").drop(columns=['matrix']),left_on='sample',right_on='column',validate='one_to_one')
assert design['exact_filename_prefix_matches'].eq(1).all()
design.to_csv(O/'GSE310202_processed/macrophage_design_SRA_verified.tsv',sep='\t',index=False)
print(json.dumps({'mapping':mapping},ensure_ascii=True))
