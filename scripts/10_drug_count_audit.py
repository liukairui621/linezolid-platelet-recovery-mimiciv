"""Download and audit human linezolid perturbation tables, without rounding or relabeling."""
from pathlib import Path
import json,urllib.request,gzip,hashlib,datetime,re
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[1];O=R/'public_data'
meta=json.loads((O/'GSE310202_metadata.json').read_text())
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'datasets':[]}
for source in meta['series']['!Series_supplementary_file']:
    if '_mice_' in source:continue
    url=source.replace('ftp://','https://');p=O/url.rsplit('/',1)[-1]
    if not p.exists():
        with urllib.request.urlopen(url,timeout=90) as r:data=r.read()
        gzip.decompress(data);p.write_bytes(data)
    x=pd.read_csv(p,sep='\t',index_col=0);v=x.to_numpy(dtype=float)
    clean=x.index.str.replace(r'\.\d+$','',regex=True)
    entry={'file':p.name,'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
           'features':len(x),'samples':len(x.columns),'sample_columns':list(x.columns),
           'nonfinite_cells':int((~np.isfinite(v)).sum()),'negative_cells':int((v<0).sum()),
           'noninteger_cells':int((np.isfinite(v)&(v!=np.floor(v))).sum()),
           'duplicate_original_features':int(x.index.duplicated().sum()),
           'duplicate_features_after_version_strip':int(clean.duplicated().sum()),
           'library_sums':dict(zip(x.columns,v.sum(axis=0).tolist()))}
    if 'macrophage' in p.name:
        assert entry['nonfinite_cells']==0 and entry['negative_cells']==0
        assert not clean.duplicated().any()
        output=O/'GSE310202_processed';output.mkdir(exist_ok=True)
        x.index=clean;x.index.name='ensembl_id';x.to_csv(output/'macrophage_counts_as_supplied.tsv.gz',sep='\t',compression='gzip')
        samples=[]
        for name in x.columns:
            treated=name.startswith('TGFL')
            samples.append({'sample':name,'condition':'TGFbeta_linezolid' if treated else 'TGFbeta',
                'block_from_column_suffix':name[4:] if treated else name[3:],
                'block_interpretation':'Matching suffix as supplied; donor/clone meaning not independently confirmed',
                'GSM_mapping':'unresolved','duration_hours':48,'model':'human iPSC-derived macrophage'})
        pd.DataFrame(samples).to_csv(output/'macrophage_design_provisional.tsv',sep='\t',index=False)
        entry['status']='Finite nonnegative fractional values preserved. Treatment prefixes explicit. Verify quantification provenance, suffix pairing and GSM correspondence. No integer conversion.'
    elif 'pcss' in p.name:
        entry['status']='Hold for provenance clarification: fractional values and EN1_Mut_185 among presumed controls. Do not rename this column to DMSO_185.'
    else:
        entry['status']='Fractional counts require quantification provenance and appropriate count import or log-count model; do not round silently for DESeq2.'
    report['datasets'].append(entry)
report['scope']='Cross-cell drug perturbation evidence only; no human platelet drug intervention is established by these datasets. Fractional estimated/weighted counts can be legitimate; their presence alone does not prove normalization or invalid data. Establish provenance and choose a compatible model.'
(R/'reports/GSE310202_HUMAN_COUNT_QC_2026-09-11.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
