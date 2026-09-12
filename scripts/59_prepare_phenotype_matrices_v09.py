"""Sample/probe identity and matrix integrity, no case-control expression tests."""
from pathlib import Path
import gzip,json,csv,re,datetime,hashlib
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];O=R/'public_data/phenotype_evidence_v0.9'
sm=json.loads((O/'GSE65682_samples.json').read_text());meta=[]
for s in sm:
 d={v.split(':',1)[0].strip():v.split(':',1)[1].strip() for v in s['!Sample_characteristics_ch1']}
 if d['thrombocytopenia']=='NA':continue
 d.update(GSM=s['GSM'],title=s['!Sample_title'][0],sample_suffix=s['!Sample_title'][0].split('_')[-1].rstrip(']'))
 meta.append(d)
meta=pd.DataFrame(meta);assert len(meta)==99 and meta.GSM.is_unique and meta.sample_suffix.is_unique
assert not meta.icu_acquired_infection_paired.str.contains('Post|Follow|ICUAI_',case=False).any()
assert set(meta.thrombocytopenia)=={'A_very_low','B_medium_low','C_low','D_normal'}
meta.to_csv(O/'GSE65682_phenotype_samples.tsv',sep='\t',index=False)
want=set(meta.GSM);arrays={};annotations=[];inplatform=False;insample=False;current=None;rows=[];header=None
with gzip.open(O/'GSE65682_family.soft.gz','rt') as f:
 for ln in f:
  ln=ln.rstrip('\r\n')
  if ln=='!platform_table_begin':inplatform=True;header=next(f).rstrip().split('\t');continue
  if ln=='!platform_table_end':inplatform=False;continue
  if inplatform:
   a=ln.split('\t');annotations.append(dict(zip(header,a)));continue
  if ln.startswith('^SAMPLE = '):current=ln.split(' = ')[1]
  elif ln=='!sample_table_begin':insample=current in want;rows=[];hdr=next(f).rstrip().split('\t')
  elif ln=='!sample_table_end':
   if insample:
    a=pd.DataFrame(rows,columns=hdr);assert set(['ID_REF','VALUE']).issubset(a.columns)
    arrays[current]=pd.Series(pd.to_numeric(a.VALUE).values,index=a.ID_REF,name=current)
   insample=False
  elif insample:rows.append(ln.split('\t'))
assert len(arrays)==99
mat=pd.concat(arrays,axis=1);assert not mat.index.duplicated().any() and np.isfinite(mat.values).all()
mat.index.name='probe_id';mat.to_csv(O/'GSE65682_selected_expression.tsv.gz',sep='\t',compression='gzip')
ann=pd.DataFrame(annotations);ann=ann[['ID','Gene Symbol','Ensembl','Entrez Gene']].rename(columns={'ID':'probe_id','Gene Symbol':'symbol'})
ann.to_csv(O/'GPL13667_annotation.tsv.gz',sep='\t',index=False,compression='gzip')
src=R/'public_data/reference_sources_v0.8/GSE273700_family.soft.gz.json';ss=json.loads(src.read_text());rows=[]
for s in ss:
 d={v.split(':',1)[0].strip():v.split(':',1)[1].strip() for v in s['!Sample_characteristics_ch1']}
 rows.append(dict(GSM=s['GSM'],sample=s['!Sample_title'][0].split(' [')[0],patient=d['patient'],time=d['time'],condition=d['condition'],declared_tissue=d['tissue']))
ss=pd.DataFrame(rows);assert len(ss)==104 and ss.GSM.is_unique and ss['sample'].is_unique and not ss.duplicated(['patient','time']).any()
with gzip.open(O/'GSE273700_All_read_count_matrix.txt.gz','rt') as f:
 hdr=next(f).strip().split('\t');a=pd.read_csv(f,sep='\t',header=None)
assert a.shape[1]==len(hdr)==105 and hdr[-1]=='Geneid'
hdr=hdr[:-1]
ids=a.iloc[:,-1].astype(str);cc=a.iloc[:,:-1];cc.columns=hdr;cc.index=ids;cc.index.name='ensembl_id'
assert set(hdr)==set(ss['sample']) and ids.is_unique and np.isfinite(cc.values).all() and (cc.values>=0).all() and (cc.values==np.floor(cc.values)).all()
cc.to_csv(O/'GSE273700_counts_mapped.tsv.gz',sep='\t',compression='gzip');ss.to_csv(O/'GSE273700_samples.tsv',sep='\t',index=False)
day1=ss[ss.time=='Day1'];assert len(day1)==52 and day1.patient.is_unique
report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),GSE65682=dict(samples=len(meta),group_counts=meta.thrombocytopenia.value_counts().to_dict(),unique_suffixes=meta.sample_suffix.nunique(),probes=mat.shape[0],value_range=[float(mat.min().min()),float(mat.max().max())],age_missing=int(pd.to_numeric(meta.age,errors='coerce').isna().sum()),pneumonia_counts=meta['pneumonia diagnoses'].value_counts().to_dict(),preprocessed_as_deposited=True),
 GSE273700=dict(total_samples=len(ss),Day1_samples=len(day1),Day1_group_counts=day1.condition.value_counts().to_dict(),genes=len(ids),all_columns_mapped=True,nonnegative_integer_counts=True,tissue_conflict_retained=True,serial_recovery_labels_available=False))
(R/'reports/INDEPENDENT_PHENOTYPE_INPUT_AUDIT_v0.9.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
