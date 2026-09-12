"""Build all-sample matrices with stable feature identity and explicit count-column QC."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import zipfile,json,re,hashlib,csv,datetime
from lxml import etree as ET
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1];O=R/'public_data';F=O/'GSE252275_raw_count_data_GeneLevelExpression_.xlsx'
N='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
STRINGS=[]
def initialize():
    global STRINGS
    with zipfile.ZipFile(F) as z:
        x=ET.fromstring(z.read('xl/sharedStrings.xml'))
        STRINGS=[''.join(si.itertext()) for si in x]
def decode(cell):
    v=cell.find(N+'v')
    if cell.get('t')=='inlineStr':return ''.join(cell.itertext())
    if v is None:return ''
    return STRINGS[int(v.text)] if cell.get('t')=='s' else v.text
def sheet_task(spec):
    name,path=spec;ident=[];values=[];header=None
    cols=['D','L','M','Q','R']; wanted=set(['A','B','H','J','K']+cols)
    with zipfile.ZipFile(F) as z, z.open(path) as f:
        for _,row in ET.iterparse(f,events=['end'],tag=N+'row'):
            if header is None:header={re.sub('[0-9]','',c.get('r')):decode(c) for c in row}
            else:
                d={re.sub('[0-9]','',c.get('r')):decode(c) for c in row if re.sub('[0-9]','',c.get('r')) in wanted}
                if any(d.values()):
                    ident.append([d.get(c,'') for c in ['A','B','H','J','K']])
                    values.append([float(d[c]) if d.get(c,'')!='' else np.nan for c in cols])
            row.clear()
            while row.getprevious() is not None:del row.getparent()[0]
    a=np.asarray(values,dtype=np.float64);ids=np.asarray(ident)
    out=O/'GSE252275_sample_arrays';out.mkdir(exist_ok=True)
    np.savez_compressed(out/f'{name}.npz',identifiers=ids,values=a)
    qc={'sample':name,'rows':len(a),'header':header,'count_column_names':[header[c] for c in cols],
        'nonfinite_per_column':(~np.isfinite(a)).sum(axis=0).tolist(),
        'negative_per_column':(a<0).sum(axis=0).tolist(),
        'noninteger_per_column':(np.isfinite(a)&(a!=np.floor(a))).sum(axis=0).tolist(),
        'expression_vs_total_exon_differences':int(np.count_nonzero(a[:,0]!=a[:,4])),
        'expression_vs_total_gene_differences':int(np.count_nonzero(a[:,0]!=a[:,2])),
        'expression_vs_unique_exon_differences':int(np.count_nonzero(a[:,0]!=a[:,3])),
        'column_sums':a.sum(axis=0).tolist(),'ensembl_unique':len(set(ids[:,3])),
        'gene_id_unique':len(set(ids[:,2])),'missing_ensembl':int(np.count_nonzero(ids[:,3]==''))}
    return qc
def main():
    (R/'manifests').mkdir(exist_ok=True)
    with zipfile.ZipFile(F) as z:
        w=ET.fromstring(z.read('xl/workbook.xml'))
        rel=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        mapping={r.get('Id'):r.get('Target').lstrip('/') for r in rel}
        specs=[]
        for s in w.find(N+'sheets'):
            target=mapping[s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')]
            specs.append((s.get('name'),target if target.startswith('xl/') else 'xl/'+target))
    meta=json.loads((O/'GSE252275_metadata.json').read_text())
    m={s['!Sample_title'][0]:s for s in meta['samples']}
    assert set(n for n,p in specs)==set(m), 'Workbook/GEO sample mismatch'
    q=[]
    with ProcessPoolExecutor(max_workers=4,initializer=initialize) as ex:
        for f in as_completed([ex.submit(sheet_task,s) for s in specs]):
            x=f.result();q.append(x);print(json.dumps({'sample':x['sample'],'rows':x['rows'],'completed':len(q)}),flush=True)
    q.sort(key=lambda x:x['sample']);names=[n for n,p in specs]
    arrays=[];expressions=[];ref=None;sample_rows=[]
    for name in names:
        with np.load(O/'GSE252275_sample_arrays'/f'{name}.npz') as a:
            ids=a['identifiers'];v=a['values']
            if ref is None:ref=ids
            else:assert np.array_equal(ref,ids),'Gene annotations/order differ; explicit join needed'
            arrays.append(v[:,4]);expressions.append(v[:,0])
        group='control' if name.startswith('control_') else 'sepsis_decrease' if name.startswith('sepsis_decrease_') else 'sepsis_preserve'
        sample_rows.append({'sample':name,'GSM':m[name]['accession'],'group':group,'patient_id':'not_reported','paired_sampling':'not_reported'})
    count=np.column_stack(arrays);expression=np.column_stack(expressions)
    assert np.isfinite(count).all() and (count>=0).all() and (count==np.floor(count)).all(),'Invalid raw total-exon counts'
    assert len(set(ref[:,3]))==len(ref) and (ref[:,3]!='').all(),'Unresolved Ensembl identity'
    out=O/'GSE252275_processed';out.mkdir(exist_ok=True)
    count=count.astype(np.int64)
    pd.DataFrame(count,index=pd.Index(ref[:,3],name='ensembl_id'),columns=names).to_csv(out/'total_exon_counts.tsv.gz',sep='\t',compression='gzip')
    pd.DataFrame(expression,index=pd.Index(ref[:,3],name='ensembl_id'),columns=names).to_csv(out/'expression_value_as_supplied.tsv.gz',sep='\t',compression='gzip')
    annotations=pd.DataFrame(ref,columns=['symbol','chromosome','gene_id_as_supplied','ensembl_id','biotype'])
    annotations.to_csv(out/'gene_annotations.tsv',sep='\t',index=False)
    sampledf=pd.DataFrame(sample_rows);sampledf['library_total_exon_counts']=count.sum(axis=0)
    sampledf['detected_features']=(count>0).sum(axis=0)
    sampledf['mitochondrial_count_fraction']=count[ref[:,1]=='MT'].sum(axis=0)/count.sum(axis=0)
    sampledf.to_csv(out/'samples_and_qc.tsv',sep='\t',index=False)
    same_expression=bool(np.array_equal(count,expression))
    cpm=count/count.sum(axis=0)*1e6
    sepsis=np.array([s['group']!='control' for s in sample_rows])
    keep=(cpm[:,sepsis]>=1).sum(axis=1)>=min(sum(sampledf.group=='sepsis_decrease'),sum(sampledf.group=='sepsis_preserve'))
    corr=np.corrcoef(np.log2(cpm[keep]+0.5).T)
    pd.DataFrame(corr,index=names,columns=names).to_csv(out/'logcpm_sample_correlation.tsv',sep='\t')
    duplicate_pairs=[(names[i],names[j]) for i in range(len(names)) for j in range(i) if np.array_equal(count[:,i],count[:,j])]
    report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':len(names),'features':len(ref),
        'groups':sampledf.group.value_counts().to_dict(),'workbook_to_GEO_mapping_exact':True,
        'all_samples_same_ordered_feature_annotations':True,'duplicate_feature_ensembl_ids':0,
        'total_exon_counts_nonnegative_integers':True,'expression_equals_total_exon_counts_all_cells':same_expression,
        'exact_duplicate_sample_vectors':duplicate_pairs,'all_zero_features':int((count.sum(axis=1)==0).sum()),
        'cpm_at_least_1_in_at_least_9_sepsis_samples':int(keep.sum()),
        'library_size_min':int(count.sum(axis=0).min()),'library_size_max':int(count.sum(axis=0).max()),
        'gene_id_equals_ensembl_all_rows':bool(np.array_equal(ref[:,2],ref[:,3])),
        'sample_qc':sampledf.to_dict('records'),'per_sheet_count_audit':q,
        'count_choice':'Total exon reads; Expression value is equivalent only if all-cell equality check passes. TPM/RPKM never used as raw counts.',
        'limits':['Group labels are as deposited; low-platelet threshold and clinical covariates not supplied.',
                  'Forty samples are not yet confirmed as forty independent donors.',
                  'Processed-count QC does not validate upstream mapping or rule out leukocyte contamination.',
                  'No differential-expression test, selected candidate gene, or drug-effect conclusion in this step.']}
    (R/'reports/GSE252275_FULL_COUNT_QC_2026-09-11.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    manifest={p.name:{'size':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
    (R/'manifests/GSE252275_processed.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k not in ['sample_qc','per_sheet_count_audit']}),flush=True)
if __name__=='__main__':main()
