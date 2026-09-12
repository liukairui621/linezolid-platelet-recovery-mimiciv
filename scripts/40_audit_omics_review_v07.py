"""Audit public saved differential-expression outputs; no new feature selection."""
from pathlib import Path
import sys,json,datetime,hashlib
import pandas as pd,numpy as np
R=Path(sys.argv[1]);root=R/'public_data/GSE252275_DE_v0.5'
df=pd.read_csv(root/'all20_group.tsv.gz',sep='\t')
loo=pd.read_csv(root/'loo_BH_q.tsv.gz',sep='\t')
samples=pd.read_csv(root/'samples_used.tsv',sep='\t')
allqc=pd.read_csv(R/'public_data/GSE252275_processed/samples_and_qc.tsv',sep='\t')
p=df['P.Value'].to_numpy();assert len(p)==11681 and np.isfinite(p).all()
qcounts=[]
for k in loo.columns[1:]:
    qcounts.append({'sample_removed':k,'n_q_lt_005':int((loo[k]<.05).sum()),'n_q_lt_020':int((loo[k]<.20).sum())})
qt=pd.DataFrame(qcounts).sort_values('n_q_lt_020',ascending=False)
qt=qt.merge(samples,left_on='sample_removed',right_on='sample',how='left',validate='1:1')
assert qt['sample'].notna().all()
qt.to_csv(R/'reports/OMICS_LOO_AUDIT_v0.7.tsv',sep='\t',index=False)
nominal=[]
for alpha in [.05,.0001]:
    n=int((p<alpha).sum());nominal.append({'alpha':alpha,'n':n,'all_null_expectation_if_calibrated':float(len(p)*alpha),'ratio':float(n/(len(p)*alpha))})
families=[]
for name,regex in [('RPL_RPS_symbol_prefix','^(RPL|RPS)'),('MRPL_MRPS_symbol_prefix','^(MRPL|MRPS)')]:
    z=df[df.symbol.fillna('').str.match(regex)]
    families.append({'label':name,'regex':regex,'n':len(z),'mean_t':float(z.t.mean()),'median_t':float(z.t.median()),'mean_logFC':float(z.logFC.mean())})
    z[['gene_id','symbol','logFC','t','P.Value','adj.P.Val']].to_csv(R/f'reports/{name}_AUDIT_v0.7.tsv',sep='\t',index=False)
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'tested_genes':len(p),'nominal_p':nominal,
 'Storey_raw_lambda_diagnostics':[{'lambda':v,'pi0_raw':float((p>v).mean()/(1-v))} for v in [.5,.6,.7,.8,.9]],
 'reference_BH_lt_005':int((df['adj.P.Val']<.05).sum()),'reference_BH_lt_020':int((df['adj.P.Val']<.2).sum()),
 'LOO_q020_min':int(qt.n_q_lt_020.min()),'LOO_q020_max':int(qt.n_q_lt_020.max()),
 'LOO_top':qt[['sample_removed','n_q_lt_020','detected_features','globin_fraction']].head(3).to_dict('records'),
 'all40_max_detected':allqc.loc[allqc.detected_features.idxmax(),['sample','detected_features']].to_dict(),
 'individual_high_globin_LOO':qt[qt.flagged_high_globin][['sample_removed','n_q_lt_020']].to_dict('records'),
 'combined_exclusion3_q020':int((pd.read_csv(root/'exclude3_group.tsv.gz',sep='\t')['adj.P.Val']<.2).sum()),
 'prefix_families':families,
 'interpretation':[
  'Small-p excess is compatible with signal but also with misspecified variance, correlated genes, confounding or sample composition; it does not prove disease biology.',
  'Raw pi0(lambda) is a diagnostic under valid null p-values, not a measured proportion of truly different genes.',
  'The existing full-rank design is estimable, not mathematically underdetermined. Gene-level FDR findings are unstable and no reference BH<0.05 candidates qualify.',
  'Individual leave-one-out high-globin exclusions are distinct from the combined exclusion-of-three branch.',
  'Two influential samples are not automatically low-quality; high detected-feature count alone is not a deletion criterion.',
  'Symbol-prefix averages are descriptive, include whatever annotated features pass the fixed filter, and are not formal enrichment or evidence of reticulocyte contamination.',
  'Review supplied no exact alpha-granule or nuclear-OXPHOS gene membership or test; those reported means cannot be independently reproduced from a named gene set yet.',
  'Transcript direction cannot directly measure mitochondrial translation; keep the candidate mechanism provisional rather than assert its presence or absence.'
 ],
 'input_sha256':{f.relative_to(R).as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in [root/'all20_group.tsv.gz',root/'loo_BH_q.tsv.gz',root/'samples_used.tsv']}}
(R/'reports/OMICS_REVIEW_AUDIT_v0.7.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
