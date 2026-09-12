"""Baseline-only medication/context extension; no endpoints loaded."""
from pathlib import Path
import json,datetime,hashlib,os,duckdb,pandas as pd
os.umask(0o077);R=Path('/root/projects/linezolid_platelet_recovery');B=Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
cfg=json.loads((R/'config/analysis_plan_v0.9.json').read_text())['clinical'];assert hashlib.sha256((R/'config/analysis_plan_v0.9.json').read_bytes()).hexdigest()=='5a510f4ddc03d7ab91f816aff7fa0cae208ae9856c877d0697027092ae5ae5cc'
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'));c.execute('SET threads=6');c.execute("SET memory_limit='12GB'")
def cached(name,q):
 p=R/'cache'/f'{name}.parquet'
 if not p.exists():c.execute(f"COPY ({q}) TO '{p}.partial' (FORMAT PARQUET,COMPRESSION ZSTD)");Path(str(p)+'.partial').replace(p)
 c.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{p}')");print(name,'ready',flush=True)
rx=cfg['beta_lactam']['emar_regex']
cached('phase9_beta_emar',f'''SELECT m.emar_id,e.hadm_id,m.event_txt,m.medication,try_cast(m.charttime AS TIMESTAMP) charttime,try_cast(m.storetime AS TIMESTAMP) storetime
 FROM read_csv('{B/'hosp/emar.csv.gz'}',all_varchar=true) m JOIN phase61_eligible_points e USING(hadm_id)
 WHERE regexp_matches(lower(m.medication),'{rx}') AND try_cast(m.charttime AS TIMESTAMP)<e.t0
 AND try_cast(m.charttime AS TIMESTAMP)>=e.t0-INTERVAL 24 HOUR AND try_cast(m.storetime AS TIMESTAMP)<=e.t0''')
cached('phase9_beta_detail',f'''SELECT d.emar_id,d.route,try_cast(replace(d.dose_given,',','') AS DOUBLE) dose_given
 FROM read_csv('{B/'hosp/emar_detail.csv.gz'}',all_varchar=true) d SEMI JOIN phase9_beta_emar e USING(emar_id)''')
ids=','.join("'"+x+"'" for x in cfg['beta_lactam']['icu_itemids'])
cached('phase9_beta_icu',f'''SELECT e.hadm_id,i.itemid,try_cast(i.starttime AS TIMESTAMP) starttime,try_cast(i.storetime AS TIMESTAMP) storetime
 FROM read_csv('{B/'icu/inputevents.csv.gz'}',all_varchar=true) i JOIN phase61_eligible_points e USING(hadm_id)
 WHERE i.itemid IN ({ids}) AND try_cast(i.amount AS DOUBLE)>0 AND coalesce(i.statusdescription,'')<>'Rewritten'
 AND try_cast(i.starttime AS TIMESTAMP)<e.t0 AND try_cast(i.starttime AS TIMESTAMP)>=e.t0-INTERVAL 24 HOUR AND try_cast(i.storetime AS TIMESTAMP)<=e.t0''')
c.execute('''CREATE OR REPLACE TABLE phase9_medication_context AS
 SELECT e.hadm_id,
 coalesce(m.ufh,false) OR coalesce(i.ufh,false) ufh_active_available7d,
 coalesce(m.lmwh,false) lmwh_available7d,
 coalesce(m.other,false) marrowdrug_available7d,
 EXISTS(SELECT 1 FROM phase9_beta_icu b WHERE b.hadm_id=e.hadm_id) OR
 EXISTS(SELECT 1 FROM phase9_beta_emar b JOIN
 (SELECT emar_id,bool_or(regexp_matches(upper(coalesce(route,'')),'^IV|INTRAVENOUS')) is_iv,
 bool_or(dose_given>0) positive,bool_or(dose_given=0) zero FROM phase9_beta_detail GROUP BY emar_id) d USING(emar_id)
 WHERE b.hadm_id=e.hadm_id AND d.is_iv AND NOT(coalesce(d.zero,false) AND NOT coalesce(d.positive,false))
 AND b.event_txt IN ('Administered','Delayed Administered','Administered in Other Location','Started','Started in Other Location','Delayed Started','Restarted','Partial Administered')) beta_lactam_available24h
 FROM phase61_eligible_points e
 LEFT JOIN LATERAL (SELECT bool_or(drug_group='ufh_including_flush' AND NOT explicit_flush_label) ufh,
 bool_or(drug_group='lmwh') lmwh,bool_or(drug_group='selected_other_drugs') other
 FROM phase3_admin_drugs d WHERE d.hadm_id=e.hadm_id AND d.charttime<e.t0 AND d.charttime>=e.t0-INTERVAL 7 DAY AND d.storetime<=e.t0) m ON true
 LEFT JOIN LATERAL (SELECT bool_or(itemid='225152') ufh FROM phase3_support_inputs d WHERE d.hadm_id=e.hadm_id
 AND d.starttime<e.t0 AND coalesce(d.endtime,d.starttime)>=e.t0-INTERVAL 7 DAY AND d.storetime<=e.t0
 AND amount>0 AND coalesce(statusdescription,'')<>'Rewritten') i ON true''')
ct=pd.read_csv(R/'cache/phase8_context_with_chart.csv',dtype={'hadm_id':str,'subject_id':str});md=c.execute('SELECT * FROM phase9_medication_context').fetchdf()
ct=ct.merge(md,on='hadm_id',validate='1:1');assert len(ct)==1569
audit=[]
for pop,base in [('physical_shared_units','physical_icu'),('culture_shared_units','physical_icu_recent_culture')]:
 counts=pd.crosstab(ct.loc[ct[base],'icu_units'],ct.loc[ct[base],'drug']);shared=counts.index[(counts>0).all(axis=1)]
 ct[pop]=ct[base]&ct.icu_units.isin(shared)
 for drug,g in ct.groupby('drug'):audit.append(dict(population=pop,drug=drug,n=int(g[pop].sum()),shared_units=shared.tolist()))
ct.to_csv(R/'cache/phase9_context.csv',index=False)
medcounts=[]
for drug,g in ct.groupby('drug'):
 for col in ['ufh_active_available7d','lmwh_available7d','marrowdrug_available7d','beta_lactam_available24h']:medcounts.append(dict(drug=drug,variable=col,n=int(g[col].sum())))
# Feasibility only: no discharge sepsis criterion, no post-treatment variables summarized.
feas=c.execute('''WITH e AS (SELECT subject_id,hadm_id,drug,t0,coded_sepsis FROM initiation_audit
 WHERE age>=18 AND initiation_type='no_prior_other_drug' AND started_in_icu
 AND t0>=admittime AND t0<dischtime AND (deathtime IS NULL OR t0<deathtime)),
 known AS (SELECT e.*,l.charttime, l.valuenum,
 (SELECT count(DISTINCT p.valuenum) FROM platelet_measurements p WHERE p.hadm_id=e.hadm_id AND p.charttime=l.charttime AND p.storetime<=e.t0) nvalues
 FROM e JOIN LATERAL (SELECT charttime,valuenum FROM platelet_measurements p WHERE p.hadm_id=e.hadm_id
 AND p.charttime>=e.t0-INTERVAL 24 HOUR AND p.charttime<e.t0 AND p.storetime<=e.t0
 ORDER BY charttime DESC,storetime DESC LIMIT 1) l ON true)
 SELECT drug,count(*) n_baseline_known_low,count(*) FILTER(WHERE coded_sepsis) n_with_discharge_sepsis_code
 FROM known WHERE valuenum<100 AND nvalues=1 GROUP BY drug''').fetchdf().to_dict('records')
def hide(x):
 if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
 if isinstance(x,list):return [hide(v) for v in x]
 if isinstance(x,int) and not isinstance(x,bool) and 0<x<10:return '<10'
 return x
out=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),no_outcomes_loaded=True,populations=audit,medication_counts=medcounts,
 future_code_free_feasibility=feas,feasibility_limits='Counts only; physical location, known VRE, culture timing and source completeness not yet verified for additional admissions. Not a new validated cohort or independent dataset.')
(R/'reports/CLINICAL_CONTEXT_v0.9.json').write_text(json.dumps(hide(out),indent=2));print(json.dumps(hide(out),indent=2))
