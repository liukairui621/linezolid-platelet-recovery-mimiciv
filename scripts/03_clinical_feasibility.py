"""Aggregate feasibility only. No treatment-effect estimates or recovery comparison."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import os
import numpy as np
import pandas as pd
import duckdb

os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery')
C=R/'cache'
variant=os.environ.get('AUDIT_VARIANT','emar')
if variant not in ('emar','combined'): raise ValueError('Unknown audit variant')
con=duckdb.connect(str(C/('audit.duckdb' if variant=='emar' else 'audit_combined.duckdb')))
con.execute("SET threads=6")
con.execute("SET memory_limit='12GB'")

con.execute('''CREATE OR REPLACE TABLE administration_audit AS
 WITH detail AS (
 SELECT emar_id,first(route) FILTER (WHERE route IS NOT NULL) detail_route,
        count(DISTINCT route) detail_route_n,
        bool_or(try_cast(replace(dose_given,',','') AS DOUBLE)>0) dose_positive,
        bool_or(try_cast(replace(dose_given,',','') AS DOUBLE)=0) dose_zero,
        count(*) detail_rows
 FROM target_emar_detail GROUP BY emar_id),
 rx AS (SELECT hadm_id,pharmacy_id,first(route) rx_route,first(drug_class) rx_drug,
               count(DISTINCT route) rx_route_n
        FROM target_prescriptions GROUP BY hadm_id,pharmacy_id)
 SELECT e.*,d.detail_route,d.detail_route_n,d.dose_positive,d.dose_zero,
        r.rx_route,r.rx_drug,r.rx_route_n,
        coalesce(d.detail_route,r.rx_route) route,
        e.event_txt IN ('Administered','Delayed Administered','Administered in Other Location',
          'Started','Started in Other Location','Delayed Started','Restarted','Partial Administered') documented_given
 FROM target_emar e LEFT JOIN detail d USING(emar_id)
 LEFT JOIN rx r USING(hadm_id,pharmacy_id)
''')
con.execute('''CREATE OR REPLACE TABLE systemic_administrations AS
 SELECT * FROM administration_audit
 WHERE documented_given AND hadm_id IS NOT NULL AND charttime IS NOT NULL
 AND coalesce(detail_route_n,0)<=1 AND coalesce(rx_route_n,0)<=1
 AND (rx_drug IS NULL OR rx_drug=drug)
 AND ((drug='vancomycin' AND upper(route) IN ('IV','PB','IV DRIP','IV BOLUS'))
   OR (drug='linezolid' AND upper(route) IN ('IV','PB','IV DRIP','IV BOLUS','PO','NG','PO/NG','ORAL','G TUBE','J TUBE')))
 AND NOT (coalesce(dose_zero,false) AND NOT coalesce(dose_positive,false))
''')

# Combined-source presence and initiation times; duplicate sources are not counted as independent doses.
if variant=='combined':
 con.execute('''CREATE OR REPLACE TABLE emar_systemic_administrations AS SELECT * FROM systemic_administrations''')
 con.execute('''CREATE OR REPLACE TABLE systemic_administrations AS
  SELECT subject_id,hadm_id,drug,charttime,emar_id,'emar' AS data_source FROM emar_systemic_administrations
  UNION ALL
  SELECT subject_id,hadm_id,drug,starttime charttime,
         'icu:'||orderid||':'||cast(starttime AS VARCHAR) emar_id,'inputevents' AS data_source
  FROM target_icu_antibiotics
  WHERE amount>0 AND statusdescription IS DISTINCT FROM 'Rewritten' AND starttime IS NOT NULL
 ''')

# One row per first in-hospital initiation of each drug, preserving switches as a separate stratum.
con.execute('''CREATE OR REPLACE TABLE drug_initiations AS
 SELECT subject_id,hadm_id,drug,min(charttime) t0,count(*) administration_records,
        count(DISTINCT cast(charttime AS DATE)) observed_treatment_days
 FROM systemic_administrations GROUP BY subject_id,hadm_id,drug
''')

# Deduplicate identical platelet results from the same specimen/time; do not average incompatible results.
con.execute('''CREATE OR REPLACE TABLE platelet_measurements AS
 SELECT subject_id,hadm_id,specimen_id,charttime,valuenum,first(storetime) storetime,
        count(*) source_rows
 FROM target_labs
 WHERE itemid IN ('51265','53189') AND valueuom='K/uL' AND valuenum>=0 AND valuenum<5000
 GROUP BY subject_id,hadm_id,specimen_id,charttime,valuenum
''')

con.execute('''CREATE OR REPLACE TABLE initiation_audit AS
 SELECT e.*, p.gender,
        try_cast(p.anchor_age AS INTEGER)+year(e.t0)-try_cast(p.anchor_year AS INTEGER) age,
        p.anchor_year_group,
        try_cast(a.admittime AS TIMESTAMP) admittime,
        try_cast(a.dischtime AS TIMESTAMP) dischtime,
        try_cast(a.deathtime AS TIMESTAMP) deathtime,
        d.coded_sepsis,d.coded_pneumonia,d.coded_marrow_disease,
        o.other_first_time,
        CASE WHEN o.other_first_time<e.t0 THEN 'prior_other_drug'
             WHEN o.other_first_time=e.t0 THEN 'simultaneous'
             ELSE 'no_prior_other_drug' END initiation_type,
        (SELECT count(*) FROM icustays i WHERE i.hadm_id=e.hadm_id)>0 icu_admission,
        (SELECT count(*) FROM icustays i WHERE i.hadm_id=e.hadm_id
          AND e.t0>=try_cast(i.intime AS TIMESTAMP) AND e.t0<try_cast(i.outtime AS TIMESTAMP))>0 started_in_icu,
        (SELECT count(*) FROM icustays i WHERE i.hadm_id=e.hadm_id
          AND e.t0>=try_cast(i.intime AS TIMESTAMP)-INTERVAL 24 HOUR AND e.t0<try_cast(i.outtime AS TIMESTAMP))>0 started_icu_or_24h_before,
        bl.platelets_baseline,bl.platelets_time,bl.baseline_available_at_t0,
        bp.platelets_prior48h,
        cr.creatinine_baseline,
        labs.n_platelet_14d,labs.days_platelet_14d,labs.last_platelet_hours,
        labs.labs_after_day3,labs.labs_after_day7,labs.labs_after_day10,
        tx.icu_transfusion_before48h,tx.icu_transfusion_after14d,
        least(14.0,date_diff('second',e.t0,try_cast(a.dischtime AS TIMESTAMP))/86400.0) hospital_observation_days
 FROM drug_initiations e
 JOIN patients p USING(subject_id)
 JOIN admissions a USING(hadm_id)
 LEFT JOIN diagnosis_flags d USING(hadm_id)
 LEFT JOIN LATERAL (
  SELECT min(x.t0) other_first_time FROM drug_initiations x
  WHERE x.hadm_id=e.hadm_id AND x.drug<>e.drug
 ) o ON true
 LEFT JOIN LATERAL (
  SELECT valuenum platelets_baseline,charttime platelets_time,storetime<=e.t0 baseline_available_at_t0
  FROM platelet_measurements x WHERE x.hadm_id=e.hadm_id
   AND x.charttime<e.t0 AND x.charttime>=e.t0-INTERVAL 24 HOUR
  ORDER BY charttime DESC,storetime DESC LIMIT 1
 ) bl ON true
 LEFT JOIN LATERAL (
  SELECT valuenum platelets_prior48h FROM platelet_measurements x WHERE x.hadm_id=e.hadm_id
   AND x.charttime<e.t0-INTERVAL 24 HOUR AND x.charttime>=e.t0-INTERVAL 72 HOUR
  ORDER BY charttime DESC LIMIT 1
 ) bp ON true
 LEFT JOIN LATERAL (
  SELECT valuenum creatinine_baseline FROM target_labs x WHERE x.hadm_id=e.hadm_id
   AND x.itemid IN ('50912','52546') AND x.valueuom='mg/dL' AND valuenum>0
   AND x.charttime<e.t0 AND x.charttime>=e.t0-INTERVAL 24 HOUR
  ORDER BY charttime DESC LIMIT 1
 ) cr ON true
 LEFT JOIN LATERAL (
  SELECT count(*) n_platelet_14d,count(DISTINCT floor(date_diff('second',e.t0,x.charttime)/86400.0)) days_platelet_14d,
    max(date_diff('second',e.t0,x.charttime)/3600.0) last_platelet_hours,
    count(*) FILTER (WHERE x.charttime>=e.t0+INTERVAL 3 DAY) labs_after_day3,
    count(*) FILTER (WHERE x.charttime>=e.t0+INTERVAL 7 DAY) labs_after_day7,
    count(*) FILTER (WHERE x.charttime>=e.t0+INTERVAL 10 DAY) labs_after_day10
  FROM platelet_measurements x WHERE x.hadm_id=e.hadm_id
   AND x.charttime>=e.t0 AND x.charttime<e.t0+INTERVAL 14 DAY
   AND x.charttime<=try_cast(a.dischtime AS TIMESTAMP)
 ) labs ON true
 LEFT JOIN LATERAL (
  SELECT count(*) FILTER(WHERE x.starttime<e.t0 AND x.starttime>=e.t0-INTERVAL 48 HOUR) icu_transfusion_before48h,
         count(*) FILTER(WHERE x.starttime>=e.t0 AND x.starttime<e.t0+INTERVAL 14 DAY) icu_transfusion_after14d
  FROM platelet_transfusions_icu x WHERE x.hadm_id=e.hadm_id AND x.amount>0
   AND x.statusdescription IS DISTINCT FROM 'Rewritten'
 ) tx ON true
''')

df=con.execute('SELECT * FROM initiation_audit').fetchdf()
assert not df.duplicated(['hadm_id','drug']).any()
assert (df.loc[df.platelets_time.notna(),'platelets_time']<df.loc[df.platelets_time.notna(),'t0']).all()
if variant=='emar':
 assert not con.execute('SELECT emar_id FROM systemic_administrations GROUP BY emar_id HAVING count(*)>1').fetchall()
df=df[df.age>=18].copy()

stages={
 '01_adult_systemic_drug_initiation':pd.Series(True,index=df.index),
 '02_hospitalization_with_icu':df.icu_admission,
 '03_icu_and_coded_sepsis':df.icu_admission & df.coded_sepsis.fillna(False),
 '04_above_and_baseline_platelets24h':df.icu_admission & df.coded_sepsis.fillna(False) & df.platelets_baseline.notna(),
 '05_above_and_baseline_platelets_lt100':df.icu_admission & df.coded_sepsis.fillna(False) & (df.platelets_baseline<100),
 '06_above_and_no_prior_other_drug':df.icu_admission & df.coded_sepsis.fillna(False) & (df.platelets_baseline<100) & (df.initiation_type=='no_prior_other_drug'),
 '07_above_and_no_marrow_code':df.icu_admission & df.coded_sepsis.fillna(False) & (df.platelets_baseline<100) & (df.initiation_type=='no_prior_other_drug') & ~df.coded_marrow_disease.fillna(False),
}
rows=[]
for name,mask in stages.items():
 for drug in ['linezolid','vancomycin']:
  g=df[mask & (df.drug==drug)]
  rows.append({'stage':name,'drug':drug,'admissions':len(g),'patients':g.subject_id.nunique()})
flow=pd.DataFrame(rows)

key=df[stages['05_above_and_baseline_platelets_lt100']].copy()
summ=[]
for (drug,kind),g in key.groupby(['drug','initiation_type']):
 row={'drug':drug,'initiation_type':kind,'admissions':len(g),'patients':g.subject_id.nunique(),
      'pneumonia_code':int(g.coded_pneumonia.fillna(False).sum()),
      'started_in_icu':int(g.started_in_icu.sum()),'started_icu_or_24h_before':int(g.started_icu_or_24h_before.sum()),
      'baseline_result_available':int(g.baseline_available_at_t0.fillna(False).sum()),
      'baseline_creatinine_available':int(g.creatinine_baseline.notna().sum()),
      'earlier_platelet_measurement_available':int(g.platelets_prior48h.notna().sum()),
      'followup_at_least_two_platelets':int((g.n_platelet_14d>=2).sum()),
      'icu_transfusion_before48h':int((g.icu_transfusion_before48h>0).sum()),
      'icu_transfusion_after14d':int((g.icu_transfusion_after14d>0).sum()),
      'marrow_code':int(g.coded_marrow_disease.fillna(False).sum())}
 for col in ['age','platelets_baseline','creatinine_baseline','hospital_observation_days','days_platelet_14d']:
  for quant,label in [(0.25,'q25'),(0.5,'median'),(0.75,'q75')]:
   row[col+'_'+label]=float(g[col].quantile(quant)) if g[col].notna().any() else None
 for day in [3,7,10]:
  observed=g[(g.dischtime>=g.t0+pd.to_timedelta(day+1,unit='D')) & (g.deathtime.isna() | (g.deathtime>=g.t0+pd.to_timedelta(day+1,unit='D')))]
  row[f'in_hospital_at_day{day}_plus1']=len(observed)
  row[f'has_platelets_after_day{day}_among_observed']=int((observed[f'labs_after_day{day}']>0).sum())
 summ.append(row)

def suppress_records(records):
 for row in records:
  for k,v in list(row.items()):
   if isinstance(v,(int,np.integer)) and not isinstance(v,bool) and 0<int(v)<10:row[k]='<10'
   elif isinstance(v,(float,np.floating)):
    row[k]=round(float(v),2) if np.isfinite(v) else None
 return records

flow_safe=suppress_records(flow.to_dict('records'))
summary_safe=suppress_records(summ)
qc={}
for name,sql in {
 'route_coverage': '''SELECT drug,documented_given,route,count(*) n FROM administration_audit GROUP BY ALL ORDER BY drug,n DESC''',
 'given_dose_status': '''SELECT drug,dose_positive,dose_zero,count(*) n FROM administration_audit WHERE documented_given GROUP BY ALL ORDER BY n DESC''',
 'route_conflicts': '''SELECT drug,count(*) FILTER(WHERE detail_route IS NOT NULL AND rx_route IS NOT NULL AND detail_route IN ('PO','NG') AND rx_route IN ('IV','PB')) n FROM administration_audit GROUP BY drug''',
 'unit_or_value_exclusions': '''SELECT itemid,valueuom,count(*) n FROM target_labs WHERE itemid IN ('51265','53189') AND (valueuom IS DISTINCT FROM 'K/uL' OR valuenum<0 OR valuenum>=5000 OR valuenum IS NULL) GROUP BY ALL''',
}.items():
 qc[name]=suppress_records(con.execute(sql).fetchdf().to_dict('records'))

result={'generated_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'source_variant':variant,
 'stage':'Feasibility screening, not final eligibility or treatment-effect analysis',
 'definitions':{'baseline_platelets':'latest measured value in [-24h,0); K/uL; storetime availability audited separately',
  'sepsis':'discharge ICD screening proxy; not Sepsis-3; timing not established',
  'initial_treatment':'no earlier documented systemic use of the other target drug in the same hospitalization; not lifetime new use',
  'population':'adult hospitalizations containing ICU care; exact start-in-ICU strata shown separately',
  'transfusion':'ICU/OR/PACU inputevents only; completeness outside ICU is not established'},
 'flow':flow_safe,'stratum_coverage':summary_safe,'qc':qc,
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
suffix='' if variant=='emar' else '_combined'
(R/f'reports/clinical_feasibility{suffix}.json').write_text(json.dumps(result,indent=2,allow_nan=False))
pd.DataFrame(flow_safe).to_csv(R/f'reports/cohort_flow{suffix}.csv',index=False)
pd.DataFrame(summary_safe).to_csv(R/f'reports/stratum_coverage{suffix}.csv',index=False)
print(json.dumps({'flow':flow_safe,'stratum_coverage':summary_safe},indent=2,allow_nan=False))
print('CLINICAL_FEASIBILITY_COMPLETE',flush=True)
