"""Fresh complete-admission outcome sources for frozen v0.11 membership, server only."""
from pathlib import Path
import os,json,hashlib,datetime,platform,shutil
import duckdb,pandas as pd
os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery');B=Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1');C=R/'cache'
plan=json.loads((R/'config/analysis_plan_v0.11.json').read_text())
for p,h in {**plan['frozen_inputs'],**plan['frozen_definitions']}.items():assert hashlib.sha256((R/p).read_bytes()).hexdigest()==h,p
d=pd.read_csv(C/'phase10_baseline.csv',dtype={'hadm_id':str,'subject_id':str})
s=d[d.routine_bacterial_culture72h_shared].copy();assert len(s)==2593 and s.hadm_id.nunique()==2593
primary=s[~s.lmwh_available7d];ct=pd.crosstab(primary.icu_unit,primary.drug);assert (ct>0).all().all()
assert len(primary)==2557
c=duckdb.connect(str(C/'phase11.duckdb'));c.execute("SET threads=4");c.execute("SET memory_limit='12GB'")
c.register('selected_frame',s);c.execute('CREATE OR REPLACE TABLE phase11_selected AS SELECT subject_id,hadm_id,drug,try_cast(t0 AS TIMESTAMP) t0,known_platelets,try_cast(known_platelets_time AS TIMESTAMP) known_platelets_time,try_cast(known_platelets_store AS TIMESTAMP) known_platelets_store,NOT lmwh_available7d primary_population FROM selected_frame')
def csv(p):return f"read_csv('{B/p}',header=true,all_varchar=true)"
queries={}
def extract(name,q,source):
 p=C/f'{name}.parquet';queries[name]=dict(query=q,raw_files=[dict(path=str(B/x),bytes=(B/x).stat().st_size,mtime_ns=(B/x).stat().st_mtime_ns) for x in source])
 if not p.exists():
  print('EXTRACT_START',name,flush=True);tmp=Path(str(p)+'.partial');c.execute(f"COPY ({q}) TO '{tmp}' (FORMAT PARQUET,COMPRESSION ZSTD)");tmp.replace(p)
 c.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{p}')")
 print('EXTRACT_READY',name,c.execute(f'SELECT count(*) FROM {name}').fetchone()[0],flush=True)
extract('phase11_admissions',f"SELECT a.subject_id,a.hadm_id,try_cast(admittime AS TIMESTAMP) admittime,try_cast(dischtime AS TIMESTAMP) dischtime,try_cast(deathtime AS TIMESTAMP) deathtime,hospital_expire_flag FROM {csv('hosp/admissions.csv.gz')} a SEMI JOIN phase11_selected USING(hadm_id)",['hosp/admissions.csv.gz'])
extract('phase11_transfers',f"SELECT t.hadm_id,t.careunit,try_cast(intime AS TIMESTAMP) intime,try_cast(outtime AS TIMESTAMP) outtime FROM {csv('hosp/transfers.csv.gz')} t SEMI JOIN phase11_selected USING(hadm_id)",['hosp/transfers.csv.gz'])
extract('phase11_icustays',f"SELECT i.hadm_id,try_cast(intime AS TIMESTAMP) intime,try_cast(outtime AS TIMESTAMP) outtime FROM {csv('icu/icustays.csv.gz')} i SEMI JOIN phase11_selected USING(hadm_id)",['icu/icustays.csv.gz'])
extract('phase11_platelet_inputs',f"SELECT t.hadm_id,t.itemid,try_cast(starttime AS TIMESTAMP) starttime,try_cast(endtime AS TIMESTAMP) endtime,try_cast(amount AS DOUBLE) amount,amountuom,statusdescription FROM {csv('icu/inputevents.csv.gz')} t SEMI JOIN phase11_selected USING(hadm_id) WHERE itemid IN ('225170','226369','227071')",['icu/inputevents.csv.gz'])
extract('phase11_platelet_procedures',f"SELECT p.hadm_id,try_cast(p.chartdate AS DATE) chartdate,p.icd_code,p.icd_version FROM {csv('hosp/procedures_icd.csv.gz')} p JOIN {csv('hosp/d_icd_procedures.csv.gz')} d USING(icd_code,icd_version) SEMI JOIN phase11_selected USING(hadm_id) WHERE lower(d.long_title) LIKE 'transfusion of%platelet%'",['hosp/procedures_icd.csv.gz','hosp/d_icd_procedures.csv.gz'])
extract('phase11_platelet_emar',f"SELECT m.hadm_id,m.emar_id,m.medication,m.event_txt,try_cast(charttime AS TIMESTAMP) charttime,try_cast(storetime AS TIMESTAMP) storetime FROM {csv('hosp/emar.csv.gz')} m SEMI JOIN phase11_selected USING(hadm_id) WHERE regexp_matches(lower(medication),'platelet|thrombocyte')",['hosp/emar.csv.gz'])
extract('phase11_platelet_labs',f"SELECT l.subject_id,l.hadm_id,l.labevent_id,l.specimen_id,l.itemid,try_cast(charttime AS TIMESTAMP) charttime,try_cast(storetime AS TIMESTAMP) storetime,try_cast(valuenum AS DOUBLE) valuenum,l.valueuom FROM {csv('hosp/labevents.csv.gz')} l SEMI JOIN phase11_selected USING(hadm_id) WHERE itemid IN ('51265','53189')",['hosp/labevents.csv.gz'])
c.execute('CREATE OR REPLACE TABLE phase11_platelet_canonical AS SELECT subject_id,hadm_id,specimen_id,charttime,valuenum,min(storetime) storetime,count(*) source_rows FROM phase11_platelet_labs WHERE valueuom=\'K/uL\' AND valuenum>=0 AND valuenum<5000 GROUP BY subject_id,hadm_id,specimen_id,charttime,valuenum')
drift=c.execute('''SELECT count(*) FROM phase11_selected e LEFT JOIN LATERAL(SELECT charttime,min(valuenum) valuenum,count(DISTINCT valuenum) nvals FROM phase11_platelet_canonical x WHERE x.hadm_id=e.hadm_id AND x.charttime<e.t0 AND x.charttime>=e.t0-INTERVAL 24 HOUR AND x.storetime<=e.t0 GROUP BY charttime ORDER BY charttime DESC LIMIT 1)b ON true WHERE b.nvals IS DISTINCT FROM 1 OR b.valuenum IS DISTINCT FROM e.known_platelets OR b.charttime IS DISTINCT FROM e.known_platelets_time''').fetchone()[0]
assert drift==0
c.execute(f"ATTACH '{C/'audit_combined.duckdb'}' AS prior (READ_ONLY)")
exit_drift=c.execute('''SELECT count(*) FROM phase11_admissions a JOIN prior.phase10_baseline e USING(hadm_id) WHERE a.admittime IS DISTINCT FROM e.admittime OR a.dischtime IS DISTINCT FROM e.dischtime OR a.deathtime IS DISTINCT FROM e.deathtime''').fetchone()[0];assert exit_drift==0
dose_drift=c.execute('''SELECT count(*) FROM phase11_selected e LEFT JOIN LATERAL(SELECT min(charttime) firsttime FROM prior.systemic_administrations x WHERE x.hadm_id=e.hadm_id AND x.drug=e.drug)a ON true WHERE a.firsttime IS DISTINCT FROM e.t0''').fetchone()[0];assert dose_drift==0
rows=[dict(source=name,rows=c.execute(f'SELECT count(*) FROM {name}').fetchone()[0]) for name in queries]
qc=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=platform.python_version(),duckdb=duckdb.__version__,admissions=2593,primary_admissions=2557,baseline_value_time_drift=drift,admission_interval_drift=exit_drift,first_dose_time_drift=dose_drift,sources=rows,limits='Fresh full-admission scans. ICU input and dated procedures do not guarantee hospital-wide transfusion ascertainment; medication-name search is not a complete blood-bank registry.',source_queries=queries)
(R/'reports/OUTCOME_SOURCE_AUDIT_v0.11.json').write_text(json.dumps(qc,indent=2))
print(json.dumps({k:v for k,v in qc.items() if k!='source_queries'},indent=2))
