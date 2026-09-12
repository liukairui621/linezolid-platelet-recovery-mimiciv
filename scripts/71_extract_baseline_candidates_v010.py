"""New baseline candidates and complete infection/location sources. No recovery endpoints loaded."""
from pathlib import Path
import os,json,datetime,hashlib
import duckdb
os.umask(0o077);R=Path('/root/projects/linezolid_platelet_recovery');B=Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
plan=R/'config/analysis_plan_v0.10.json';assert hashlib.sha256(plan.read_bytes()).hexdigest()=='b14a735f7c58771faa8b9a766d8cf329f1e0193e02011d1e6f5898803fd361b6'
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'));c.execute('SET threads=6');c.execute("SET memory_limit='12GB'")
def rows(q):return c.execute(q).fetchdf().to_dict('records')
def hide(x):
 if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
 if isinstance(x,list):return [hide(v) for v in x]
 if isinstance(x,int) and not isinstance(x,bool) and 0<x<10:return '<10'
 return x
def cache(name,q):
 p=R/'cache'/f'{name}.parquet'
 if not p.exists():
  c.execute(f"COPY ({q}) TO '{p}.partial' (FORMAT PARQUET,COMPRESSION ZSTD)");Path(str(p)+'.partial').replace(p)
 c.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{p}')");print(name,'ready',flush=True)
c.execute('''CREATE OR REPLACE TABLE phase10_base_pool AS
 SELECT subject_id,hadm_id,drug,t0,age,gender,anchor_year_group,admittime,dischtime,deathtime,coalesce(coded_sepsis,false) coded_sepsis
 FROM initiation_audit WHERE age>=18 AND initiation_type='no_prior_other_drug' AND started_in_icu
 AND t0>=admittime AND t0<dischtime AND (deathtime IS NULL OR t0<deathtime)''')
c.execute('''CREATE OR REPLACE TABLE phase10_base_audit AS
 SELECT e.*,l.charttime known_platelets_time,CASE WHEN l.nval=1 THEN l.val ELSE NULL END known_platelets,l.storetime known_platelets_store,
 l.nval>1 baseline_ambiguous,p.charttime legacy_known_time,CASE WHEN pc.nval=1 THEN p.valuenum ELSE NULL END legacy_known_platelets
 FROM phase10_base_pool e
 LEFT JOIN LATERAL(SELECT charttime,min(valuenum) val,count(DISTINCT valuenum) nval,min(storetime) storetime FROM target_labs
 WHERE hadm_id=e.hadm_id AND itemid IN ('51265','53189') AND valueuom='K/uL' AND valuenum>=0 AND valuenum<5000
 AND charttime>=e.t0-INTERVAL 24 HOUR AND charttime<e.t0 AND storetime<=e.t0 GROUP BY charttime ORDER BY charttime DESC LIMIT 1) l ON true
 LEFT JOIN LATERAL(SELECT charttime,valuenum FROM platelet_measurements WHERE hadm_id=e.hadm_id
 AND charttime>=e.t0-INTERVAL 24 HOUR AND charttime<e.t0 AND storetime<=e.t0 ORDER BY charttime DESC,storetime DESC LIMIT 1) p ON true
 LEFT JOIN LATERAL(SELECT count(DISTINCT valuenum) nval FROM platelet_measurements WHERE hadm_id=e.hadm_id AND charttime=p.charttime AND storetime<=e.t0) pc ON true''')
c.execute('''CREATE OR REPLACE TABLE phase10_candidates AS SELECT subject_id,hadm_id,drug,t0,age,gender,anchor_year_group,admittime,dischtime,deathtime,coded_sepsis,known_platelets,known_platelets_time,known_platelets_store
 FROM phase10_base_audit WHERE known_platelets<100''')
assert c.execute('SELECT count(*)=count(DISTINCT hadm_id) FROM phase10_candidates').fetchone()[0]
assert c.execute('SELECT count(*) FROM phase10_candidates WHERE known_platelets_time>=t0 OR known_platelets_store>t0').fetchone()[0]==0
cache('phase10_transfers',f'''SELECT t.hadm_id,t.careunit,try_cast(t.intime AS TIMESTAMP) intime,try_cast(t.outtime AS TIMESTAMP) outtime
 FROM read_csv('{B/'hosp/transfers.csv.gz'}',all_varchar=true) t SEMI JOIN phase10_candidates e USING(hadm_id)''')
cache('phase10_microbiology',f'''WITH joined AS (
 SELECT e.hadm_id,m.microevent_id,m.micro_specimen_id,
 coalesce(try_cast(m.charttime AS TIMESTAMP),try_cast(m.chartdate AS TIMESTAMP)+INTERVAL 1 DAY) sample_available_by,
 coalesce(try_cast(m.storetime AS TIMESTAMP),try_cast(m.storedate AS TIMESTAMP)+INTERVAL 1 DAY) result_available_by,
 try_cast(m.charttime AS TIMESTAMP) IS NULL sample_date_only,try_cast(m.storetime AS TIMESTAMP) IS NULL result_date_only,
 m.spec_type_desc,m.test_name,
 CASE WHEN coalesce(try_cast(m.storetime AS TIMESTAMP),try_cast(m.storedate AS TIMESTAMP)+INTERVAL 1 DAY)<=e.t0 THEN m.org_name END known_org_name,
 CASE WHEN coalesce(try_cast(m.storetime AS TIMESTAMP),try_cast(m.storedate AS TIMESTAMP)+INTERVAL 1 DAY)<=e.t0 THEN m.ab_name END known_ab_name,
 CASE WHEN coalesce(try_cast(m.storetime AS TIMESTAMP),try_cast(m.storedate AS TIMESTAMP)+INTERVAL 1 DAY)<=e.t0 THEN m.interpretation END known_interpretation,
 m.hadm_id IS NULL linked_missing_hadm,count(*) OVER(PARTITION BY m.microevent_id) linked_admissions
 FROM read_csv('{B/'hosp/microbiologyevents.csv.gz'}',all_varchar=true) m JOIN phase10_candidates e
 ON m.subject_id=e.subject_id AND (m.hadm_id=e.hadm_id OR (m.hadm_id IS NULL AND
 coalesce(try_cast(m.charttime AS TIMESTAMP),try_cast(m.chartdate AS TIMESTAMP))>=e.admittime AND
 coalesce(try_cast(m.charttime AS TIMESTAMP),try_cast(m.chartdate AS TIMESTAMP))<e.dischtime))
 WHERE coalesce(try_cast(m.charttime AS TIMESTAMP),try_cast(m.chartdate AS TIMESTAMP)+INTERVAL 1 DAY)>=e.t0-INTERVAL 7 DAY
 AND coalesce(try_cast(m.charttime AS TIMESTAMP),try_cast(m.chartdate AS TIMESTAMP)+INTERVAL 1 DAY)<=e.t0)
 SELECT * FROM joined''')
qc=rows('''SELECT drug,count(*) baseline_pool,count(*) FILTER(WHERE known_platelets<100) raw_known_low,count(*) FILTER(WHERE legacy_known_platelets<100) legacy_known_low,
 count(*) FILTER(WHERE (known_platelets<100) IS DISTINCT FROM (legacy_known_platelets<100)) threshold_or_missing_disagreements,
 count(*) FILTER(WHERE known_platelets<100 AND coded_sepsis) raw_low_with_discharge_code FROM phase10_base_audit GROUP BY drug''')
counts=rows('''SELECT e.drug,count(*) candidates,count(DISTINCT subject_id) unique_patients,count(*) FILTER(WHERE coalesce(loc.in_old,false)) old_cohort,
 count(*) FILTER(WHERE EXISTS(SELECT 1 FROM phase10_transfers t WHERE t.hadm_id=e.hadm_id)) transfers_covered,
 count(*) FILTER(WHERE EXISTS(SELECT 1 FROM phase10_microbiology t WHERE t.hadm_id=e.hadm_id)) prior7d_microbiology_recorded
 FROM phase10_candidates e LEFT JOIN LATERAL(SELECT true in_old FROM phase61_eligible_points o WHERE o.hadm_id=e.hadm_id) loc ON true GROUP BY drug''')
tests=rows('''SELECT coalesce(test_name,'[MISSING]') test_name,count(DISTINCT hadm_id) admissions,count(*) source_rows FROM phase10_microbiology WHERE linked_admissions=1 GROUP BY test_name ORDER BY admissions DESC''')
out=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),no_recovery_outcomes_loaded=True,plan_sha256=hashlib.sha256(plan.read_bytes()).hexdigest(),baseline_reproduction=qc,source_counts=counts,ambiguous_microbiology_rows=c.execute('SELECT count(*) FROM phase10_microbiology WHERE linked_admissions>1').fetchone()[0],test_name_inventory=tests)
(R/'reports/BASELINE_EXTENSION_INPUTS_v0.10.json').write_text(json.dumps(hide(out),indent=2));(R/'cache/internal_reports/BASELINE_EXTENSION_INPUTS_v0.10_exact.json').write_text(json.dumps(out,indent=2))
c.execute(f"COPY phase10_candidates TO '{R/'cache/phase10_candidates.csv'}' (HEADER,DELIMITER ',')")
print(json.dumps(hide(out),indent=2))
