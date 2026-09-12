"""Decision-time eligibility and switch risk-set feasibility; no outcome/effect analysis.
Patient-level tables stay on the controlled server. Prior discharge codes are history
proxies; current-admission discharge codes remain screening labels only.
"""
from pathlib import Path
import os,json,datetime,hashlib
import duckdb
os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery')
B=Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'))
c.execute("SET threads=6");c.execute("SET memory_limit='12GB'")
def cached(name,query):
    p=R/'cache'/f'{name}.parquet'
    if not p.exists():c.execute(f"COPY ({query}) TO '{p}' (FORMAT PARQUET,COMPRESSION ZSTD)")
    c.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{p}')")
def rows(sql):
    cur=c.execute(sql);cols=[x[0] for x in cur.description]
    return [dict(zip(cols,r)) for r in cur.fetchall()]
def publish(data):
    if isinstance(data,dict):return {k:publish(v) for k,v in data.items()}
    if isinstance(data,list):return [publish(v) for v in data]
    if isinstance(data,bool):return data
    if isinstance(data,int):return '<10' if 0<data<10 else data
    if isinstance(data,float):return round(data,3)
    return data
c.execute('''CREATE OR REPLACE TABLE phase2_admissions AS
 SELECT DISTINCT subject_id,hadm_id,admittime,dischtime FROM initiation_audit
 WHERE age>=18 AND icu_admission AND coded_sepsis''')
cached('phase2_microbiology',f"""
 SELECT e.hadm_id,m.micro_specimen_id,
 coalesce(try_cast(m.charttime AS TIMESTAMP),try_cast(m.chartdate AS TIMESTAMP)+INTERVAL 1 DAY) sample_available_by,
 coalesce(try_cast(m.storetime AS TIMESTAMP),try_cast(m.storedate AS TIMESTAMP)+INTERVAL 1 DAY) result_available_by,
 m.charttime IS NULL sample_date_only,m.storetime IS NULL result_date_only,
 m.spec_type_desc,m.test_name,m.org_name,m.ab_name,m.interpretation
 FROM read_csv('{B/'hosp/microbiologyevents.csv.gz'}',all_varchar=true) m
 JOIN phase2_admissions e ON m.subject_id=e.subject_id AND
 (m.hadm_id=e.hadm_id OR (m.hadm_id IS NULL AND
 coalesce(try_cast(m.charttime AS TIMESTAMP),try_cast(m.chartdate AS TIMESTAMP))>=e.admittime AND
 coalesce(try_cast(m.charttime AS TIMESTAMP),try_cast(m.chartdate AS TIMESTAMP))<e.dischtime))
""")
cached('phase2_transfers',f"""SELECT t.hadm_id,t.careunit,try_cast(t.intime AS TIMESTAMP) intime,
 try_cast(t.outtime AS TIMESTAMP) outtime FROM read_csv('{B/'hosp/transfers.csv.gz'}',all_varchar=true) t
 SEMI JOIN phase2_admissions e USING(hadm_id)""")
cached('phase2_platelet_procedures',f"""SELECT p.hadm_id,try_cast(p.chartdate AS DATE) chartdate,p.icd_code,p.icd_version
 FROM read_csv('{B/'hosp/procedures_icd.csv.gz'}',all_varchar=true) p
 JOIN read_csv('{B/'hosp/d_icd_procedures.csv.gz'}',all_varchar=true) d USING(icd_code,icd_version)
 SEMI JOIN phase2_admissions e USING(hadm_id)
 WHERE lower(d.long_title) LIKE 'transfusion of%platelet%'""")
cached('phase2_platelet_emar',f"""SELECT m.hadm_id,m.medication,m.event_txt,try_cast(m.charttime AS TIMESTAMP) charttime
 FROM read_csv('{B/'hosp/emar.csv.gz'}',all_varchar=true) m SEMI JOIN phase2_admissions e USING(hadm_id)
 WHERE regexp_matches(lower(m.medication),'platelet|thrombocyte')""")
c.execute('''CREATE OR REPLACE TABLE phase2_prior_marrow AS
 SELECT e.hadm_id,bool_or(d.coded_marrow_disease) prior_recorded_marrow
 FROM phase2_admissions e JOIN admissions a ON e.subject_id=a.subject_id
 AND try_cast(a.dischtime AS TIMESTAMP)<e.admittime
 JOIN diagnosis_flags d ON d.hadm_id=a.hadm_id GROUP BY e.hadm_id''')
c.execute('''CREATE OR REPLACE TEMP TABLE phase2_icu AS
 SELECT hadm_id,stay_id,try_cast(intime AS TIMESTAMP) intime,try_cast(outtime AS TIMESTAMP) outtime
 FROM icustays''')

def enrich(source,destination):
    # source must expose hadm_id,t0. Strictly historical specimen AND stored result.
    c.execute(f'''CREATE OR REPLACE TABLE {destination} AS
    SELECT e.*,CASE WHEN pc.n_distinct=1 THEN p.valuenum ELSE NULL END known_platelets,
      pc.n_distinct>1 ambiguous_latest_platelets,p.charttime known_platelets_time,p.storetime known_platelets_store,
      pold.valuenum earlier_platelets,cr.valuenum known_creatinine,
      coalesce(m.vre_before,false) vre_before,coalesce(m.culture72h,false) culture72h,
      coalesce(m.respiratory7d,false) respiratory7d,
      coalesce(h.prior_recorded_marrow,false) prior_recorded_marrow,
      v.charttime last_vancomycin,
      EXISTS(SELECT 1 FROM phase2_icu i WHERE i.hadm_id=e.hadm_id AND e.t0>=i.intime AND e.t0<i.outtime) decision_in_icu
    FROM {source} e
    LEFT JOIN LATERAL (SELECT * FROM platelet_measurements p WHERE p.hadm_id=e.hadm_id
       AND p.charttime>=e.t0-INTERVAL 24 HOUR AND p.charttime<e.t0 AND p.storetime<=e.t0
       ORDER BY p.charttime DESC,p.storetime DESC,p.valuenum LIMIT 1) p ON true
    LEFT JOIN LATERAL (SELECT count(DISTINCT x.valuenum) n_distinct FROM platelet_measurements x
       WHERE x.hadm_id=e.hadm_id AND x.charttime=p.charttime AND x.storetime<=e.t0) pc ON true
    LEFT JOIN LATERAL (SELECT * FROM platelet_measurements p WHERE p.hadm_id=e.hadm_id
       AND p.charttime>=e.t0-INTERVAL 72 HOUR AND p.charttime<e.t0-INTERVAL 24 HOUR AND p.storetime<=e.t0
       ORDER BY p.charttime DESC,p.storetime DESC,p.valuenum LIMIT 1) pold ON true
    LEFT JOIN LATERAL (SELECT * FROM target_labs p WHERE p.hadm_id=e.hadm_id
       AND p.itemid IN ('50912','52546') AND p.valuenum>0
       AND p.charttime>=e.t0-INTERVAL 24 HOUR AND p.charttime<e.t0 AND p.storetime<=e.t0
       ORDER BY p.charttime DESC,p.storetime DESC,p.valuenum LIMIT 1) cr ON true
    LEFT JOIN LATERAL (SELECT bool_or(upper(m.org_name) LIKE '%ENTEROCOCCUS%' AND upper(m.ab_name)='VANCOMYCIN'
       AND upper(m.interpretation)='R' AND m.result_available_by<=e.t0) vre_before,
       bool_or(m.sample_available_by>=e.t0-INTERVAL 72 HOUR) culture72h,
       bool_or(regexp_matches(upper(m.spec_type_desc),'SPUTUM|BRONCH|TRACHE|RESPIRATORY|LUNG')) respiratory7d
       FROM phase2_microbiology m WHERE m.hadm_id=e.hadm_id
       AND m.sample_available_by>=e.t0-INTERVAL 7 DAY AND m.sample_available_by<=e.t0) m ON true
    LEFT JOIN phase2_prior_marrow h ON h.hadm_id=e.hadm_id
    LEFT JOIN LATERAL (SELECT charttime FROM systemic_administrations v WHERE v.hadm_id=e.hadm_id
       AND v.drug='vancomycin' AND v.charttime<e.t0 ORDER BY charttime DESC LIMIT 1) v ON true
    ''')
c.execute('''CREATE OR REPLACE TABLE phase2_initial_source AS
 SELECT * FROM initiation_audit WHERE age>=18 AND icu_admission AND coded_sepsis''')
enrich('phase2_initial_source','phase2_decisions')
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'definitions':{'baseline':'Latest specimen within prior 24h with storetime <= decision; no future lab results. Conflicting values at latest specimen time are set missing, not averaged or selected for eligibility.',
 'infection':'Current discharge sepsis code remains a screening proxy, not validated Sepsis-3. Respiratory specimen is not proof of pneumonia.',
 'microbiology':'Last known update <= decision; date-only records conservatively available by next midnight. Missing hadm linked by subject and hospital interval.',
 'switch':'At least 24h since first observed vancomycin; most recent vancomycin in prior 48h; in ICU; known platelets <100; no recorded prior VRE in 7d.',
 'controls':'At same elapsed time since first vancomycin as index switch, same anchor-year group, admission elapsed time within 48h; alive/in hospital/in ICU, no linezolid yet, vancomycin in previous 48h, known platelets <100, no known VRE.',
 'interpretation':'Candidate risk sets only. Not matched/weighted causal comparators; antibiotic continuation from past dose is provisional, esp renal dosing. No future survival or treatment-duration restriction.'}}
out['baseline_reclassification']=rows('''SELECT drug,initiation_type,count(*) screened,
 count(*) FILTER(WHERE platelets_baseline<100) previous_specimen_based_low,
 count(*) FILTER(WHERE known_platelets<100) decision_known_low,
 count(*) FILTER(WHERE platelets_baseline<100 AND (known_platelets>=100 OR known_platelets IS NULL)) previous_low_not_confirmed,
 count(*) FILTER(WHERE known_platelets<100 AND (platelets_baseline>=100 OR platelets_baseline IS NULL)) newly_low_by_known_result
 FROM phase2_decisions GROUP BY ALL''')
out['decision_cohorts']=rows('''SELECT drug,initiation_type,count(*) admissions,count(DISTINCT subject_id) patients,
 count(*) FILTER(WHERE decision_in_icu) decision_in_icu,
 count(*) FILTER(WHERE vre_before) known_vre,
 count(*) FILTER(WHERE culture72h) recent_culture,
 count(*) FILTER(WHERE respiratory7d) respiratory_specimen,
 count(*) FILTER(WHERE prior_recorded_marrow) historical_marrow_code,
 count(*) FILTER(WHERE coded_marrow_disease) current_discharge_marrow_code,
 median(age) median_age,median(known_platelets) median_platelets,median(known_creatinine) median_creatinine,
 count(*) FILTER(WHERE earlier_platelets IS NOT NULL) prior_platelet_trend_available,
 median(known_platelets-earlier_platelets) median_platelet_change,
 median(date_diff('second',admittime,t0)/86400.0) median_days_from_admission
 FROM phase2_decisions WHERE known_platelets<100 GROUP BY ALL''')
c.execute('''CREATE OR REPLACE TABLE phase2_switches AS SELECT *, t0-other_first_time switch_delay
 FROM phase2_decisions WHERE drug='linezolid' AND initiation_type='prior_other_drug'
 AND t0>=other_first_time+INTERVAL 24 HOUR AND last_vancomycin>=t0-INTERVAL 48 HOUR
 AND decision_in_icu AND known_platelets<100 AND NOT vre_before''')
out['switch_flow']=rows('''SELECT count(*) all_prior_van_known_low,
 count(*) FILTER(WHERE t0>=other_first_time+INTERVAL 24 HOUR) at_least_24h,
 count(*) FILTER(WHERE t0>=other_first_time+INTERVAL 24 HOUR AND last_vancomycin>=t0-INTERVAL 48 HOUR) recent_van,
 count(*) FILTER(WHERE t0>=other_first_time+INTERVAL 24 HOUR AND last_vancomycin>=t0-INTERVAL 48 HOUR AND decision_in_icu) in_icu,
 count(*) FILTER(WHERE t0>=other_first_time+INTERVAL 24 HOUR AND last_vancomycin>=t0-INTERVAL 48 HOUR AND decision_in_icu AND NOT vre_before) without_vre
 FROM phase2_decisions WHERE drug='linezolid' AND initiation_type='prior_other_drug' AND known_platelets<100''')
c.execute('''CREATE OR REPLACE TABLE phase2_risk_candidates AS
 SELECT s.hadm_id index_hadm,s.subject_id index_subject,v.subject_id,v.hadm_id,
 v.t0+s.switch_delay t0,v.t0 vancomycin_start,v.admittime,v.dischtime,v.deathtime,v.age,v.gender,
 v.anchor_year_group,s.respiratory7d index_respiratory7d
 FROM phase2_switches s JOIN phase2_decisions v ON v.drug='vancomycin'
 AND v.initiation_type='no_prior_other_drug' AND v.subject_id<>s.subject_id
 AND v.anchor_year_group=s.anchor_year_group
 AND abs(date_diff('second',v.admittime,v.t0)-date_diff('second',s.admittime,s.other_first_time))<=172800
 WHERE v.t0+s.switch_delay<v.dischtime AND (v.deathtime IS NULL OR v.t0+s.switch_delay<v.deathtime)
 AND (v.other_first_time IS NULL OR v.other_first_time>v.t0+s.switch_delay)
 AND EXISTS(SELECT 1 FROM phase2_icu i WHERE i.hadm_id=v.hadm_id
 AND v.t0+s.switch_delay>=i.intime AND v.t0+s.switch_delay<i.outtime)
 AND EXISTS(SELECT 1 FROM systemic_administrations a WHERE a.hadm_id=v.hadm_id AND a.drug='vancomycin'
 AND a.charttime<v.t0+s.switch_delay AND a.charttime>=v.t0+s.switch_delay-INTERVAL 48 HOUR)''')
print(json.dumps({'stage':'risk_candidates','rows':c.execute('select count(*) from phase2_risk_candidates').fetchone()[0]}),flush=True)
enrich('phase2_risk_candidates','phase2_risk_enriched')
c.execute('''CREATE OR REPLACE TABLE phase2_risk_eligible AS SELECT * FROM phase2_risk_enriched
 WHERE known_platelets<100 AND NOT vre_before AND decision_in_icu''')
out['risk_sets']=rows('''WITH n AS (SELECT s.hadm_id,count(r.hadm_id) n,
 count(r.hadm_id) FILTER(WHERE r.respiratory7d=s.respiratory7d) n_resp,
 count(r.hadm_id) FILTER(WHERE abs(r.known_platelets-s.known_platelets)<=20
 AND abs(r.age-s.age)<=15 AND ((r.known_creatinine<2 AND s.known_creatinine<2)
 OR (r.known_creatinine>=2 AND s.known_creatinine>=2))) n_simple
 FROM phase2_switches s LEFT JOIN phase2_risk_eligible r ON s.hadm_id=r.index_hadm
 GROUP BY s.hadm_id)
 SELECT count(*) switch_cases,count(*) FILTER(WHERE n>0) with_any_control,
 count(*) FILTER(WHERE n>=5) with_at_least_five_controls,median(n)::DOUBLE median_controls,
 count(*) FILTER(WHERE n_resp>0) with_same_respiratory_specimen_status,
 count(*) FILTER(WHERE n_simple>0) with_basic_clinical_calipers,
 median(n_simple)::DOUBLE median_controls_basic_calipers FROM n''')
out['risk_pool']=rows('''SELECT count(*) eligible_case_control_pairs,count(DISTINCT hadm_id) unique_control_admissions,
 count(DISTINCT subject_id) unique_control_patients FROM phase2_risk_eligible''')
# Follow-up ascertainment audit, not an endpoint count.
c.execute('''CREATE OR REPLACE TABLE phase2_low_decisions AS
 SELECT * FROM phase2_decisions WHERE known_platelets<100''')
out['followup_location']=rows('''WITH f AS (SELECT e.hadm_id,e.drug,e.initiation_type,
 count(p.charttime) n, count(p.charttime) FILTER(WHERE NOT EXISTS(SELECT 1 FROM phase2_icu i
 WHERE i.hadm_id=e.hadm_id AND p.charttime>=i.intime AND p.charttime<i.outtime)) outside_icu
 FROM phase2_low_decisions e LEFT JOIN platelet_measurements p ON p.hadm_id=e.hadm_id
 AND p.charttime>=e.t0 AND p.charttime<e.t0+INTERVAL 14 DAY AND p.charttime<e.dischtime GROUP BY ALL)
 SELECT drug,initiation_type,count(*) admissions,count(*) FILTER(WHERE outside_icu>0) with_non_icu_platelet_labs,
 count(*) FILTER(WHERE outside_icu=n AND n>0) all_followup_labs_non_icu FROM f GROUP BY ALL''')
out['transfusion_source_coverage']=rows('''WITH h AS (SELECT DISTINCT hadm_id FROM phase2_low_decisions),
 p AS (SELECT x.*,EXISTS(SELECT 1 FROM platelet_transfusions_icu t WHERE t.hadm_id=x.hadm_id
 AND t.amount>0 AND coalesce(t.statusdescription,'')<>'Rewritten'
 AND cast(t.starttime AS DATE)=x.chartdate) matching_icu_day,
 EXISTS(SELECT 1 FROM phase2_icu i WHERE i.hadm_id=x.hadm_id
 AND i.intime<cast(x.chartdate AS TIMESTAMP)+INTERVAL 1 DAY AND i.outtime>cast(x.chartdate AS TIMESTAMP)) overlaps_icu_day
 FROM phase2_platelet_procedures x SEMI JOIN h USING(hadm_id))
 SELECT count(*) procedure_records,count(DISTINCT hadm_id) admissions_with_platelet_procedure,
 count(*) FILTER(WHERE NOT matching_icu_day) procedure_records_without_same_day_icu_transfusion,
 count(*) FILTER(WHERE NOT overlaps_icu_day) procedure_records_on_non_icu_days FROM p''')
out['platelet_emar_terms']=rows('''SELECT medication,event_txt,count(*) records FROM phase2_platelet_emar
 GROUP BY ALL ORDER BY records DESC LIMIT 15''')
out['checks']={
 'known_baseline_has_no_future_specimen_or_result':c.execute('''SELECT count(*)=0 FROM phase2_decisions
 WHERE known_platelets_time>=t0 OR known_platelets_store>t0''').fetchone()[0],
 'risk_set_no_same_subject':c.execute('SELECT count(*)=0 FROM phase2_risk_eligible WHERE subject_id=index_subject').fetchone()[0],
 'one_control_admission_per_risk_set':c.execute('''SELECT count(*)=count(DISTINCT (index_hadm,hadm_id)) FROM phase2_risk_eligible''').fetchone()[0]}
assert all(out['checks'].values()),out['checks']
out['limitations']=['No clinical treatment effect or outcome counts computed.',
 'Current discharge sepsis/pneumonia/marrow labels cannot establish baseline onset; prior hospital records are incomplete history.',
 'Past vancomycin dose within 48h does not establish a clinician decision to continue; renal dosing/temporary holds require prescription-course adjudication.',
 'ICD procedure dates cannot prove exact timing or complete transfusion capture; medication terms require semantic review.',
 'Detailed severity, ventilation, renal replacement, contemporaneous antibiotics and treatment indication remain to be added before weighting.',
 'A patient may appear as a control before later switching; valid for a risk-set design, requires correct subsequent-treatment handling and clustered uncertainty.']
out=publish(out)
(R/'reports/DECISION_TIME_AUDIT_2026-09-11.json').write_text(json.dumps(out,indent=2))
(R/'manifests/phase2_clinical.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
print(json.dumps(out,indent=2),flush=True)
