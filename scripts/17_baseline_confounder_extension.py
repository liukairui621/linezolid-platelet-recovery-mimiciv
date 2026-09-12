"""First clinical step under v0.3: time-anchored renal, support, heparin and drug history.
No treatment effect or post-treatment endpoint analysis. Patient data stay on server.
"""
from pathlib import Path
import os,json,datetime,hashlib
import duckdb
os.umask(0o077);R=Path('/root/projects/linezolid_platelet_recovery');B=Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'));c.execute('SET threads=6');c.execute("SET memory_limit='12GB'")
def csv(name):return f"read_csv('{B/name}',all_varchar=true,header=true)"
def cached(name,q):
    p=R/'cache'/f'{name}.parquet'
    if not p.exists():
        temp=str(p)+'.partial';c.execute(f"COPY ({q}) TO '{temp}' (FORMAT PARQUET,COMPRESSION ZSTD)");Path(temp).replace(p)
    c.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{p}')")
    print(json.dumps({'asset':name,'ready':True}),flush=True)
def rows(q):
    cur=c.execute(q);cols=[x[0] for x in cur.description];return [dict(zip(cols,r)) for r in cur.fetchall()]
def hide(x):
    if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
    if isinstance(x,list):return [hide(v) for v in x]
    if isinstance(x,bool):return x
    if isinstance(x,int) and 0<x<10:return '<10'
    if isinstance(x,float):return round(x,3)
    return x
c.execute('''CREATE OR REPLACE TABLE phase3_baseline_points AS SELECT * FROM phase2_decisions
 WHERE known_platelets<100 AND initiation_type='no_prior_other_drug' AND decision_in_icu AND NOT vre_before''')
cohort='(SELECT DISTINCT hadm_id FROM phase3_baseline_points)'
medrx='heparin|enoxaparin|dalteparin|fondaparinux|argatroban|bivalirudin|sulfamethoxazole|trimethoprim|bactrim|ganciclovir|valpro|divalpro|mycophen|azathiopr|methotrex|cyclophosphamide|cytarabine|daunorubicin|doxorubicin|etoposide|fludarabine|tacrolimus|sirolimus|eptifibatide|tirofiban'
cached('phase3_med_emar',f'''SELECT m.hadm_id,m.emar_id,m.medication,m.event_txt,
 try_cast(m.charttime AS TIMESTAMP) charttime,try_cast(m.storetime AS TIMESTAMP) storetime
 FROM {csv('hosp/emar.csv.gz')} m SEMI JOIN {cohort} h USING(hadm_id)
 WHERE regexp_matches(lower(m.medication),'{medrx}')''')
cached('phase3_med_detail',f'''SELECT d.emar_id,d.route,try_cast(replace(d.dose_given,',','') AS DOUBLE) dose_given
 FROM {csv('hosp/emar_detail.csv.gz')} d SEMI JOIN phase3_med_emar e USING(emar_id)''')
cached('phase3_selected_rx',f'''SELECT m.hadm_id,m.drug,m.route,try_cast(m.starttime AS TIMESTAMP) starttime,
 try_cast(m.stoptime AS TIMESTAMP) stoptime FROM {csv('hosp/prescriptions.csv.gz')} m
 SEMI JOIN {cohort} h USING(hadm_id) WHERE regexp_matches(lower(m.drug),'{medrx}')''')
cached('phase3_support_inputs',f'''SELECT i.hadm_id,i.itemid,try_cast(i.starttime AS TIMESTAMP) starttime,
 try_cast(i.endtime AS TIMESTAMP) endtime,try_cast(i.storetime AS TIMESTAMP) storetime,
 try_cast(i.amount AS DOUBLE) amount,i.amountuom,i.statusdescription
 FROM {csv('icu/inputevents.csv.gz')} i SEMI JOIN {cohort} h USING(hadm_id)
 WHERE i.itemid IN ('225152','221906','221289','222315','221749','221662')''')
cached('phase3_support_procedures',f'''SELECT i.hadm_id,i.itemid,try_cast(i.starttime AS TIMESTAMP) starttime,
 try_cast(i.endtime AS TIMESTAMP) endtime,try_cast(i.storetime AS TIMESTAMP) storetime,i.statusdescription
 FROM {csv('icu/procedureevents.csv.gz')} i SEMI JOIN {cohort} h USING(hadm_id)
 WHERE i.itemid IN ('225441','225802','225803','225809','225955','225805','225792','225794','224272','229529','229530')''')
cached('phase3_coag_labs',f'''SELECT l.hadm_id,l.itemid,try_cast(l.charttime AS TIMESTAMP) charttime,
 try_cast(l.storetime AS TIMESTAMP) storetime,try_cast(l.valuenum AS DOUBLE) valuenum,l.valueuom
 FROM {csv('hosp/labevents.csv.gz')} l SEMI JOIN {cohort} h USING(hadm_id)
 WHERE l.itemid IN ('50915','51196','52551','51214','51623','52116','52117','51274','51275','51237','50885','50862')''')
c.execute('''CREATE OR REPLACE TABLE phase3_admin_drugs AS
 SELECT e.*,d.route,
 CASE WHEN regexp_matches(lower(e.medication),'enoxaparin|dalteparin') THEN 'lmwh'
 WHEN regexp_matches(lower(e.medication),'heparin') THEN 'ufh_including_flush'
 WHEN regexp_matches(lower(e.medication),'fondaparinux|argatroban|bivalirudin') THEN 'non_heparin_anticoagulant'
 ELSE 'selected_other_drugs' END drug_group,
 regexp_matches(lower(e.medication),'flush|lock|hep-lock') explicit_flush_label
 FROM phase3_med_emar e LEFT JOIN (SELECT emar_id,first(route) FILTER(WHERE route IS NOT NULL) route,
 bool_or(dose_given>0) positive_dose,bool_or(dose_given=0) zero_dose FROM phase3_med_detail GROUP BY emar_id) d USING(emar_id)
 WHERE e.event_txt IN ('Administered','Delayed Administered','Administered in Other Location','Started',
 'Started in Other Location','Delayed Started','Restarted','Partial Administered')
 AND NOT (coalesce(d.zero_dose,false) AND NOT coalesce(d.positive_dose,false))''')
c.execute('''CREATE OR REPLACE TABLE phase3_baseline_covariates AS
 SELECT e.*,renal.n_creat_7d,renal.min_creat_48h,renal.min_creat_7d,
 e.known_creatinine-renal.min_creat_48h creatinine_rise_vs_prior48h,
 e.known_creatinine/nullif(renal.min_creat_7d,0) creatinine_ratio_vs_prior7d,
 coalesce(proc.rrt_before24h,false) rrt_documented_prior24h,
 coalesce(proc.rrt_stored_before24h,false) rrt_update_available_prior24h,
 coalesce(proc.crrt_before24h,false) crrt_documented_prior24h,
 coalesce(proc.ihd_before24h,false) ihd_documented_prior24h,
 coalesce(proc.invasive_vent24h,false) invasive_vent_documented_prior24h,
 coalesce(vaso.pressors6h,false) pressors_documented_prior6h,
 coalesce(meds.ufh7d,false) OR coalesce(vaso.heparin7d,false) ufh_documented_prior7d,
 coalesce(meds.lmwh7d,false) lmwh_documented_prior7d,
 coalesce(meds.flush7d,false) explicit_heparin_flush_prior7d,
 coalesce(meds.ufh_available7d,false) OR coalesce(vaso.heparin_available7d,false) ufh_update_available_prior7d,
 coalesce(meds.lmwh_available7d,false) lmwh_update_available_prior7d,
 coalesce(meds.non_heparin7d,false) non_heparin_anticoagulant_prior7d,
 coalesce(meds.other7d,false) selected_other_drug_prior7d,
 coalesce(rx.any_heparin_rx7d,false) prescribed_heparin_prior7d,
 labs.fibrinogen_known24h,labs.ddimer_known24h,labs.inr_known24h,labs.bilirubin_known24h,
 (SELECT count(*) FROM emar_systemic_administrations a WHERE a.hadm_id=e.hadm_id)>0 has_target_emar_source
 FROM phase3_baseline_points e
 LEFT JOIN LATERAL (SELECT count(DISTINCT l.charttime) n_creat_7d,
 min(l.valuenum) FILTER(WHERE l.charttime>=e.t0-INTERVAL 48 HOUR) min_creat_48h,min(l.valuenum) min_creat_7d
 FROM target_labs l WHERE l.hadm_id=e.hadm_id AND l.itemid IN ('50912','52546') AND l.valuenum>0
 AND l.charttime<e.t0 AND l.charttime>=e.t0-INTERVAL 7 DAY AND l.storetime<=e.t0) renal ON true
 LEFT JOIN LATERAL (SELECT
 bool_or(p.itemid IN ('225441','225802','225803','225809','225955','225805')) rrt_before24h,
 bool_or(p.itemid IN ('225441','225802','225803','225809','225955','225805') AND p.storetime<=e.t0) rrt_stored_before24h,
 bool_or(p.itemid IN ('225802','225803','225809','225955')) crrt_before24h,
 bool_or(p.itemid='225441') ihd_before24h,bool_or(p.itemid='225792') invasive_vent24h
 FROM phase3_support_procedures p WHERE p.hadm_id=e.hadm_id AND p.starttime<e.t0
 AND coalesce(p.endtime,p.starttime)>=e.t0-INTERVAL 24 HOUR AND coalesce(p.statusdescription,'')<>'Rewritten') proc ON true
 LEFT JOIN LATERAL (SELECT
 bool_or(v.itemid<>'225152' AND coalesce(v.endtime,v.starttime)>=e.t0-INTERVAL 6 HOUR) pressors6h,
 bool_or(v.itemid='225152') heparin7d,bool_or(v.itemid='225152' AND v.storetime<=e.t0) heparin_available7d
 FROM phase3_support_inputs v WHERE v.hadm_id=e.hadm_id AND v.starttime<e.t0
 AND coalesce(v.endtime,v.starttime)>=e.t0-INTERVAL 7 DAY AND v.amount>0 AND coalesce(v.statusdescription,'')<>'Rewritten') vaso ON true
 LEFT JOIN LATERAL (SELECT bool_or(m.drug_group='ufh_including_flush') ufh7d,bool_or(m.drug_group='lmwh') lmwh7d,
 bool_or(m.explicit_flush_label) flush7d,
 bool_or(m.drug_group='ufh_including_flush' AND m.storetime<=e.t0) ufh_available7d,
 bool_or(m.drug_group='lmwh' AND m.storetime<=e.t0) lmwh_available7d,
 bool_or(m.drug_group='non_heparin_anticoagulant') non_heparin7d,
 bool_or(m.drug_group='selected_other_drugs') other7d
 FROM phase3_admin_drugs m WHERE m.hadm_id=e.hadm_id AND m.charttime<e.t0 AND m.charttime>=e.t0-INTERVAL 7 DAY) meds ON true
 LEFT JOIN LATERAL (SELECT bool_or(regexp_matches(lower(r.drug),'heparin|enoxaparin|dalteparin')) any_heparin_rx7d
 FROM phase3_selected_rx r WHERE r.hadm_id=e.hadm_id AND r.starttime<e.t0 AND r.starttime>=e.t0-INTERVAL 7 DAY) rx ON true
 LEFT JOIN LATERAL (SELECT
 count(*) FILTER(WHERE l.itemid IN ('51214','51623','52116','52117'))>0 fibrinogen_known24h,
 count(*) FILTER(WHERE l.itemid IN ('50915','51196','52551'))>0 ddimer_known24h,
 count(*) FILTER(WHERE l.itemid='51237')>0 inr_known24h,
 count(*) FILTER(WHERE l.itemid='50885')>0 bilirubin_known24h
 FROM phase3_coag_labs l WHERE l.hadm_id=e.hadm_id AND l.charttime<e.t0
 AND l.charttime>=e.t0-INTERVAL 24 HOUR AND l.storetime<=e.t0 AND l.valuenum IS NOT NULL) labs ON true
''')
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'cohort':'Initial systemic drug during admission, at decision already in ICU, known platelets<100, no prior available VRE. Current sepsis code still only a screen.',
 'coverage':rows('''SELECT drug,count(*) admissions,count(DISTINCT subject_id) patients,
 count(*) FILTER(WHERE known_creatinine IS NOT NULL) known_creatinine,
 median(known_creatinine) median_creatinine,quantile_cont(known_creatinine,.75) q75_creatinine,
 count(*) FILTER(WHERE n_creat_7d>=2) at_least_two_prior_creatinines,
 count(*) FILTER(WHERE creatinine_rise_vs_prior48h>=.3 OR creatinine_ratio_vs_prior7d>=1.5) observed_creatinine_increase_criterion,
 count(*) FILTER(WHERE rrt_documented_prior24h) prior_rrt,
 count(*) FILTER(WHERE rrt_update_available_prior24h) prior_rrt_update_available,
 count(*) FILTER(WHERE crrt_documented_prior24h) prior_crrt,
 count(*) FILTER(WHERE ihd_documented_prior24h) prior_ihd,
 count(*) FILTER(WHERE invasive_vent_documented_prior24h) prior_invasive_vent,
 count(*) FILTER(WHERE pressors_documented_prior6h) prior_pressors,
 count(*) FILTER(WHERE ufh_documented_prior7d) prior_ufh_any,
 count(*) FILTER(WHERE lmwh_documented_prior7d) prior_lmwh,
 count(*) FILTER(WHERE explicit_heparin_flush_prior7d) explicit_flush,
 count(*) FILTER(WHERE ufh_update_available_prior7d) prior_ufh_update_available,
 count(*) FILTER(WHERE prescribed_heparin_prior7d) prior_heparin_prescription,
 count(*) FILTER(WHERE selected_other_drug_prior7d) selected_other_drugs,
 count(*) FILTER(WHERE fibrinogen_known24h) known_fibrinogen,
 count(*) FILTER(WHERE ddimer_known24h) known_ddimer,
 count(*) FILTER(WHERE inr_known24h) known_inr,
 count(*) FILTER(WHERE bilirubin_known24h) known_bilirubin,
 count(*) FILTER(WHERE has_target_emar_source) target_emar_present
 FROM phase3_baseline_covariates GROUP BY drug'''),
 'limits':['No full KDIGO AKI stage, eGFR diagnosis, SOFA score, HIT diagnosis or validated 4Ts score constructed.',
 'Renal minima are from observed admission labs, not premorbid renal baseline; chronic dialysis and AKI need separation.',
 'Procedure/input start and end times are retrospective clinical timing; final storetime can follow the decision despite clinicians knowing the ongoing treatment. Clinical status and last-update availability are therefore separate.',
 'No pretreatment drug record is not proven nonexposure. eMAR era coverage, flush labeling, doses/routes and prescription-only records still require reconciliation.',
 'HIT-specific PF4/SRA tests were not identified by the inspected lab dictionary. Generic heparin level or antiplatelet antibody is not a HIT test.',
 'Reported coagulation results are availability counts only; units and validated DIC definitions remain to be assessed.',
 'No propensity model, overlap weights, endpoint rates or clinical effect estimated in this step.']}
report=hide(report);(R/'reports/BASELINE_CONFOUNDER_EXTENSION_v0.3.json').write_text(json.dumps(report,indent=2))
(R/'manifests/phase3_confounder_script.json').write_text(json.dumps({'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
print(json.dumps(report,indent=2),flush=True)
