"""Outcome-free v0.8 baseline extraction. Individual data never leave server cache."""
from pathlib import Path
import json,datetime,hashlib,os
import duckdb
os.umask(0o077);R=Path('/root/projects/linezolid_platelet_recovery')
p=R/'config/analysis_plan_v0.8.json';plan=json.loads(p.read_text())
assert hashlib.sha256(p.read_bytes()).hexdigest()=='7ec85c9067c25bc65b2f643f1fca0dd82d0cf01151d74087571d6eeef3ec3d83'
units=json.loads((R/'config/analysis_plan_v0.7.json').read_text())['phase7']['ICU']['careunits']
unit_sql=','.join("'"+x.replace("'","''")+"'" for x in units)
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'));c.execute('SET threads=4')
c.execute(f'''CREATE OR REPLACE TABLE phase8_context AS
SELECT e.subject_id,e.hadm_id,e.drug,
 coalesce(loc.physical_icu,false) physical_icu,coalesce(loc.conflicting_non_icu,false) conflicting_non_icu,
 loc.icu_units,
 coalesce(m.nonscreen_culture72h,false) nonscreen_culture72h,
 coalesce(m.nonscreen_respiratory72h,false) nonscreen_respiratory72h,
 coalesce(m.nonscreen_blood72h,false) nonscreen_blood72h,
 coalesce(m.nonscreen_urine72h,false) nonscreen_urine72h,
 coalesce(m.screen72h,false) screen72h,
 coalesce(m.date_only_culture72h,false) date_only_culture72h,
 coalesce(m.known_named_nonscreen7d,false) known_named_nonscreen7d,
 coalesce(m.known_mrsa7d,false) known_mrsa7d,
 coalesce(m.known_enterococcus7d,false) known_enterococcus7d,
 coalesce(m.known_staph_aureus7d,false) known_staph_aureus7d,
 coalesce(proc.rrt_available24h,false) rrt_available24h,
 coalesce(proc.vent_available24h,false) vent_available24h,
 coalesce(vaso.pressor_available6h,false) pressor_available6h,
 cr.valuenum current_creatinine,prior.min48h prior_creatinine_min48h,prior.min7d prior_creatinine_min7d,
 coalesce(cr.valuenum-prior.min48h>=0.3 OR cr.valuenum/prior.min7d>=1.5,false) observed_creatinine_aki,
 prior.min7d IS NOT NULL AND cr.valuenum IS NOT NULL creatinine_change_evaluable,
 b.rrt_documented_prior24h AND NOT coalesce(proc.rrt_available24h,false) rrt_late_record,
 b.invasive_vent_documented_prior24h AND NOT coalesce(proc.vent_available24h,false) vent_late_record,
 b.pressors_documented_prior6h AND NOT coalesce(vaso.pressor_available6h,false) pressor_late_record
FROM phase61_eligible_points e JOIN phase4_baseline_only b USING(hadm_id)
LEFT JOIN LATERAL (SELECT bool_or(t.careunit IN ({unit_sql})) physical_icu,
 bool_or(t.careunit NOT IN ({unit_sql})) conflicting_non_icu,
 string_agg(DISTINCT t.careunit, ';' ORDER BY t.careunit) FILTER(WHERE t.careunit IN ({unit_sql})) icu_units
 FROM phase2_transfers t WHERE t.hadm_id=e.hadm_id AND t.intime<=e.t0 AND t.outtime>e.t0) loc ON true
LEFT JOIN LATERAL (SELECT
 bool_or(NOT screening AND sample_available_by>=e.t0-INTERVAL 72 HOUR) nonscreen_culture72h,
 bool_or(NOT screening AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND regexp_matches(upper(spec_type_desc),'SPUTUM|BRONCH|TRACHE|RESPIRATORY|LUNG')) nonscreen_respiratory72h,
 bool_or(NOT screening AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND regexp_matches(upper(spec_type_desc),'BLOOD')) nonscreen_blood72h,
 bool_or(NOT screening AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND regexp_matches(upper(spec_type_desc),'URINE')) nonscreen_urine72h,
 bool_or(screening AND sample_available_by>=e.t0-INTERVAL 72 HOUR) screen72h,
 bool_or(NOT screening AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND sample_date_only) date_only_culture72h,
 bool_or(NOT screening AND result_available_by<=e.t0 AND trim(coalesce(org_name,'')) NOT IN ('','CANCELLED','POSITIVE','NEGATIVE','NO GROWTH')) known_named_nonscreen7d,
 bool_or(NOT screening AND result_available_by<=e.t0 AND
 (regexp_matches(upper(org_name),'MRSA|METHICILLIN.?RESISTANT.*STAPH') OR
 (upper(org_name) LIKE '%STAPH% AUREUS%' AND upper(ab_name) IN ('OXACILLIN','METHICILLIN','CEFOXITIN') AND upper(interpretation)='R'))) known_mrsa7d,
 bool_or(NOT screening AND result_available_by<=e.t0 AND upper(org_name) LIKE '%ENTEROCOCCUS%') known_enterococcus7d,
 bool_or(NOT screening AND result_available_by<=e.t0 AND upper(org_name) LIKE '%STAPH% AUREUS%') known_staph_aureus7d
 FROM (SELECT *,regexp_matches(upper(coalesce(spec_type_desc,'')||' '||coalesce(test_name,'')),'SCREEN|SURVEILLANCE') screening FROM phase2_microbiology) mm
 WHERE mm.hadm_id=e.hadm_id AND sample_available_by>=e.t0-INTERVAL 7 DAY AND sample_available_by<=e.t0) m ON true
LEFT JOIN LATERAL (SELECT
 bool_or(itemid IN ('225441','225802','225803','225809','225955','225805')) rrt_available24h,
 bool_or(itemid='225792') vent_available24h
 FROM phase3_support_procedures s WHERE s.hadm_id=e.hadm_id AND s.starttime<e.t0 AND s.storetime<=e.t0
 AND coalesce(s.endtime,s.starttime)>=e.t0-INTERVAL 24 HOUR AND coalesce(statusdescription,'')<>'Rewritten') proc ON true
LEFT JOIN LATERAL (SELECT bool_or(itemid<>'225152') pressor_available6h
 FROM phase3_support_inputs s WHERE s.hadm_id=e.hadm_id AND s.starttime<e.t0 AND s.storetime<=e.t0
 AND coalesce(s.endtime,s.starttime)>=e.t0-INTERVAL 6 HOUR AND s.amount>0 AND coalesce(statusdescription,'')<>'Rewritten') vaso ON true
LEFT JOIN LATERAL (SELECT valuenum,charttime FROM target_labs l WHERE l.hadm_id=e.hadm_id
 AND itemid IN ('50912','52546') AND valuenum>0 AND valueuom='mg/dL'
 AND charttime<e.t0 AND charttime>=e.t0-INTERVAL 24 HOUR AND storetime<=e.t0
 ORDER BY charttime DESC,storetime DESC,valuenum LIMIT 1) cr ON true
LEFT JOIN LATERAL (SELECT min(valuenum) FILTER(WHERE charttime>=e.t0-INTERVAL 48 HOUR) min48h,min(valuenum) min7d
 FROM target_labs l WHERE l.hadm_id=e.hadm_id AND itemid IN ('50912','52546') AND valuenum>0 AND valueuom='mg/dL'
 AND charttime<cr.charttime AND charttime>=e.t0-INTERVAL 7 DAY AND storetime<=e.t0) prior ON true
''')
d=c.execute('SELECT * FROM phase8_context ORDER BY hadm_id').fetchdf();assert len(d)==1569 and d.hadm_id.nunique()==1569
# Ambiguous concurrent ICU/non-ICU locations are not declared verified physical ICU.
d['physical_icu']=d.physical_icu & ~d.conflicting_non_icu
d['availability_corrected_all']=True
d['physical_icu_recent_culture']=d.physical_icu & d.nonscreen_culture72h
d['physical_icu_respiratory_culture']=d.physical_icu & d.nonscreen_respiratory72h
d['physical_icu_recent_culture_aki']=d.physical_icu_recent_culture & d.observed_creatinine_aki
d.to_csv(R/'cache/phase8_context.csv',index=False)
counts=[]
for drug,z in d.groupby('drug'):
 for col in d.columns:
  if str(d[col].dtype)=='bool':counts.append(dict(drug=drug,definition=col,n=int(z[col].sum()),denominator=len(z)))
units=d.groupby(['drug','icu_units'],dropna=False).size().reset_index(name='n').fillna('not_explicit_icu').to_dict('records')
report=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),no_outcomes_loaded=True,plan_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),counts=counts,unit_counts=units,
 limitations=plan['clinical']['sepsis_limit'],indication_conclusion='Source-level proxies only. Notes or adjudicated baseline indication unavailable in current raw modules; common anti-MRSA indication not established by this audit.')
def hide(x):
 if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
 if isinstance(x,list):return [hide(v) for v in x]
 if isinstance(x,int) and not isinstance(x,bool) and 0<x<10:return '<10'
 return x
(R/'cache/internal_reports/CLINICAL_CONTEXT_v0.8_exact.json').write_text(json.dumps(report,indent=2))
(R/'reports/CLINICAL_CONTEXT_v0.8.json').write_text(json.dumps(hide(report),indent=2))
print(json.dumps(hide(report),indent=2))
