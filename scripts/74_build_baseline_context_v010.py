"""Build baseline-only infection, care context and source readiness for the new candidate pool."""
from pathlib import Path
import json,os,datetime
import duckdb,pandas as pd,numpy as np
os.umask(0o077);R=Path('/root/projects/linezolid_platelet_recovery');c=duckdb.connect(str(R/'cache/audit_combined.duckdb'));c.execute('SET threads=4');c.execute("SET memory_limit='12GB'")
quote=lambda seq:','.join("'"+x.replace("'","''")+"'" for x in seq)
units=json.loads((R/'config/analysis_plan_v0.7.json').read_text())['phase7']['ICU']['careunits'];us=quote(units)
vent=json.loads((R/'config/source_addendum_v0.8.json').read_text())['ventilation_modes'];vm=quote(vent['ventilator_mode']);hm=quote(vent['ventilator_mode_hamilton'])
beta=json.loads((R/'config/analysis_plan_v0.9.json').read_text())['clinical']['beta_lactam'];bids=quote(beta['icu_itemids']);brx=beta['emar_regex']
otherrx='sulfamethoxazole|trimethoprim|bactrim|ganciclovir|valpro|divalpro|mycophen|azathiopr|methotrex|cyclophosphamide|cytarabine|daunorubicin|doxorubicin|etoposide|fludarabine|tacrolimus|sirolimus|eptifibatide|tirofiban'
c.execute(f'''CREATE OR REPLACE TEMP VIEW phase10_micro_class AS
 SELECT m.*,coalesce(t.routine_bacterial,false) AND NOT regexp_matches(upper(coalesce(m.spec_type_desc,'')),'SCREEN|SURVEILLANCE') routine,
 coalesce(t.colonization_screen,false) OR regexp_matches(upper(coalesce(m.spec_type_desc,'')||' '||coalesce(m.test_name,'')),'SCREEN|SURVEILLANCE') screening
 FROM phase10_microbiology m LEFT JOIN read_csv('{R/'config/microbiology_test_taxonomy_v0.10.tsv'}',delim='\t',header=true) t ON coalesce(m.test_name,'[MISSING]')=t.test_name
 WHERE m.linked_admissions=1''')
c.execute('''CREATE OR REPLACE TEMP VIEW phase10_detail_flags AS SELECT emar_id,bool_or(regexp_matches(upper(coalesce(route,'')),'^IV|INTRAVENOUS')) is_iv,bool_or(dose_given>0) positive,bool_or(dose_given=0) zero FROM phase10_med_detail GROUP BY emar_id''')
c.execute('''CREATE OR REPLACE TEMP TABLE phase10_available_labs AS SELECT l.* FROM target_labs l JOIN phase10_candidates e ON l.hadm_id=e.hadm_id WHERE l.charttime>=e.t0-INTERVAL 7 DAY AND l.charttime<e.t0 AND l.storetime<=e.t0''')
print('Available pre-index labs ready',flush=True)
query=f'''CREATE OR REPLACE TABLE phase10_baseline AS SELECT e.*,
 coalesce(loc.icu_n=1 AND NOT loc.non_icu,false) physical_verified,loc.icu_unit,
 coalesce(mi.diagnostic72,false) diagnostic_microbiology72h,coalesce(mi.routine72,false) routine_bacterial_culture72h,coalesce(mi.routine_dateonly72,false) routine_dateonly72h,
 coalesce(mi.blood72,false) routine_blood72h,coalesce(mi.urine72,false) routine_urine72h,coalesce(mi.resp72,false) routine_respiratory72h,coalesce(mi.nongut72,false) routine_nongut72h,
 coalesce(mi.screen72,false) screening72h,coalesce(mi.vre,false) known_vre7d,coalesce(mi.vre_diagnostic,false) known_vre_diagnostic7d,
 coalesce(mi.named,false) known_named_routine7d,coalesce(mi.entero,false) known_enterococcus7d,coalesce(mi.mrsa,false) known_mrsa_routine7d,coalesce(mi.blood_named,false) known_blood_organism7d,
 coalesce(md.ufh,false) OR coalesce(inp.ufh,false) ufh_active_available7d,coalesce(md.lmwh,false) lmwh_available7d,coalesce(md.otherdrug,false) marrowdrug_available7d,
 coalesce(md.beta,false) OR coalesce(inp.beta,false) beta_lactam_available24h,coalesce(md.any_before,false) any_emar_available_before,
 coalesce(inp.pressor,false) pressors_documented_prior6h,coalesce(pr.rrt,false) OR coalesce(ch.rrt,false) rrt_documented_prior24h,
 coalesce(pr.vent,false) OR coalesce(ch.vent,false) invasive_vent_documented_prior24h,
 cr.valuenum known_creatinine,coalesce(cr.valuenum-prior.min48>=0.3 OR cr.valuenum/prior.min7>=1.5,false) observed_creatinine_aki,
 inr.valuenum baseline_inr,bili.valuenum baseline_bilirubin,e.known_platelets-oldplt.val platelet_change,
 coalesce(hist.marrow,false) prior_recorded_marrow,date_diff('second',e.admittime,e.t0)/86400.0 days_since_admission,
 CASE WHEN try_cast(right(e.anchor_year_group,4) AS INTEGER)+year(e.t0)-try_cast(pt.anchor_year AS INTEGER)<=2013 THEN 'mapped_early'
 WHEN try_cast(substr(e.anchor_year_group,1,4) AS INTEGER)+year(e.t0)-try_cast(pt.anchor_year AS INTEGER)>=2014 THEN 'mapped_late' ELSE 'boundary_uncertain' END mapped_era,
 EXISTS(SELECT 1 FROM phase10_inputs i WHERE i.hadm_id=e.hadm_id AND i.starttime=e.t0 AND i.itemid=CASE WHEN e.drug='linezolid' THEN '225881' ELSE '225798' END AND i.amount>0 AND coalesce(i.statusdescription,'')<>'Rewritten') index_input_verified,
 EXISTS(SELECT 1 FROM emar_systemic_administrations i JOIN phase10_emar v USING(emar_id) WHERE i.hadm_id=e.hadm_id AND i.drug=e.drug AND i.charttime=e.t0) index_emar_verified,
 EXISTS(SELECT 1 FROM phase61_eligible_points o WHERE o.hadm_id=e.hadm_id) in_original_cohort
 FROM phase10_candidates e JOIN patients pt USING(subject_id)
 LEFT JOIN LATERAL(SELECT count(DISTINCT careunit) FILTER(WHERE careunit IN ({us})) icu_n,bool_or(careunit NOT IN ({us})) non_icu,min(careunit) FILTER(WHERE careunit IN ({us})) icu_unit
 FROM phase10_transfers t WHERE t.hadm_id=e.hadm_id AND t.intime<=e.t0 AND t.outtime>e.t0) loc ON true
 LEFT JOIN LATERAL(SELECT bool_or(NOT screening AND sample_available_by>=e.t0-INTERVAL 72 HOUR) diagnostic72,
 bool_or(routine AND sample_available_by>=e.t0-INTERVAL 72 HOUR) routine72,
 bool_or(routine AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND sample_date_only) routine_dateonly72,
 bool_or(routine AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND upper(spec_type_desc) LIKE '%BLOOD%') blood72,
 bool_or(routine AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND upper(spec_type_desc) LIKE '%URINE%') urine72,
 bool_or(routine AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND regexp_matches(upper(spec_type_desc),'SPUTUM|BRONCH|TRACHE|RESPIRATORY|LUNG')) resp72,
 bool_or(routine AND sample_available_by>=e.t0-INTERVAL 72 HOUR AND NOT regexp_matches(upper(spec_type_desc),'STOOL|FECAL|RECTAL')) nongut72,
 bool_or(screening AND sample_available_by>=e.t0-INTERVAL 72 HOUR) screen72,
 bool_or((upper(known_org_name) LIKE '%ENTEROCOCCUS%' AND upper(known_ab_name)='VANCOMYCIN' AND upper(known_interpretation)='R') OR regexp_matches(upper(known_org_name),'VANCOMYCIN.?RESISTANT.*ENTEROCOCC|^VRE$')) vre,
 bool_or(NOT screening AND ((upper(known_org_name) LIKE '%ENTEROCOCCUS%' AND upper(known_ab_name)='VANCOMYCIN' AND upper(known_interpretation)='R') OR regexp_matches(upper(known_org_name),'VANCOMYCIN.?RESISTANT.*ENTEROCOCC|^VRE$'))) vre_diagnostic,
 bool_or(routine AND trim(coalesce(known_org_name,'')) NOT IN ('','CANCELLED','POSITIVE','NEGATIVE','NO GROWTH')) named,
 bool_or(routine AND upper(known_org_name) LIKE '%ENTEROCOCCUS%') entero,
 bool_or(routine AND (regexp_matches(upper(known_org_name),'MRSA|METHICILLIN.?RESISTANT.*STAPH') OR (upper(known_org_name) LIKE '%STAPH% AUREUS%' AND upper(known_ab_name) IN ('OXACILLIN','METHICILLIN','CEFOXITIN') AND upper(known_interpretation)='R'))) mrsa,
 bool_or(routine AND upper(spec_type_desc) LIKE '%BLOOD%' AND trim(coalesce(known_org_name,'')) NOT IN ('','CANCELLED','POSITIVE','NEGATIVE','NO GROWTH')) blood_named
 FROM phase10_micro_class mm WHERE mm.hadm_id=e.hadm_id) mi ON true
 LEFT JOIN LATERAL(SELECT bool_or(m.storetime<=e.t0) any_before,
 bool_or(m.storetime<=e.t0 AND m.charttime>=e.t0-INTERVAL 7 DAY AND regexp_matches(lower(m.medication),'heparin') AND NOT regexp_matches(lower(m.medication),'flush|lock|hep-lock') AND m.event_txt IN ('Administered','Delayed Administered','Administered in Other Location','Started','Started in Other Location','Delayed Started','Restarted','Partial Administered')) ufh,
 bool_or(m.storetime<=e.t0 AND m.charttime>=e.t0-INTERVAL 7 DAY AND regexp_matches(lower(m.medication),'enoxaparin|dalteparin') AND m.event_txt IN ('Administered','Delayed Administered','Administered in Other Location','Started','Started in Other Location','Delayed Started','Restarted','Partial Administered')) lmwh,
 bool_or(m.storetime<=e.t0 AND m.charttime>=e.t0-INTERVAL 7 DAY AND regexp_matches(lower(m.medication),'{otherrx}') AND m.event_txt IN ('Administered','Delayed Administered','Administered in Other Location','Started','Started in Other Location','Delayed Started','Restarted','Partial Administered')) otherdrug,
 bool_or(m.storetime<=e.t0 AND m.charttime>=e.t0-INTERVAL 24 HOUR AND regexp_matches(lower(m.medication),'{brx}') AND df.is_iv AND NOT(coalesce(df.zero,false) AND NOT coalesce(df.positive,false)) AND m.event_txt IN ('Administered','Delayed Administered','Administered in Other Location','Started','Started in Other Location','Delayed Started','Restarted','Partial Administered')) beta
 FROM phase10_emar m LEFT JOIN phase10_detail_flags df USING(emar_id) WHERE m.hadm_id=e.hadm_id AND m.charttime<e.t0) md ON true
 LEFT JOIN LATERAL(SELECT bool_or(itemid='225152') ufh,bool_or(itemid IN ({bids}) AND starttime>=e.t0-INTERVAL 24 HOUR) beta,
 bool_or(itemid IN ('221906','221289','222315','221749','221662') AND coalesce(endtime,starttime)>=e.t0-INTERVAL 6 HOUR) pressor
 FROM phase10_inputs i WHERE i.hadm_id=e.hadm_id AND i.starttime<e.t0 AND i.storetime<=e.t0 AND i.amount>0 AND coalesce(statusdescription,'')<>'Rewritten') inp ON true
 LEFT JOIN LATERAL(SELECT bool_or(itemid IN ('225441','225802','225803','225809','225955','225805')) rrt,bool_or(itemid='225792') vent FROM phase10_procedures p WHERE p.hadm_id=e.hadm_id AND p.storetime<=e.t0 AND coalesce(statusdescription,'')<>'Rewritten') pr ON true
 LEFT JOIN LATERAL(SELECT bool_or((itemid='223849' AND value IN ({vm})) OR (itemid='229314' AND value IN ({hm})) OR (itemid='226732' AND value='Endotracheal tube')) vent,
 bool_or((itemid IN ('224144','224154','226499','225806') AND valuenum>0) OR (itemid='227290' AND upper(value) IN ('CVVH','CVVHD','CVVHDF','SCUF'))) rrt FROM phase10_charted_support s WHERE s.hadm_id=e.hadm_id) ch ON true
 LEFT JOIN LATERAL(SELECT valuenum,charttime FROM phase10_available_labs l WHERE l.hadm_id=e.hadm_id AND itemid IN ('50912','52546') AND valueuom='mg/dL' AND valuenum>0 AND charttime>=e.t0-INTERVAL 24 HOUR AND charttime<e.t0 AND storetime<=e.t0 ORDER BY charttime DESC,storetime DESC,valuenum LIMIT 1) cr ON true
 LEFT JOIN LATERAL(SELECT min(valuenum) FILTER(WHERE charttime>=e.t0-INTERVAL 48 HOUR) min48,min(valuenum) min7 FROM phase10_available_labs l WHERE l.hadm_id=e.hadm_id AND itemid IN ('50912','52546') AND valueuom='mg/dL' AND valuenum>0 AND charttime>=e.t0-INTERVAL 7 DAY AND charttime<cr.charttime AND storetime<=e.t0) prior ON true
 LEFT JOIN LATERAL(SELECT valuenum FROM phase10_available_labs l WHERE l.hadm_id=e.hadm_id AND itemid='51237' AND valuenum>0 AND charttime>=e.t0-INTERVAL 24 HOUR AND charttime<e.t0 AND storetime<=e.t0 ORDER BY charttime DESC,storetime DESC,valuenum LIMIT 1) inr ON true
 LEFT JOIN LATERAL(SELECT valuenum FROM phase10_available_labs l WHERE l.hadm_id=e.hadm_id AND itemid='50885' AND valueuom='mg/dL' AND valuenum>=0 AND charttime>=e.t0-INTERVAL 24 HOUR AND charttime<e.t0 AND storetime<=e.t0 ORDER BY charttime DESC,storetime DESC,valuenum LIMIT 1) bili ON true
 LEFT JOIN LATERAL(SELECT charttime,min(valuenum) val FROM phase10_available_labs l WHERE l.hadm_id=e.hadm_id AND itemid IN ('51265','53189') AND valueuom='K/uL' AND valuenum>=0 AND valuenum<5000 AND charttime>=e.t0-INTERVAL 72 HOUR AND charttime<e.t0-INTERVAL 24 HOUR AND storetime<=e.t0 GROUP BY charttime HAVING count(DISTINCT valuenum)=1 ORDER BY charttime DESC LIMIT 1) oldplt ON true
 LEFT JOIN LATERAL(SELECT bool_or(d.coded_marrow_disease) marrow FROM admissions a JOIN diagnosis_flags d USING(hadm_id) WHERE a.subject_id=e.subject_id AND try_cast(a.dischtime AS TIMESTAMP)<e.admittime) hist ON true'''
# Materialize each independent covariate block to avoid a giant multi-lateral join plan.
import re
pat=re.compile(r'LEFT JOIN LATERAL\((.*?)\) (\w+) ON true(?=\n LEFT JOIN LATERAL|$)',re.S)
parts=list(pat.finditer(query));assert len(parts)==12,len(parts)
for match in parts:
 body,alias=match.groups();extra='LEFT JOIN phase10_cov_cr cr USING(hadm_id)' if alias=='prior' else ''
 c.execute(f'CREATE OR REPLACE TEMP TABLE phase10_cov_{alias} AS SELECT e.hadm_id,{alias}.* FROM phase10_candidates e {extra} LEFT JOIN LATERAL({body}) {alias} ON true')
 assert c.execute(f'SELECT count(*)=3674 AND count(DISTINCT hadm_id)=3674 FROM phase10_cov_{alias}').fetchone()[0]
 print('Covariate block',alias,'ready',flush=True)
query=pat.sub(lambda m:f'LEFT JOIN phase10_cov_{m.group(2)} {m.group(2)} ON {m.group(2)}.hadm_id=e.hadm_id',query)
c.execute(query)
d=c.execute('SELECT * FROM phase10_baseline ORDER BY hadm_id').fetchdf();assert len(d)==3674 and d.hadm_id.nunique()==3674;assert (d.index_input_verified|d.index_emar_verified).all()
d['physical_no_known_vre']=d.physical_verified & ~d.known_vre7d
base=d.physical_no_known_vre
spec={'physical_no_known_vre':base,'diagnostic_microbiology72h':base&d.diagnostic_microbiology72h,'routine_bacterial_culture72h':base&d.routine_bacterial_culture72h,'routine_bacterial_culture72h_plus_IV_beta24h':base&d.routine_bacterial_culture72h&d.beta_lactam_available24h,'respiratory_routine_culture72h':base&d.routine_respiratory72h}
flow=[];context=[];unitrows=[]
for pop,mask in spec.items():
 xx=pd.crosstab(d.loc[mask,'icu_unit'],d.loc[mask,'drug']);shared=xx.index[(xx>0).all(axis=1)];d[pop+'_shared']=mask&d.icu_unit.isin(shared)
 for drug,g in d.groupby('drug'):
  kept=g[g[pop+'_shared']];flow.append(dict(population=pop,drug=drug,before_shared_units=int(mask.loc[g.index].sum()),n=len(kept),patients=int(kept.subject_id.nunique()),in_original=int(kept.in_original_cohort.sum()),without_discharge_sepsis_code=int((~kept.coded_sepsis).sum()),shared_units=list(shared)))
  for col in ['routine_dateonly72h','routine_blood72h','routine_urine72h','routine_respiratory72h','known_named_routine7d','known_enterococcus7d','known_mrsa_routine7d','known_blood_organism7d','any_emar_available_before','index_input_verified','index_emar_verified','beta_lactam_available24h','ufh_active_available7d','lmwh_available7d','marrowdrug_available7d','rrt_documented_prior24h','invasive_vent_documented_prior24h','pressors_documented_prior6h','observed_creatinine_aki']:
   context.append(dict(population=pop,drug=drug,variable=col,n=int(kept[col].sum()),denominator=len(kept)))
  for col in ['known_creatinine','baseline_inr','baseline_bilirubin','platelet_change']:
   context.append(dict(population=pop,drug=drug,variable=col+'_missing',n=int(kept[col].isna().sum()),denominator=len(kept)))
  for unit,u in kept.groupby('icu_unit'):unitrows.append(dict(population=pop,drug=drug,unit=unit,n=len(u)))
d.to_csv(R/'cache/phase10_baseline.csv',index=False)
# Source-equivalence check among original admissions, excluding intentionally revised microbiology definitions.
old=pd.read_csv(R/'cache/phase9_context.csv',dtype={'hadm_id':str});v=d[d.in_original_cohort].merge(old,on='hadm_id',suffixes=('_new','_old'),validate='1:1');comp={}
for key in ['ufh_active_available7d','lmwh_available7d','marrowdrug_available7d','beta_lactam_available24h','observed_creatinine_aki']:
 comp[key]=int((v[key+'_new']!=v[key+'_old']).sum())
for new,prev in [('rrt_documented_prior24h','combined_rrt24h'),('invasive_vent_documented_prior24h','combined_vent24h'),('pressors_documented_prior6h','pressor_available6h')]:comp[new]=int((v[new]!=v[prev]).sum())
report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),no_recovery_outcomes_loaded=True,flow=flow,context=context,unit_counts=unitrows,old_cohort_source_disagreements=comp,physical_failures=d.loc[~d.physical_verified].groupby('drug').size().to_dict(),known_vre=d[d.known_vre7d].groupby('drug').size().to_dict(),stool_only_routine=d[d.routine_bacterial_culture72h&~d.routine_nongut72h].groupby('drug').size().to_dict())
def hide(x):
 if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
 if isinstance(x,list):return [hide(v) for v in x]
 if isinstance(x,int) and not isinstance(x,bool) and 0<x<10:return '<10'
 return x
(R/'cache/internal_reports/BASELINE_TREATMENT_CONTEXT_v0.10_exact.json').write_text(json.dumps(report,indent=2));(R/'reports/BASELINE_TREATMENT_CONTEXT_v0.10.json').write_text(json.dumps(hide(report),indent=2))
print(json.dumps(hide({k:report[k] for k in ['flow','old_cohort_source_disagreements','physical_failures','known_vre','stool_only_routine']}),indent=2))
