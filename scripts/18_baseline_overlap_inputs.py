"""v0.4 baseline-only snapshot and medication/indication provenance audit.
Patient rows and unsuppressed small aggregates remain on the controlled server.
No endpoint, future lab frequency, survival, or treatment-effect output is loaded.
"""
from pathlib import Path
import os, json, datetime, hashlib
import duckdb

os.umask(0o077)
R = Path('/root/projects/linezolid_platelet_recovery')
B = Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
c = duckdb.connect(str(R/'cache/audit_combined.duckdb'))
c.execute("SET threads=6"); c.execute("SET memory_limit='12GB'")
(R/'cache/internal_reports').mkdir(exist_ok=True)
def rows(q):
    cur=c.execute(q)
    return [dict(zip([x[0] for x in cur.description],r)) for r in cur.fetchall()]
def suppressed(x):
    if isinstance(x,dict): return {k:suppressed(v) for k,v in x.items()}
    if isinstance(x,list): return [suppressed(v) for v in x]
    if isinstance(x,bool): return x
    if isinstance(x,int) and 0<x<10: return '<10'
    if isinstance(x,float): return round(x,4)
    return x

# Any-drug coverage is distinct from the presence of the two target antibiotics.
# Whole-admission coverage below is an audit descriptor, NEVER a baseline PS term.
emar_cache=R/'cache/phase4_emar_coverage.parquet'
if not emar_cache.exists():
    c.execute(f"""COPY (SELECT e.hadm_id,
       count(*)>0 any_emar_during_admission,
       count(*) FILTER(WHERE try_cast(m.charttime AS TIMESTAMP)<e.t0)>0 any_emar_before,
       count(*) FILTER(WHERE try_cast(m.charttime AS TIMESTAMP)<e.t0
          AND try_cast(m.storetime AS TIMESTAMP)<=e.t0)>0 any_emar_available_before,
       count(*) FILTER(WHERE try_cast(m.charttime AS TIMESTAMP)<e.t0
          AND try_cast(m.charttime AS TIMESTAMP)>=e.t0-INTERVAL 7 DAY)>0 any_emar_prior7d,
       count(*) FILTER(WHERE try_cast(m.charttime AS TIMESTAMP)<e.t0
          AND try_cast(m.charttime AS TIMESTAMP)>=e.t0-INTERVAL 7 DAY
          AND try_cast(m.storetime AS TIMESTAMP)<=e.t0)>0 any_emar_available_prior7d
       FROM read_csv('{B/'hosp/emar.csv.gz'}',all_varchar=true) m
       JOIN phase3_baseline_points e USING(hadm_id) GROUP BY e.hadm_id
       ) TO '{emar_cache}' (FORMAT PARQUET,COMPRESSION ZSTD)""")
c.execute(f"CREATE OR REPLACE VIEW phase4_emar_coverage AS SELECT * FROM read_parquet('{emar_cache}')")
print('any-drug eMAR coverage scan ready',flush=True)

c.execute('''CREATE OR REPLACE TABLE phase4_baseline_only AS
 SELECT e.subject_id,e.hadm_id,e.drug,e.age,e.gender,e.anchor_year_group,
 year(e.t0)-try_cast(p.anchor_year AS INTEGER) decision_year_minus_anchor_year,
 date_diff('second',e.admittime,e.t0)/86400.0 days_since_admission,
 e.known_platelets,e.known_creatinine,
 oldplt.valuenum earlier_platelets_unique,
 e.known_platelets-oldplt.valuenum platelet_change,
 date_diff('second',oldplt.charttime,e.known_platelets_time)/86400.0 platelet_change_interval_days,
 e.rrt_documented_prior24h,e.crrt_documented_prior24h,e.ihd_documented_prior24h,
 e.rrt_update_available_prior24h,
 e.invasive_vent_documented_prior24h,e.pressors_documented_prior6h,
 e.n_creat_7d,e.creatinine_rise_vs_prior48h,e.creatinine_ratio_vs_prior7d,
 e.ufh_documented_prior7d,e.lmwh_documented_prior7d,e.explicit_heparin_flush_prior7d,
 e.prescribed_heparin_prior7d,e.ufh_update_available_prior7d,e.lmwh_update_available_prior7d,
 e.selected_other_drug_prior7d,e.prior_recorded_marrow,
 e.has_target_emar_source,
 coalesce(em.any_emar_during_admission,false) any_emar_during_admission,
 coalesce(em.any_emar_before,false) any_emar_before,
 coalesce(em.any_emar_available_before,false) any_emar_available_before,
 coalesce(em.any_emar_prior7d,false) any_emar_prior7d,
 coalesce(em.any_emar_available_prior7d,false) any_emar_available_prior7d,
 e.culture72h,e.respiratory7d,
 coalesce(mi.respiratory72h,false) respiratory72h,
 coalesce(mi.blood72h,false) blood72h,
 coalesce(mi.known_gpc_nonscreen,false) known_gpc_nonscreen,
 coalesce(mi.known_mrsa_nonscreen,false) known_mrsa_nonscreen,
 coalesce(mi.known_enterococcus_nonscreen,false) known_enterococcus_nonscreen,
 coalesce(mi.known_any_organism_nonscreen,false) known_any_organism_nonscreen,
 coalesce(ad.emar_at_decision,false) target_emar_at_decision,
 coalesce(ad.input_at_decision,false) target_input_at_decision,
 inr.valuenum baseline_inr,inr.valueuom baseline_inr_unit,
 bili.valuenum baseline_bilirubin,bili.valueuom baseline_bilirubin_unit
 FROM phase3_baseline_covariates e LEFT JOIN patients p USING(subject_id)
 LEFT JOIN phase4_emar_coverage em USING(hadm_id)
 LEFT JOIN LATERAL (SELECT charttime,min(valuenum) valuenum
    FROM platelet_measurements x WHERE x.hadm_id=e.hadm_id
    AND x.charttime>=e.t0-INTERVAL 72 HOUR AND x.charttime<e.t0-INTERVAL 24 HOUR
    AND x.storetime<=e.t0 GROUP BY charttime HAVING count(DISTINCT valuenum)=1
    ORDER BY charttime DESC LIMIT 1) oldplt ON true
 LEFT JOIN LATERAL (SELECT
    bool_or(regexp_matches(upper(m.spec_type_desc),'SPUTUM|BRONCH|TRACHE|RESPIRATORY|LUNG')
       AND m.sample_available_by>=e.t0-INTERVAL 72 HOUR) respiratory72h,
    bool_or(upper(m.spec_type_desc) LIKE '%BLOOD%' AND m.sample_available_by>=e.t0-INTERVAL 72 HOUR) blood72h,
    bool_or(m.result_available_by<=e.t0 AND regexp_matches(upper(m.org_name),'STAPH|STREP|ENTEROCOCCUS')
       AND NOT regexp_matches(upper(coalesce(m.test_name,'')||' '||coalesce(m.spec_type_desc,'')),'SCREEN|SURVEILLANCE')) known_gpc_nonscreen,
    bool_or(m.result_available_by<=e.t0 AND (upper(m.org_name) LIKE '%STAPH AUREUS%' OR upper(m.org_name) LIKE '%STAPHYLOCOCCUS AUREUS%')
       AND upper(m.ab_name) IN ('OXACILLIN','METHICILLIN','CEFOXITIN') AND upper(m.interpretation)='R'
       AND NOT regexp_matches(upper(coalesce(m.test_name,'')||' '||coalesce(m.spec_type_desc,'')),'SCREEN|SURVEILLANCE')) known_mrsa_nonscreen,
    bool_or(m.result_available_by<=e.t0 AND upper(m.org_name) LIKE '%ENTEROCOCCUS%'
       AND NOT regexp_matches(upper(coalesce(m.test_name,'')||' '||coalesce(m.spec_type_desc,'')),'SCREEN|SURVEILLANCE')) known_enterococcus_nonscreen,
    bool_or(m.result_available_by<=e.t0 AND m.org_name IS NOT NULL
       AND NOT regexp_matches(upper(coalesce(m.test_name,'')||' '||coalesce(m.spec_type_desc,'')),'SCREEN|SURVEILLANCE')) known_any_organism_nonscreen
    FROM phase2_microbiology m WHERE m.hadm_id=e.hadm_id
       AND m.sample_available_by>=e.t0-INTERVAL 7 DAY AND m.sample_available_by<=e.t0) mi ON true
 LEFT JOIN LATERAL (SELECT bool_or(a.data_source='emar') emar_at_decision,
     bool_or(a.data_source<>'emar') input_at_decision
     FROM systemic_administrations a WHERE a.hadm_id=e.hadm_id AND a.drug=e.drug AND a.charttime=e.t0) ad ON true
 LEFT JOIN LATERAL (SELECT charttime,min(valuenum) valuenum,first(valueuom) valueuom
    FROM phase3_coag_labs x WHERE x.hadm_id=e.hadm_id AND x.itemid='51237'
    AND x.charttime<e.t0 AND x.charttime>=e.t0-INTERVAL 24 HOUR AND x.storetime<=e.t0 AND x.valuenum>0
    GROUP BY charttime HAVING count(DISTINCT valuenum)=1 ORDER BY charttime DESC LIMIT 1) inr ON true
 LEFT JOIN LATERAL (SELECT charttime,min(valuenum) valuenum,first(valueuom) valueuom
    FROM phase3_coag_labs x WHERE x.hadm_id=e.hadm_id AND x.itemid='50885'
    AND x.charttime<e.t0 AND x.charttime>=e.t0-INTERVAL 24 HOUR AND x.storetime<=e.t0 AND x.valuenum>=0
    GROUP BY charttime HAVING count(DISTINCT valuenum)=1 ORDER BY charttime DESC LIMIT 1) bili ON true
''')

cols=[x[0] for x in c.execute('describe phase4_baseline_only').fetchall()]
for bad in ['deathtime','dischtime','n_platelet_14d','observed_treatment_days','hospital_observation_days']:
    assert bad not in cols
assert c.execute('select count(*)=count(distinct hadm_id) from phase4_baseline_only').fetchone()[0]
assert c.execute('select count(*) from phase4_baseline_only').fetchone()[0] == c.execute('select count(*) from phase3_baseline_points').fetchone()[0]
c.execute(f"COPY phase4_baseline_only TO '{R/'cache/phase4_baseline_only.csv'}' (HEADER,DELIMITER ',')")
coverage_cols=['has_target_emar_source','any_emar_during_admission','any_emar_before','any_emar_available_before',
  'any_emar_prior7d','any_emar_available_prior7d','target_emar_at_decision','target_input_at_decision',
  'ufh_documented_prior7d','lmwh_documented_prior7d','prescribed_heparin_prior7d',
  'rrt_documented_prior24h','culture72h','respiratory72h','blood72h',
  'known_gpc_nonscreen','known_mrsa_nonscreen','known_enterococcus_nonscreen']
qcols=','.join(f'count(*) FILTER(WHERE {x}) {x}' for x in coverage_cols)
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'scope':'baseline-only audit; whole-admission eMAR flags are retrospective coverage descriptors and excluded from PS',
 'coverage_by_drug':rows(f'SELECT drug,count(*) n,{qcols} FROM phase4_baseline_only GROUP BY drug'),
 'emar_by_anchor':rows('''SELECT anchor_year_group,drug,count(*) n,
    count(*) FILTER(WHERE any_emar_available_before) baseline_any_emar,
    count(*) FILTER(WHERE has_target_emar_source) target_emar_ever
    FROM phase4_baseline_only GROUP BY ALL ORDER BY anchor_year_group,drug'''),
 'heparin_by_baseline_coverage':rows('''SELECT drug,any_emar_available_before,count(*) n,
    count(*) FILTER(WHERE ufh_documented_prior7d OR lmwh_documented_prior7d) recorded_ufh_or_lmwh,
    count(*) FILTER(WHERE prescribed_heparin_prior7d) prescription
    FROM phase4_baseline_only GROUP BY ALL'''),
 'indication_screen':rows('''SELECT drug,count(*) n,
    count(*) FILTER(WHERE culture72h) recent_culture,
    count(*) FILTER(WHERE respiratory72h) recent_respiratory_sample,
    median(days_since_admission) median_decision_day,
    count(*) FILTER(WHERE prior_recorded_marrow) prior_marrow,
    count(*) FILTER(WHERE platelet_change IS NOT NULL) platelet_trend_available
    FROM phase4_baseline_only GROUP BY drug'''),
 'units':rows('SELECT baseline_inr_unit,baseline_bilirubin_unit,count(*) n FROM phase4_baseline_only GROUP BY ALL'),
 'limitations':['Any eMAR record is not proof of complete administration capture. Target-drug eMAR presence is not a hospital eMAR interface measure.',
   'Anchor-year groups are deidentification anchor ranges, not exact current-admission calendar years; admission-to-anchor year displacement also recorded.',
   'Recorded respiratory/blood specimens and nonscreen organisms do not establish pneumonia, bloodstream infection or clinical treatment indication.',
   'INR, bilirubin and platelet trend remain separate measured variables, not a newly validated score. No DIC score will be constructed.',
   'Public suppression is an output rule; exact counts and patient records remain available for server-side modeling.'],
 'checks':{'one_row_per_admission':True,'no_outcome_columns_exported':True,'cohort_size_preserved':True}}
(R/'cache/internal_reports/BASELINE_PROVENANCE_v0.4_exact.json').write_text(json.dumps(report,indent=2))
(R/'reports/BASELINE_PROVENANCE_v0.4.json').write_text(json.dumps(suppressed(report),indent=2))
manifest={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'snapshot_sha256':hashlib.sha256((R/'cache/phase4_baseline_only.csv').read_bytes()).hexdigest(),
 'columns':cols,'created_utc':report['created_utc']}
(R/'manifests/phase4_baseline_inputs.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(suppressed(report),indent=2))
