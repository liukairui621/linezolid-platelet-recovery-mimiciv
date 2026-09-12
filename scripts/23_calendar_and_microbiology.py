"""v0.5: map decision-time calendar intervals; re-audit known organisms.
No clinical endpoints; individual-level outputs stay on server.
"""
from pathlib import Path
import duckdb,json,os,datetime,hashlib
os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery')
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'))
def rows(q):
 cur=c.execute(q);return [dict(zip([x[0] for x in cur.description],r)) for r in cur.fetchall()]
def hide(x):
 if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
 if isinstance(x,list):return [hide(v) for v in x]
 if isinstance(x,bool):return x
 if isinstance(x,int) and 0<x<10:return '<10'
 return x
c.execute('''CREATE OR REPLACE TABLE phase5_calendar AS SELECT e.hadm_id,e.subject_id,e.drug,
 e.anchor_year_group,
 try_cast(substr(p.anchor_year_group,1,4) AS INTEGER)+year(e.t0)-try_cast(p.anchor_year AS INTEGER) mapped_year_lower,
 try_cast(substr(p.anchor_year_group,8,4) AS INTEGER)+year(e.t0)-try_cast(p.anchor_year AS INTEGER) mapped_year_upper,
 CASE WHEN mapped_year_upper<=2013 THEN 'mapped_early'
 WHEN mapped_year_lower>=2014 THEN 'mapped_late' ELSE 'boundary_uncertain' END mapped_era,
 try_cast(substr(p.anchor_year_group,1,4) AS INTEGER)<=2013 early_anchor,
 e.t0 FROM phase3_baseline_points e JOIN patients p USING(subject_id)''')
assert c.execute('SELECT count(*)=1576 AND count(*)=count(DISTINCT hadm_id) FROM phase5_calendar').fetchone()[0]
c.execute(f"COPY (SELECT * EXCLUDE(t0) FROM phase5_calendar) TO '{R/'cache/phase5_calendar.csv'}' (HEADER,DELIMITER ',')")
c.execute('''CREATE OR REPLACE TABLE phase5_known_organisms AS SELECT e.hadm_id,e.drug,
 m.spec_type_desc,m.test_name,m.org_name,m.ab_name,m.interpretation,
 regexp_matches(upper(coalesce(m.test_name,'')||' '||coalesce(m.spec_type_desc,'')),'SCREEN|SURVEILLANCE') screen_label,
 (regexp_matches(upper(coalesce(m.org_name,'')),'MRSA|METHICILLIN.?RESISTANT.*STAPH') OR
 ((upper(m.org_name) LIKE '%STAPH AUREUS%' OR upper(m.org_name) LIKE '%STAPHYLOCOCCUS AUREUS%')
 AND upper(m.ab_name) IN ('OXACILLIN','METHICILLIN','CEFOXITIN') AND upper(m.interpretation)='R')) expanded_mrsa_flag
 FROM phase3_baseline_points e JOIN phase2_microbiology m ON m.hadm_id=e.hadm_id
 WHERE m.sample_available_by>=e.t0-INTERVAL 7 DAY AND m.sample_available_by<=e.t0
 AND m.result_available_by<=e.t0 AND m.org_name IS NOT NULL''')
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'calendar_mapping':'Decision-year range = anchor-year range + (year(t0)-anchor_year). Approximate calendar intervals, not exact reidentified dates.',
 'calendar_counts':rows('''SELECT drug,mapped_era,count(*) n FROM phase5_calendar GROUP BY ALL ORDER BY drug,mapped_era'''),
 'anchor_vs_calendar':rows('''SELECT drug,early_anchor,mapped_era,count(*) n FROM phase5_calendar GROUP BY ALL ORDER BY drug,early_anchor,mapped_era'''),
 'organism_coverage':rows('''SELECT e.drug,count(*) n,
 count(*) FILTER(WHERE x.any_known) known_any_nonempty_org_field,
 count(*) FILTER(WHERE x.nonscreen_known) known_nonscreen_nonempty_org_field,
 count(*) FILTER(WHERE x.nonscreen_named_result) known_nonscreen_named_result,
 count(*) FILTER(WHERE x.mrsa_known) expanded_mrsa_including_screen,
 count(*) FILTER(WHERE x.mrsa_nonscreen) expanded_mrsa_nonscreen
 FROM phase3_baseline_points e LEFT JOIN LATERAL(SELECT count(*)>0 any_known,
 count(*) FILTER(WHERE NOT screen_label)>0 nonscreen_known,
 count(*) FILTER(WHERE NOT screen_label AND trim(upper(org_name)) NOT IN ('','CANCELLED','POSITIVE','NEGATIVE','NO GROWTH'))>0 nonscreen_named_result,
 bool_or(expanded_mrsa_flag) mrsa_known,
 bool_or(expanded_mrsa_flag AND NOT screen_label) mrsa_nonscreen
 FROM phase5_known_organisms x WHERE x.hadm_id=e.hadm_id) x ON true GROUP BY e.drug'''),
 'known_organism_terms':rows('''SELECT drug,org_name,screen_label,count(DISTINCT hadm_id) admissions FROM phase5_known_organisms GROUP BY ALL ORDER BY drug,org_name'''),
 'limits':['No practice trend can be inferred from patient anchor groups. Case-mix and denominator selection remain even after correct calendar mapping.',
 'Nonempty org_name fields include generic POSITIVE/CANCELLED tokens. The named-result count excludes these fixed tokens but remains a record-label diagnostic, not adjudicated infection.',
 'Absence of nonscreen MRSA by this rule is not absence of screening evidence, other pathogens, outside-hospital history, prior colonization or clinician rationale.',
 'Culture timing and blood-culture proportions support a similar infection-investigation setting, not therapeutic equivalence or proven equipoise.']}
(R/'cache/internal_reports/CALENDAR_MICROBIOLOGY_v0.5_exact.json').write_text(json.dumps(report,indent=2))
(R/'reports/CALENDAR_MICROBIOLOGY_v0.5.json').write_text(json.dumps(hide(report),indent=2))
(R/'manifests/phase5_calendar.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'snapshot_sha256':hashlib.sha256((R/'cache/phase5_calendar.csv').read_bytes()).hexdigest()},indent=2))
print(json.dumps(hide(report),indent=2))
