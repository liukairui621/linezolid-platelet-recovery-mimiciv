"""Correct an impossible baseline clock; retain v0.6 failed run and provenance."""
from pathlib import Path
import duckdb,json,hashlib,datetime
R=Path('/root/projects/linezolid_platelet_recovery');c=duckdb.connect(str(R/'cache/audit_combined.duckdb'))
c.execute('''CREATE OR REPLACE TABLE phase61_eligible_points AS SELECT * FROM phase3_baseline_points
 WHERE t0>=admittime AND t0<dischtime AND (deathtime IS NULL OR t0<deathtime)''')
counts=dict(c.execute('SELECT drug,count(*) FROM phase61_eligible_points GROUP BY drug').fetchall())
reason=c.execute('''SELECT drug,count(*) FILTER(WHERE t0>=dischtime) at_or_after_discharge,
 count(*) FILTER(WHERE deathtime<=t0) at_or_after_death,
 count(*) FILTER(WHERE t0<admittime) before_admission,
 count(*)-count(*) FILTER(WHERE t0>=admittime AND t0<dischtime AND (deathtime IS NULL OR t0<deathtime)) excluded_union
 FROM phase3_baseline_points GROUP BY drug''')
rows=[dict(zip([d[0] for d in reason.description],r)) for r in reason.fetchall()]
c.execute(f"COPY (SELECT hadm_id FROM phase61_eligible_points ORDER BY hadm_id) TO '{R/'cache/phase61_eligible_ids.csv'}' (HEADER,DELIMITER ',')")
p=json.loads((R/'config/analysis_plan_v0.6.json').read_text());p['version']='0.6.1'
p['amendment_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
p['replaces_for_future_analysis']='0.6 after endpoint precondition failure; historical model and failure retained'
p['registration_status']='Internal baseline-clock correction after Firth diagnostic results, before any group recovery summary or treatment contrast was read. Not external preregistration.'
p['clinical']['phase6_execution']['population']='phase61_eligible_points: require t0 within index admission and before any recorded death. Corrects impossible timestamps; no requirement for future treatment duration, testing, survival or recovery.'
p['clinical']['phase6_execution']['cohort_n']=sum(counts.values());p['clinical']['phase6_execution']['cohort_counts']=counts
p['clinical']['phase6_execution']['clock_amendment']='Original endpoint constructor stopped on t0>=discharge; aggregate audit found VAN records after discharge/death. They have no valid index-admission baseline risk interval. Rebuild cohort and all weights; never code them as0followup recovery failures. Original datasets unchanged.'
q=R/'config/analysis_plan_v0.6.1.json';assert not q.exists();q.write_text(json.dumps(p,indent=2,ensure_ascii=False))
def hide(x):
 if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
 if isinstance(x,list):return [hide(v) for v in x]
 if isinstance(x,int) and not isinstance(x,bool) and 0<x<10:return '<10'
 return x
report={'created_utc':p['amendment_utc'],'counts_after_clock_check':counts,'reason_counts_overlap':rows,
 'plan_sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'eligibility_ids_sha256':hashlib.sha256((R/'cache/phase61_eligible_ids.csv').read_bytes()).hexdigest(),
 'reason':'Baseline time must be during index hospitalization and before death. Original13_ verification covered switch risk sets, not all first-initiation entries; this check is now explicit.'}
(R/'reports/ELIGIBILITY_CLOCK_CORRECTION_v0.6.1.json').write_text(json.dumps(hide(report),indent=2))
(R/'cache/internal_reports/ELIGIBILITY_CLOCK_CORRECTION_v0.6.1_exact.json').write_text(json.dumps(report,indent=2))
print(json.dumps(hide(report),indent=2))
