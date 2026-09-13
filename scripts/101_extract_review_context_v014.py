"""Extract additional exposure context inside the controlled server.
Patient-level rows never leave cache/internal_review_v014. Public outputs suppress small and linked cells.
"""
from pathlib import Path
import datetime, hashlib, json, os, platform
import duckdb
import pandas as pd

os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery')
A=R/'analysis_v0.14'
P=A/'config/analysis_plan_v0.14.json'
plan=json.loads(P.read_text(encoding='utf-8-sig'))
expected=plan['inputs']
for rel,h in expected.items():
    got=hashlib.sha256((R/rel).read_bytes()).hexdigest()
    if got!=h: raise RuntimeError(f'input hash mismatch: {rel}: {got}')

controlled=R/'cache/internal_review_v014'
controlled.mkdir(parents=True,exist_ok=True)
os.chmod(controlled,0o700)
out=A/'output'; out.mkdir(parents=True,exist_ok=True)

c=duckdb.connect(str(R/'cache/audit_combined.duckdb'),read_only=True)
c.execute("ATTACH '/root/projects/linezolid_platelet_recovery/cache/phase11.duckdb' AS p11 (READ_ONLY)")
c.execute("CREATE TEMP TABLE cohort AS SELECT * FROM read_csv_auto('/root/projects/linezolid_platelet_recovery/cache/phase10_baseline.csv') WHERE routine_bacterial_culture72h_shared AND NOT lmwh_available7d")
assert c.execute('select count(*) from cohort').fetchone()[0]==2557
assert c.execute("select count(*) from cohort where drug='linezolid'").fetchone()[0]==65

q="""
WITH idx AS (
 SELECT c.hadm_id,c.subject_id,c.drug,c.t0,c.dischtime,c.deathtime,
        bool_or(a.data_source='inputevents') AS has_input,
        bool_or(upper(coalesce(ad.route,''))='IV') AS has_emar_iv,
        bool_or(upper(coalesce(ad.route,'')) IN ('PO','NG','ORAL','ENTERAL')) AS has_enteral,
        string_agg(DISTINCT coalesce(ad.route,CASE WHEN a.data_source='inputevents' THEN 'IV-input' ELSE '[missing]' END),';' ORDER BY coalesce(ad.route,CASE WHEN a.data_source='inputevents' THEN 'IV-input' ELSE '[missing]' END)) AS route_sources
 FROM cohort c
 LEFT JOIN systemic_administrations a ON a.hadm_id=c.hadm_id AND a.drug=c.drug AND a.charttime=c.t0
 LEFT JOIN administration_audit ad ON ad.hadm_id=a.hadm_id AND ad.drug=a.drug AND ad.emar_id=a.emar_id
 GROUP BY c.hadm_id,c.subject_id,c.drug,c.t0,c.dischtime,c.deathtime
), classified AS (
 SELECT *,CASE WHEN coalesce(has_input,false) OR coalesce(has_emar_iv,false) THEN 'parenteral'
               WHEN coalesce(has_enteral,false) THEN 'enteral'
               ELSE 'unknown' END AS initiation_route,
        CASE WHEN deathtime IS NOT NULL AND deathtime<dischtime THEN deathtime ELSE dischtime END AS observed_exit
 FROM idx
), traj AS (
 SELECT s.hadm_id,s.subject_id,s.drug,s.t0,s.initiation_route,s.route_sources,
        count(a.charttime) FILTER(WHERE a.drug=s.drug) AS index_records_14d,
        count(DISTINCT cast(a.charttime AS DATE)) FILTER(WHERE a.drug=s.drug) AS index_calendar_days_14d,
        min(a.charttime) FILTER(WHERE a.drug<>s.drug AND a.charttime>s.t0) AS other_first_time
 FROM classified s
 LEFT JOIN systemic_administrations a ON a.hadm_id=s.hadm_id AND a.charttime>=s.t0
      AND a.charttime<=least(s.t0+INTERVAL 14 DAY,s.observed_exit)
 GROUP BY s.hadm_id,s.subject_id,s.drug,s.t0,s.initiation_route,s.route_sources
)
SELECT hadm_id,subject_id,drug,initiation_route,route_sources,index_records_14d,index_calendar_days_14d,
       other_first_time IS NOT NULL AS switched_14d,
       CASE WHEN other_first_time IS NULL THEN NULL ELSE date_diff('minute',t0,other_first_time)/60.0 END AS switch_delay_hours
FROM traj ORDER BY hadm_id
"""
d=c.execute(q).fetchdf()
assert len(d)==2557 and d.hadm_id.nunique()==2557
assert set(d.initiation_route)<= {'parenteral','enteral','unknown'}
exact=controlled/'EXPOSURE_CONTEXT_EXACT_v0.14.csv'
d.to_csv(exact,index=False)
os.chmod(exact,0o600)

# Public aggregate display. Values 1-9 and linked/complementary cells are suppressed.
def cell_count(n):
    n=int(n)
    return '<10' if 0<n<10 else str(n)
def qfmt(x,digits=2):
    return '' if pd.isna(x) else f'{float(x):.{digits}f}'
rows=[]
for drug,label in [('linezolid','LZD'),('vancomycin','VAN')]:
    z=d[d.drug==drug]
    rows.extend([
      dict(metric='Admissions',arm=label,value=cell_count(len(z)),unit='n'),
      dict(metric='Parenteral initiation',arm=label,value=('Suppressed linked value' if ((z.initiation_route=='enteral').sum()+(z.initiation_route=='unknown').sum()) in range(1,10) else cell_count((z.initiation_route=='parenteral').sum())),unit='n'),
      dict(metric='Enteral or unknown initiation',arm=label,value=cell_count((z.initiation_route!='parenteral').sum()),unit='n'),
      dict(metric='Index-drug administration records through observed day 14',arm=label,value=qfmt(z.index_records_14d.mean()),unit='mean'),
      dict(metric='Index-drug administration records through observed day 14',arm=label,value=qfmt(z.index_records_14d.median(),1),unit='median'),
      dict(metric='Distinct index-drug calendar days through observed day 14',arm=label,value=qfmt(z.index_calendar_days_14d.mean()),unit='mean'),
      dict(metric='Distinct index-drug calendar days through observed day 14',arm=label,value=qfmt(z.index_calendar_days_14d.median(),1),unit='median'),
      dict(metric='Opposite target drug initiated through observed day 14',arm=label,value=cell_count(z.switched_14d.sum()),unit='n'),
      dict(metric='Time to opposite target drug among switchers',arm=label,value=qfmt(z.loc[z.switched_14d,'switch_delay_hours'].median(),1),unit='median hours')
    ])
pub=pd.DataFrame(rows)
pub.to_csv(out/'TREATMENT_EXPOSURE_CONTEXT_v0.14.csv',index=False)

route_counts=d.groupby(['drug','initiation_route']).size().reset_index(name='n')
route_safe=[]
for drug in route_counts.drug.unique():
    g=route_counts[route_counts.drug==drug]
    linked_small=any(0<int(n)<10 for n in g.n)
    for _,r in g.iterrows():
        n=int(r.n)
        route_safe.append({'drug':r.drug,'route':r.initiation_route,'n':('Suppressed linked value' if linked_small else cell_count(n))})
summary={
 'status':'COMPLETE',
 'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'python':platform.python_version(),
 'duckdb':duckdb.__version__,
 'plan_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),
 'primary_admissions':len(d),
 'route_counts_suppressed':route_safe,
 'exact_output_server_only':str(exact),
 'public_output':'analysis_v0.14/output/TREATMENT_EXPOSURE_CONTEXT_v0.14.csv',
 'privacy':'No identifiers, shifted dates, patient-level trajectories or small exact cells are included in public outputs.',
 'interpretation':'Observed administration trajectory through death/discharge/day14; descriptive only and not an as-treated causal estimand.'
}
(A/'reports/SOURCE_AND_EXPOSURE_AUDIT_v0.14.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
c.close()
print(json.dumps({k:v for k,v in summary.items() if k!='exact_output_server_only'},indent=2))
