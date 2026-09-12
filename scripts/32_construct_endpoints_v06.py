"""Fixed first-event endpoint; individual data remain on server."""
from pathlib import Path
import os,json,datetime,hashlib,random
import duckdb,pandas as pd,numpy as np
os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery');H=14*86400
plan=R/'config/analysis_plan_v0.6.json'
assert hashlib.sha256(plan.read_bytes()).hexdigest()=='3a149ee3f5b15441bd7ca7d3f6cfd4f7ec4f924378a01d0323d66e5329d3fb18'
def sequential(t,high,maxgap=None):
    candidates=[]
    for j,(s,yes) in enumerate(zip(t,high)):
        if not yes:candidates=[];continue
        eligible=[i for i in candidates if s-t[i]>=86400 and (maxgap is None or s-t[i]<=maxgap)]
        if eligible:return (eligible[0],j)
        candidates.append(j)
    return None
def pair_search(t,high,maxgap=None):
    pairs=[(i,j) for j in range(len(t)) for i in range(j) if t[j]-t[i]>=86400
           and (maxgap is None or t[j]-t[i]<=maxgap) and all(high[i:j+1])]
    return min(pairs,key=lambda ij:(t[ij[1]],ij[0])) if pairs else None
def classify(rec,exit_s,dead):
    options=[(H,2,0),(exit_s,0,2 if dead else 3)]
    if rec is not None:options.append((rec,1,1))
    t,_,k=min(options);return t/86400,k
fixtures=[([3600,90000],[1,1],(0,1)),([3600,90000],[1,0],None),
          ([3600,45000,100000],[1,0,1],None),([3600,350000],[1,1],(0,1)),
          ([],[],None),([1000,87400,90000],[1,1,1],(0,1))]
for t,h,expect in fixtures:assert sequential(t,h)==pair_search(t,h)==expect
assert sequential([3600,350000],[1,1],72*3600) is None
assert classify(86400,86400,True)==(1,2) and classify(86400,86400,False)==(1,3)
assert classify(H,H+1000,False)==(14,1) and classify(None,H,False)==(14,3)
random.seed(42)
for _ in range(500):
    t=sorted(random.sample(range(1,1000),random.randrange(0,22)));t=[x*3600 for x in t]
    h=[random.choice([False,True]) for _ in t]
    for maxgap in [None,72*3600]:assert sequential(t,h,maxgap)==pair_search(t,h,maxgap)
print('Synthetic endpoint tests passed before real endpoint extraction',flush=True)
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'))
c.execute('''CREATE OR REPLACE TABLE phase6_platelet_canonical AS SELECT p.subject_id,p.hadm_id,p.specimen_id,
 p.charttime,p.valuenum,min(p.storetime) storetime,count(*) source_rows FROM target_labs p
 SEMI JOIN phase3_baseline_points e USING(hadm_id)
 WHERE p.itemid IN ('51265','53189') AND p.valueuom='K/uL' AND p.valuenum>=0 AND p.valuenum<5000
 GROUP BY p.subject_id,p.hadm_id,p.specimen_id,p.charttime,p.valuenum''')
drift=c.execute('''SELECT count(*) FROM phase3_baseline_points e LEFT JOIN LATERAL (
 SELECT charttime,min(valuenum) valuenum,count(DISTINCT valuenum) nvals FROM phase6_platelet_canonical x
 WHERE x.hadm_id=e.hadm_id AND x.charttime<e.t0 AND x.charttime>=e.t0-INTERVAL 24 HOUR AND x.storetime<=e.t0
 GROUP BY charttime ORDER BY charttime DESC LIMIT 1) b ON true
 WHERE b.nvals IS DISTINCT FROM 1 OR b.valuenum IS DISTINCT FROM e.known_platelets OR b.charttime IS DISTINCT FROM e.known_platelets_time''').fetchone()[0]
assert drift==0,'Canonical dedup changed baseline; no outcome analysis permitted'
admin_drift=c.execute('''SELECT count(*) FROM phase3_baseline_points e JOIN LATERAL (
 SELECT min(charttime) firsttime FROM systemic_administrations x WHERE x.hadm_id=e.hadm_id AND x.drug=e.drug) a ON true WHERE a.firsttime IS DISTINCT FROM e.t0''').fetchone()[0]
assert admin_drift==0
co=c.execute('''SELECT e.subject_id,e.hadm_id,e.drug,e.t0,e.admittime,e.dischtime,e.deathtime,
 e.known_platelets,e.known_platelets_time,e.known_platelets_store,e.age,e.decision_in_icu,e.ambiguous_latest_platelets,
 a.hospital_expire_flag FROM phase3_baseline_points e JOIN admissions a USING(hadm_id) ORDER BY e.hadm_id''').fetchdf()
assert len(co)==1576 and co.hadm_id.nunique()==1576
assert ((co.known_platelets<100)&(co.age>=18)&co.decision_in_icu&~co.ambiguous_latest_platelets).all()
assert (co.known_platelets_time<co.t0).all() and (co.known_platelets_store<=co.t0).all()
la=c.execute('SELECT hadm_id,charttime,valuenum,storetime FROM phase6_platelet_canonical ORDER BY hadm_id,charttime').fetchdf()
groups={h:g for h,g in la.groupby('hadm_id',sort=False)}
output=[];conflicts=0;proxy=0;stored_late=0;stored_missing=0;store_lags=[]
for e in co.itertuples():
    assert pd.notna(e.dischtime) and e.dischtime>e.t0
    if pd.notna(e.deathtime):assert e.deathtime>e.t0
    valid_death=pd.notna(e.deathtime) and e.deathtime<=e.dischtime
    died=valid_death or str(e.hospital_expire_flag)=='1'
    proxy_death=died and not valid_death;proxy+=int(proxy_death)
    exit_time=e.deathtime if valid_death else e.dischtime
    exit_s=(exit_time-e.t0).total_seconds();assert exit_s>0
    labs=groups[e.hadm_id];labs=labs[(labs.charttime>=e.admittime)&(labs.charttime<exit_time)]
    follow=labs[labs.charttime>e.t0]
    times=[];high=[];stores=[]
    for ct,gg in follow.groupby('charttime',sort=True):
        vals=gg.valuenum.unique();times.append((ct-e.t0).total_seconds())
        high.append(bool((vals>=100).all()))
        conflicts+=int((vals>=100).any() and (vals<100).any())
        stores.append(gg.storetime.min())
    pair=sequential(times,high);pair72=sequential(times,high,72*3600)
    assert pair==pair_search(times,high) and pair72==pair_search(times,high,72*3600)
    rec=None if pair is None else times[pair[1]]
    rec72=None if pair72 is None else times[pair72[1]]
    day,status=classify(rec,exit_s,died);day72,status72=classify(rec72,exit_s,died)
    lag=np.nan
    if status==1:
        store=stores[pair[1]]
        if pd.isna(store):stored_missing+=1
        else:
            lag=(store-(e.t0+pd.Timedelta(seconds=rec))).total_seconds()/3600
            store_lags.append(lag);stored_late+=int(store>=exit_time)
    last_value=np.nan;last_lag=np.nan;last_ambiguous=False;last_post=False
    if not died and len(labs):
        lasttime=labs.charttime.max();vs=labs.loc[labs.charttime==lasttime,'valuenum'].unique()
        last_ambiguous=len(vs)!=1
        if not last_ambiguous:last_value=float(vs[0])
        last_lag=(e.dischtime-lasttime).total_seconds()/3600;last_post=lasttime>e.t0
    output.append({'subject_id':e.subject_id,'hadm_id':e.hadm_id,'drug':e.drug,'event_day':day,'status':status,
      'event_day_72h':day72,'status_72h':status72,'live_discharge':not died,'discharge_day':(e.dischtime-e.t0).total_seconds()/86400,
      'confirmed_before_discharge':bool(not died and rec is not None and rec<exit_s),
      'last_platelets':last_value,'last_draw_to_discharge_hours':last_lag,'last_ambiguous':last_ambiguous,'last_post_t0':last_post,
      'confirmation_store_lag_hours':lag,'death_time_proxy':proxy_death})
out=pd.DataFrame(output);assert (out.event_day.between(0,14)).all() and (out.event_day_72h.between(0,14)).all()
out.to_csv(R/'cache/phase6_endpoints.csv',index=False)
def hide(x):
    if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
    if isinstance(x,list):return [hide(v) for v in x]
    if isinstance(x,bool):return x
    if isinstance(x,int) and 0<x<10:return '<10'
    return x
summary=out.groupby(['drug','status']).size().rename('n').reset_index().to_dict('records')
qc={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cohort_n':len(out),
 'synthetic_fixtures_passed':True,'randomized_pair_checks':1000,'real_admission_pair_checks':len(out)*2,
 'baseline_value_or_time_drift':drift,'first_dose_time_drift':admin_drift,
 'ambiguous_threshold_followup_timepoints':conflicts,'death_time_proxies':proxy,
 'confirmation_storetime_missing':stored_missing,'confirmation_storetime_at_or_after_exit':stored_late,
 'confirmation_store_lag_hours_quantiles':np.quantile(store_lags,[0,.25,.5,.75,1]).tolist() if store_lags else None,
 'first_event_counts':summary,'status_labels':{'0':'no first event by14d','1':'confirmed recovery','2':'death before recovery','3':'live discharge before recovery'},
 'limits':['Endpoint is retrospectively determined by specimen time, not online availability of its result.',
 'No early censoring at last lab. Primary recovery may be transfusion-supported; it is not spontaneous recovery.',
 'Group-specific recovery counts first opened in this run; no subsequent weighting choice may use their differences.']}
(R/'reports/ENDPOINT_CONSTRUCTION_v0.6.json').write_text(json.dumps(hide(qc),indent=2))
(R/'cache/internal_reports/ENDPOINT_CONSTRUCTION_v0.6_exact.json').write_text(json.dumps(qc,indent=2))
(R/'manifests/endpoints_v0.6.json').write_text(json.dumps({'plan_sha256':hashlib.sha256(plan.read_bytes()).hexdigest(),
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'endpoint_cache_sha256':hashlib.sha256((R/'cache/phase6_endpoints.csv').read_bytes()).hexdigest()},indent=2))
print(json.dumps(hide(qc),indent=2))
