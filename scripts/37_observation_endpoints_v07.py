"""Independent reconstruction and pre-frozen observation/transfusion sensitivities."""
from pathlib import Path
import os,json,datetime,hashlib,random
import duckdb,pandas as pd,numpy as np
os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery');D=86400;H=14*D
plan=R/'config/analysis_plan_v0.7.json';cfg=json.loads(plan.read_text())['phase7']
def merge(intervals):
    out=[]
    for s,e in sorted(intervals):
        if not out or s>out[-1][1]:out.append([s,e])
        else:out[-1][1]=max(out[-1][1],e)
    return out
def blocked(lo,hi,blocks):
    return any(s<=hi and (e>=lo if closed else e>lo) for s,e,closed in blocks)
def pair(t,high,blocks=()):
    run=[]
    for j in range(len(t)):
        if not high[j]:run=[];continue
        for i in run:
            if t[j]-t[i]>=D and not blocked(t[i]-2*D,t[j],blocks):return i,j
        run.append(j)
    return None
def brute(t,high,blocks=()):
    return next(((i,j) for j in range(len(t)) for i in range(j) if t[j]-t[i]>=D and all(high[i:j+1]) and not blocked(t[i]-2*D,t[j],blocks)),None)
def event(rec,end,dead):
    # Exit wins at same timestamp; a specimen at exactly14d can qualify.
    opts=[(H,2,0),(end,0,2 if dead else 3)]
    if rec is not None:opts.append((rec,1,1))
    tt,_,status=min(opts);return tt/D,status
def reason(t,high,end):
    pp=pair(t,high)
    if pp:return 'confirmed'
    hh=np.flatnonzero(high)
    if not len(hh):return 'never_threshold'
    i=hh[0]
    if not all(high[i:]):return 'unconfirmed_observed_low'
    eligible_time=t[i]+D
    if eligible_time>end or eligible_time==end and end<=H or eligible_time>H:
        # Caller separately supplies actual exit, so end==H remains an exit only when actual exit==H.
        return 'unconfirmed_insufficient_clock'
    assert not any(t[j]>=eligible_time for j in range(i+1,len(t)))
    return 'unconfirmed_no_late_measurement'
random.seed(77)
for _ in range(500):
    t=sorted(random.sample(range(1,500),random.randrange(0,25)));t=[v*3600 for v in t]
    h=[random.choice([True,False]) for _ in t];b=[(30*3600,31*3600,True)]
    assert pair(t,h)==brute(t,h) and pair(t,h,b)==brute(t,h,b)
assert pair([D,2*D,4*D,5*D],[1,1,1,1],[(D,D,True)])==(2,3)
assert event(D,D,True)==(1,2) and event(H,H+1,False)==(14,1)
assert blocked(0,10,[(10,11,True)]) and not blocked(10,20,[(0,10,False)])
assert merge([(0,1),(1,2),(3,4)])==[[0,2],[3,4]]
print('Endpoint and interval fixtures passed',flush=True)
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'),read_only=True)
for name in ['phase2_transfers','phase2_platelet_emar','phase2_platelet_procedures']:
    c.execute(f"CREATE TEMP VIEW {name} AS SELECT * FROM read_parquet('{R/'cache'/name}.parquet')")
co=c.execute('''SELECT e.subject_id,e.hadm_id,e.drug,e.t0,e.dischtime,e.deathtime,a.hospital_expire_flag
 FROM phase61_eligible_points e JOIN admissions a USING(hadm_id) ORDER BY hadm_id''').fetchdf()
assert len(co)==1569 and co.hadm_id.nunique()==1569
labs=c.execute('SELECT * FROM phase61_platelet_canonical ORDER BY hadm_id,charttime').fetchdf()
icu=c.execute('''SELECT hadm_id,try_cast(intime AS TIMESTAMP) intime,try_cast(outtime AS TIMESTAMP) outtime
 FROM icustays SEMI JOIN phase61_eligible_points USING(hadm_id)''').fetchdf()
loc=c.execute('SELECT * FROM phase2_transfers SEMI JOIN phase61_eligible_points USING(hadm_id)').fetchdf()
loc=loc[loc.careunit.isin(cfg['ICU']['careunits'])]
tx=c.execute('''SELECT t.* FROM platelet_transfusions_icu t SEMI JOIN phase61_eligible_points USING(hadm_id)
 WHERE amount>0 AND coalesce(statusdescription,'')<>'Rewritten' ''').fetchdf()
dated=c.execute('SELECT * FROM phase2_platelet_procedures SEMI JOIN phase61_eligible_points USING(hadm_id)').fetchdf()
emar=c.execute('SELECT count(*) FROM phase2_platelet_emar').fetchone()[0];assert emar==0
def groups(df):return {h:g for h,g in df.groupby('hadm_id',sort=False)}
lg,ig,pg,tg,dg=map(groups,[labs,icu,loc,tx,dated])
out=[];no_physical=0;icu_mismatch=0;invalid_icu=0;tx_bad_end=0;tx_missing_end=0
qc={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'synthetic_tests_passed':True}
for e in co.itertuples():
    sec=lambda ts:(ts-e.t0).total_seconds()
    valid_death=pd.notna(e.deathtime) and e.deathtime<=e.dischtime
    died=valid_death or str(e.hospital_expire_flag)=='1'
    end=sec(e.deathtime if valid_death else e.dischtime)
    assert end>0
    def intervals(group):
        global invalid_icu
        vals=[]
        if group is not None:
            for x in group.itertuples():
                if pd.isna(x.intime) or pd.isna(x.outtime) or x.outtime<=x.intime:
                    invalid_icu+=1;continue
                vals.append((sec(x.intime),sec(x.outtime)))
        return merge(vals)
    all_icu=intervals(ig.get(e.hadm_id));containing=[x for x in all_icu if x[0]<=0<x[1]]
    assert len(containing)==1,'Baseline ICU coverage mismatch; stop branch'
    icu_end=containing[0][1]
    physical=intervals(pg.get(e.hadm_id));physical0=[x for x in physical if x[0]<=0<x[1]]
    if len(physical0)!=1:no_physical+=1;physical_end=np.nan
    else:physical_end=physical0[0][1];icu_mismatch+=int(physical_end!=icu_end)
    ff=lg[e.hadm_id];ff=ff[(ff.charttime>e.t0)&(ff.charttime<e.t0+pd.Timedelta(seconds=end))&(ff.charttime<=e.t0+pd.Timedelta(days=14))]
    gg=ff.groupby('charttime',sort=True).valuenum.min()
    t=np.array([sec(x) for x in gg.index]);high=(gg.values>=100)
    blocks=[];dateblocks=[]
    if e.hadm_id in tg:
        for x in tg[e.hadm_id].itertuples():
            assert pd.notna(x.starttime)
            ee=x.endtime
            if pd.isna(ee):ee=x.starttime;tx_missing_end+=1
            if ee<x.starttime:ee=x.starttime;tx_bad_end+=1
            blocks.append((sec(x.starttime),sec(ee),True))
    if e.hadm_id in dg:
        for x in dg[e.hadm_id].itertuples():
            if pd.notna(x.chartdate):
                s=sec(pd.Timestamp(x.chartdate));dateblocks.append((s,s+D,False))
    r={'subject_id':e.subject_id,'hadm_id':e.hadm_id,'drug':e.drug}
    for prefix,limit,dd in [('hospital',end,died),('icu',min(end,icu_end),died and end<=icu_end),('physical',min(end,physical_end) if np.isfinite(physical_end) else np.nan,died and end<=physical_end)]:
        if not np.isfinite(limit):continue
        sel=t<limit;tt=t[sel];hh=high[sel]
        for suffix,bb in [('',()),('_timed',blocks),('_dated',blocks+dateblocks)]:
            pp=pair(tt,hh,bb);assert pp==brute(tt,hh,bb)
            rec=None if pp is None else tt[pp[1]]
            day,status=event(rec,limit,dd)
            r[prefix+suffix+'_day']=day;r[prefix+suffix+'_status']=status
        first=next((x for x,y in zip(tt,hh) if y),None)
        day,status=event(first,limit,dd)
        r[prefix+'_single_day']=day;r[prefix+'_single_status']=status
        r[prefix+'_tests']=len(tt);r[prefix+'_days']=min(limit,H)/D
        r[prefix+'_tests_per_day']=len(tt)/r[prefix+'_days']
    pp=pair(t,high);r['reason']=reason(t,high,end)
    first=next((x for x,y in zip(t,high) if y),None)
    r['first_high_day']=np.nan if first is None else first/D
    r['specimen_groups']=len(ff[['specimen_id','charttime']].drop_duplicates())
    r['canonical_rows']=len(ff)
    r['hospital_icu_days']=sum(max(0,min(b,H,end)-max(0,a)) for a,b in all_icu)/D
    r['hospital_outside_icu_days']=r['hospital_days']-r['hospital_icu_days']
    r['hospital_icu_tests']=sum(any(a<=v<b for a,b in all_icu) for v in t)
    r['hospital_outside_icu_tests']=len(t)-r['hospital_icu_tests']
    r['any_observed_low_after_first_high']=False if first is None else any(not y for v,y in zip(t,high) if v>first)
    r['any_test_24h_after_first_high']=False if first is None else any(v>=first+D for v in t)
    r['physical_exit_difference_hours']=(physical_end-icu_end)/3600
    out.append(r)
out=pd.DataFrame(out)
old=pd.read_csv(R/'cache/phase61_endpoints.csv',dtype={'hadm_id':str})
j=out.merge(old,on='hadm_id',validate='1:1')
assert (j.hospital_status==j.status).all() and np.max(abs(j.hospital_day-j.event_day))<1e-10
assert ((out.hospital_status==1)==(out.reason=='confirmed')).all()
assert ((out.hospital_single_status!=1)==(out.reason=='never_threshold')).all()
assert (out.hospital_icu_days<=out.hospital_days+1e-9).all()
for prefix in ['hospital','icu']:
    assert ((out[prefix+'_dated_status']!=1)|(out[prefix+'_timed_status']==1)).all()
    assert ((out[prefix+'_timed_status']!=1)|(out[prefix+'_status']==1)).all()
out.to_csv(R/'cache/phase7_endpoints.csv',index=False)
qc.update({'n_admissions':len(out),'primary_endpoint_mismatches':0,'physical_missing_t0':no_physical,'physical_vs_episode_exit_mismatches':icu_mismatch,
 'invalid_interval_rows':invalid_icu,'transfusion_missing_end_rows':tx_missing_end,'transfusion_negative_interval_rows':tx_bad_end,
 'platelet_emar_rows':emar,'physical_branch_eligible':no_physical==0,
 'reason_counts':out.groupby(['drug','reason']).size().rename('n').reset_index().to_dict('records'),
 'source_rows':{'timed_icu':len(tx),'dated_hospital_procedures':len(dated)},
 'limits':['Reason categories are hierarchical descriptions, not causal mechanisms. Observed lows can coexist with inadequate confirmation opportunity.',
 'Hospital transfusion-free status is not known from an ICU input source and day-resolution procedure codes.',
 'ICU sensitivity changes the competing-exit estimand and does not standardize monitoring intensity.']})
def hide(v):
    if isinstance(v,dict):return {k:hide(x) for k,x in v.items()}
    if isinstance(v,list):return [hide(x) for x in v]
    if isinstance(v,int) and not isinstance(v,bool) and 0<v<10:return '<10'
    return v
(R/'cache/internal_reports/OBSERVATION_ENDPOINTS_v0.7_exact.json').write_text(json.dumps(qc,indent=2))
(R/'reports/OBSERVATION_ENDPOINTS_v0.7.json').write_text(json.dumps(hide(qc),indent=2))
(R/'manifests/endpoints_v0.7.json').write_text(json.dumps({'plan_sha256':hashlib.sha256(plan.read_bytes()).hexdigest(),
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'endpoint_sha256':hashlib.sha256((R/'cache/phase7_endpoints.csv').read_bytes()).hexdigest()},indent=2))
print(json.dumps(hide(qc),indent=2))
