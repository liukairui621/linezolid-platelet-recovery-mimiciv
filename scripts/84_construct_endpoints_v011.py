"""Rebuild fixed endpoints from refreshed sources and compare overlapping historical admissions."""
from pathlib import Path
import ast,os,json,random,datetime,hashlib
import duckdb,pandas as pd,numpy as np
os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery');D=86400;H=14*D
cfg=json.loads((R/'config/analysis_plan_v0.7.json').read_text())['phase7']
for source,names in [('32b_construct_endpoints_v061.py',['sequential','pair_search','classify']),('37_observation_endpoints_v07.py',['merge','blocked','pair','brute','event','reason'])]:
 tree=ast.parse((R/'scripts'/source).read_text())
 for node in tree.body:
  if isinstance(node,ast.FunctionDef) and node.name in names:exec(compile(ast.Module(body=[node],type_ignores=[]),source,'exec'))
random.seed(20260915)
for _ in range(1000):
 t=sorted(random.sample(range(1,500),random.randrange(0,25)));t=[x*3600 for x in t];hi=[random.choice([False,True]) for _ in t]
 for gap in [None,72*3600]:assert sequential(t,hi,gap)==pair_search(t,hi,gap)
 assert pair(t,hi)==brute(t,hi)==sequential(t,hi)
 assert pair(t,hi,[(30*3600,31*3600,True)])==brute(t,hi,[(30*3600,31*3600,True)])
assert sequential([],[]) is None
assert sequential([1,1+D],[True,True])==(0,1)
assert sequential([1,1+D,1+2*D],[True,False,True]) is None
assert sequential([1,1+4*D],[True,True],72*3600) is None
assert classify(D,D,True)==(1,2) and classify(D,D,False)==(1,3)
assert classify(H,H+1,False)==(14,1) and classify(None,H,False)==(14,3)
assert pair([D,2*D,4*D,5*D],[1,1,1,1],[(D,D,True)])==(2,3)
assert not blocked(10,20,[(0,10,False)]) and blocked(0,10,[(10,11,True)])
assert merge([(0,1),(1,2),(3,4)])==[[0,2],[3,4]]
print('Synthetic endpoint tests passed before production summaries',flush=True)
c=duckdb.connect(str(R/'cache/phase11.duckdb'),read_only=True)
co=c.execute('SELECT e.*,a.admittime,a.dischtime,a.deathtime,a.hospital_expire_flag FROM phase11_selected e JOIN phase11_admissions a USING(hadm_id,subject_id) ORDER BY hadm_id').fetchdf()
assert len(co)==2593 and co.hadm_id.nunique()==2593
labs=c.execute('SELECT * FROM phase11_platelet_canonical ORDER BY hadm_id,charttime').fetchdf()
icu=c.execute('SELECT * FROM phase11_icustays').fetchdf();loc=c.execute('SELECT * FROM phase11_transfers').fetchdf();loc=loc[loc.careunit.isin(cfg['ICU']['careunits'])]
tx=c.execute("SELECT * FROM phase11_platelet_inputs WHERE amount>0 AND coalesce(statusdescription,'')<>'Rewritten'").fetchdf()
dated=c.execute('SELECT * FROM phase11_platelet_procedures').fetchdf();emar=c.execute('SELECT count(*) FROM phase11_platelet_emar').fetchone()[0]
assert emar==0,'Nonempty platelet eMAR requires administration/route audit before transfusion branch'
groups=lambda df:{h:g for h,g in df.groupby('hadm_id',sort=False)}
lg,ig,pg,tg,dg=map(groups,[labs,icu,loc,tx,dated])
out=[];qc=dict(synthetic_tests_passed=True,randomized_case_sets=1000,physical_missing_t0=0,invalid_interval_rows=0,transfusion_missing_end_rows=0,transfusion_negative_interval_rows=0,threshold_ambiguous_timepoints=0)
for e in co.itertuples():
 sec=lambda ts:(ts-e.t0).total_seconds()
 assert pd.notna(e.dischtime) and e.dischtime>e.t0 and (pd.isna(e.deathtime) or e.deathtime>e.t0)
 valid_death=pd.notna(e.deathtime) and e.deathtime<=e.dischtime;died=valid_death or str(e.hospital_expire_flag)=='1';exit_time=e.deathtime if valid_death else e.dischtime;end=sec(exit_time)
 def intervals(g):
  vals=[]
  if g is not None:
   for x in g.itertuples():
    if pd.isna(x.intime) or pd.isna(x.outtime) or x.outtime<=x.intime:qc['invalid_interval_rows']+=1;continue
    vals.append((sec(x.intime),sec(x.outtime)))
  return merge(vals)
 all_icu=intervals(ig.get(e.hadm_id));ep0=[x for x in all_icu if x[0]<=0<x[1]];assert len(ep0)==1;icu_end=ep0[0][1]
 physical=intervals(pg.get(e.hadm_id));ph0=[x for x in physical if x[0]<=0<x[1]];assert len(ph0)==1,'Frozen cohort lacks physical ICU episode';physical_end=ph0[0][1]
 all_labs=lg[e.hadm_id];all_labs=all_labs[(all_labs.charttime>=e.admittime)&(all_labs.charttime<exit_time)]
 follow=all_labs[all_labs.charttime>e.t0];gt=follow.groupby('charttime',sort=True)
 times=np.array([sec(x) for x in gt.groups]);high=np.array([bool((g.valuenum>=100).all()) for _,g in gt]);stores=[g.storetime.min() for _,g in gt]
 qc['threshold_ambiguous_timepoints']+=sum((g.valuenum>=100).any() and (g.valuenum<100).any() for _,g in gt)
 fullpair=sequential(times,high);assert fullpair==pair_search(times,high)
 fullrec=None if fullpair is None else times[fullpair[1]]
 t=times[times<=H];hi=high[times<=H];blocks=[];dateblocks=[]
 if e.hadm_id in tg:
  for x in tg[e.hadm_id].itertuples():
   assert pd.notna(x.starttime);ee=x.endtime
   if pd.isna(ee):ee=x.starttime;qc['transfusion_missing_end_rows']+=1
   if ee<x.starttime:ee=x.starttime;qc['transfusion_negative_interval_rows']+=1
   blocks.append((sec(x.starttime),sec(ee),True))
 if e.hadm_id in dg:
  for x in dg[e.hadm_id].itertuples():
   if pd.notna(x.chartdate):a=sec(pd.Timestamp(x.chartdate));dateblocks.append((a,a+D,False))
 r=dict(subject_id=e.subject_id,hadm_id=e.hadm_id,drug=e.drug,primary_population=e.primary_population)
 for prefix,limit,dead in [('hospital',end,died),('icu',min(end,icu_end),died and end<=icu_end),('physical',min(end,physical_end),died and end<=physical_end)]:
  sel=t<limit;tt=t[sel];hh=hi[sel]
  for suffix,bb in [('',()),('_timed',blocks),('_dated',blocks+dateblocks)]:
   pp=pair(tt,hh,bb);assert pp==brute(tt,hh,bb);rec=None if pp is None else tt[pp[1]]
   r[prefix+suffix+'_day'],r[prefix+suffix+'_status']=event(rec,limit,dead)
  first=next((v for v,h in zip(tt,hh) if h),None)
  r[prefix+'_single_day'],r[prefix+'_single_status']=event(first,limit,dead)
  r[prefix+'_tests']=len(tt);r[prefix+'_days']=min(limit,H)/D;r[prefix+'_tests_per_day']=len(tt)/r[prefix+'_days']
 pp72=sequential(t,hi,72*3600);assert pp72==pair_search(t,hi,72*3600)
 r['hospital72_day'],r['hospital72_status']=event(None if pp72 is None else t[pp72[1]],end,died)
 r['reason']=reason(t,hi,end);first=next((v for v,h in zip(t,hi) if h),None);r['first_high_day']=np.nan if first is None else first/D
 ff=follow[follow.charttime<=e.t0+pd.Timedelta(days=14)]
 r['specimen_groups']=len(ff[['specimen_id','charttime']].drop_duplicates());r['canonical_rows']=len(ff)
 r['hospital_icu_days']=sum(max(0,min(b,H,end)-max(0,a)) for a,b in all_icu)/D;r['hospital_outside_icu_days']=r['hospital_days']-r['hospital_icu_days']
 r['hospital_icu_tests']=sum(any(a<=v<b for a,b in all_icu) for v in t);r['hospital_outside_icu_tests']=len(t)-r['hospital_icu_tests']
 r['any_observed_low_after_first_high']=False if first is None else any(not h for v,h in zip(t,hi) if v>first)
 r['any_test_24h_after_first_high']=False if first is None else any(v>=first+D for v in t)
 r['physical_exit_difference_hours']=(physical_end-icu_end)/3600
 r.update(live_discharge=not died,discharge_day=sec(e.dischtime)/D,confirmed_before_discharge=bool(not died and fullrec is not None and fullrec<end),last_platelets=np.nan,last_draw_to_discharge_hours=np.nan,last_ambiguous=False,last_post_t0=False,confirmation_store_lag_hours=np.nan,confirmation_stored_at_or_after_exit=False,death_time_proxy=bool(died and not valid_death))
 if r['hospital_status']==1:
  store=stores[fullpair[1]]
  if pd.notna(store):r['confirmation_store_lag_hours']=(store-(e.t0+pd.Timedelta(seconds=fullrec))).total_seconds()/3600;r['confirmation_stored_at_or_after_exit']=store>=exit_time
 if not died and len(all_labs):
  lt=all_labs.charttime.max();vs=all_labs.loc[all_labs.charttime==lt,'valuenum'].unique();r['last_ambiguous']=len(vs)!=1
  if len(vs)==1:r['last_platelets']=float(vs[0])
  r['last_draw_to_discharge_hours']=sec(e.dischtime)/3600-sec(lt)/3600;r['last_post_t0']=lt>e.t0
 out.append(r)
out=pd.DataFrame(out)
old7=pd.read_csv(R/'cache/phase7_endpoints.csv',dtype={'hadm_id':str,'subject_id':str});j=out.merge(old7,on='hadm_id',suffixes=('_new','_old'),validate='1:1');assert len(j)==1185
cmp={}
for col in set(out)&set(old7)-{'hadm_id'}:
 a,b=j[col+'_new'],j[col+'_old']
 if pd.api.types.is_numeric_dtype(a):bad=~np.isclose(a.astype(float),b.astype(float),atol=1e-9,rtol=0,equal_nan=True)
 else:bad=~((a==b)|(a.isna()&b.isna()))
 cmp[col]=int(bad.sum())
assert not any(cmp.values()),json.dumps(cmp)
old61=pd.read_csv(R/'cache/phase61_endpoints.csv',dtype={'hadm_id':str,'subject_id':str})
j=out.merge(old61,on='hadm_id',suffixes=('_new','_old'),validate='1:1');cmp61={}
mapping={'hospital_day':'event_day','hospital_status':'status','hospital72_day':'event_day_72h','hospital72_status':'status_72h'}
for nc,oc in mapping.items():cmp61[nc]=int((~np.isclose(j[nc],j[oc],atol=1e-9,rtol=0,equal_nan=True)).sum())
for col in ['live_discharge','discharge_day','confirmed_before_discharge','last_platelets','last_draw_to_discharge_hours','last_ambiguous','last_post_t0','confirmation_store_lag_hours','death_time_proxy']:
 cmp61[col]=int((~np.isclose(j[col+'_new'].astype(float),j[col+'_old'].astype(float),atol=1e-8,rtol=0,equal_nan=True)).sum())
assert not any(cmp61.values()),json.dumps(cmp61)
assert ((out.hospital_status==1)==(out.reason=='confirmed')).all()
assert ((out.hospital_single_status!=1)==(out.reason=='never_threshold')).all()
for prefix in ['hospital','icu','physical']:
 for restrictive,liberal in [('_dated','_timed'),('_timed',''),('','_single')]:assert ((out[prefix+restrictive+'_status']!=1)|(out[prefix+liberal+'_status']==1)).all()
 assert out[prefix+'_day'].between(0,14).all()
assert ((out.physical_status!=1)|(out.hospital_status==1)).all()
out.to_csv(R/'cache/phase11_endpoints.csv',index=False)
qc.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),admissions=len(out),primary_admissions=int(out.primary_population.sum()),overlapping_historical_admissions=len(j),old_v07_field_disagreements=cmp,old_v061_field_disagreements=cmp61,platelet_emar_rows=emar,timed_transfusion_rows=len(tx),dated_transfusion_rows=len(dated),death_time_proxies=int(out.death_time_proxy.sum()),confirmation_storetime_missing=int(out.loc[out.hospital_status==1,'confirmation_store_lag_hours'].isna().sum()),confirmation_stored_at_or_after_exit=int(out.confirmation_stored_at_or_after_exit.sum()),endpoint_sha256=hashlib.sha256((R/'cache/phase11_endpoints.csv').read_bytes()).hexdigest(),limits=['No record of transfusion does not establish spontaneous recovery.','Historical overlap is an implementation check, not an independent cohort.','No outcome-conditioned changes to cohort or weights.'])
def hide(x):
 if isinstance(x,dict):return {k:hide(v) for k,v in x.items()}
 if isinstance(x,list):return [hide(v) for v in x]
 if isinstance(x,(int,np.integer)) and not isinstance(x,bool):return '<10' if 0<x<10 else int(x)
 return x
(R/'cache/internal_reports/ENDPOINT_CHECKS_v0.11_exact.json').write_text(json.dumps(qc,indent=2,default=lambda x:x.item()))
(R/'reports/ENDPOINT_CHECKS_v0.11.json').write_text(json.dumps(hide(qc),indent=2))
print(json.dumps(hide(qc),indent=2))
