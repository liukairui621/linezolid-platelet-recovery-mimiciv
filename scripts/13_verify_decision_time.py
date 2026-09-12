"""Independent pandas/numpy checks against source administration/lab rows on server."""
from pathlib import Path
import json,os
import duckdb,pandas as pd,numpy as np
os.umask(0o077);R=Path('/root/projects/linezolid_platelet_recovery')
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'),read_only=True)
r=c.execute('select * from phase2_risk_eligible').fetchdf()
s=c.execute('select hadm_id,t0,other_first_time from phase2_switches').fetchdf().set_index('hadm_id')
admins=c.execute("select hadm_id,drug,charttime from systemic_administrations").fetchdf()
drug_times={(h,d):np.sort(g.charttime.to_numpy()) for (h,d),g in admins.groupby(['hadm_id','drug'])}
labs=c.execute('select hadm_id,charttime,storetime,valuenum from platelet_measurements').fetchdf()
lab_groups={h:(g.charttime.to_numpy(),g.storetime.to_numpy(),g.valuenum.to_numpy()) for h,g in labs.groupby('hadm_id')}
icu=c.execute('select hadm_id,intime,outtime from icustays').fetchdf()
icu['intime']=pd.to_datetime(icu.intime);icu['outtime']=pd.to_datetime(icu.outtime)
icu_groups={h:(g.intime.to_numpy(),g.outtime.to_numpy()) for h,g in icu.groupby('hadm_id')}
checks={k:True for k in ['same_elapsed_vancomycin_time','no_prior_linezolid','vancomycin_in_previous48h','latest_known_platelet_value','decision_during_icu','no_future_outcome_eligibility']}
for row in r.itertuples():
    t=np.datetime64(row.t0);cs=s.loc[row.index_hadm]
    checks['same_elapsed_vancomycin_time'] &= bool(row.t0-row.vancomycin_start==cs.t0-cs.other_first_time)
    l=drug_times.get((row.hadm_id,'linezolid'),np.array([],dtype='datetime64[ns]'))
    checks['no_prior_linezolid'] &= not np.any(l<=t)
    v=drug_times[(row.hadm_id,'vancomycin')]
    checks['vancomycin_in_previous48h'] &= bool(np.any((v<t)&(v>=t-np.timedelta64(48,'h'))))
    times,stored,values=lab_groups[row.hadm_id]
    valid=(times<t)&(times>=t-np.timedelta64(24,'h'))&(stored<=t)
    if not valid.any():checks['latest_known_platelet_value']=False
    else:
        latest=times[valid].max(); vals=np.unique(values[valid&(times==latest)])
        checks['latest_known_platelet_value'] &= bool(len(vals)==1 and vals[0]==row.known_platelets)
    ins,outs=icu_groups[row.hadm_id]
    checks['decision_during_icu'] &= bool(np.any((ins<=t)&(outs>t)))
    checks['no_future_outcome_eligibility'] &= bool(row.t0<row.dischtime and (pd.isna(row.deathtime) or row.t0<row.deathtime))
out={'pairs_checked':len(r),'implementation':'pandas/numpy from source drug and lab views; no SQL eligibility predicates reused','checks':{k:bool(v) for k,v in checks.items()},
     'limit':'Verifies implementation of temporal rules, not causal exchangeability or true clinical indication.'}
assert all(checks.values()),out
(R/'reports/phase2_independent_verification.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
