"""Audit the selected baseline cohort and hash controlled inputs; no recovery outcomes."""
from pathlib import Path
import json, hashlib, datetime
import pandas as pd
R=Path(__file__).resolve().parents[1]
d=pd.read_csv(R/'cache/phase10_baseline.csv',dtype={'hadm_id':str,'subject_id':str})
keep=d.routine_bacterial_culture72h_shared & ~d.lmwh_available7d
ct=pd.crosstab(d.loc[keep,'icu_unit'],d.loc[keep,'drug'])
keep &= d.icu_unit.isin(ct.index[(ct>0).all(axis=1)])
s=d[keep].copy()
assert s.groupby('drug').size().to_dict()=={'linezolid':65,'vancomycin':2492}
assert s.hadm_id.nunique()==len(s)
assert (s.physical_verified & ~s.known_vre7d & s.routine_bacterial_culture72h & ~s.lmwh_available7d).all()
assert (s.index_input_verified|s.index_emar_verified).all()
assert s.known_platelets.lt(100).all()
rows=[]
for drug,g in s.groupby('drug'):
 for col in ['index_input_verified','index_emar_verified','any_emar_available_before','routine_blood72h','routine_urine72h','routine_respiratory72h','known_mrsa_routine7d','known_named_routine7d','beta_lactam_available24h','ufh_active_available7d','lmwh_available7d','rrt_documented_prior24h','invasive_vent_documented_prior24h','pressors_documented_prior6h']:
  n=int(g[col].sum());rows.append(dict(drug=drug,variable=col,n='<10' if 0<n<10 else n,denominator=len(g)))
paths=['cache/phase10_baseline.csv','cache/phase10_readiness_routine_culture_no_recorded_lmwh_shared.rds','cache/phase10_readiness_routine_bacterial_culture72h_shared.rds']
sha=lambda p:hashlib.sha256((R/p).read_bytes()).hexdigest()
out=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),passed=True,no_recovery_outcomes_loaded=True,patient_rows_exported=False,membership='Routine culture shared ICU population, no available prior7d LMWH, shared ICU units recomputed',index_actual_administration_verified_all=True,counts=s.groupby('drug').size().to_dict(),patients=s.groupby('drug').subject_id.nunique().to_dict(),context=rows,controlled_input_sha256={p:sha(p) for p in paths})
(R/'reports/SELECTED_BASELINE_AUDIT_v0.10.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='context'},indent=2))
