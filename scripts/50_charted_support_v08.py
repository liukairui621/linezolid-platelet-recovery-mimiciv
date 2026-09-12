"""Source repair with baseline chartevents; no outcomes. Raw data read-only."""
from pathlib import Path
import os,json,duckdb,datetime
os.umask(0o077);R=Path('/root/projects/linezolid_platelet_recovery');B=Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
cfg=json.loads((R/'config/source_addendum_v0.8.json').read_text())
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'));c.execute('SET threads=6');c.execute("SET memory_limit='12GB'")
p=R/'cache/phase8_charted_support.parquet'
if not p.exists():
 c.execute(f'''COPY (SELECT e.hadm_id,ce.itemid,ce.value,try_cast(ce.valuenum AS DOUBLE) valuenum,
 try_cast(ce.charttime AS TIMESTAMP) charttime,try_cast(ce.storetime AS TIMESTAMP) storetime
 FROM read_csv('{B/'icu/chartevents.csv.gz'}',all_varchar=true) ce
 JOIN phase61_eligible_points e ON ce.hadm_id=e.hadm_id
 WHERE ce.itemid IN ('223849','229314','226732','224144','224154','226499','225806','227290')
 AND try_cast(ce.charttime AS TIMESTAMP)>=e.t0-INTERVAL 24 HOUR AND try_cast(ce.charttime AS TIMESTAMP)<e.t0
 AND try_cast(ce.storetime AS TIMESTAMP)<=e.t0 AND ce.value IS NOT NULL)
 TO '{p}.partial' (FORMAT PARQUET,COMPRESSION ZSTD)''');Path(str(p)+'.partial').replace(p)
def quote(v):return ','.join("'"+x.replace("'","''")+"'" for x in v)
vm=quote(cfg['ventilation_modes']['ventilator_mode']);hm=quote(cfg['ventilation_modes']['ventilator_mode_hamilton'])
q=f'''SELECT hadm_id,bool_or((itemid='223849' AND value IN ({vm})) OR (itemid='229314' AND value IN ({hm})) OR (itemid='226732' AND value='Endotracheal tube')) charted_vent24h,
 bool_or((itemid IN ('224144','224154','226499','225806') AND valuenum>0) OR (itemid='227290' AND upper(value) IN ('CVVH','CVVHD','CVVHDF','SCUF'))) charted_rrt24h
 FROM read_parquet('{p}') GROUP BY hadm_id'''
z=c.execute(q).fetchdf();import pandas as pd
d=pd.read_csv(R/'cache/phase8_context.csv',dtype={'hadm_id':str,'subject_id':str});z.hadm_id=z.hadm_id.astype(str)
d=d.merge(z,on='hadm_id',how='left',validate='1:1');d[['charted_vent24h','charted_rrt24h']]=d[['charted_vent24h','charted_rrt24h']].fillna(False).astype(bool)
d['combined_vent24h']=d.vent_available24h|d.charted_vent24h;d['combined_rrt24h']=d.rrt_available24h|d.charted_rrt24h
d.to_csv(R/'cache/phase8_context_with_chart.csv',index=False)
cols=['charted_vent24h','charted_rrt24h','combined_vent24h','combined_rrt24h'];out=[]
for drug,g in d.groupby('drug'):
 for col in cols:
  n=int(g[col].sum());out.append(dict(drug=drug,source=col,n='<10' if 0<n<10 else n))
report=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),no_outcomes_loaded=True,counts=out,limitations='Charted support at baseline improves evidence coverage; absence of an available row is not clinical absence. These are history indicators, not full organ dysfunction scores.')
(R/'reports/CHARTED_SUPPORT_v0.8.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
