"""Add a second documented-administration source without replacing eMAR audit outputs."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import os
import duckdb

os.umask(0o077)
R=Path('/root/projects/linezolid_platelet_recovery')
B=Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
C=R/'cache'
c=duckdb.connect(str(C/'audit_combined.duckdb'))
c.execute("SET threads=6")
c.execute("SET memory_limit='12GB'")
for name in ['target_emar','target_emar_detail','target_prescriptions','patients','admissions','icustays','diagnosis_flags','platelet_transfusions_icu']:
 c.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{C/name}.parquet')")

def extract(name,q):
 path=C/f'{name}.parquet'
 if not path.exists():
  temp=str(path)+'.partial'
  c.execute(f"COPY ({q}) TO '{temp}' (FORMAT PARQUET,COMPRESSION ZSTD)")
  Path(temp).replace(path)
 c.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{path}')")
 print(name,c.execute(f'SELECT count(*) FROM {name}').fetchone()[0],flush=True)

extract('target_icu_antibiotics',f"""
 SELECT subject_id,hadm_id,stay_id,itemid,try_cast(starttime AS TIMESTAMP) starttime,
 try_cast(endtime AS TIMESTAMP) endtime,try_cast(amount AS DOUBLE) amount,amountuom,
 orderid,linkorderid,statusdescription,
 CASE WHEN itemid='225881' THEN 'linezolid' ELSE 'vancomycin' END drug
 FROM read_csv('{B/'icu/inputevents.csv.gz'}',header=true,all_varchar=true)
 WHERE itemid IN ('225881','225798')
""")
extract('extra_target_labs',f"""
 SELECT l.subject_id,l.hadm_id,l.labevent_id,l.specimen_id,l.itemid,
        try_cast(l.charttime AS TIMESTAMP) charttime,
        try_cast(l.storetime AS TIMESTAMP) storetime,
        try_cast(l.valuenum AS DOUBLE) valuenum,l.valueuom
 FROM read_csv('{B/'hosp/labevents.csv.gz'}',header=true,all_varchar=true) l
 SEMI JOIN (SELECT DISTINCT hadm_id FROM target_icu_antibiotics
            EXCEPT SELECT DISTINCT hadm_id FROM target_emar) e USING(hadm_id)
 WHERE l.itemid IN ('51265','53189','50912','52546','51006','52647','50885','51237','51301','51222','50813')
""")
c.execute(f"CREATE OR REPLACE VIEW target_labs AS SELECT * FROM read_parquet('{C/'target_labs.parquet'}') UNION ALL SELECT * FROM extra_target_labs")
qc=c.execute('SELECT drug,statusdescription,amountuom,count(*) n FROM target_icu_antibiotics GROUP BY ALL ORDER BY drug,n DESC').fetchdf().to_dict('records')
for row in qc:
 if 0<row['n']<10:row['n']='<10'
(R/'reports/icu_source_qc.json').write_text(json.dumps(qc,indent=2))
print(json.dumps(qc),flush=True)
print('ICU_SOURCE_EXTENSION_COMPLETE',flush=True)
