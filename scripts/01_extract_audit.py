"""Read-only MIMIC-IV asset extraction; individual records stay on bpa-server."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import os
import platform
import time

import duckdb

os.umask(0o077)
ROOT = Path('/root/projects/linezolid_platelet_recovery')
SOURCE = Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
CACHE = ROOT / 'cache'
started = time.time()
con = duckdb.connect(str(CACHE / 'audit.duckdb'))
con.execute("SET threads=6")
con.execute("SET memory_limit='12GB'")
con.execute(f"SET temp_directory='{CACHE / 'tmp'}'")


def log(message):
    print(dt.datetime.now(dt.timezone.utc).isoformat(), message, flush=True)


def csv(name):
    return f"read_csv('{SOURCE / name}', header=true, all_varchar=true, sample_size=20000)"


def extract(name, query):
    target = CACHE / f'{name}.parquet'
    if not target.exists():
        log(f'EXTRACT_START {name}')
        temp = str(target) + '.partial'
        con.execute(f"COPY ({query}) TO '{temp}' (FORMAT PARQUET, COMPRESSION ZSTD)")
        Path(temp).replace(target)
    con.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{target}')")
    log(f'EXTRACT_READY {name}: {con.execute(f"SELECT count(*) FROM {name}").fetchone()[0]} records')


# No result-based exclusions, treatment effect estimation, or gene screening in this run.
extract('target_emar', f"""
 SELECT subject_id, hadm_id, emar_id, pharmacy_id,
        try_cast(charttime AS TIMESTAMP) charttime, medication, event_txt,
        CASE WHEN regexp_matches(lower(medication),'linezolid|zyvox') THEN 'linezolid'
             ELSE 'vancomycin' END drug
 FROM {csv('hosp/emar.csv.gz')}
 WHERE regexp_matches(lower(medication),'linezolid|zyvox|vanco')
""")

extract('target_emar_detail', f"""
 SELECT d.subject_id, d.emar_id, d.parent_field_ordinal, d.administration_type,
        d.complete_dose_not_given, d.dose_given, d.dose_given_unit,
        d.product_amount_given, d.product_unit, d.product_description,
        d.route, d.infusion_complete, d.new_iv_bag_hung
 FROM {csv('hosp/emar_detail.csv.gz')} d
 SEMI JOIN target_emar e ON d.emar_id=e.emar_id
""")

extract('target_prescriptions', f"""
 SELECT subject_id,hadm_id,pharmacy_id,try_cast(starttime AS TIMESTAMP) starttime,
        try_cast(stoptime AS TIMESTAMP) stoptime,drug,route,dose_val_rx,dose_unit_rx,
        CASE WHEN regexp_matches(lower(drug),'linezolid|zyvox') THEN 'linezolid'
             ELSE 'vancomycin' END drug_class
 FROM {csv('hosp/prescriptions.csv.gz')}
 WHERE regexp_matches(lower(drug),'linezolid|zyvox|vanco')
""")

extract('admissions', f"SELECT * FROM {csv('hosp/admissions.csv.gz')}")
extract('patients', f"SELECT * FROM {csv('hosp/patients.csv.gz')}")
extract('icustays', f"SELECT * FROM {csv('icu/icustays.csv.gz')}")

extract('diagnosis_flags', f"""
 SELECT hadm_id,
 bool_or((icd_version='9' AND (icd_code LIKE '038%' OR icd_code IN ('99591','99592','78552')))
      OR (icd_version='10' AND (icd_code LIKE 'A40%' OR icd_code LIKE 'A41%' OR icd_code IN ('R6520','R6521')))) coded_sepsis,
 bool_or((icd_version='9' AND (substr(icd_code,1,3) IN ('480','481','482','483','484','485','486') OR icd_code='99731'))
      OR (icd_version='10' AND (substr(icd_code,1,3) IN ('J12','J13','J14','J15','J16','J17','J18') OR icd_code='J95851'))) coded_pneumonia,
 bool_or((icd_version='9' AND (substr(icd_code,1,3) IN ('204','205','206','207','208') OR icd_code='2841'))
      OR (icd_version='10' AND (substr(icd_code,1,3) IN ('C91','C92','C93','C94','C95') OR icd_code LIKE 'D61%'))) coded_marrow_disease
 FROM {csv('hosp/diagnoses_icd.csv.gz')}
 GROUP BY hadm_id
""")

# Include both platelet-count item IDs found in the current version's dictionary.
extract('target_labs', f"""
 SELECT l.subject_id,l.hadm_id,l.labevent_id,l.specimen_id,l.itemid,
        try_cast(l.charttime AS TIMESTAMP) charttime,
        try_cast(l.storetime AS TIMESTAMP) storetime,
        try_cast(l.valuenum AS DOUBLE) valuenum,l.valueuom
 FROM {csv('hosp/labevents.csv.gz')} l
 SEMI JOIN (SELECT DISTINCT hadm_id FROM target_emar WHERE hadm_id IS NOT NULL) e
 ON l.hadm_id=e.hadm_id
 WHERE l.itemid IN ('51265','53189','50912','52546','51006','52647','50885','51237','51301','51222','50813')
""")

extract('platelet_transfusions_icu', f"""
 SELECT subject_id,hadm_id,stay_id,itemid,try_cast(starttime AS TIMESTAMP) starttime,
        try_cast(endtime AS TIMESTAMP) endtime,try_cast(amount AS DOUBLE) amount,
        amountuom,statusdescription
 FROM {csv('icu/inputevents.csv.gz')}
 WHERE itemid IN ('225170','226369','227071')
""")

summary = {}
for name, sql in {
 'medication_event_types': 'SELECT drug,event_txt,count(*) n FROM target_emar GROUP BY ALL ORDER BY drug,n DESC',
 'detail_routes': '''SELECT e.drug,d.route,d.administration_type,count(*) n
                    FROM target_emar e JOIN target_emar_detail d USING(emar_id)
                    GROUP BY ALL ORDER BY e.drug,n DESC''',
 'detail_ordinals': 'SELECT parent_field_ordinal,count(*) n FROM target_emar_detail GROUP BY ALL ORDER BY n DESC',
 'lab_units': 'SELECT itemid,valueuom,count(*) n FROM target_labs GROUP BY ALL ORDER BY itemid,n DESC',
}.items():
    data=con.execute(sql).fetchdf().to_dict('records')
    # No patient identifiers; conservative display suppression for small cells.
    for row in data:
        if 0 < row['n'] < 10:
            row['n']='<10'
    summary[name]=data
(ROOT/'reports/extraction_qc.json').write_text(json.dumps(summary,indent=2))
inputs=[]
for p in SOURCE.rglob('*.csv.gz'):
    inputs.append({'path':str(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns})
manifest={'stage':'asset_extraction','completed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
          'elapsed_seconds':round(time.time()-started,2),'python':platform.python_version(),
          'duckdb':duckdb.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'source_version':'MIMIC-IV 3.1','inputs':inputs,
          'restrictions':'Patient-level cache remains on restricted server; diagnosis flags are screening proxies, not validated Sepsis-3 or baseline disease status.'}
(ROOT/'manifests/extraction.json').write_text(json.dumps(manifest,indent=2))
log('EXTRACTION_COMPLETE')
