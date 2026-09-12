"""Fresh pre-index medication and support caches for extended candidates; no recovery outcomes."""
from pathlib import Path
import duckdb,json,os,datetime
os.umask(0o077);R=Path('/root/projects/linezolid_platelet_recovery');B=Path('/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1')
c=duckdb.connect(str(R/'cache/audit_combined.duckdb'));c.execute('SET threads=6');c.execute("SET memory_limit='12GB'")
def cached(name,q):
 p=R/'cache'/f'{name}.parquet'
 if not p.exists():c.execute(f"COPY ({q}) TO '{p}.partial' (FORMAT PARQUET,COMPRESSION ZSTD)");Path(str(p)+'.partial').replace(p)
 c.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{p}')");print(name,'ready',flush=True)
# Includes prior records plus exactly t0 for source verification; covariates later require strict pre-t0 and store<=t0.
cached('phase10_emar',f'''SELECT m.emar_id,m.hadm_id,m.medication,m.event_txt,try_cast(m.charttime AS TIMESTAMP) charttime,try_cast(m.storetime AS TIMESTAMP) storetime
 FROM read_csv('{B/'hosp/emar.csv.gz'}',all_varchar=true) m JOIN phase10_candidates e USING(hadm_id)
 WHERE try_cast(m.charttime AS TIMESTAMP)<=e.t0''')
medrx='heparin|enoxaparin|dalteparin|fondaparinux|argatroban|bivalirudin|sulfamethoxazole|trimethoprim|bactrim|ganciclovir|valpro|divalpro|mycophen|azathiopr|methotrex|cyclophosphamide|cytarabine|daunorubicin|doxorubicin|etoposide|fludarabine|tacrolimus|sirolimus|eptifibatide|tirofiban'
beta=json.loads((R/'config/analysis_plan_v0.9.json').read_text())['clinical']['beta_lactam'];rx=medrx+'|'+beta['emar_regex']+'|linezolid|vancomycin'
cached('phase10_med_detail',f'''SELECT d.emar_id,d.route,try_cast(replace(d.dose_given,',','') AS DOUBLE) dose_given FROM read_csv('{B/'hosp/emar_detail.csv.gz'}',all_varchar=true) d
 SEMI JOIN (SELECT m.emar_id FROM phase10_emar m JOIN phase10_candidates e USING(hadm_id) WHERE m.charttime>=e.t0-INTERVAL 7 DAY AND regexp_matches(lower(m.medication),'{rx}')) x USING(emar_id)''')
ids=beta['icu_itemids']+['225152','221906','221289','222315','221749','221662','225881','225798'];ids=','.join("'"+x+"'" for x in ids)
cached('phase10_inputs',f'''SELECT i.hadm_id,i.itemid,try_cast(i.starttime AS TIMESTAMP) starttime,try_cast(i.endtime AS TIMESTAMP) endtime,try_cast(i.storetime AS TIMESTAMP) storetime,try_cast(i.amount AS DOUBLE) amount,i.statusdescription
 FROM read_csv('{B/'icu/inputevents.csv.gz'}',all_varchar=true) i JOIN phase10_candidates e USING(hadm_id)
 WHERE i.itemid IN ({ids}) AND try_cast(i.starttime AS TIMESTAMP)<=e.t0 AND coalesce(try_cast(i.endtime AS TIMESTAMP),try_cast(i.starttime AS TIMESTAMP))>=e.t0-INTERVAL 7 DAY''')
cached('phase10_procedures',f'''SELECT i.hadm_id,i.itemid,try_cast(i.starttime AS TIMESTAMP) starttime,try_cast(i.endtime AS TIMESTAMP) endtime,try_cast(i.storetime AS TIMESTAMP) storetime,i.statusdescription
 FROM read_csv('{B/'icu/procedureevents.csv.gz'}',all_varchar=true) i JOIN phase10_candidates e USING(hadm_id)
 WHERE i.itemid IN ('225441','225802','225803','225809','225955','225805','225792') AND try_cast(i.starttime AS TIMESTAMP)<e.t0
 AND coalesce(try_cast(i.endtime AS TIMESTAMP),try_cast(i.starttime AS TIMESTAMP))>=e.t0-INTERVAL 24 HOUR''')
cached('phase10_charted_support',f'''SELECT ce.hadm_id,ce.itemid,ce.value,try_cast(ce.valuenum AS DOUBLE) valuenum,try_cast(ce.charttime AS TIMESTAMP) charttime,try_cast(ce.storetime AS TIMESTAMP) storetime
 FROM read_csv('{B/'icu/chartevents.csv.gz'}',all_varchar=true) ce JOIN phase10_candidates e USING(hadm_id)
 WHERE ce.itemid IN ('223849','229314','226732','224144','224154','226499','225806','227290') AND try_cast(ce.charttime AS TIMESTAMP)>=e.t0-INTERVAL 24 HOUR
 AND try_cast(ce.charttime AS TIMESTAMP)<e.t0 AND try_cast(ce.storetime AS TIMESTAMP)<=e.t0 AND ce.value IS NOT NULL''')
# Non-patient result: query scope and source-file metadata only.
report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),no_recovery_outcomes_loaded=True,selection='All3674 baseline candidates. Raw scans do not depend on discharge diagnosis codes.',caches=[dict(name=n,rows=c.execute(f'SELECT count(*) FROM {n}').fetchone()[0]) for n in ['phase10_emar','phase10_med_detail','phase10_inputs','phase10_procedures','phase10_charted_support']])
(R/'reports/EXTENDED_SOURCE_COVERAGE_v0.10.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
