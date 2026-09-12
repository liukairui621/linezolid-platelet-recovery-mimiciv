from pathlib import Path
import json,csv,datetime
R=Path(r'G:\Research\MIMIC_IV\projects\linezolid_platelet_recovery');j=json.loads((R/'reports/BASELINE_EXTENSION_INPUTS_v0.10.json').read_text())
routine={'Blood Culture, Routine','URINE CULTURE','REFLEX URINE CULTURE','RESPIRATORY CULTURE','ANAEROBIC CULTURE','FLUID CULTURE','Fluid Culture in Bottles','WOUND CULTURE','FECAL CULTURE','Sonication culture, prosthetic joint','TISSUE CULTURE-TISSUE','TISSUE CULTURE-FLUID'}
screen={'MRSA SCREEN','Staph aureus Screen','Staph aureus Preop PCR','R/O VANCOMYCIN RESISTANT ENTEROCOCCUS'}
uncertain={'TISSUE','AEROBIC BOTTLE','Tissue Culture-Bone Marrow','Stem Cell Aer/Ana Culture'}
rows=[]
for x in j['test_name_inventory']:
 n=x['test_name'];u=n.upper()
 cat='routine_bacterial' if n in routine else 'colonization_screen' if n in screen else 'uncertain_or_product' if n in uncertain else 'specialized_culture' if 'CULTURE' in u else 'nonculture_or_other'
 rows.append({'test_name':n,'category':cat,'routine_bacterial':n in routine,'colonization_screen':n in screen,'rationale':'Exact named routine culture; specimen screening still overrides' if n in routine else 'Screening or uncertain/nonroutine tests do not establish primary routine bacterial-culture evidence'})
assert routine.issubset({x['test_name'] for x in rows})
p=R/'config/microbiology_test_taxonomy_v0.10.tsv'
with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
(R/'config/microbiology_taxonomy_provenance_v0.10.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selection_basis':'Test names only, no treatment-specific outcomes or organism findings used. Twelve explicit routine culture names; generic TISSUE and AEROBIC BOTTLE retained as uncertain. Stool-only context audited separately.','freeze_stage':'Before any extended-cohort outcome comparison','routine_test_names':sorted(routine),'screen_test_names':sorted(screen)},indent=2))
print(len(rows),'test names mapped;',len(routine),'routine culture names')
