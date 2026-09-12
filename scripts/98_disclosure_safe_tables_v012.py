"""Create disclosure-safe distribution copies; never modify frozen analysis inputs.

Suppress rare binary cells, their complementary proportions and derived SMDs.
Categorical families receive secondary suppression to prevent subtraction from
the known group total. Suppression is propagated across overlapping populations.
"""
from pathlib import Path
import csv, json, hashlib, math, re
R = Path(__file__).resolve().parents[1]
D = R / 'submission_v0.12/disclosure_safe'
CONTINUOUS = {'age', 'known_platelets', 'log_creatinine',
              'platelet_change_filled', 'log_inr', 'log1p_bilirubin',
              'log1p_admission_days'}
MARK = 'Suppressed'
def read(rel):
    with (R / rel).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))
def write(rel, rows):
    p = D / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def rare(k, n):
    return 0 < k < 10 or 0 < n-k < 10
sel = {r['population']: r for r in read('reports/submission_v0.12/READINESS_SELECTION_v0.12.csv')}
raw = read('reports/READINESS_BALANCE_v0.10.csv')
assert set(r['population'] for r in raw) == set(sel)
blocked = set()
for r in raw:
    if r['term'] in CONTINUOUS: continue
    for arm in ['LZD', 'VAN']:
        n = int(sel[r['population']]['n_'+arm]); p = float(r['mean_'+arm+'_before'])
        assert 0 <= p <= 1, r['term']
        k = round(n*p)
        assert math.isclose(n*p, k, abs_tol=1e-7), r['term']
        if rare(k, n): blocked.add(r['term'])
# Protect categorical totals, including the reference category shown in Table 1.
families = [prefix for prefix in ['icu_unit', 'mapped_era'] if any(t.startswith(prefix) for t in blocked)]
blocked.update(r['term'] for r in raw if any(r['term'].startswith(p) for p in families))
def is_blocked(term):
    return term in blocked or any(term.startswith(p) for p in families)
numeric = ['mean_LZD_before', 'mean_VAN_before', 'mean_LZD_after',
           'mean_VAN_after', 'SMD_before', 'SMD_after']
sources = ['reports/READINESS_BALANCE_v0.10.csv',
           'reports/submission_v0.12/TABLE_S2_MODEL_BALANCE_v0.12.csv']
checks = []
for rel in sources:
    rr = read(rel)
    for r in rr:
        if is_blocked(r['term']):
            for key in numeric: r[key] = MARK
        elif r['term'] not in CONTINUOUS:
            for arm in ['LZD', 'VAN']:
                n = int(sel[r['population']]['n_'+arm])
                assert not rare(round(n*float(r['mean_'+arm+'_before'])), n)
    write(rel, rr)
    checks.append({'file':rel, 'suppressed_rows':sum(is_blocked(r['term']) for r in rr)})
rel = 'reports/submission_v0.12/TABLE1_BASELINE_v0.12.csv'
rr = read(rel)
for r in rr:
    b = is_blocked(r['model_term'])
    if r['model_term'] not in CONTINUOUS:
        for arm in ['LZD', 'VAN']:
            value = r[arm+'_before']
            if 'suppressed' in value.lower() or '<10' in value: b = True
            else:
                m = re.match(r'^(\d+) \(', value)
                assert m, (r['key'], 'invalid count format')
                if rare(int(m[1]), int(sel[r['population']]['n_'+arm])): b = True
    if b:
        for key in ['LZD_before', 'LZD_after', 'VAN_before', 'VAN_after', 'SMD_before', 'SMD_after']:
            r[key] = MARK
    r['abs_SMD_before'] = r.pop('SMD_before')
    r['abs_SMD_after'] = r.pop('SMD_after')
write(rel, rr)
checks.append({'file':rel, 'suppressed_rows':sum(r['abs_SMD_before']==MARK for r in rr)})
audit = {'source_changes':False, 'outcome_changes':False,
         'policy':'Suppress binary counts or complements of 1-9; suppress all numeric cells and SMDs for affected terms across populations; categorical-family secondary suppression.',
         'checks':checks, 'direct_binary_count_inversion_checked':True,
         'all_linked_derived_SMDs_suppressed':True,
         'scope':'Baseline distribution tables. This is a check of direct count/complement/SMD inversion, not a guarantee against every possible inference across all published statistics.',
         'source_sha256':{p:hashlib.sha256((R/p).read_bytes()).hexdigest() for p in sources+[rel]},
         'output_sha256':{p:hashlib.sha256((D/p).read_bytes()).hexdigest() for p in sources+[rel]}}
(D/'DISCLOSURE_AUDIT_v0.12.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps({'distribution_tables':len(checks), 'direct_inversion_check':'passed', 'frozen_inputs_modified':False}))
