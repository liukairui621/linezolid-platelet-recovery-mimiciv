"""Independent aggregate-level checks for the v0.14 revision artifacts."""
from pathlib import Path
import hashlib
import json
import math
import pandas as pd
from scipy.stats import t as student_t

R = Path(__file__).resolve().parents[1]
A = R / "reports/additional_analyses_v0.14"
checks = {}

plan = R / "config/analysis_plan_v0.14.json"
checks["plan_sha256"] = hashlib.sha256(plan.read_bytes()).hexdigest()
checks["plan_matches_freeze"] = checks["plan_sha256"] == json.loads((A / "PLAN_FREEZE_v0.14.json").read_text(encoding="utf-8-sig"))["sha256"]

route = pd.read_csv(A / "PARENTERAL_ROUTE_SENSITIVITY_v0.14.csv").iloc[0]
checks["route_valid_draws"] = int(route.valid_draws) == int(route.requested_draws) == 2000
checks["route_balance_below_0_10"] = float(route.max_abs_SMD) < 0.10
checks["route_interval_finite"] = all(math.isfinite(float(x)) for x in [route.RD, route.lower, route.upper]) and route.lower < route.RD < route.upper

mi = pd.read_csv(A / "MULTIPLE_IMPUTATION_SENSITIVITY_v0.14.csv")
parts = mi[mi.imputation.notna()].copy()
pooled = mi[mi.imputation.isna()].iloc[0]
checks["mi_completed_datasets"] = len(parts) == 20
checks["mi_valid_bootstraps"] = int(pooled.bootstrap_valid) == int(pooled.bootstrap_requested) == 10000
qbar = parts.RD.mean()
ubar = parts.within_variance.mean()
bvar = parts.RD.var(ddof=1)
tvar = ubar + (1 + 1/len(parts)) * bvar
df = (len(parts)-1) * (1 + ubar / ((1 + 1/len(parts)) * bvar))**2 if bvar > 0 else math.inf
crit = student_t.ppf(.975, df) if math.isfinite(df) else 1.959963984540054
calc = (qbar, qbar-crit*math.sqrt(tvar), qbar+crit*math.sqrt(tvar))
reported = (float(pooled.RD), float(pooled.lower), float(pooled.upper))
checks["rubin_recalculation_max_error"] = max(abs(a-b) for a,b in zip(calc,reported))
checks["rubin_recalculation_matches"] = checks["rubin_recalculation_max_error"] < 1e-10
checks["mi_balance_below_0_10"] = float(pooled.max_abs_SMD) < 0.10

traj = pd.read_csv(A / "TREATMENT_EXPOSURE_CONTEXT_v0.14.csv", dtype=str)
small = traj[(traj.arm == "LZD") & traj.metric.isin(["Parenteral initiation", "Enteral or unknown initiation"])]
checks["route_linked_cells_suppressed"] = len(small) == 2 and all("Suppress" in x or x == "<10" for x in small.value)
for name in ["PROPENSITY_WEIGHT_DENSITY_v0.14.csv", "PROPENSITY_WEIGHT_QUANTILES_v0.14.csv", "PARENTERAL_ROUTE_SENSITIVITY_v0.14.csv", "MULTIPLE_IMPUTATION_SENSITIVITY_v0.14.csv", "OBSERVATION_OPPORTUNITY_v0.14.csv", "CLINICAL_TREATMENT_CONTEXT_v0.14.csv", "TREATMENT_EXPOSURE_CONTEXT_v0.14.csv"]:
    cols = [x.lower() for x in pd.read_csv(A / name, nrows=1).columns]
    checks[f"no_identifier_columns_{name}"] = not any(x in cols for x in ["subject_id", "hadm_id", "charttime", "t0", "dischtime", "deathtime"])

# Pandas and SciPy comparisons can return numpy scalar booleans. Normalize them
# before both the blocking decision and JSON serialization.
checks = {k: (v.item() if hasattr(v, "item") else v) for k, v in checks.items()}
blocking = [k for k,v in checks.items() if isinstance(v,bool) and not v]
report = {
    "status": "REVIEWED" if not blocking else "BLOCKING",
    "classification": "public aggregate-level independent implementation check",
    "checks": checks,
    "blocking": blocking,
    "scope": "Recalculated Rubin pooling from per-imputation outputs and checked route diagnostics and privacy columns. Manuscript-only checks are recorded in the local submission audit. The frozen v0.11 primary estimate already had independent reproduction and was not re-estimated here.",
}
(A / "QUALITY_GATE_v0.14.json").write_bytes((json.dumps(report, indent=2) + "\n").encode("utf-8"))
print(json.dumps(report, indent=2))
if blocking:
    raise SystemExit(1)
