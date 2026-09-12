"""Freeze the next exploratory clinical analysis before endpoint extraction."""
from pathlib import Path
import json,hashlib,datetime
R=Path('/root/projects/linezolid_platelet_recovery')
p=json.loads((R/'config/analysis_plan_v0.5.json').read_text())
p['version']='0.6';p['replaces_for_future_analysis']='0.5; historical artifacts unchanged'
p['amendment_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
p['registration_status']='Internal amendment after baseline QC and disease DE, before clinical recovery counts or effect contrasts. Not external preregistration.'
p['selected_strategy']='Keep v0.5 corrections; test mean-bias-reduced logistic overlap weights and execute first exploratory clinical association analysis using the fixed first-initiation cohort and endpoint.'
p['clinical']['phase6_execution']={
 'claim':'Exploratory standardized association in the currently screened cohort; not a confirmed causal treatment effect or validated shared-indication target trial.',
 'population':'Frozen1576admissions from phase3_baseline_points. No outcome-based exclusion. Discharge-coded sepsis remains a screening proxy. Baseline exposure, ICU and lab availability checks must pass.',
 'preprocessing':'Same raw phase4 baseline inputs and v0.5 covariates. Positive creatinine/INR, nonnegative bilirubin and platelet trend: pooled median fills plus missingness indicators; medians recomputed within each cluster bootstrap. This does not assert MAR or solve missing-data bias.',
 'weight_model':'Same21df calendar_source formula, glm(method=brglm2::brglmFit,type=AS_mean), binomial logit. Weights LZD=1-e,VAN=e. Firth/mean bias reduction addresses infinite estimates, not unmeasured confounding; exact mean balance is not guaranteed.',
 'bootstrap':{'B':500,'seed':20260911,'sampling':'Stratify patient clusters by observed treatment pattern; resample whole patients with all admissions, including cross-arm patients. Recompute fills, refit propensity and both arm summaries in each draw.',
  'uncertainty':'Percentile95% intervals; counts of failed/nonconverged/nonfinite draws explicit. The same stored outcome-free resamples are used for diagnostics and outcome contrasts.',
  'interpretation':'Algorithmic sampling uncertainty conditional on measured-data assumptions; not bias bounds or equivalence tests.'},
 'reporting_gate_before_outcomes':'Permit weighted exploratory association output if original fit converged with finite interior weights, maximum absolute SMD over included terms<=0.10, and >=98% bootstrap refits valid. If not, still build and report unweighted endpoint descriptions; retain failed model report. Never tune the model using outcome differences.',
 'original_model_sensitivity':'v0.5 unpenalized21df model: descriptive point summaries only; no substitution of its separated bootstrap as valid primary intervals.',
 'time_rules':{
  'primary':'Specimen-clock retrospective outcome, preserving v0.5 second-confirming-specimen time. Both samples strictly aftert0, before death/discharge; second sample<=t0+14d. Samples>=100 separated>=24h with no intervening observed low result.',
  'followup_storetime':'Outcome uses retrospectively available lab values indexed by specimen time. Unlike baseline, no claim the clinician knew every outcome value at specimen time. Report confirmation-result storetime lag and late/missing storetime; no future lab count enters weights.',
  'duplicates':'Build a new canonical platelet view from raw target_labs using existing item/unit/range rules. Identical specimen/time/value duplicates use earliest nonmissing storetime. Do not overwrite prior views. Verify no baseline eligibility/value drift.',
  'ambiguous_times':'If simultaneous measurements straddle100, do not confirm and break the high-count run. If all>=100, the threshold status is known even when values differ; use one timepoint. Numeric last-count summaries require a unique value.',
  'ties':'Death or live discharge at exactly confirmation time precedes confirmation; endpoint events exactly14d count before administrative censoring. Death takes precedence over discharge at a tied exit.',
  'death':'In-hospital deathtime when valid; hospital_expire_flag=1 with missing/after-discharge deathtime uses dischtime as explicitly flagged proxy. Invalid discharge or death<=t0 is a construction error, not silent outcome exclusion.',
  'exit':'Death-before-recovery and live-discharge-before-recovery compete. Otherwise follow to14d. No early censoring simply for last blood test; this is a documented outcome, not latent recovery.'},
 'summaries':['Three first-event CIFs plus no-first-event fraction at0..14d','Recovery risk difference LZD-VAN at14d','Integral0..14of recovery CIF, and LZD-VAN area difference','Risk tables at0,3,7,14d','Last pre-live-discharge count, draw-to-discharge lag, missing/ambiguous status, and prior confirmation; split discharge<=14d versus later'],
 'estimator':'Weighted Aalen-Johansen. With only fixed14d administrative censoring and death/live-discharge as competing exits, verify its algebraic equality to weighted first-event empirical CDF; area is weighted mean(14-T)*I(recovery). Patient bootstrap jointly captures all summaries.',
 'sensitivity_in_this_run':'Original24..72h confirmation gap, without reselecting population or weights. Transfusion-restricted and ICU-restricted endpoints remain later named analyses; do not label the present recovery spontaneous or transfusion-free.',
 'test_before_extract':'Synthetic sequential and independent pair-search implementations must agree for threshold equality,24h equality,intervening low,>72h gap, simultaneous conflicting result, death/discharge tie, no test and14d boundary.',
 'no_primary_claim_upgrade':'Any observed drug contrast may reflect indication, severity, missingness, discharge and monitoring differences. Equal counts or nonsignificance cannot establish safety/equivalence.'}
# Resolve inherited active-text conflicts; keep original historical versions on disk.
p['clinical']['weighting']['candidate']='v0.6 mean-bias-reduced logistic overlap weights for exploratory association, subject to explicit outcome-free reporting gate; v0.5 MLE retained as point sensitivity.'
p['clinical']['medication_coverage']['era']='Use mapped decision-year intervals plus baseline source availability; anchor categories remain provenance diagnostics, not mandatory substitutes for admission era.'
p['clinical']['baseline_covariate_priority'][-1]='Age, sex, mapped decision-year interval, baseline eMAR availability, pretreatment platelets and trend, admission-to-treatment time; retain prior recorded marrow disease.'
out=R/'config/analysis_plan_v0.6.json';assert not out.exists(),'Preserve existing freeze; do not overwrite'
out.write_text(json.dumps(p,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'created_utc':p['amendment_utc'],'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'path':str(out)}))
