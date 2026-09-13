"""Patch the copied v0.12 DOCX builder for the v0.14 clinical-only supplement."""
from pathlib import Path
import re

p = Path(__file__).with_name("105_build_submission_documents_v014.py")
s = p.read_text(encoding="utf-8-sig")
start = "para(d,'S4. Exploratory transcriptomic context','Heading 1')"
end = "d.save(O/'Additional_File_1_Supplement.docx')"
if start not in s or end not in s:
    raise RuntimeError("supplement replacement anchors missing")

replacement = r'''para(d,'S4. Additional sensitivity analyses','Heading 1')
para(d,'These post-primary analyses examined exposure route, missing baseline covariates, propensity-score behavior, treatment trajectories, and platelet observation opportunities.')

pq=rows('analysis_v0.14/output/PROPENSITY_WEIGHT_QUANTILES_v0.14.csv')
pqkeep=[]
for x in pq:
 if round(float(x['quantile']),2) in [.05,.50,.95,.99]:
  pqkeep.append([x['metric'].replace('_',' '),x['arm'],f"{100*float(x['quantile']):.0f}th",f"{float(x['value']):.4f}",f"{float(x['ESS']):.1f}"])
table(d,'Table S11. Propensity-score and overlap-weight distributions',['Measure','Arm','Quantile','Value','ESS'],pqkeep,[2.25,.65,.8,1.25,1.15],'Values are aggregate distribution summaries from the frozen primary propensity model. ESS is repeated within each arm for readability. Individual propensity scores and weights are not redistributed.')

rs=rows('analysis_v0.14/output/PARENTERAL_ROUTE_SENSITIVITY_v0.14.csv')[0]
mir=rows('analysis_v0.14/output/MULTIPLE_IMPUTATION_SENSITIVITY_v0.14.csv')
mp=next(x for x in mir if not x['imputation'])
sens=[
 ['Parenteral linezolid initiation',f"{100*float(rs['RD']):.2f} ({100*float(rs['lower']):.2f}, {100*float(rs['upper']):.2f})",f"{float(rs['ESS_LZD']):.1f} / {float(rs['ESS_VAN']):.1f}",f"{float(rs['max_abs_SMD']):.3f}",f"{rs['valid_draws']} / {rs['requested_draws']}"],
 ['Multiple imputation, pooled',f"{100*float(mp['RD']):.2f} ({100*float(mp['lower']):.2f}, {100*float(mp['upper']):.2f})",f"{float(mp['ESS_LZD']):.1f} / {float(mp['ESS_VAN']):.1f}",f"{float(mp['max_abs_SMD']):.3f}",f"{mp['bootstrap_valid']} / {mp['bootstrap_requested']}"]]
table(d,'Table S12. Route and missing-data sensitivity analyses',['Analysis','RD (95% CI), pp','ESS LZD / VAN','Maximum |SMD|','Valid / requested'],sens,[1.85,1.85,1.15,.9,.75],'The parenteral analysis used 2000 new patient-cluster refits. The multiple-imputation analysis used 20 completed datasets and 500 patient-cluster refits per dataset; variances were pooled with Rubin rules. Both were planned after the primary result was known and are exploratory.')

tr=rows('analysis_v0.14/output/TREATMENT_EXPOSURE_CONTEXT_v0.14.csv')
keys=[]
for x in tr:
 k=(x['metric'],x['unit'])
 if k not in keys:keys.append(k)
td=[]
for metric,unit in keys:
 z=[x for x in tr if x['metric']==metric and x['unit']==unit]
 td.append([metric,unit,next(x['value'] for x in z if x['arm']=='LZD'),next(x['value'] for x in z if x['arm']=='VAN')])
table(d,'Table S13. Observed treatment route and exposure trajectory',['Measure','Summary','LZD','VAN'],td,[3.25,1.0,1.0,1.0],'The observation window ended at day 14, hospital discharge or recorded death, whichever came first. Small and complementary initiation-route cells are suppressed. These are descriptive records and do not estimate a sustained-exposure or switching effect.')

ob=rows('analysis_v0.14/output/OBSERVATION_OPPORTUNITY_v0.14.csv')
om=[]
for metric in dict.fromkeys(x['metric'] for x in ob):
 z=[x for x in ob if x['metric']==metric]
 def fmt(arm):
  x=next(v for v in z if v['arm']==arm)
  if metric=='at_least_two_postinitiation_tests':return f"{100*float(x['weighted_mean']):.1f}%"
  return f"{float(x['weighted_mean']):.2f}; {float(x['weighted_median']):.2f} ({float(x['weighted_q25']):.2f}–{float(x['weighted_q75']):.2f})"
 om.append([metric.replace('_',' '),'Weighted mean; median (IQR)' if metric!='at_least_two_postinitiation_tests' else 'Weighted percentage',fmt('LZD'),fmt('VAN')])
table(d,'Table S14. Platelet observation opportunities',['Measure','Summary','LZD','VAN'],om,[2.6,1.7,1.1,1.1],'Testing was summarized through each hospital observation exit and includes measurements after first recovery. Similar aggregate opportunities do not establish absence of informative measurement timing.')

para(d,'Supplementary figure legends','Heading 1')
para(d,'Figure S1. Propensity-score, overlap-weight and baseline-balance diagnostics. Panels A and B show aggregate kernel-density curves for the frozen primary propensity scores and overlap weights. Panel C shows the 20 largest baseline absolute standardized mean differences among disclosure-safe fitted terms before and after weighting. The dashed line marks 0.10. Individual scores and weights are not displayed.')
para(d,'Figure S2. Platelet counts before live discharge and recorded monitoring. Panels A and B describe all live discharges within and after 14 days, including discharge after confirmed recovery. Panel C presents pooled tests per observation-day and the weighted mean of individual test rates. These summaries describe the recorded platelet observation process.')
d.save(O/'Additional_File_1_Supplement.docx')'''

s = s[: s.index(start)] + replacement + s[s.index(end) + len(end) :]
s = s.replace(
    "('17','Other analyses','Results: Sensitivity; Monitoring; Transcriptomic context','Secondary/exploratory labels; no selective pathway display')",
    "('17','Other analyses','Results: Sensitivity; Monitoring; Supplement S4','Post-result route, missing-data and observation-process analyses labeled exploratory')",
)
s = s.replace(
    "Figure_1_flow.pdf; Figure_2_first_events.pdf; Figure_3_recovery.pdf; Figure_4_observation_context.pdf\\nAdditional_File_1_Supplement.docx",
    "Figure_1_flow.pdf; Figure_2_first_events.pdf; Figure_3_recovery.pdf\\nAdditional_File_1_Supplement.docx; Figure_S1_propensity_diagnostics.pdf; Figure_S2_observation_context.pdf",
)
p.write_text(s, encoding="utf-8")
print(p)
