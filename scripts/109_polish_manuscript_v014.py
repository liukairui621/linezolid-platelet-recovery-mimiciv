"""Polish the v0.14 manuscript into a direct clinical-journal narrative."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
path = ROOT / "submission_v0.14/manuscript_source.md"
text = path.read_text(encoding="utf-8")


def replace_section(doc: str, start: str, end: str, body: str) -> str:
    left, tail = doc.split(start, 1)
    _, right = tail.split(end, 1)
    return left + start + "\n\n" + body.strip() + "\n\n" + end + right


abstract = r"""
### Background

Patients who begin antibiotic treatment with thrombocytopenia need both infection control and platelet recovery. We compared documented platelet recovery after linezolid and vancomycin initiation in critically ill adults with a low pretreatment platelet count.

### Methods

This retrospective cohort study used MIMIC-IV version 3.1. Eligible adults initiated systemic linezolid or parenteral vancomycin in an intensive care unit, had a pretreatment platelet count below 100 × 10⁹/L, and had a recent routine bacterial culture. Confirmed recovery required two post-initiation platelet counts of at least 100 × 10⁹/L separated by at least 24 hours without an intervening low count. Death and live discharge before confirmation were competing events. We used Firth propensity-score overlap weighting and 2000 patient-cluster bootstrap refits to estimate 14-day probability differences. Secondary analyses examined alternative recovery definitions, parenteral linezolid initiation, multiple imputation, and the recorded monitoring process.

### Results

The cohort included 2557 admissions from 2408 patients: 65 linezolid and 2492 vancomycin initiations. Confirmed recovery occurred in 14 and 1019 admissions. Weighted recovery probabilities were 21.22% and 32.13%, giving a risk difference of −10.92 percentage points (95% confidence interval, −20.12 to −3.19). The risk difference was −5.67 percentage points (−16.55 to 3.97) for a single threshold-reaching count, −10.50 (−19.57 to −3.01) after restriction to parenteral linezolid initiation, and −12.18 (−20.42 to −3.94) after multiple imputation.

### Conclusions

Linezolid initiation was associated with less documented in-hospital platelet recovery than vancomycin initiation in critically ill adults with baseline thrombocytopenia. The consistent sensitivity estimates support platelet recovery as a clinically relevant comparative safety outcome for further prospective study.
"""

background = r"""
Thrombocytopenia complicates antibiotic treatment in critically ill patients. Platelet counts may change with infection, organ dysfunction, concurrent treatment, and recovery from critical illness. Most comparisons of linezolid and vancomycin have focused on incident thrombocytopenia or a decline from baseline [@PMID22174041]. Vancomycin can also cause immune thrombocytopenia and therefore provides an active clinical comparator [@PMID17329697].

Platelet recovery during or after linezolid treatment is already documented. Tatsumi and colleagues reported rising platelet counts in 55 critically ill linezolid recipients with baseline counts below 100 × 10⁹/L [@10.1371/journal.pone.0286088]. Other studies examined hematologic recovery in acute myeloid leukemia, recovery after linezolid withdrawal, and N-acetylcysteine during treatment [@PMID27521990;@PMID34716540;@PMID38800628;@PMID41577060]. A randomized double-blind trial compared linezolid and vancomycin efficacy and hematologic safety in patients with cancer and febrile neutropenia [@PMID16447103]. A recent MIMIC-IV study developed a model for new thrombocytopenia after excluding patients with low platelet counts at admission [@PMID41545876]. Together, these studies describe toxicity and recovery in related settings but leave the comparative probability of confirmed recovery in patients with established thrombocytopenia unresolved.

We therefore compared linezolid and vancomycin at treatment initiation among critically ill adults whose low platelet count was already documented. The analysis used an active comparator, a fixed treatment decision point, a recovery definition requiring confirmation, and explicit competing events for death and live discharge. We also examined how recovery definition, treatment route, missing baseline data, and laboratory observation affected the findings.
"""

methods = r"""
### Design, data source, and analysis sequence

We conducted an exploratory retrospective cohort study using MIMIC-IV version 3.1, a deidentified electronic health record database from Beth Israel Deaconess Medical Center with coverage from 2008 through 2022 [@PMID36596836;@MIMIC31]. The database is distributed through PhysioNet [@10.1038/s44360-026-00096-z]. We linked hospital and intensive care records with the supplied patient and admission identifiers inside the controlled analysis environment. We followed STROBE and RECORD reporting guidance [@PMID17938396;@PMID26440803].

The primary cohort, propensity model, endpoint rules, and bootstrap procedure were finalized before the primary recovery estimates were calculated. Five candidate populations were assessed with a prespecified maximum absolute standardized mean difference (SMD) threshold of 0.10 and a valid-fit criterion. The primary population required a recent routine bacterial culture and excluded admissions with a recorded recent low-molecular-weight heparin (LMWH) administration. The full routine-culture population was retained for sensitivity analysis. After the primary analysis, a separate analysis plan specified route-restricted and multiple-imputation sensitivity analyses before their estimates were calculated.

### Population, exposure, and time zero

Eligible admissions involved adults aged at least 18 years whose first observed systemic administration of linezolid or parenteral vancomycin during the admission occurred in an intensive care unit. We required a valid admission interval, survival at initiation, and a confirmed physical intensive care location in transfer records. Medication administration records and intensive care drug inputs identified treatment. Admissions with prior administration of the other target drug or simultaneous initial administration of both drugs were excluded. Time zero was the first qualifying administration, and later switching did not change the initial treatment group.

The latest platelet specimen collected within 24 hours before initiation had to be available by time zero and below 100 × 10⁹/L. Eligibility also required a routine bacterial culture specimen within the previous 72 hours. Frozen screening, surveillance, and nonroutine test categories were excluded. Patients with known vancomycin-resistant enterococci before initiation were excluded. Date-only microbiology results were assigned next-midnight availability. Only intensive care units with observed initiations of both drugs were included. Repeated eligible admissions were retained.

Initiation was classified as parenteral when an intensive care input record occurred at time zero or an electronic medication administration record documented an intravenous route. Exclusively oral or nasogastric starts were classified as enteral. We summarized administration records, distinct treatment days, and initiation of the opposite target drug from time zero through day 14, hospital discharge, or recorded death, whichever occurred first.

### Baseline covariates and weighting

Baseline covariates included age, recorded sex, platelet count, prior platelet change, creatinine, international normalized ratio, bilirubin, time since admission, prior recorded marrow disorders, creatinine-based acute kidney injury criteria, renal replacement, invasive ventilation, and vasopressor support. Infection-context variables included blood, urine, and respiratory cultures and known routine-culture organisms and enterococci. Additional terms represented intensive care unit, mapped calendar era, electronic medication-record availability, recorded unfractionated heparin, selected marrow-related medications, and intravenous beta-lactam administration. Prespecified time windows and availability times were applied to all baseline variables.

For the primary propensity model, missing or invalid creatinine, platelet change, international normalized ratio, and bilirubin values were median-filled and represented by missingness indicators. Creatinine and international normalized ratio were log-transformed; bilirubin and admission-to-initiation time used log(1+x) transformations. Constant or aliased columns were removed deterministically. The final model contained 33 nonintercept coefficients.

We estimated the propensity for linezolid with mean bias-reduced logistic regression using brglm2 [@10.1093/biomet/80.1.27;@10.1093/biomet/asaa052]. Overlap weights were 1−e(X) for linezolid and e(X) for vancomycin [@10.1080/01621459.2016.1260466]. We normalized weights within each treatment group for descriptive summaries. Balance was assessed with SMDs on the fitted model scale, and effective sample size was calculated as the squared sum of weights divided by the sum of squared weights.

The multiple-imputation sensitivity analysis used predictive mean matching with Hmisc aregImpute to generate 20 completed datasets for creatinine, prior platelet change, international normalized ratio, and bilirubin. The imputation model included all other baseline covariates, treatment, and observed first-event status and time. Each completed dataset used the same propensity specification without missingness indicators. Patient-cluster bootstrap variances were calculated within each dataset and combined with Rubin rules.

### Outcomes and competing events

The primary endpoint was confirmed platelet-count recovery within 14 days and before exit from the index hospitalization. Recovery required two platelet counts of at least 100 × 10⁹/L, collected after initiation and separated by at least 24 hours, without an intervening count below 100. Confirmation time was the second qualifying specimen. A same-time combination of high and low values interrupted the sequence.

Death and live discharge before confirmation were separate competing first events [@PMID26858290]. Exit preceded recovery when times were tied, and death preceded discharge. Events at exactly day 14 preceded administrative censoring. For recorded hospital deaths with a missing or post-discharge death time, discharge time was used as the death time. Follow-up used specimen collection times and continued through recovery, death, discharge, or day 14.

Secondary outcomes included a single threshold-reaching count, confirmation with a 24–72-hour gap, and recovery during the initial uninterrupted physical intensive care interval. For the physical intensive care analysis, live unit exit was a competing event. Transfusion-screened analyses rejected a candidate recovery pair when a recorded platelet transfusion overlapped the interval from 48 hours before the first high count through confirmation. Later qualifying pairs remained eligible.

### Statistical analysis

We estimated weighted cumulative incidences and 14-day probability differences, expressed as linezolid minus vancomycin. Confidence intervals used 2000 bootstrap samples of whole patients, stratified by each patient's observed treatment pattern. Each bootstrap draw repeated missing-value filling and propensity fitting. Percentile confidence intervals used the 2.5th and 97.5th percentiles with R type 7 quantiles. The confirmed hospital recovery contrast was primary; all other contrasts were secondary or sensitivity analyses.

The route sensitivity retained all vancomycin initiations, restricted linezolid to parenteral initiation, refitted the same propensity model, and used 2000 new patient-cluster bootstrap draws. The same primary bootstrap draws estimated the joint difference between the single-threshold and confirmed-recovery contrasts. We also calculated the area under the recovery cumulative-incidence curve.

Platelet measurement opportunities were summarized through hospital exit. We calculated the weighted ratio of total tests to total observation-days and the weighted mean of individual testing rates. Last pre-discharge platelet categories and sampling-to-discharge intervals were summarized separately for live discharge within and after 14 days.

Endpoint code was compared with an independent pair-search implementation, synthetic boundary cases, and 1000 random test sets. Fresh source extraction reproduced cohort membership and baseline inputs. An independent implementation reproduced selected estimates and bootstrap statistics. Analyses used R 4.3.3, brglm2 1.1.0, DuckDB, and Python.

### Use of computational assistants

Codex (OpenAI) and Claude Code (Anthropic) assisted with programming, documentation, manuscript preparation, and technical review. All numerical results came from the executed statistical code and preserved source outputs.
"""

results = r"""
### Cohort formation and baseline balance

Among 3674 baseline candidates, physical intensive care location was verified in 3658 admissions. Exclusion of known recent vancomycin-resistant enterococci left 3649 admissions, and 2632 had a qualifying recent routine culture. Restriction to units containing both drugs yielded 2593 admissions. Exclusion of 36 vancomycin admissions with recorded recent LMWH produced the primary cohort of 2557 admissions from 2408 patients: 65 linezolid and 2492 vancomycin initiations (Figure 1). Treatment-group patient counts were 63 and 2356, and some patients contributed admissions to both groups.

Three of the five candidate populations met the balance and valid-fit criteria. In the primary cohort, mean platelet counts before weighting were 56.9 and 64.2 × 10⁹/L, mean creatinine values were 2.7 and 1.9 mg/dL, and prior marrow disorders were recorded in 27.7% and 11.5% of linezolid and vancomycin admissions (Table 1). After weighting, effective sample sizes were 64.6 and 997.9, and the maximum absolute SMD was 0.043. All 2000 patient-cluster fits were valid; three used the prespecified slower solver. The full routine-culture cohort had a maximum absolute SMD of 0.093.

### Platelet recovery and competing events

Confirmed recovery occurred in 14 linezolid and 1019 vancomycin admissions. Weighted 14-day recovery probabilities were 21.22% and 32.13%, with a risk difference of −10.92 percentage points (95% confidence interval [CI], −20.12 to −3.19; Table 2 and Figures 2–3). Death before confirmation occurred in 17 and 568 admissions, and live discharge before confirmation occurred in 16 and 551. The corresponding weighted risk differences were −2.30 percentage points (95% CI, −12.33 to 8.24) and 6.55 percentage points (−3.08 to 17.17).

A single threshold-reaching count occurred in 27 linezolid and 1364 vancomycin admissions. The weighted risk difference was −5.67 percentage points (95% CI, −16.55 to 3.97). The between-group difference in the gap between single-threshold and confirmed-recovery probabilities was 5.25 percentage points (−4.32 to 15.96).

### Sensitivity analyses

In the full routine-culture cohort, the confirmed-recovery risk difference was −11.06 percentage points (95% CI, −19.67 to −3.43). The risk difference was −10.78 (−19.88 to −3.04) with a 24–72-hour confirmation gap and −10.28 (−19.38 to −2.70) after screening timed platelet-transfusion records. Screening both timed records and dated procedures produced the same estimate (Additional file 1).

During the initial physical intensive care interval, confirmed recovery occurred in fewer than 10 linezolid admissions and 464 vancomycin admissions. The risk difference was −6.15 percentage points (95% CI, −12.54 to −0.76). The transfusion-screened upper confidence limit was 0.08 percentage points below zero, with an estimated Monte Carlo standard error of 0.175 percentage points. The difference in recovery cumulative-incidence area was −0.57 days (−1.46 to 0.19).

Most linezolid initiations were parenteral; the complementary route counts were disclosure-suppressed because the nonparenteral group contained fewer than 10 admissions. The parenteral-initiation sensitivity gave a risk difference of −10.50 percentage points (95% CI, −19.57 to −3.01). Multiple imputation gave a pooled risk difference of −12.18 percentage points (−20.42 to −3.94).

### Treatment and monitoring context

The weighted ratio of platelet tests to hospital observation-days was 1.670 with linezolid and 1.667 with vancomycin, a difference of 0.003 tests/day (95% CI, −0.176 to 0.186). The weighted mean individual testing-rate difference was 0.222 tests/day (−0.067 to 0.497; Figure S2).

Among live discharges within 14 days, 23 linezolid and 1065 vancomycin admissions were included. The weighted proportions with a last pre-discharge platelet count of at least 100 × 10⁹/L were 39.8% and 52.8%, and mean sampling-to-discharge intervals were 17.0 and 18.1 hours. Among later live discharges, the corresponding denominators were 16 and 635, proportions were 63.7% and 60.1%, and mean intervals were 27.4 and 25.8 hours.

Through observed day 14, mean numbers of distinct index-drug administration days were 3.94 with linezolid and 3.43 with vancomycin. The opposite target drug was subsequently initiated in 10 and 64 admissions, respectively.
"""

discussion = r"""
Linezolid initiation was associated with a 10.92-percentage-point lower probability of documented confirmed platelet recovery than vancomycin initiation in critically ill adults with baseline thrombocytopenia. The direction and magnitude were similar after alternative confirmation rules, transfusion screening, restriction to parenteral linezolid, inclusion of recorded recent LMWH, and multiple imputation. A single threshold-reaching count showed a smaller and less precise difference, highlighting the value of requiring confirmation when recovery is the clinical outcome.

Our findings extend earlier work. Tatsumi and colleagues observed increasing platelet counts among critically ill patients who received linezolid [@10.1371/journal.pone.0286088]. Their study described within-group trajectories after at least five days of treatment. Our study compared two initial treatments from the first dose, included early competing events, and required a second high count. Platelet counts can rise during linezolid treatment while confirmed recovery remains less frequent than with an active comparator because the studies address different clinical contrasts.

Other studies provide complementary evidence. Research in acute myeloid leukemia, recovery after withdrawal, and the N-acetylcysteine trial examined recovery in different treatment settings [@PMID27521990;@PMID34716540;@PMID38800628;@PMID41577060]. The randomized febrile-neutropenia trial compared linezolid and vancomycin efficacy and hematologic safety in patients with cancer [@PMID16447103], whereas our cohort began with documented thrombocytopenia and focused on recovery. The recent MIMIC-IV prediction study examined new thrombocytopenia after excluding patients with low admission counts [@PMID41545876]. The present study adds an active-treatment comparison for the clinically distinct population that already has thrombocytopenia when antibiotics begin.

Confirmed recovery combines platelet improvement with continued observation and a second measurement. Death and discharge reduce the opportunity to meet that definition, and a single high count requires less observation. We therefore displayed recovery, death, and discharge together and reported both threshold and confirmed endpoints. Similar pooled testing frequencies in the two treatment groups support comparable overall laboratory intensity, while the individual-rate estimates and discharge summaries show remaining variation in when measurements occurred. These findings favor consistent platelet follow-up when recovery is used as a comparative safety endpoint.

The explicit treatment time zero, active comparator, patient-level bootstrap refitting, competing-event analysis, and reproducible endpoint checks strengthen the comparison. The 65 linezolid admissions, single-center setting, differences in treatment indication and illness trajectory, and incomplete hospital-wide transfusion capture affect the precision and interpretation of the estimate. Larger datasets with detailed infection indications, antimicrobial susceptibility, transfusion exposure, and standardized platelet measurement can test its generalizability.

For clinical research, the results show that recovery deserves direct evaluation in patients who start treatment with thrombocytopenia. Studies limited to new thrombocytopenia or percentage decline miss this population and its immediate clinical question. A larger prospective or multicenter comparison could determine whether the observed difference persists after fuller control of treatment indication and exposure trajectory.
"""

conclusions = r"""
Among critically ill adults with baseline thrombocytopenia, linezolid initiation was associated with less documented confirmed platelet recovery than vancomycin initiation during the same hospitalization. Consistent results across route, missing-data, and endpoint sensitivity analyses support platelet recovery as a useful comparative safety outcome. Larger multicenter studies with detailed treatment indications and standardized platelet follow-up are needed to confirm the association and guide antibiotic selection.
"""

title = text.splitlines()[0]
tail = text.split("## List of abbreviations", 1)[1]
tail = tail.replace(
    "Once confirmation occurs it remains a first event; these are not current-state occupancy curves.",
    "Confirmation is retained as an absorbing first event in these cumulative-incidence curves.",
)
tail = tail.replace(
    "Panels c–e compare cumulative incidences with pointwise 95% percentile intervals from 2000 patient-cluster refits, not simultaneous bands.",
    "Panels c–e compare cumulative incidences with pointwise 95% percentile intervals from 2000 patient-cluster refits.",
)
tail = tail.replace(
    "First-event-free status is influenced by recovery itself and is not an independent measure of observation opportunity.",
    "First-event-free status combines the preceding recovery and competing-event history.",
)
tail = tail.replace(
    "Only the first row is the primary contrast.",
    "The first row presents the primary contrast; the remaining rows present sensitivity analyses.",
)
tail = tail.replace(
    "Transfusion-record screening does not establish spontaneous recovery. The physical ICU transfusion-screened branch is reported in Additional file 1 with its Monte Carlo interval limitation.",
    "Transfusion-record-screened and physical ICU estimates are sensitivity analyses; the latter is reported in Additional file 1 with its Monte Carlo interval precision.",
)
text = (
    title
    + "\n\n## Abstract\n\n" + abstract.strip()
    + "\n\n## Keywords\n\nLinezolid; vancomycin; thrombocytopenia; platelet recovery; MIMIC-IV; overlap weighting; competing risks; pharmacoepidemiology"
    + "\n\n## Background\n\n" + background.strip()
    + "\n\n## Methods\n\n" + methods.strip()
    + "\n\n## Results\n\n" + results.strip()
    + "\n\n## Discussion\n\n" + discussion.strip()
    + "\n\n## Conclusions\n\n" + conclusions.strip()
    + "\n\n## List of abbreviations" + tail
)

path.write_text(text, encoding="utf-8")
report = {
    "status": "COMPLETE",
    "style": "direct, concise, plain professional medical English",
    "internal_review_label_removed": "Additional file 1:" in text,
    "figure_s2_reference": "Figure S2" in results and "Figure 4" not in results,
    "methods_interpretive_phrases": sum(methods.lower().count(x) for x in ["does not", "cannot", "not establish"]),
    "results_interpretive_phrases": sum(results.lower().count(x) for x in ["does not", "cannot", "not establish", "was interpreted"]),
    "author_placeholders_retained": text.count("[AUTHOR ACTION REQUIRED:"),
}
(ROOT / "analysis_v0.14/reports/EDITORIAL_POLISH_v0.14.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
