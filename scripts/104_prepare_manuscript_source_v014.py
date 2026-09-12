"""Create the v0.14 manuscript source from v0.12 and frozen reviewer-response outputs."""
from pathlib import Path
import json
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
A = ROOT / "analysis_v0.14"
OUT = ROOT / "submission_v0.14"
source = (ROOT / "submission_v0.12/manuscript_source.md").read_text(encoding="utf-8")
route = pd.read_csv(A / "output/PARENTERAL_ROUTE_SENSITIVITY_v0.14.csv").iloc[0]
mi = pd.read_csv(A / "output/MULTIPLE_IMPUTATION_SENSITIVITY_v0.14.csv")
mi = mi[mi.imputation.isna()].iloc[0]
traj = pd.read_csv(A / "output/TREATMENT_EXPOSURE_CONTEXT_v0.14.csv", dtype=str)

def pp(x):
    return f"{100 * float(x):.2f}".replace("-", "−")

def num(metric, arm, unit):
    z = traj[(traj.metric == metric) & (traj.arm == arm) & (traj.unit == unit)]
    if len(z) != 1:
        raise RuntimeError((metric, arm, unit, len(z)))
    return z.iloc[0].value

route_text = (
    f"The post-result parenteral-initiation sensitivity produced a risk difference of "
    f"{pp(route.RD)} percentage points (95% CI, {pp(route.lower)} to {pp(route.upper)}). "
    f"Multiple imputation of the four incomplete baseline covariates produced a pooled difference of "
    f"{pp(mi.RD)} percentage points ({pp(mi.lower)} to {pp(mi.upper)})."
)
switch_lzd = num("Opposite target drug initiated through observed day 14", "LZD", "n")
switch_van = num("Opposite target drug initiated through observed day 14", "VAN", "n")
mean_days_lzd = num("Distinct index-drug calendar days through observed day 14", "LZD", "mean")
mean_days_van = num("Distinct index-drug calendar days through observed day 14", "VAN", "mean")

# Abstract.
source = source.replace(
    "Firth propensity-score overlap weighting and 2000 patient-cluster bootstrap refits estimated 14-day probability differences.",
    "Firth propensity-score overlap weighting and 2000 patient-cluster bootstrap refits estimated 14-day probability differences. Post-result sensitivity analyses examined parenteral initiation and multiple imputation of incomplete baseline covariates.",
)
source = source.replace(
    "The primary association was similar in the full routine-culture cohort retaining recorded heparin exposure. Pooled weighted platelet-testing frequencies were 1.670 and 1.667 tests per observation-day.",
    f"The primary association was similar in the full routine-culture cohort retaining recorded heparin exposure. {route_text}",
)
# Keep the structured abstract within the journal's 350-word limit while the
# full Results section retains these secondary contrasts.
source = source.replace(
    "The joint difference between single-threshold and confirmed-recovery contrasts was 5.25 percentage points (−4.32 to 15.96). ",
    "",
    1,
)
source = source.replace(
    "The primary association was similar in the full routine-culture cohort retaining recorded heparin exposure. ",
    "",
    1,
)
source = source.replace(
    "Sparse linezolid events, residual confounding, incomplete exposure capture, and monitoring and discharge processes limit interpretation.",
    "The single-center result was based on sparse linezolid events; residual confounding, incomplete exposure capture, and monitoring and discharge processes limit interpretation.",
)

# Direct clinical precedent.
needle = "Thus, neither recovery as a topic nor linezolid-associated platelet toxicity is a new observation."
source = source.replace(
    needle,
    needle
    + "\n\nRandomized comparisons of linezolid and vancomycin also predate the present study. In a multicenter double-blind trial in patients with cancer and febrile neutropenia, Jaksic and colleagues evaluated clinical efficacy and hematologic safety [@PMID16447103]. That population and treatment context differ from the present analysis of patients with thrombocytopenia already documented before treatment initiation. Our contribution is therefore a population- and endpoint-specific observational comparison of documented recovery, rather than the first comparison of these antibiotics.",
)

# Development status and exposure trajectory definitions.
needle = "Two narrower candidates failed balance criteria and were not promoted to outcome analyses. All candidate diagnostics are disclosed in Additional file 1."
source = source.replace(
    needle,
    needle
    + "\n\nAfter the primary result was known, a separate reviewer-requested extension plan was frozen before calculating the route-restricted and alternative missing-data estimates. These additions are explicitly post-result exploratory analyses. They did not change the primary population, endpoint, weighting model, bootstrap result, or interpretation boundary.",
)
needle = "The design concerns initial observed inpatient treatment choice; it does not establish outpatient new use, adherence, or a sustained regimen."
source = source.replace(
    needle,
    needle
    + " For the reviewer-requested exposure description, initiation was classified as parenteral when an intensive-care input record occurred at time zero or an audited electronic administration route was intravenous; an exclusively oral or nasogastric start was classified as enteral. Administration records, distinct calendar treatment days, and initiation of the opposite target drug were summarized from time zero through the earliest of day 14, hospital discharge, or recorded death. These post-initiation descriptions did not redefine treatment or create an as-treated analysis.",
)
needle = "The primary fit contained 33 nonintercept coefficients."
source = source.replace(
    needle,
    needle
    + " An exploratory missing-data sensitivity used predictive mean matching with Hmisc aregImpute to generate 20 completed datasets for creatinine, prior platelet change, international normalized ratio, and bilirubin. The imputation model included the other baseline covariates, treatment, and observed first-event status and time. Each completed dataset used the same substantive propensity specification without missing-indicator terms. Patient-cluster bootstrap variances were estimated within each dataset and pooled with Rubin rules.",
)
needle = "One primary contrast was designated; other confidence intervals were exploratory and not adjusted for multiple comparisons."
source = source.replace(
    needle,
    needle
    + " We additionally reported propensity-score and overlap-weight distributions, all disclosure-safe fitted-term balance diagnostics, and effective sample sizes. The parenteral-initiation sensitivity retained all vancomycin initiations, excluded nonparenteral or unclassified linezolid starts, refitted the same propensity specification, and used 2000 new patient-cluster bootstrap draws. These analyses were reported regardless of direction.",
)

# Remove transcriptomics from the clinical submission.
source = re.sub(
    r"\n### Exploratory public transcriptomic context\n\n.*?(?=\n### Use of computational assistants)",
    "",
    source,
    flags=re.S,
)
source = re.sub(r"\n### Transcriptomic context\n\n.*?(?=\n## Discussion)", "", source, flags=re.S)

# Results additions.
needle = "The recovery cumulative-incidence area difference was −0.57 days (−1.46 to 0.19), which did not establish an overall difference in the timing of first confirmation."
source = source.replace(
    needle,
    needle
    + f"\n\nMost linezolid initiations were parenteral; the exact complementary route counts were disclosure-suppressed because the nonparenteral group contained fewer than 10 admissions. {route_text} Both analyses were specified after the primary result was known and were treated as exploratory sensitivity evidence rather than independent confirmation.",
)
needle = "These future-selected descriptions differ from the competing event of discharge before recovery and cannot establish the unobserved recovery status at day 14."
source = source.replace(
    needle,
    needle
    + f" Through observed day 14, the mean numbers of distinct index-drug administration days were {mean_days_lzd} and {mean_days_van}; the opposite target drug was subsequently initiated in {switch_lzd} and {switch_van} admissions, respectively. These summaries describe recorded treatment trajectories and do not estimate a sustained-exposure effect.",
)

# Discussion precedent and sensitivity interpretation.
needle = "The contribution is a comparative estimate for an already thrombocytopenic population using an explicit first-dose decision point and documented recovery endpoint."
source = source.replace(
    needle,
    "The contribution is a population- and endpoint-specific comparative estimate for an already thrombocytopenic population using an explicit first-dose decision point and documented recovery endpoint. A prior randomized double-blind comparison in febrile neutropenic patients with cancer already evaluated linezolid and vancomycin efficacy and hematologic safety [@PMID16447103]; the present study does not claim to be the first clinical comparison of these antibiotics.",
)
needle = "This comparison therefore cannot identify the absolute hematologic effect of either drug."
source = source.replace(
    needle,
    needle
    + " The parenteral-initiation and multiple-imputation sensitivities retained the direction of the primary estimate, but neither analysis supplies the unmeasured indication information required for causal exchangeability.",
)
needle = "Median filling with missingness indicators is a transparent modeling choice, not a guarantee that missingness is ignorable."
source = source.replace(
    needle,
    needle
    + " The multiple-imputation sensitivity addressed dependence on that single handling rule, but its assumptions also cannot be verified from the observed records.",
)
needle = "No sustained exposure contrast, dose-response analysis, adjudicated drug-induced thrombocytopenia, or treatment-switch effect was estimated."
source = source.replace(
    needle,
    f"Observed treatment trajectories showed initiation of the opposite target drug in {switch_lzd} linezolid and {switch_van} vancomycin admissions, but no sustained exposure contrast, dose-response analysis, adjudicated drug-induced thrombocytopenia, or treatment-switch effect was estimated.",
)
source = source.replace(
    " Finally, independent public transcriptomic cohorts supplied heterogeneous disease-context results and had unresolved differences in phenotype and sample provenance. They were not treatment or recovery experiments and cannot supply a mechanistic explanation for the clinical association.",
    "",
)
source = source.replace(
    "It should be interpreted as comparative observational evidence rather than a causal drug effect or validated molecular mechanism.",
    "It should be interpreted as comparative observational evidence rather than a causal drug effect or treatment recommendation.",
)

# Data availability and supplementary figure placement.
source = source.replace(" Public transcriptomic accession numbers and provenance are described in Additional file 1.", "")
source = source.replace(
    "We acknowledge the investigators and data contributors who made MIMIC-IV and the public transcriptomic datasets available.",
    "We acknowledge the investigators and data contributors who made MIMIC-IV available.",
)
source = re.sub(
    r"\n### Figure 4\. Platelet counts before live discharge and recorded monitoring\n\n.*?(?=\n## Additional file information)",
    "",
    source,
    flags=re.S,
)
source = source.replace(
    "Additional file 1: Supplementary methods, baseline completeness and balance, all frozen outcome contrasts, observation-process descriptions, and exploratory transcriptomic context. Editable Word document with three-line tables.",
    "Additional file 1: Supplementary methods, baseline completeness and balance, all frozen outcome contrasts, reviewer-requested propensity, route, missing-data and observation-process analyses, and supplementary figure legends. Editable Word document with three-line tables.",
)
source = source.replace(
    "Additional file 2: Machine-readable aggregate clinical and transcriptomic tables, analysis definitions, and source provenance. ZIP archive; no individual clinical records.",
    "Additional file 2: Machine-readable aggregate clinical tables, analysis definitions, diagnostics, and source provenance. ZIP archive; no individual clinical records.",
)

# Verified reference set v0.14.
refdir = ROOT / "literature/submission_v0.14"
refdir.mkdir(parents=True, exist_ok=True)
refs = json.loads((ROOT / "literature/submission_v0.12/verified_references.json").read_text(encoding="utf-8"))
if not any(x["key"] == "PMID16447103" for x in refs["records"]):
    refs["records"].append(
        {
            "key": "PMID16447103",
            "pmid": "16447103",
            "doi": "10.1086/500139",
            "title": "Efficacy and safety of linezolid compared with vancomycin in a randomized, double-blind study of febrile neutropenic patients with cancer.",
            "authors": ["Jaksic B", "Martinelli G", "Perez-Oteyza J", "Hartman CS", "Leonard LB", "Tack KJ"],
            "journal": "Clin Infect Dis",
            "year": "2006",
            "volume": "42",
            "issue": "5",
            "pages": "597-607",
            "publication_types": ["Randomized Controlled Trial", "Multicenter Study", "Journal Article"],
            "source": "https://pubmed.ncbi.nlm.nih.gov/16447103/",
            "verification": "PubMed and Europe PMC metadata archived in analysis_v0.13/literature/lzd_2006_rct.json; wording limited to abstract-supported design and hematologic safety.",
        }
    )
refs["verified_utc"] = "2026-09-12T07:00:00+00:00"
(refdir / "verified_references.json").write_text(json.dumps(refs, ensure_ascii=False, indent=2), encoding="utf-8")
OUT.mkdir(exist_ok=True)
(OUT / "manuscript_source.md").write_text(source, encoding="utf-8")
audit = {
    "status": "COMPLETE",
    "source": "submission_v0.12/manuscript_source.md",
    "output": "submission_v0.14/manuscript_source.md",
    "route_sensitivity_rd_pp": 100 * float(route.RD),
    "mi_sensitivity_rd_pp": 100 * float(mi.RD),
    "jaksic_added": True,
    "transcriptomic_main_sections_removed": "### Transcriptomic context" not in source and "### Exploratory public transcriptomic context" not in source,
    "author_fact_placeholders_retained": source.count("[AUTHOR ACTION REQUIRED:"),
    "claim_boundary_check": not any(x in source.lower() for x in ["validated molecular mechanism", "causal drug effect or validated"]),
}
(A / "reports/MANUSCRIPT_SOURCE_UPDATE_v0.14.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
print(json.dumps(audit, indent=2))
