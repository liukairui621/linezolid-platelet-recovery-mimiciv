"""Build the disclosure-safe v0.14 local submission package and manifests."""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
SUB = ROOT / "submission_v0.14"
FINAL = SUB / "final_upload_polished"
AGG = SUB / "aggregate_data_v0.14_polished_final"
REPO = ROOT / "submission_v0.12/repository"


def copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


FINAL.mkdir(parents=True, exist_ok=True)
AGG.mkdir(parents=True, exist_ok=True)

# Clinical-only aggregate supplement. Historical omics files remain in the
# public repository for provenance but are outside the revised clinical paper.
aggregate_sources = [
    REPO / "config/analysis_plan_v0.10.json",
    REPO / "config/analysis_plan_v0.11.json",
    REPO / "reports/ASSOCIATION_CHECKS_v0.11.json",
    REPO / "reports/CLINICAL_ASSOCIATIONS_v0.11.csv",
    REPO / "reports/CONFIRMATION_GAP_v0.11.csv",
    REPO / "reports/ENDPOINT_CHECKS_v0.11.json",
    REPO / "reports/FIRST_EVENT_COUNTS_v0.11.csv",
    REPO / "reports/FIRST_EVENT_CURVES_v0.11.csv",
    REPO / "reports/INFERENCE_REPRODUCTION_v0.11.json",
    REPO / "reports/INTERVAL_MONTE_CARLO_v0.11.csv",
    REPO / "reports/LAST_DISCHARGE_COUNT_BINS_v0.11.csv",
    REPO / "reports/LAST_DISCHARGE_LAGS_v0.11.csv",
    REPO / "reports/OUTCOME_SOURCE_AUDIT_v0.11.json",
    REPO / "reports/READINESS_BALANCE_v0.10.csv",
    REPO / "reports/RISK_TABLE_v0.11.csv",
    REPO / "reports/R_sessionInfo_v0.11.txt",
    REPO / "reports/submission_v0.12/DISCLOSURE_AUDIT_v0.12.json",
    REPO / "reports/submission_v0.12/FLOW_v0.12.csv",
    REPO / "reports/submission_v0.12/READINESS_SELECTION_v0.12.csv",
    REPO / "reports/submission_v0.12/TABLE1_BASELINE_v0.12.csv",
    REPO / "reports/submission_v0.12/TABLE_S1_MISSINGNESS_v0.12.csv",
    REPO / "reports/submission_v0.12/TABLE_S2_MODEL_BALANCE_v0.12.csv",
]
aggregate_sources += sorted((ROOT / "analysis_v0.14/output").glob("*.csv"))
aggregate_sources += [
    ROOT / "analysis_v0.14/reports/SOURCE_AND_EXPOSURE_AUDIT_v0.14.json",
    ROOT / "analysis_v0.14/reports/REVIEW_SENSITIVITY_RESULTS_v0.14.json",
    ROOT / "analysis_v0.14/reports/R_sessionInfo_v0.14.txt",
]
for src in aggregate_sources:
    if not src.exists():
        raise FileNotFoundError(src)
    if "analysis_v0.14" in src.parts:
        rel = Path("additional_analyses_v0.14") / src.relative_to(ROOT / "analysis_v0.14")
    else:
        rel = src.relative_to(REPO)
    dst = AGG / rel
    copy(src, dst)

(AGG / "README.txt").write_text(
    "Additional file 2: disclosure-safe aggregate clinical data\n\n"
    "This archive contains the frozen primary clinical aggregate outputs and the "
    "post-result v0.14 route, missing-data, propensity and "
    "observation-process summaries. No patient-level records, dates, identifiers, "
    "propensity scores, weights, imputed rows or bootstrap memberships are included.\n\n"
    "The v0.14 additions were planned after the primary result was known. They are "
    "exploratory sensitivity evidence and did not change the frozen primary estimand.\n"
    "Historical transcriptomic analyses are excluded because they are outside the "
    "revised clinical manuscript. Source code and complete provenance are available at "
    "https://github.com/liukairui621/linezolid-platelet-recovery-mimiciv.\n",
    encoding="utf-8",
)

agg_zip = SUB / "Additional_File_2_Aggregate_Data_v0.14.zip"
with zipfile.ZipFile(agg_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for p in sorted(AGG.rglob("*")):
        if p.is_file():
            z.write(p, p.relative_to(AGG).as_posix())

start_here = SUB / "START_HERE_v0.14.md"
start_here.write_text(
    """# BMC Pharmacology and Toxicology submission package v0.14

The frozen primary clinical result is unchanged. This revision adds sensitivity analyses of parenteral linezolid initiation and multiple imputation, together with propensity, treatment-trajectory, and observation summaries. The clinical manuscript focuses on the comparative recovery result and its clinical meaning.

Use `final_upload_polished/` for portal files. `BMCPT_Submission_Package_v0.14.zip` contains the same set. Tables in the Word files are editable three-line tables. Main and supplementary figures are separate final-size vector PDFs.

## Author facts still required before submission

1. Confirm the authorized MIMIC investigator and the applicable local ethics approval, exemption, or waiver, including committee and reference number when applicable.
2. Provide the actual CRediT contributions of all four authors and confirm that all authors approved the final manuscript.
3. Confirm originality and no simultaneous submission before replacing the bracketed paragraph in the cover letter.

These facts were not inferred or fabricated. Search for `AUTHOR ACTION REQUIRED` in the manuscript and `BEFORE SENDING` in the cover letter.
""",
    encoding="utf-8",
)

revision_memo = SUB / "INTERNAL_EDIT_LOG_v0.14.md"
revision_memo.write_text(
    """# Internal edit log v0.14

## Changes accepted and implemented

- Added the prior randomized linezolid-versus-vancomycin febrile-neutropenia trial and narrowed the novelty claim to this population, time zero, and documented recovery endpoint.
- Added a route-restricted sensitivity analysis, exposure trajectory summaries, complete propensity/weight diagnostics, effective sample sizes, and an alternative multiple-imputation analysis.
- Kept the primary v0.11 analysis frozen and labeled all new estimates as post-result exploratory extensions.
- Removed transcriptomic sections and related supplementary material from the clinical submission because heterogeneous public expression cohorts could not validate a treatment mechanism.
- Rebuilt Figure 3, added supplementary diagnostic Figure S1, moved observation-context Figure 4 to Figure S2, and checked final-size typography and panel spacing.
- Rebuilt every manuscript and supplementary table as an editable three-line table and rendered all six Word documents for visual inspection.

## Reporting boundaries retained

- The estimates remain observational associations. Route restriction and imputation do not solve unmeasured indication confounding.
- Exposure trajectory and platelet-testing summaries remain descriptive; they do not estimate sustained treatment or an observation-process causal effect.
- Ethics/access facts, author contributions, final author approval, and exclusivity were not available and remain clearly marked for author completion.
""",
    encoding="utf-8",
)

qa_root = SUB / "qa/rendered_polished"
page_counts = {}
for d in sorted(qa_root.iterdir()) if qa_root.exists() else []:
    if d.is_dir():
        page_counts[d.name] = len(list(d.glob("page-*.png")))
quality = {
    "status": "PASS_WITH_AUTHOR_FACTS_PENDING",
    "statistics": {
        "primary_unchanged": True,
        "route_draws": "2000/2000 valid",
        "mi_draws": "10000/10000 valid across 20 completed datasets",
        "rubin_max_recalculation_error": 9.992007221626409e-16,
        "max_abs_smd_route": 0.0443076292851,
        "max_abs_smd_mi": 0.04955233357525,
    },
    "privacy": "No identifiers or patient-level records in public outputs; small linked route cells suppressed.",
    "documents": {"rendered_page_counts": page_counts, "three_line_tables": True, "visual_review": "PASS"},
    "figures": {"font": "Arial", "width_mm": 170, "visual_review": "PASS", "main": 3, "supplementary": 2},
    "references": {"numbered": 19, "jaksic_pmid": "16447103", "jaksic_doi": "10.1086/500139"},
    "pending": ["ethics/access facts", "actual CRediT roles and all-author approval", "originality/exclusivity confirmation"],
}
(SUB / "QUALITY_CHECKS_v0.14.json").write_text(json.dumps(quality, indent=2), encoding="utf-8")

upload_sources = [
    SUB / "Main_Manuscript.docx",
    SUB / "Title_Page.docx",
    SUB / "Cover_Letter.docx",
    SUB / "Additional_File_1_Supplement.docx",
    agg_zip,
    SUB / "STROBE_RECORD_Checklist.docx",
    SUB / "Submission_Copy_Paste_Sheet.docx",
    SUB / "Submission_Copy_Paste_Sheet.txt",
    start_here,
    SUB / "QUALITY_CHECKS_v0.14.json",
    SUB / "figures/Figure_1_flow.pdf",
    SUB / "figures/Figure_2_first_events.pdf",
    SUB / "figures/Figure_3_recovery.pdf",
    SUB / "figures/Figure_S1_propensity_diagnostics.pdf",
    SUB / "figures/Figure_S2_observation_context.pdf",
]
for src in upload_sources:
    copy(src, FINAL / src.name)

manifest = {
    "version": "0.14",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "target": "BMC Pharmacology and Toxicology",
    "status": "TECHNICALLY READY; AUTHOR FACTS PENDING",
    "files": [
        {"name": p.name, "bytes": p.stat().st_size, "sha256": sha(p)}
        for p in sorted(FINAL.iterdir()) if p.is_file()
    ],
    "remaining_author_actions": [
        "Confirm authorized MIMIC investigator and applicable local ethics approval, exemption, or waiver.",
        "Provide actual CRediT contributions and confirm final approval by all authors.",
        "Confirm originality and no simultaneous submission before using the cover-letter declaration.",
    ],
}
(FINAL / "UPLOAD_MANIFEST_v0.14.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

package_zip = SUB / "BMCPT_Submission_Package_v0.14.zip"
with zipfile.ZipFile(package_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for p in sorted(FINAL.iterdir()):
        if p.is_file():
            z.write(p, p.name)

print(json.dumps({
    "status": "COMPLETE",
    "aggregate_zip": str(agg_zip),
    "submission_zip": str(package_zip),
    "upload_files": len(list(FINAL.iterdir())),
    "submission_sha256": sha(package_zip),
}, indent=2))
