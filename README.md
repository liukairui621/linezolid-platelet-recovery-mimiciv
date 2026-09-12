# Linezolid initiation and documented platelet recovery in MIMIC-IV

Source code and aggregate research outputs for an exploratory observational comparison in adults with baseline thrombocytopenia. Manuscript preparation release v0.12, 11 September 2026.

## Main result and interpretation

The primary population comprised 65 linezolid and 2492 vancomycin admissions from 2408 unique patients. Unweighted confirmed-recovery events were 14 and 1019. Overlap-weighted 14-day probabilities were 21.22% and 32.13%; the difference was −10.92 percentage points (patient-cluster percentile 95% CI −20.12 to −3.19). This is a selected observational association. It does not establish a causal drug effect, platelet-production impairment, or a treatment recommendation. Sparse exposed events, residual indication differences, recording completeness, death/discharge and monitoring limit interpretation.

The primary outcome requires two post-initiation counts ≥100 ×10⁹/L, ≥24 hours apart without an intervening observed low count. Death and live hospital discharge before confirmation compete with the event. The first-event curves do not represent current platelet-state occupancy. The eight-program transcriptomic analyses were heterogeneous and did not validate a treatment mechanism.

## Data access

This repository contains source code and aggregate statistics, not a patient dataset or a fitted patient prediction model. Individual clinical records, identifiers, shifted dates/times, cohort membership, propensity scores, weights, bootstrap membership/draw objects, caches, local access credentials, and internal reviewer records are not included. Small clinical cells are suppressed in public count tables.

Disclosure correction, 12 September 2026 (UTC): the initial commit `6b0149e15e56382f6bea4484dd74787c90463079` included exact binary baseline proportions and derived SMDs that could reconstruct small aggregate counts despite suppression in the displayed tables. The current distribution suppresses those linked values across overlapping populations, with secondary suppression of categorical families. No clinical outcome estimates changed. Earlier Git history retains the initial aggregate disclosure; this correction does not retract copies already obtained. Use the corrected current files for analysis or redistribution. The filename `READINESS_BALANCE_v0.10.csv` identifies its analysis source; its public copy is now disclosure-redacted, while the original frozen analysis output remains unchanged internally.

In `TABLE1_BASELINE_v0.12.csv`, `abs_SMD_before` and `abs_SMD_after` are absolute values. Model-balance CSV columns `SMD_before` and `SMD_after` are signed. `Suppressed` can denote a linked or complementary value, not necessarily a small cell itself. Reported maximum imbalance uses the complete matrix before disclosure suppression. `reports/submission_v0.12/DISCLOSURE_AUDIT_v0.12.json` records the transformation and its checks. Script 98 requires the original internal aggregate inputs; it is not rerunnable from already redacted copies.

MIMIC-IV 3.1 requires PhysioNet credentialing, training and its data-use agreement: https://physionet.org/content/mimiciv/3.1/ ; dataset DOI https://doi.org/10.13026/kpb9-mt58 . Source data remain governed by their own terms. Public transcriptomic resources include GSE252275, GSE65682, GSE273700 and GSE310202. Prior broader context GSE210797 is retained in historical pathway outputs and is not presented as validation of the final recovery association. GEO accessions can be opened at https://www.ncbi.nlm.nih.gov/geo/ .

## Reproduce the public figures

From this directory, install numpy, pandas and matplotlib in an appropriate Python environment, then run:

```text
python scripts/94_submission_figures_v012.py
```

The script reads only included aggregate tables. Outputs are written to `submission_v0.12/figures/`. Exact typeface reproduction requires Arial; the Windows reference rendering used Arial at a final width of 170 mm. Fonts may substitute on other systems. PDF figures embed fonts; PNGs are 450 dpi.

## Reproduce the controlled clinical analysis

These are preserved research scripts, not a packaged single-command pipeline. They retain the original project root `/root/projects/linezolid_platelet_recovery` and controlled source root `/root/controlled_data/MIMIC_IV/3.1/raw/mimic-iv-3.1`. A credentialed analyst must provision the raw source files, historical extraction products and dependencies, and adapt paths in a separate working copy. Never place recreated patient-level outputs in this public repository.

The historical extraction code begins with scripts 01, 03 and 04. The extended baseline extraction and selection use 71, 73, 74, 76, 77 and 78. The final clinical extension uses 82 (2000 patient-cluster refits per population), 83 (source refresh), 84 (endpoint construction), 85 (estimation), 86 (independent inference checks), and 87 (joint contrast and simulation uncertainty). Scripts 91 and 92 create submission baseline and flow summaries. Earlier scripts provide source definitions and reusable functions; in particular 76 reads selected functions from 52 and 84 loads audited endpoint functions from earlier code. Do not execute every historical script indiscriminately or treat older results as final.

The frozen v0.11 analysis plan SHA-256 is `523f35a63946e3e3b19f60ec88dd7da8ed09a2620843a2d1224cef1b231abd12`. The plan was frozen locally after earlier analyses in overlapping data; it is not independent preregistration. Its historical PLANNED status is preserved. Subsequent execution generated the included v0.11 results. The v0.12 release changes presentation and adds audit summaries; it does not refit or replace those results.

R session information and dependency snapshots are supplied under `reports/` and `environment/`. R 4.3.3 and brglm2 1.1.0 were used for the final clinical analysis. Patient resampling, median filling and propensity refitting are all repeated within draws. Physical-ICU and transfusion-record-screened outcomes are sensitivity estimands. Intervals beyond the primary contrast are exploratory. A near-zero screened-ICU interval endpoint is smaller than its Monte Carlo uncertainty and is not independent positive evidence.

## Transcriptomic context

Scripts 27, 51, 59 and 63 document expression and gene-set processing. Frozen eight-program membership is in `config/pathway_sets_v0.8.json`; full pathway test tables and selected gene-level results are included. These are disease-context and cross-cell-type perturbation analyses, not matched treatment/recovery experiments. Tissue metadata conflicts, threshold differences, globin composition, limited covariate availability and discovery donor uncertainty are retained in reporting.

## Integrity and licensing

Run `python verify_release.py` to verify the file hashes in `SHA256SUMS.json`. The source-code license is MIT; licenses/terms of the original scientific datasets and third-party resources remain unchanged. The repository is a reproducibility resource, not clinical software. No manuscript has been submitted through this repository.
