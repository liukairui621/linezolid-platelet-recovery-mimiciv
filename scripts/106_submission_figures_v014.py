"""Build revised main and supplementary figures from aggregate outputs only."""
from pathlib import Path
import json
import shutil
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

R = Path(__file__).resolve().parents[1]
O = R / "submission_v0.14/figures"
O.mkdir(parents=True, exist_ok=True)
old = R / "submission_v0.12/figures"
for stem in ["Figure_1_flow", "Figure_2_first_events"]:
    for ext in ["pdf", "png"]:
        shutil.copy2(old / f"{stem}.{ext}", O / f"{stem}.{ext}")

# Observation-process graphic moves from the main article to supplementary material.
for ext in ["pdf", "png"]:
    shutil.copy2(old / f"Figure_4_observation_context.{ext}", O / f"Figure_S2_observation_context.{ext}")

if Path("C:/Windows/Fonts/arial.ttf").exists():
    font_manager.fontManager.addfont("C:/Windows/Fonts/arial.ttf")
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8.3, "axes.labelsize": 8.5,
    "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": .65,
})
P = "routine_culture_no_recorded_lmwh_shared"
S = "routine_bacterial_culture72h_shared"
a = pd.read_csv(R / "reports/CLINICAL_ASSOCIATIONS_v0.11.csv")
ev = pd.read_csv(R / "reports/FIRST_EVENT_COUNTS_v0.11.csv", dtype={"n": str})
route = pd.read_csv(R / "reports/additional_analyses_v0.14/PARENTERAL_ROUTE_SENSITIVITY_v0.14.csv").iloc[0]
mi = pd.read_csv(R / "reports/additional_analyses_v0.14/MULTIPLE_IMPUTATION_SENSITIVITY_v0.14.csv")
mi = mi[mi.imputation.isna()].iloc[0]

def frozen(pop, metric, label):
    z = a[(a.population == pop) & (a.metric == metric)].iloc[0]
    endpoint, state = metric.split("__")
    def n(drug):
        q = ev[(ev.population == pop) & (ev.endpoint == endpoint) & (ev.state == state) & (ev.drug == drug)]
        return q.iloc[0].n
    return dict(label=label, events=f"{n('linezolid')} / {n('vancomycin')}", rd=100*z.RD, lo=100*z.lower, hi=100*z.upper, primary=(metric == "hospital__recovery" and pop == P))

rows = [
    frozen(P, "hospital__recovery", "Hospital confirmed\n(primary)"),
    frozen(P, "hospital_single__recovery", "Hospital single threshold"),
    frozen(P, "hospital72__recovery", "Confirmation gap 24–72 h"),
    frozen(P, "hospital_timed__recovery", "Transfusion-record screened"),
    frozen(P, "physical__recovery", "Initial physical ICU\nconfirmed"),
    frozen(S, "hospital__recovery", "Full routine-culture cohort"),
    dict(label="Parenteral initiation\npost-result sensitivity", events=f"{route.LZD_events} / {route.VAN_events}", rd=100*route.RD, lo=100*route.lower, hi=100*route.upper, primary=False),
    dict(label="Multiple imputation\npost-result sensitivity", events="14 / 1019", rd=100*mi.RD, lo=100*mi.lower, hi=100*mi.upper, primary=False),
]

W = 170 / 25.4
fig = plt.figure(figsize=(W, 4.9))
gs = fig.add_gridspec(1, 4, width_ratios=[2.5, 1.25, 2.05, 2.55], left=.02, right=.985, bottom=.18, top=.87, wspace=.04)
axs = [fig.add_subplot(gs[0, j]) for j in range(4)]
for ax in axs:
    ax.set_ylim(len(rows)-.45, -.95)
for ax in [axs[0], axs[1], axs[3]]:
    ax.set_axis_off(); ax.set_xlim(0, 1)
for ax, title in zip(axs, ["Outcome / population", "Events\nLZD / VAN", "Risk difference", "RD (95% CI), pp"]):
    ax.set_title(title, loc="left", fontsize=8.2, fontweight="bold", pad=9)
for i, row in enumerate(rows):
    axs[0].text(0, i, row["label"], va="center", fontsize=8)
    axs[1].text(.02, i, row["events"], va="center", fontsize=8)
    color = "#087f8c" if row["primary"] else "#535e69"
    axs[2].errorbar(row["rd"], i, xerr=[[row["rd"]-row["lo"]], [row["hi"]-row["rd"]]], fmt="o", ms=3.7, capsize=2, color=color)
    val = f"{row['rd']:.2f} ({row['lo']:.2f}, {row['hi']:.2f})".replace("-", "−")
    axs[3].text(.02, i, val, va="center", fontsize=8)
axs[2].axvline(0, color="#7b8790", ls=":", lw=.8)
axs[2].set(xlim=(-25, 5), xticks=[-20, -10, 0], yticks=[], xlabel="LZD − VAN (pp)")
axs[2].spines["left"].set_visible(False)
fig.text(.025, .035, "Only the first row is primary. The final two analyses were specified after the primary result was known.\nAll other rows are secondary or sensitivity analyses; negative values indicate less documented recovery with LZD.", fontsize=7.8)
fig.savefig(O / "Figure_3_recovery.pdf", facecolor="white")
fig.savefig(O / "Figure_3_recovery.png", dpi=450, facecolor="white")
plt.close(fig)

meta = {
    "status": "COMPLETE", "version": "0.14", "width_mm": 170,
    "main_figures": ["Figure_1_flow", "Figure_2_first_events", "Figure_3_recovery"],
    "supplementary_figures": ["Figure_S1_propensity_diagnostics", "Figure_S2_observation_context"],
    "font": "Arial", "minimum_planned_font_pt": 7.8,
    "main_figure_3_change": "Added post-result parenteral-route and multiple-imputation sensitivity rows.",
    "observation_context_change": "Moved unchanged v0.12 figure to Supplementary Figure S2.",
}
(O / "FIGURE_LAYOUT_SPEC_v0.14.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
print(json.dumps(meta))
