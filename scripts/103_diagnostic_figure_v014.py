"""Render disclosure-safe propensity, weight, and balance diagnostics."""
from pathlib import Path
import json
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[1]
A=ROOT/'reports/additional_analyses_v0.14'; O=ROOT/'submission_v0.14/figures'; O.mkdir(parents=True,exist_ok=True)
D=pd.read_csv(A/'PROPENSITY_WEIGHT_DENSITY_v0.14.csv')
B=pd.read_csv(ROOT/'reports/submission_v0.12/TABLE_S2_MODEL_BALANCE_v0.12.csv',na_values=['Suppressed'])
P='routine_culture_no_recorded_lmwh_shared'
B=B[B.population.eq(P)].copy()
for c in ['SMD_before','SMD_after']:B[c]=pd.to_numeric(B[c],errors='coerce').abs()
B=B.dropna(subset=['SMD_before','SMD_after']).sort_values('SMD_before',ascending=False).head(20).sort_values('SMD_before')
if Path('C:/Windows/Fonts/arial.ttf').exists(): font_manager.fontManager.addfont('C:/Windows/Fonts/arial.ttf')
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8.5,'axes.titlesize':9,'legend.fontsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.2,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'lines.linewidth':1.2})
COL={'LZD':'#087f8c','VAN':'#b4562a'}
def label(t):
 t=t.replace('known_platelets','Baseline platelet count').replace('log_creatinine','Creatinine, log').replace('platelet_change_filled','Prior platelet change').replace('log1p_admission_days','Admission-to-initiation time, log').replace('prior_recorded_marrow','Recorded marrow disorder').replace('beta_lactam_available24h','IV beta-lactam, 24 h').replace('nonscreen_blood72h','Routine blood culture, 72 h').replace('nonscreen_urine72h','Routine urine culture, 72 h').replace('respiratory72h','Routine respiratory culture, 72 h').replace('known_named_nonscreen7d','Known named organism, 7 d').replace('known_enterococcus7d','Known enterococcus, 7 d').replace('observed_creatinine_aki','Observed creatinine AKI').replace('rrt_documented_prior24h','Renal replacement, 24 h').replace('invasive_vent_documented_prior24h','Invasive ventilation, 24 h').replace('pressors_documented_prior6h','Vasopressor, 6 h').replace('any_emar_available_before','Prior eMAR availability').replace('marrowdrug_available7d','Marrow-related medication, 7 d').replace('ufh_active_available7d','Unfractionated heparin, 7 d').replace('log1p_bilirubin','Bilirubin, log').replace('log_inr','INR, log').replace('genderM','Male sex').replace('age','Age')
 t=t.replace('mapped_era','Calendar era: ').replace('boundary_uncertain','boundary').replace('mapped_early','early').replace('mapped_late','late')
 return t
fig=plt.figure(figsize=(170/25.4,185/25.4))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.55],hspace=.50,wspace=.36,left=.235,right=.97,bottom=.09,top=.95)
for j,(metric,title,xlab) in enumerate([('propensity_score','Propensity-score overlap','Estimated probability of linezolid'),('overlap_weight','Overlap-weight distribution','Unnormalized overlap weight')]):
 ax=fig.add_subplot(gs[0,j])
 for arm in ['LZD','VAN']:
  z=D[(D.metric==metric)&(D.arm==arm)]
  ax.plot(z.x,z.density,color=COL[arm],label=arm)
  ax.fill_between(z.x,0,z.density,color=COL[arm],alpha=.10)
 ax.set_xlim(0,1);ax.set_xlabel(xlab);ax.set_ylabel('Kernel density');ax.grid(axis='y',alpha=.15)
 ax.set_title(f'{chr(65+j)}  {title}',loc='left',fontweight='bold',pad=7)
 ax.legend(frameon=False)
ax=fig.add_subplot(gs[1,:]);y=np.arange(len(B))
ax.scatter(B.SMD_before,y,color='#6b7780',s=20,label='Before weighting',zorder=3)
ax.scatter(B.SMD_after,y,color='#087f8c',s=20,label='After weighting',zorder=3)
for yy,a,b in zip(y,B.SMD_before,B.SMD_after):ax.plot([a,b],[yy,yy],color='#b8c0c5',lw=.7,zorder=1)
ax.axvline(.1,color='#b4562a',ls='--',lw=.9,label='0.10 threshold')
ax.set_yticks(y);ax.set_yticklabels(['\n'.join(textwrap.wrap(label(x), width=27)) for x in B.term]);ax.set_xlabel('Absolute standardized mean difference');ax.grid(axis='x',alpha=.15);ax.set_title('C  Largest disclosure-safe baseline imbalances',loc='left',fontweight='bold',pad=7);ax.legend(frameon=False,ncol=3,loc='lower right')
fig.savefig(O/'Figure_S1_propensity_diagnostics.pdf',facecolor='white')
fig.savefig(O/'Figure_S1_propensity_diagnostics.png',dpi=450,facecolor='white')
plt.close(fig)
meta={'status':'COMPLETE','source_files':['reports/additional_analyses_v0.14/PROPENSITY_WEIGHT_DENSITY_v0.14.csv','reports/submission_v0.12/TABLE_S2_MODEL_BALANCE_v0.12.csv'],'width_mm':170,'height_mm':185,'font':'Arial','minimum_font_pt':7.2,'panels':3,'privacy':'Love plot uses disclosure-safe rows only; density curves do not publish individual propensity scores or weights.'}
(A/'FIGURE_S1_BUILD_v0.14.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print(json.dumps(meta))
