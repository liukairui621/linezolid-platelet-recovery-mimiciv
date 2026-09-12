"""Final-size figures using frozen aggregate data; no statistical re-estimation."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import font_manager
R=Path(__file__).resolve().parents[1];O=R/'submission_v0.12/figures';O.mkdir(parents=True,exist_ok=True)
P='routine_culture_no_recorded_lmwh_shared';S='routine_bacterial_culture72h_shared';W=170/25.4
if Path('C:/Windows/Fonts/arial.ttf').exists():font_manager.fontManager.addfont('C:/Windows/Fonts/arial.ttf')
plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'lines.linewidth':1.1})
a=pd.read_csv(R/'reports/CLINICAL_ASSOCIATIONS_v0.11.csv');cur=pd.read_csv(R/'reports/FIRST_EVENT_CURVES_v0.11.csv');b=pd.read_csv(R/'reports/LAST_DISCHARGE_COUNT_BINS_v0.11.csv');risk=pd.read_csv(R/'reports/RISK_TABLE_v0.11.csv');ev=pd.read_csv(R/'reports/FIRST_EVENT_COUNTS_v0.11.csv',dtype={'n':str})
COL={'linezolid':'#087f8c','vancomycin':'#b4562a'};states=['recovery','death','discharge','no_first_event'];sc=['#39988e','#c47d83','#d9b360','#e4e9ed'];slabels=['Confirmed recovery','Death before recovery','Live discharge before recovery','No first event']
audit=[]
def save(fig,name):
 fig.canvas.draw()
 for ext in ['pdf','png']:fig.savefig(O/f'{name}.{ext}',dpi=450,facecolor='white')
 audit.append({'name':name,'width_mm':round(fig.get_figwidth()*25.4,2),'height_mm':round(fig.get_figheight()*25.4,2),'font':'Arial','minimum_planned_font_pt':8,'dpi_png':450})
 plt.close(fig)
def tag(ax,letter,title):ax.set_title(f'{letter}  {title}',loc='left',fontweight='bold',pad=8)
def event(pop,metric,drug):
 endpoint,state=metric.split('__');z=ev[(ev.population==pop)&(ev.endpoint==endpoint)&(ev.state==state)&(ev.drug==drug)];assert len(z)==1;return z.iloc[0]['n']
# Figure 1: cumulative eligibility filters, not a future-duration-selected sample.
flow=pd.read_csv(R/'reports/submission_v0.12/FLOW_v0.12.csv');labels=['Baseline platelet count <100\nand observed ICU initiation','Physical ICU location verified','No known recent VRE','Routine bacterial culture\nwithin previous 72 h','ICU units with both drugs','No available recent LMWH record\nPrimary population']
fig,ax=plt.subplots(figsize=(W,6.2));ax.set_axis_off();ax.set(xlim=(0,1),ylim=(0,1))
ys=np.linspace(.87,.12,6)
for i,(stage,label) in enumerate(zip(flow.stage.drop_duplicates(),labels)):
 z=flow[flow.stage==stage].set_index('drug');y=ys[i]
 ax.add_patch(Rectangle((.06,y-.052),.88,.105,facecolor='#f0f5f5' if i==5 else 'white',edgecolor='#50606b',lw=.75))
 ax.text(.08,y,label,va='center',fontsize=8.5)
 for x,drug in [(.70,'linezolid'),(.86,'vancomycin')]:ax.text(x,y,str(z.loc[drug,'n']),ha='center',va='center',fontsize=9,fontweight='bold')
 if i<5:ax.annotate('',xy=(.50,ys[i+1]+.055),xytext=(.50,y-.055),arrowprops={'arrowstyle':'->','lw':.8,'color':'#50606b'})
ax.text(.70,.968,'Linezolid',ha='center',fontweight='bold');ax.text(.86,.968,'Vancomycin',ha='center',fontweight='bold')
ax.text(.08,.025,'Counts are admissions. The full routine-culture sensitivity cohort precedes\nthe LMWH restriction (65 / 2528). No minimum future treatment duration was required.',fontsize=8)
fig.subplots_adjust(left=0,right=1,top=.98,bottom=.01);save(fig,'Figure_1_flow')
# Figure 2: two arm-specific first-event displays plus three pointwise CIF panels.
fig=plt.figure(figsize=(W,6.45));gs=fig.add_gridspec(2,6,left=.095,right=.98,bottom=.085,top=.83,hspace=.86,wspace=1.15,height_ratios=[1,1])
for j,drug in enumerate(COL):
 ax=fig.add_subplot(gs[0,3*j:3*j+3]);z=cur[(cur.population==P)&(cur.drug==drug)].pivot(index='day',columns='state',values='estimate').sort_index()
 ax.stackplot(z.index,*[100*z[k] for k in states],colors=sc,step='post');ax.set(xlim=(0,14),ylim=(0,100),xticks=[0,3,7,14],xlabel='Days from initiation');tag(ax,'a' if j==0 else 'b',drug.capitalize())
 if j==0:ax.set_ylabel('First-event probability (%)')
 rr=risk[(risk.population==P)&(risk.drug==drug)]
 ax.text(0,-.35,'At risk, days 0 / 3 / 7 / 14\n'+' / '.join(rr.n_at_risk.astype(str)),transform=ax.transAxes,fontsize=8,va='top')
for j,st in enumerate(states[:3]):
 ax=fig.add_subplot(gs[1,2*j:2*j+2]);title=['Recovery','Death before\nrecovery','Live discharge\nbefore recovery'][j];tag(ax,'cde'[j],title)
 for drug,color in COL.items():
  z=cur[(cur.population==P)&(cur.drug==drug)&(cur.state==st)].sort_values('day');ax.step(z.day,100*z.estimate,where='post',color=color,ls='-' if drug=='linezolid' else '--');ax.fill_between(z.day,100*z.lower,100*z.upper,color=color,alpha=.14,step='post')
 ax.set(xlim=(0,14),ylim=(0,45),xticks=[0,7,14],xlabel='Days');ax.grid(alpha=.12)
 if j==0:ax.set_ylabel('Cumulative incidence (%)')
fig.legend([Rectangle((0,0),1,1,color=x) for x in sc],slabels,loc='upper center',bbox_to_anchor=(.53,1.005),frameon=False,ncol=2,columnspacing=1.7,handlelength=1.25)
fig.legend([plt.Line2D([0],[0],color=COL[d],ls='-' if d=='linezolid' else '--') for d in COL],['Linezolid','Vancomycin'],loc='upper center',bbox_to_anchor=(.52,.932),frameon=False,ncol=2)
save(fig,'Figure_2_first_events')
# Figure 3: forest, events and numeric interval are separate axes.
metrics=[(P,'hospital__recovery','Hospital confirmed\n(primary)'),(P,'hospital_single__recovery','Hospital single threshold'),(P,'hospital72__recovery','Confirmation gap 24–72 h'),(P,'hospital_timed__recovery','Hospital transfusion-record\nscreened'),(P,'physical__recovery','Initial physical ICU\nconfirmed'),(S,'hospital__recovery','Full routine-culture cohort\nconfirmed')]
fig=plt.figure(figsize=(W,3.75));gs=fig.add_gridspec(1,4,width_ratios=[2.45,1.28,2.0,2.5],left=.015,right=.985,bottom=.20,top=.85,wspace=.04);axs=[fig.add_subplot(gs[0,j]) for j in range(4)]
for ax in axs:ax.set_ylim(len(metrics)-.5,-.95)
for ax in [axs[0],axs[1],axs[3]]:ax.set_axis_off();ax.set_xlim(0,1)
for ax,title in zip(axs,['Outcome / population','Events\nLZD / VAN','Risk difference','RD (95% CI), pp']):ax.set_title(title,loc='left',fontsize=8,fontweight='bold',pad=9)
for i,(pop,metric,label) in enumerate(metrics):
 z=a[(a.population==pop)&(a.metric==metric)].iloc[0];axs[0].text(0,i,label,va='center',fontsize=8);axs[1].text(.02,i,event(pop,metric,'linezolid')+' / '+event(pop,metric,'vancomycin'),va='center',fontsize=8)
 color=COL['linezolid'] if i==0 else '#535e69';axs[2].errorbar(100*z.RD,i,xerr=[[100*(z.RD-z.lower)],[100*(z.upper-z.RD)]],fmt='o',ms=3.6,capsize=2,color=color)
 axs[3].text(.02,i,f'{100*z.RD:.2f} ({100*z.lower:.2f}, {100*z.upper:.2f})',va='center',fontsize=8)
axs[2].axvline(0,color='#7b8790',ls=':',lw=.8);axs[2].set(xlim=(-23,5),xticks=[-20,-10,0],yticks=[],xlabel='LZD − VAN (pp)');axs[2].spines['left'].set_visible(False)
fig.text(.025,.04,'All rows except the primary contrast are secondary or sensitivity analyses.\nLZD / VAN denominators: 65 / 2492; full routine-culture cohort: 65 / 2528.',fontsize=8);save(fig,'Figure_3_recovery')
# Figure 4: all live discharges, retaining post-recovery discharges in denominator.
fig=plt.figure(figsize=(W,5.5));gs=fig.add_gridspec(2,2,left=.13,right=.97,top=.90,bottom=.11,hspace=1.10,wspace=.42,height_ratios=[.88,1.12]);bins=['Below50','From50to99','AtLeast100'];cols=['#be7478','#d9b360','#39988e']
for j,window in enumerate(['Within14d','After14d']):
 ax=fig.add_subplot(gs[0,j]);labs=[]
 for y,drug in enumerate(COL):
  z=b[(b.population==P)&(b.drug==drug)&(b.discharge_window==window)&(b.prior_confirmation=='All')].set_index('bin');left=0
  assert abs(z.fraction.sum()-1)<1e-10
  for bn,cc in zip(bins,cols):v=100*z.loc[bn,'fraction'];ax.barh(y,v,left=left,color=cc,height=.55);left+=v
  labs.append(('LZD' if drug=='linezolid' else 'VAN')+'\nn='+str(z.iloc[0].denominator))
 ax.set(xlim=(0,100),yticks=[0,1],yticklabels=labs,ylim=(1.65,-.65),xticks=[0,50,100],xlabel='Weighted percentage');tag(ax,'ab'[j],['Discharge ≤14 days','Discharge >14 days'][j])
fig.legend([Rectangle((0,0),1,1,color=x) for x in cols],['<50','50–99','≥100 ×10⁹/L'],loc='upper center',bbox_to_anchor=(.58,1.005),ncol=3,frameon=False)
ax=fig.add_subplot(gs[1,:]);ms=['hospital_pooled_tests_per_day','hospital_mean_person_tests_per_day','physical_pooled_tests_per_day','physical_mean_person_tests_per_day'];ll=['Hospital, pooled tests / days','Hospital, mean individual rate','Physical ICU, pooled tests / days','Physical ICU, mean individual rate']
z=a[(a.population==P)&a.metric.isin(ms)].set_index('metric').loc[ms];ax.errorbar(z.RD,np.arange(4),xerr=[z.RD-z.lower,z.upper-z.RD],fmt='o',color=COL['linezolid'],ms=3.5,capsize=2);ax.axvline(0,color='#7b8790',ls=':',lw=.8);ax.set(xlim=(-.4,.6),ylim=(3.6,-.6),yticks=np.arange(4),yticklabels=ll,xlabel='LZD − VAN (platelet tests per observation-day)');tag(ax,'c','Monitoring frequency')
pos=ax.get_position();ax.set_position([.43,pos.y0,.54,pos.height]);save(fig,'Figure_4_observation_context')
(O/'FIGURE_LAYOUT_SPEC_v0.12.json').write_text(json.dumps({'figures':audit,'data_version':'0.11','numerical_refit':False,'caption_location':'Main manuscript; no figure titles embedded beyond panel descriptions'},indent=2),encoding='utf-8')
print('Four final-size composite figures saved as PDF and 450-dpi PNG.')
