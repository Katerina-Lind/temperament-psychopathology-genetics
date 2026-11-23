#!/usr/bin/env python3
"""Reproduce figures/tables from shared 23 Nov 2026 LDSC result excerpts.
Replace age CSV with --age-csv and rerun to source the 18 contrasts directly from R output.
"""
import argparse,csv,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

p=argparse.ArgumentParser(); p.add_argument('--age-csv',type=Path); p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]); args=p.parse_args()
root=args.output; fig=root/'figures'; tabs=root/'tables'; fig.mkdir(parents=True,exist_ok=True);tabs.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':12,'axes.labelsize':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
navy='#1D3247'; blue='#2B628D'; amber='#B87636'; gray='#6E7580'; faint='#E7ECF0'; red='#B5484A'
def save(name):
    plt.savefig(fig/(name+'.pdf'),bbox_inches='tight');plt.savefig(fig/(name+'.svg'),bbox_inches='tight');plt.savefig(fig/(name+'.png'),dpi=300,bbox_inches='tight');plt.close()
def write(name,rows):
    with open(tabs/name,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)

traits=['EMO','ACT','SHY','CADHD','AINT','ASUD','ARISK']
labels=['Emotionality','Activity','Shyness','Childhood ADHD','Adult internalising','Adult substance use','Risk-taking']
vals=np.array([[1,-.0165,.1513,-.0256,.0663,.0132,-.1202],[-.0165,1,-.3473,.7805,.1076,.293,.0859],[.1513,-.3473,1,-.4086,-.0368,-.0808,-.2213],[-.0256,.7805,-.4086,1,.4989,.5319,.4105],[.0663,.1076,-.0368,.4989,1,.5374,.2161],[.0132,.293,-.0808,.5319,.5374,1,.3319],[-.1202,.0859,-.2213,.4105,.2161,.3319,1]])
se={(0,1):.1018,(0,2):.0825,(0,3):.0916,(0,4):.0459,(0,5):.0651,(0,6):.0498,(1,2):.0755,(1,3):.0883,(1,4):.0449,(1,5):.0516,(1,6):.0432,(2,3):.0682,(2,4):.0339,(2,5):.0432,(2,6):.0392,(3,4):.0463,(3,5):.0597,(3,6):.0474,(4,5):.0310,(4,6):.0269,(5,6):.0335}
long=[]
for i in range(7):
 for j in range(i+1,7):long.append(dict(trait_1=traits[i],trait_2=traits[j],rg=round(float(vals[i,j]),5),SE_joint_LDSC=se[(i,j)]))
write('Table_2_joint_LDSC_genetic_correlations.csv',long)
fig1,ax=plt.subplots(figsize=(9,7.6)); masked=np.ma.masked_array(vals,mask=np.triu(np.ones_like(vals,dtype=bool),1)); im=ax.imshow(masked,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-.85,vcenter=0,vmax=.85)); ax.set_xticks(range(7),labels,rotation=40,ha='right');ax.set_yticks(range(7),labels)
for i in range(7):
 for j in range(i+1):
  v=vals[i,j]; ax.text(j,i,f'{v:.2f}' if i!=j else '1.00',va='center',ha='center',fontsize=9,fontweight='bold' if abs(v)>.5 else 'normal',color='white' if abs(v)>.55 else navy)
ax.set_title('Genetic correlation structure of early temperament and psychopathology',loc='left',pad=20,fontweight='bold',color=navy)
ax.tick_params(length=0);fig1.colorbar(im,ax=ax,fraction=.046,pad=.03,label='Genetic correlation (rg)');fig1.text(.11,.02,'Joint multivariable LDSC • European-ancestry GWAS • point estimates only',fontsize=8,color=gray);fig1.tight_layout(rect=[0,.04,1,1]);save('Figure_1_genetic_correlation_heatmap')

cross_raw=[('Activity','ADHD − substance use',.487472,.100366,1.1922e-6,2.3845e-6),('Activity','ADHD − risk-taking',.694523,.093887,1.3886e-13,5.5543e-13),('Shyness','ADHD − risk-taking',-.187277,.072083,.0093743,.012499),('Emotionality','Internalising − substance use',.053059,.070460,.45143,.45143)]
cross=[]
for t,c,d,s,pv,q in cross_raw:cross.append(dict(temperament=t,contrast=c,delta_rg=d,SE_delta=s,CI_low=d-1.96*s,CI_high=d+1.96*s,P=pv,FDR_q=q))
write('Table_3_cross_phenotype_contrasts.csv',cross)
fig2,ax=plt.subplots(figsize=(9,4.0));sorted_cross=list(reversed(cross))
for i,r in enumerate(sorted_cross):
 sig=r['FDR_q']<.05;col=blue if sig else gray
 ax.plot([r['CI_low'],r['CI_high']],[i,i],lw=2,color=col);ax.scatter([r['delta_rg']],[i],s=52,color=col,zorder=3)
 ax.text(.92,i,f"q={r['FDR_q']:.3g}",va='center',fontsize=9,color=navy)
ax.axvline(0,c=gray,lw=1,ls='--');ax.set_yticks(range(4),[r['temperament']+': '+r['contrast'] for r in sorted_cross]);ax.set_xlim(-.42,1.02);ax.set_xlabel('Difference in genetic correlation (95% CI)');ax.set_title('Cross-phenotype contrasts',loc='left',fontweight='bold',color=navy,pad=16);ax.grid(axis='x',alpha=.13);fig2.text(.13,.01,'Differences are within-temperament; BH FDR over 4 exploratory contrasts',fontsize=8,color=gray);fig2.tight_layout(rect=[0,.045,1,1]);save('Figure_2_cross_phenotype_contrasts')

age_raw=[
('EMO','CADHD',-.116701,.14080,-.257504,.11753,.0284478,.103208,.8674,.005021,1.686e-5),('EMO','AINT',.054476,.10682,-.052345,.05118,.3064353,.551584,.8674,.006836,1.730e-6),('EMO','ASUD',-.035962,.11210,-.148066,.07000,.0344026,.103208,.8674,.006296,7.640e-6),('EMO','ARISK',-.112439,-.05451,-.057928,.06426,.3673749,.601159,.8674,.006810,5.123e-6),('EMO','WOR',.200136,.21269,-.012558,.06232,.8402976,.889727,.8674,.006916,7.591e-6),('EMO','AFF',.099104,.19686,-.097755,.05589,.0803043,.206497,.8674,.006659,6.733e-6),
('ACT','CADHD',.693405,.76257,-.069163,.09084,.4464159,.656733,.9188,.006479,1.689e-5),('ACT','AINT',.044189,.11381,-.069621,.04367,.1108562,.249426,.9188,.006842,1.736e-6),('ACT','ASUD',.169359,.32871,-.159350,.05496,.0037368,.033631,.9188,.005808,8.366e-6),('ACT','ARISK',.101057,.07205,.029005,.04819,.5472774,.656733,.9188,.007054,5.280e-6),('ACT','WOR',.008133,-.02840,.036528,.05190,.4815418,.656733,.9188,.007035,8.352e-6),('ACT','AFF',.128310,.10071,.027597,.04521,.5416177,.656733,.9188,.007060,7.643e-6),
('SHY','CADHD',-.370420,-.39663,.026212,.07411,.7235634,.814009,.9720,.003177,1.625e-5),('SHY','AINT',-.028614,-.02449,-.004124,.03435,.9044595,.9044595,.9720,.003225,1.671e-6),('SHY','ASUD',-.047554,-.15104,.103490,.04790,.0307317,.103208,.9720,.002580,9.095e-6),('SHY','ARISK',-.202207,-.25208,.049876,.04339,.2504089,.500818,.9720,.003066,5.209e-6),('SHY','WOR',.048698,.18329,-.134592,.03770,.0003571,.006428,.9720,.002142,6.909e-6),('SHY','AFF',-.042786,.06709,-.109874,.04087,.0071838,.043103,.9720,.002510,8.456e-6)]
ages=[]
for r in age_raw:
 t,o,y2,y3,d,s,pv,q,rgage,es,ev=r
 ages.append(dict(temperament=t,outcome=o,rg_2Y=y2,rg_3Y=y3,rg_2Y_3Y=rgage,delta_rg_2Y_minus_3Y=d,SE_delta=s,CI_low=d-1.96*s,CI_high=d+1.96*s,P=pv,FDR_q=q,min_eigenvalue_S=es,min_eigenvalue_V=ev))
if args.age_csv:
 source=list(csv.DictReader(open(args.age_csv)))
 assert len(source)==18, f'Expected 18 rows, got {len(source)}'
 rows=[]
 for r in source:
  t=r['temperament'];o=r['outcome'];d=float(r['difference_2Y_minus_3Y']);s=float(r['SE_difference']);
  rows.append(dict(temperament=t,outcome=o,rg_2Y=float(r['rg_2Y']),rg_3Y=float(r['rg_3Y']),rg_2Y_3Y=float(r['rg_2Y_3Y']),delta_rg_2Y_minus_3Y=d,SE_delta=s,CI_low=float(r['CI_lower']),CI_high=float(r['CI_upper']),P=float(r['P']),FDR_q=float(r['FDR_BH']),min_eigenvalue_S=float(r['min_eigenvalue_S']),min_eigenvalue_V=float(r['min_eigenvalue_V'])))
 ages=rows
write('Table_4_age_specific_contrasts.csv',ages)
her=[('EMO','2Y',.0439,.0073,1.0062,.0068),('EMO','3Y',.059,.0094,1.0061,.0065),('ACT','2Y',.0781,.0085,1.0031,.0065),('ACT','3Y',.0888,.01,.9998,.0063),('SHY','2Y',.104,.009,.9946,.0062),('SHY','3Y',.1111,.0107,1.0123,.0064)]
write('Table_5_age_specific_heritability.csv',[dict(temperament=t,age=a,h2=h,SE_h2=se,h2_Z=h/se,LDSC_intercept=i,SE_intercept=ise) for t,a,h,se,i,ise in her])

outcomes=['CADHD','AINT','ASUD','ARISK','WOR','AFF']; outcome_label={'CADHD':'Childhood ADHD','AINT':'Adult internalising','ASUD':'Adult substance use','ARISK':'Risk-taking','WOR':'Worry','AFF':'Depressed affect'}
fig3,axs=plt.subplots(1,3,figsize=(13.0,6.7),sharex=True,sharey=True)
for ti,t in enumerate(['EMO','ACT','SHY']):
 ax=axs[ti];subset={r['outcome']:r for r in ages if r['temperament']==t};
 for idx,o in enumerate(outcomes):
  r=subset[o];y=5-idx
  ax.plot([r['CI_low'],r['CI_high']],[y,y],color=blue if r['FDR_q']<.05 else '#AAB5BF',lw=2.1,zorder=2)
  ax.scatter(r['delta_rg_2Y_minus_3Y'],y,s=44 if r['FDR_q']<.05 else 30,color=blue if r['FDR_q']<.05 else gray,zorder=3)
  if r['FDR_q']<.05:ax.text(.205,y,'*',va='center',color=blue,fontweight='bold')
 ax.axvline(0,color=navy,alpha=.55,lw=.9,ls='--');ax.set_title({'EMO':'Emotionality','ACT':'Activity','SHY':'Shyness'}[t],fontweight='bold',color=navy,pad=10);ax.grid(axis='x',alpha=.12);ax.set_xlim(-.47,.25);ax.set_ylim(-.8,5.8);ax.set_yticks(range(5,-1,-1),[outcome_label[o] for o in outcomes]);ax.tick_params(axis='y',length=0)
fig3.suptitle('Age-specific differences in genetic correlations',x=.08,ha='left',fontsize=14,fontweight='bold',color=navy);fig3.supxlabel('Difference in genetic correlation: second year minus third year (95% CI)',y=.055,fontsize=10);fig3.text(.08,.015,'Filled blue points: BH q < 0.05 across 18 contrasts. All estimates from trivariate joint LDSC.',color=gray,fontsize=8);fig3.tight_layout(rect=[0,.075,1,.93],w_pad=2.3);save('Figure_3_age_specific_difference_forest')

fig4,axs=plt.subplots(1,3,figsize=(13,7),sharex=True,sharey=True)
for ti,t in enumerate(['EMO','ACT','SHY']):
 ax=axs[ti];subset={r['outcome']:r for r in ages if r['temperament']==t}
 for i,o in enumerate(outcomes):
  r=subset[o];y=5-i
  ax.plot([r['rg_2Y'],r['rg_3Y']],[y,y],lw=1.4,c='#C5CCD2',zorder=1)
  ax.scatter(r['rg_2Y'],y,s=40,color=blue,marker='o',zorder=2)
  ax.scatter(r['rg_3Y'],y,s=45,color=amber,marker='s',zorder=3)
 ax.axvline(0,lw=.8,ls='--',color=gray);ax.set_title({'EMO':'Emotionality','ACT':'Activity','SHY':'Shyness'}[t],fontweight='bold',color=navy,pad=10);ax.set_xlim(-.52,.90);ax.set_yticks(range(5,-1,-1),[outcome_label[o] for o in outcomes]);ax.tick_params(axis='y',length=0);ax.grid(axis='x',alpha=.12)
fig4.suptitle('Genetic correlations by temperament assessment period',x=.08,ha='left',fontweight='bold',fontsize=14,color=navy);fig4.supxlabel('Genetic correlation (rg)',y=.055);fig4.legend([plt.Line2D([0],[0],marker='o',color='none',markerfacecolor=blue,markeredgecolor=blue),plt.Line2D([0],[0],marker='s',color='none',markerfacecolor=amber,markeredgecolor=amber)],['Second year','Third year'],loc='lower center',bbox_to_anchor=(.5,-.01),ncol=2,frameon=False);fig4.text(.08,.015,'Point estimates only. Statistical inference about age differences is in Figure 3.',fontsize=8,color=gray);fig4.tight_layout(rect=[0,.09,1,.93],w_pad=2.3);save('Figure_4_age_specific_rg_dumbbell')

print('Saved',len(list(fig.glob('*.pdf'))),'PDF figures,',len(list(tabs.glob('*.csv'))),'CSV tables to',root)
