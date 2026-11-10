from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
t=read(R/'results/paired_guild_tests_hmm.tsv');d=read(R/'results/turnover_decomposition_hmm.tsv');names=['HMM:glycine betaine transporter','HMM:ahpCF'];labels=['Betaine uptake\n(BetL)','Peroxide defence\n(AhpCF)'];fig,ax=plt.subplots(1,2,figsize=(9,3.4),gridspec_kw={'width_ratios':[1,1.15]})
for i,name in enumerate(names):
 for delta,denom,color,label in [(-.13,'known_catalogue_fraction','#246477','Mapped community'),(.13,'whole_read_bound','#bc6c25','Original-read denominator')]:
  x=t[(t.panel==name)&(t.quality=='all')&(t.cohort=='all')&(t.denominator==denom)&(t.contrast=='Surface-minus-Deep')&t.omitted_trip.isna()].iloc[0]
  ax[0].errorbar(x.mean_pp,i+delta,xerr=[[x.mean_pp-x.ci_low],[x.ci_high-x.mean_pp]],fmt='o',color=color,capsize=3,label=label if i==0 else None)
 bottom=0
 for metric,color,label in [('between_pp','#246477','Between genera'),('within_pp','#8fb4a1','Within genera'),('exclusive_pp','#bc6c25','Compartment-exclusive')]:
  x=d[(d.panel==name)&(d['rank']=='genus')&(d.metric==metric)].iloc[0]['mean'];ax[1].bar(i,x,bottom=bottom if x>=0 else 0,color=color,width=.5,label=label if i==0 else None)
  if x>=0:bottom+=x
ax[0].set_yticks([0,1],labels);ax[0].axvline(0,color='.65',lw=.7);ax[0].set_xlabel('Surface − subsurface (percentage points)');ax[0].set_title('a  Paired functional contrasts',loc='left');ax[0].legend(frameon=False,fontsize=8,loc='lower right');ax[0].set_ylim(1.65,-.65)
ax[1].set_xticks([0,1],['Betaine uptake','Peroxide defence']);ax[1].set_ylabel('Contribution (percentage points)');ax[1].axhline(0,color='.65',lw=.7);ax[1].set_title('b  Genus-level decomposition',loc='left');ax[1].legend(frameon=False,fontsize=8,loc='upper right');ax[1].set_ylim(-.7,10)
fig.tight_layout();fig.savefig(R/'results/functional_integration.pdf');fig.savefig(R/'results/functional_integration.png',dpi=180)
