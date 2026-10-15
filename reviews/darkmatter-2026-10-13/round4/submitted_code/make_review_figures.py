#!/usr/bin/env python3
"""Render current evidence figures from checksum-verified compact tables."""
import sys,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from round2_checks import read
R=Path(__file__).resolve().parents[2];out=R/'manuscript/figures';out.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.labelcolor':'#253444','text.color':'#253444'})
teal='#187f83';grey='#73849b';orange='#d18948'
def frame(p):return pd.DataFrame(read(p))
def save(fig,name):
 fig.savefig(out/(name+'.pdf'),bbox_inches='tight');fig.savefig(out/(name+'.png'),dpi=170,bbox_inches='tight');plt.close(fig)
def tag(ax,s,title):ax.set_title(title,loc='left',pad=12);ax.text(-.16,1.10,s,transform=ax.transAxes,fontweight='bold',fontsize=11)
def ecdf(ax,values,label,color):
 x=np.sort(np.asarray(values,dtype=float));ax.plot(x,np.arange(1,len(x)+1)/len(x),label=label,color=color,lw=1.4)
base='06_darkmatter/darkfam_basis_clean/results/selection/'
d=frame(base+'dnds.tsv');n=frame(base+'dnds_null_frame1.tsv');r=frame(base+'dnds_dark_subset.tsv')
fig,axes=plt.subplots(3,2,figsize=(7.2,8.8),layout='constrained')
ax=axes[0,0]
for arm,label,c in [('dark','Dark (97,200)',teal),('ctrl','Annotated (95,647)',grey)]:ecdf(ax,d.loc[d['set']==arm,'median_dNdS'],label,c)
ax.axvline(1,c='.6',ls=':',lw=.8);ax.set(xlim=(0,2),xlabel='Family median dN/dS',ylabel='Cumulative fraction');ax.legend(frameon=False,fontsize=7);tag(ax,'a','Family dN/dS')
ax=axes[0,1];ecdf(ax,r['median_dNdS'],'Annotated frame (9,182)',teal);ecdf(ax,n['median_dNdS'],'Shifted frame (9,148)',orange);ax.set(xlim=(0,3),xlabel='Family median dN/dS',ylabel='Cumulative fraction');ax.legend(frameon=False,fontsize=7);tag(ax,'b','Frame-specific calibration')
ax=axes[1,0]
stats=json.loads((Path(__file__).parent/'selection_estimands/summary.json').read_text())['strata']
y=np.array([r['difference'] for r in stats]);lo=np.array([r['bootstrap_low'] for r in stats]);hi=np.array([r['bootstrap_high'] for r in stats])
ax.errorbar(range(5),y,yerr=[y-lo,hi-y],fmt='o',c=teal,capsize=3);ax.axhline(0,c='.6',lw=.8)
ax.set(xticks=range(5),xticklabels=['0–0.1','0.1–0.25','0.25–0.5','0.5–1','1–2'],xlabel='Mean synonymous divergence stratum',ylabel='Dark − annotated median dN/dS');ax.tick_params(axis='x',labelsize=6)
for i,r in enumerate(stats):ax.text(i,.35,str(r['dark_n'])+' /\n'+str(r['control_n']),ha='center',fontsize=5.5)
ax.set_ylim(-.07,.43);tag(ax,'c','Stratum effects with 95% intervals')
p=frame('06_darkmatter/cultured_comparison_clean/results/pair_anatomy/pair_anatomy_per_pair.tsv');p=p[p['mag_rep']!='V42PRr2__single_SemiBin_33'];ax=axes[1,1]
for _,row in p.iterrows():ax.plot([0,1],[float(row['mag_annot_with_partner_pct']),float(row['mag_dark_with_partner_pct'])],color=grey,alpha=.55,lw=.7)
ax.scatter(np.zeros(len(p)),p['mag_annot_with_partner_pct'].astype(float),s=12,c=grey);ax.scatter(np.ones(len(p)),p['mag_dark_with_partner_pct'].astype(float),s=12,c=teal)
ax.set(xticks=[0,1],xticklabels=['Annotated','Dark'],ylabel='MAG proteins with reciprocal hit (%)',ylim=(-3,103),xlim=(-.25,1.25));tag(ax,'d','Protein sharing in 22 genome pairs')
calls=pd.read_csv(Path(__file__).parent/'culture_nucleotide_cigar/per_gene.tsv',sep='\t');calls=calls[calls.eligible==1]
ax=axes[2,0]
for i,dark in enumerate([0,1]):
 for j,shared in enumerate([0,1]):
  g=calls[(calls.dark==dark)&(calls.protein_shared==shared)];x=i*3+j;ax.bar(x,100*g.dna_supported.mean(),color=teal if dark else grey);ax.text(x,100*g.dna_supported.mean()+2,f'{g.dna_supported.sum():,}/{len(g):,}',ha='center',fontsize=5.5)
ax.set(xticks=[0,1,3,4],xticklabels=['Unshared','Shared','Unshared','Shared'],ylabel='DNA-supported calls (%)',ylim=(0,114),xlabel='Annotated                          Dark');ax.tick_params(axis='x',labelsize=6);tag(ax,'e','Base-level nucleotide support')
ax=axes[2,1];v=calls[calls.protein_shared==1].groupby('rep').dna_supported.mean().sort_values();ax.plot(range(1,len(v)+1),100*v,'o',c=grey,ms=3);ax.axhline(90,c=orange,ls=':',lw=.8);ax.set(xlabel='Genome pair (sorted)',ylabel='Protein-shared calls with DNA support (%)',ylim=(0,103));tag(ax,'f','Per-pair positive-control recovery');save(fig,'fig_dark_current')
lin=frame('06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv');lin['dark_fraction']=lin['dark_fraction'].astype(float)
fig,axes=plt.subplots(1,3,figsize=(7.4,3.2),layout='constrained');ax=axes[0]
groups=lin.groupby('phylum')['dark_fraction'];names=[k for k,g in groups if len(g)>=20];names=sorted(names,key=lambda k:groups.get_group(k).median())[-12:]
ax.boxplot([100*groups.get_group(k) for k in names],vert=False,showfliers=False,tick_labels=names,widths=.5);ax.tick_params(axis='y',labelsize=6);ax.set_xlabel('Dark proteins per genome (%)');tag(ax,'a','Largest phylum medians')
ax=axes[1];x=lin['first_unnamed_rank'].astype(int);levels=sorted(x.unique());ax.boxplot([100*lin.loc[x==k,'dark_fraction'] for k in levels],showfliers=False,tick_labels=[str(k) for k in levels]);ax.set(xlabel='First empty/digit-bearing rank (index)',ylabel='Dark proteins per genome (%)');tag(ax,'b','Reference naming depth')
mat=frame('06_darkmatter/lineage_darkness_clean/results/v2_matched_null.tsv');ax=axes[2]
for i,pool in enumerate(['archaea','bacteria']):
 a=mat[mat['pool']==pool];ax.boxplot([100*a['target_dark'].astype(float),100*a['control_dark'].astype(float)],positions=[4*i,4*i+1.7],showfliers=False,widths=.6)
ax.set(xticks=[0,1.7,4,5.7],xticklabels=['T','A','T','B'],ylabel='Dark proteins per genome (%)');ax.text(.24,.94,'218 pairs',transform=ax.transAxes,ha='center',fontsize=6);ax.text(.78,.94,'237 pairs',transform=ax.transAxes,ha='center',fontsize=6);tag(ax,'c','Quality/abundance matches');save(fig,'fig_lineage_current')
fig,axes=plt.subplots(1,3,figsize=(7.4,3),layout='constrained')
for arm,fn,c in [('Dark','per_rep_structures.tsv',teal),('Control','per_control_structures.tsv',grey)]:
 st=frame('04_structure/00_structure_clean/05_validate/'+fn);st=st[st['confidence_qualified']=='1'];q=st[['best_pdb_qtm','best_afdb_qtm']].astype(float).max(axis=1);a=st[['best_pdb_alntm','best_afdb_alntm']].astype(float).max(axis=1);z=(q==0)&(a==0)
 axes[0].scatter(q,a,s=3,alpha=.23,c=c,label=arm,rasterized=True)
 thresholds=np.arange(.4,.601,.025);axes[1].plot(thresholds,[100*(q<t).mean() for t in thresholds],c=c,label=arm+' query');axes[1].plot(thresholds,[100*(a<t).mean() for t in thresholds],c=c,ls='--',label=arm+' alignment')
 i=0 if arm=='Dark' else 1
 axes[2].bar(i,100*z.mean(),color=orange,label='No retrieved rows' if i==0 else None);axes[2].bar(i,100*((q<.5)&~z).mean(),bottom=100*z.mean(),color=c,label='Score below 0.5' if i==0 else None)
axes[0].set(xlabel='Best query-normalised TM-score',ylabel='Best alignment-normalised TM-score',xlim=(0,1),ylim=(0,1));axes[0].axhline(.5,c='.6',ls=':',lw=.8);axes[0].axvline(.5,c='.6',ls=':',lw=.8);axes[0].legend(frameon=False,markerscale=3,fontsize=7);tag(axes[0],'a','Metric-dependent coverage')
axes[1].set(xlabel='TM-score threshold',ylabel='Models below threshold (%)');axes[1].legend(frameon=False,fontsize=6.5);tag(axes[1],'b','Threshold sensitivity')
axes[2].set(xticks=[0,1],xticklabels=['Dark','Control'],ylabel='Models below query threshold (%)');axes[2].legend(frameon=False,fontsize=6);tag(axes[2],'c','Retrieval and score failures');save(fig,'fig_structure_current')
print('Rendered figures 2, 3 and 4 from verified compact inputs')
meta=frame('metadata/sample_metadata.tsv');meta['trip_num']=meta['trip_num'].astype(int)
fig,axes=plt.subplots(1,3,figsize=(7.4,3.0),layout='constrained')
a=meta.groupby('trip_num').size();axes[0].bar(a.index,a.values,color=teal);axes[0].set(xlabel='Expedition',ylabel='Soil samples',xticks=a.index);tag(axes[0],'a','Repeated sampling')
m=pd.crosstab(meta['compartment'],meta['trip_num']);axes[1].imshow(m.values,cmap='Blues',aspect='auto');axes[1].set(xticks=range(len(m.columns)),xticklabels=m.columns,yticks=range(len(m.index)),yticklabels=[{'Deep':'Subsurface','PlantRoot':'Root-adjacent'}.get(k,k) for k in m.index],xlabel='Expedition')
for i in range(len(m.index)):
 for j in range(len(m.columns)):axes[1].text(j,i,str(m.iloc[i,j]),ha='center',va='center',color='white' if m.iloc[i,j]>m.values.max()*.6 else '#253444')
tag(axes[1],'b','Soil positions')
axes[2].scatter(lin['completeness'].astype(float),lin['contamination'].astype(float),s=3,c=teal,alpha=.3,rasterized=True);axes[2].set(xlabel='Completeness (%)',ylabel='Contamination (%)');tag(axes[2],'c','5,438 species representatives');save(fig,'fig_study_current')

ref=frame('06_reference_bias/refbias_clean/results/per_sample_comparison.tsv');dep=frame('06_darkmatter/specificity_clean/results/content_vs_depth.tsv');occ=frame('06_darkmatter/specificity_clean/results/b6_absence_occupancy_by_length.tsv')
fig,axes=plt.subplots(1,3,figsize=(7.5,3.1),layout='constrained')
ax=axes[0];ax.boxplot([ref.recall_of_mag.astype(float),ref.spearman_shared.astype(float)],tick_labels=['Recall','Rank correlation'],showfliers=False);ax.set(ylabel='Per-sample agreement',ylim=(0,1));ax.tick_params(axis='x',labelsize=6);tag(ax,'a','448 paired profiles')
ax=axes[1]
for prefix,c,label in [('dark',teal,'Dark'),('ctrl',grey,'Annotated')]:ax.scatter(dep[prefix+'_exp_size_only'].astype(float),dep[prefix+'_obs'].astype(float),s=18,c=c,label=label)
ax.plot([0,100],[0,100],c='.6',lw=.7);ax.set(xlabel='Analytical expected absence (%)',ylabel='Observed absence (%)',xlim=(0,100),ylim=(0,100));ax.legend(frameon=False,fontsize=7);tag(ax,'b','Catalogue composition')
ax=axes[2]
for col,c,label in [('median_sites_absent',teal,'No catalogue hit'),('median_sites_has_hit',grey,'Catalogue hit')]:ax.plot(range(len(occ)),occ[col].astype(float),'o-',c=c,label=label)
ax.set(xticks=range(len(occ)),xticklabels=occ.len_stratum,xlabel='Representative length',ylabel='Median occupied locations',ylim=(0,61));ax.tick_params(axis='x',rotation=30,labelsize=5.5);ax.legend(frameon=False,fontsize=6.5);tag(ax,'c','Occupancy by length');save(fig,'fig_reference_current')
