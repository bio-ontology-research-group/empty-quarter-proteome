#!/usr/bin/env python3
"""Regenerate the manuscript figure from pinned numerical results; no inference rerun."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
N=Path(__file__).resolve().parent
R=N.parents[2]/'manuscript/figures'
C=N/'results/published_clusters'
def tab(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
allmeta=tab(N/'results/f3_comparison/taxonomy.tsv')
reads={r['id']:r for r in tab(N/'results/read_support/loci.tsv')}
coord=[r for r in tab(N/'results/expanded_analysis/coordinates.tsv')if r['seed']=='42'and r['perplexity']=='30']
plotids=[r['id']for r in coord];sources={r['id']:r['source']for r in coord}
primary=np.array([[float(r['x']),float(r['y'])]for r in coord])
assert len(coord)==1726
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
col={'EQ':'#b54232','SMAG':'#007f7b','GEM':'#456daf','curated_reference':'#555555','Published':'#b8b8b8'}
def scatter(ax,y):
    for src in ['Published','curated_reference','GEM','SMAG','EQ']:
        jj=[j for j,pid in enumerate(plotids)if sources[pid]==src]
        ax.scatter(y[jj,0],y[jj,1],s=7 if src=='Published'else 16,c=col[src],alpha=.45 if src=='Published'else.9,label=src.replace('curated_reference','Curated panel'),linewidths=0)
    ax.set(xlabel='t-SNE 1',ylabel='t-SNE 2');ax.legend(frameon=False,fontsize=6,markerscale=1.3)
fig,axes=plt.subplots(2,2,figsize=(8.5,7),layout='constrained');scatter(axes[0,0],primary);axes[0,0].set_title('a  Representation with published heliorhodopsins')
ax=axes[0,1];rr=tab(C/'rarefaction.tsv');target=[r for r in rr if r['rank']=='class'and r['taxon']=='Saccharimonadia'and r['comparison']=='EQ+GEM'and r['identity']=='50'];labels=[]
for i,(kind,src)in enumerate([(k,s)for k in ['full','core']for s in ['EQ','GEM']]):
    r=next(r for r in target if r['kind']==kind and r['source']==src);v=float(r['mean_external_clusters']);ax.bar(i,v,color=col[src],width=.65);ax.errorbar(i,v,yerr=[[v-float(r['q025'])],[float(r['q975'])-v]],fmt='none',color='#333333',capsize=3);labels.append(src+'\n'+('Full length'if kind=='full'else'Profile core'))
ax.set_xticks(range(4),labels);ax.set(title='b  Saccharimonadia: 10-sequence additions',ylabel='Clusters absent from expanded reference')
phcol={'Patescibacteriota':'#b54232','Actinomycetota':'#007f7b','Other':'#756a9e'};ph={r['id']:r['phylum']for r in allmeta};struct=tab(N/'results/structure_audit/structures.tsv');ax=axes[1,0]
for label,c in phcol.items():
    rr=[r for r in struct if r['source']=='EQ'and(ph[r['id']]if ph[r['id']]in phcol else'Other')==label];ax.scatter([float(r['mean_CA_pLDDT'])for r in rr],[float(r['query_TM'])for r in rr],s=22,c=c,label=label,alpha=.85)
ax.set(title='c  EQ predicted structural compatibility',xlabel='Mean Cα pLDDT',ylabel='Query-normalized TM-score to 6SU3',ylim=(.5,1));ax.legend(frameon=False,fontsize=6)
ax=axes[1,1];rr=sorted([r for r in reads.values()if r['track']=='F3_Helio'],key=lambda r:float(r['median_depth']))
for i,r in enumerate(rr):
    p=ph['Helio__EQ__'+r['id']];c=phcol.get(p,phcol['Other']);ax.scatter(i+1,float(r['median_depth']),s=23,c=c,marker='o'if r['coverage_support_gate']=='1'else'x')
ax.axhline(5,lw=.8,ls=':',color='#555555');ax.set(title='d  Original-read support at EQ loci',xlabel='Locus rank by depth',ylabel='Median depth (MAPQ and base quality ≥20)',yscale='log');ax.text(.03,.97,'43/45 pass depth and breadth criteria',transform=ax.transAxes,va='top',fontsize=7)
fig.savefig(R/'fig_helio_expanded.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(R/'fig_helio_expanded.png',dpi=200);plt.close(fig)
