#!/usr/bin/env python3
"""Rebuild Fig. 4 from the frozen qTM baseline and paired region census."""
from pathlib import Path
import csv, gzip, hashlib, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; HERE=Path(__file__).resolve().parent
for rel, expected in json.loads((HERE/'inputs.sha256.json').read_text()).items():
    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==expected,rel
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'text.color':'#253444','axes.labelcolor':'#253444'})
teal='#187f83';grey='#73849b';orange='#d18948'
fig,axs=plt.subplots(1,3,figsize=(8.8,3.15),layout='constrained',gridspec_kw={'width_ratios':[1,1.05,1.15]})
for name,fn,c in [('Dark (2,838)','per_rep_structures.tsv.gz',teal),('Controls (3,013)','per_control_structures.tsv.gz',grey)]:
    with gzip.open(ROOT/'reproducibility/darkmatter/inputs/04_structure/00_structure_clean/05_validate'/fn,'rt') as f:
        rr=[r for r in csv.DictReader(f,delimiter='\t') if r['confidence_qualified']=='1']
    q=np.sort([max(float(r['best_pdb_qtm']),float(r['best_afdb_qtm'])) for r in rr])
    axs[0].plot(q,np.arange(1,len(q)+1)/len(q),c=c,label=name,lw=1.6)
axs[0].axvline(.5,c='.5',ls=':',lw=.9);axs[0].set(xlabel='Best query-normalised TM-score',ylabel='Cumulative fraction',xlim=(0,1),ylim=(0,1));axs[0].legend(frameon=False,fontsize=7,loc='lower right')
with open(HERE/'v7_stratum_census.tsv') as f: rows=[r for r in csv.DictReader(f,delimiter='\t') if r['tier']=='tier_1']
groups=[['K_local_pdb','K_local_afdb_only','K_dual_none'],['Q_not_K_pdb','Q_not_K_afdb_only'],['Controls']]
def total(group,col):return sum(int(r[col]) for r in rows if r['stratum'] in group)
for i,g in enumerate(groups):
    n=total(g,'model_count')
    for j,(col,c,label) in enumerate([('wc_matched',grey,'Whole chain'),('crop_matched',teal,'Any crop')]):
        val=100*total(g,col)/n;x=i+(j-.5)*.34
        axs[1].bar(x,val,width=.31,color=c,label=label if i==0 else None)
        axs[1].text(x,val+2,f'{val:.1f}',ha='center',fontsize=6.5)
axs[1].set(xticks=range(3),xticklabels=['Initially low\nquery score','Other dark','Controls'],ylabel='Models with match (%)',ylim=(0,116));axs[1].tick_params(axis='x',labelsize=6.5);axs[1].legend(frameon=False,fontsize=7,loc='upper left')
cols=['trans_recovered','trans_persistent_low','trans_raw_zero','trans_post_empty','trans_unm_no_reg']
labels=['Crop match','Retrieved, low score','No retrieved rows','Empty after exclusion','No eligible region']
vals=[total(groups[0],c) for c in cols]
assert sum(vals)==964 and vals[0]==714
axs[2].barh(range(5),vals,color=[teal,grey,orange,'#acb8c6','#c6cdd5']);axs[2].invert_yaxis()
axs[2].set(yticks=range(5),yticklabels=labels,xlabel='Whole-chain-unmatched models',xlim=(0,830));axs[2].tick_params(axis='y',labelsize=6.5)
for i,v in enumerate(vals):axs[2].text(v+10,i,str(v),va='center',fontsize=7)
for a,letter,title in zip(axs,'abc',['Initial dual-reference screen','Paired PDB searches','Outcomes after cropping']):
    a.set_title(title,loc='left',fontsize=9,pad=14);a.text(-.12,1.10,letter,transform=a.transAxes,fontweight='bold',fontsize=12)
out=ROOT/'manuscript/figures';out.mkdir(exist_ok=True)
fig.savefig(out/'fig_structure_current.pdf',bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});fig.savefig(out/'fig_structure_current.png',dpi=180,bbox_inches='tight');plt.close(fig)
print('Figure 4 rebuilt from checksum-verified qTM inputs; pooled transition counts:',vals)
