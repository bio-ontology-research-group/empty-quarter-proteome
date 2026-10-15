from common import *
import json,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE,trustworthiness
R=N/'results/figures';R.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
meta=tab(N/'results/f3_finer/taxonomy_unique.tsv');z=np.load(D/'root_audit/f345_embeddings/embeddings.npz',allow_pickle=True);zi={str(x):i for i,x in enumerate(z['ids'])};emb=z['embeddings'];z.close();x=emb[[zi[r['id']]for r in meta]];x=x/np.linalg.norm(x,axis=1,keepdims=True);pc=PCA(n_components=50,random_state=42).fit_transform(x)
colors={'EQ':'#b54232','SMAG':'#007f7b','GEM':'#456daf','curated_reference':'#777777'};phcolors={'Patescibacteriota':'#b54232','Actinomycetota':'#007f7b','Other':'#aaa4bc'};co=[];quality=[]
for seed in [42,123,20261015]:
    for perp in [5,30,50]:
        y=TSNE(n_components=2,perplexity=perp,random_state=seed,init='pca',learning_rate='auto',max_iter=2000).fit_transform(pc)
        for r,(a,b)in zip(meta,y):co.append({'id':r['id'],'seed':seed,'perplexity':perp,'x':float(a),'y':float(b)})
        quality.append({'seed':seed,'perplexity':perp,'trustworthiness_cosine_k5':float(trustworthiness(x,y,n_neighbors=5,metric='cosine'))})
        if seed==42 and perp==30:primary=y
        fig,axes=plt.subplots(1,2,figsize=(8,3.5),layout='constrained')
        for src,c in colors.items():
            ii=[i for i,r in enumerate(meta)if r['source']==src];axes[0].scatter(y[ii,0],y[ii,1],s=17,c=c,label=src,alpha=.85)
        for ph,c in phcolors.items():
            ii=[i for i,r in enumerate(meta)if(r['phylum']if r['phylum']in phcolors else'Other')==ph];axes[1].scatter(y[ii,0],y[ii,1],s=17,c=c,label=ph,alpha=.85)
        for ax in axes:ax.set(xlabel='t-SNE 1',ylabel='t-SNE 2');ax.legend(fontsize=6,frameon=False)
        fig.savefig(R/f'helio_tsne_seed{seed}_p{perp}.pdf');plt.close(fig)
write(R/'tsne_coordinates.tsv',co);write(R/'tsne_quality.tsv',quality)
fig,axes=plt.subplots(2,2,figsize=(8.5,7.3),layout='constrained');ax=axes[0,0]
for src,c in colors.items():
    ii=[i for i,r in enumerate(meta)if r['source']==src];ax.scatter(primary[ii,0],primary[ii,1],s=19,c=c,label=src.replace('curated_reference','Reference panel'),alpha=.85,edgecolors='white',linewidths=.25)
ax.set(title='a  Protein representation by source',xlabel='t-SNE 1',ylabel='t-SNE 2');ax.legend(fontsize=7,frameon=False)
ax=axes[0,1]
for ph,c in phcolors.items():
    ii=[i for i,r in enumerate(meta)if(r['phylum']if r['phylum']in phcolors else'Other')==ph];ax.scatter(primary[ii,0],primary[ii,1],s=19,c=c,label=ph,alpha=.85,edgecolors='white',linewidths=.25)
ax.set(title='b  The same projection by host phylum',xlabel='t-SNE 1',ylabel='t-SNE 2');ax.legend(fontsize=7,frameon=False)
ax=axes[1,0];rr=tab(N/'results/f3_finer/rarefaction.tsv');target=[r for r in rr if r['rank']=='class'and r['taxon']=='Saccharimonadia'and r['comparison']=='EQ+GEM'and r['identity']=='50'];labels=[]
for i,(kind,src)in enumerate([(k,s)for k in ['full','core']for s in ['EQ','GEM']]):
    r=next(r for r in target if r['kind']==kind and r['source']==src);v=float(r['mean_external_clusters']);ax.bar(i,v,color=colors[src],width=.65);ax.errorbar(i,v,yerr=[[v-float(r['q025'])],[float(r['q975'])-v]],fmt='none',color='#333333',capsize=3);labels.append(src+'\n'+('Full length'if kind=='full'else'Profile core'))
ax.set_xticks(range(4),labels);ax.set(title='c  Matched Saccharimonadia additions',ylabel='Reference-unrepresented 50% clusters');ax.text(.02,.97,'Equal additions: '+target[0]['equal_n']+' sequences',transform=ax.transAxes,va='top',fontsize=7)
ax=axes[1,1];struct=tab(N/'results/structure_audit/structures.tsv');done=set()
for src,c in colors.items():
    rr=[]
    for r in struct:
        k=(r['source'],r['structure_id'])
        if r['source']==src and k not in done:rr.append(r);done.add(k)
    ax.scatter([float(r['mean_CA_pLDDT'])for r in rr],[float(r['query_TM'])for r in rr],c=c,s=19,alpha=.8,edgecolors='white',linewidths=.25)
ax.set(title='d  Predicted structural compatibility',xlabel='Mean Cα pLDDT',ylabel='Query-normalized TM-score to 6SU3',ylim=(0,1));ax.axhline(.5,color='#999999',lw=.7,ls=':')
fig.savefig(R/'fig_helio_followup.pdf');fig.savefig(R/'fig_helio_followup.png',dpi=200);plt.close(fig)
# Supporting ecology: same mean-based additive estimand for each decomposition component.
ec=tab(N/'results/ecology/decomposition_summary.tsv');fig,axes=plt.subplots(1,2,figsize=(8.5,3.5),layout='constrained');cc={'total_pp':'#222222','between_pp':'#456daf','within_pp':'#007f7b','exclusive_pp':'#b54232'}
for ax,trait in zip(axes,['CPF','EctABC_strict']):
    labels=[]
    for i,(rank,metric)in enumerate([(r,m)for r in ['phylum','genus']for m in cc]):
        d=next(r for r in ec if r['trait']==trait and r['rank']==rank and r['metric']==metric);v=float(d['mean']);ax.errorbar(v,i,xerr=[[v-float(d['mean_ci_low'])],[float(d['mean_ci_high'])-v]],fmt='o',color=cc[metric],capsize=2);labels.append(rank+': '+metric.removesuffix('_pp'))
    ax.set_yticks(range(8),labels);ax.invert_yaxis();ax.axvline(0,color='#bbbbbb',lw=.7);ax.set(title=trait.replace('_',' '),xlabel='Surface − deep contribution (percentage points)')
fig.savefig(R/'fig_family_ecology.pdf');fig.savefig(R/'fig_family_ecology.png',dpi=200);plt.close(fig)
(R/'provenance.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'n_source_deduplicated':len(meta),'embedding_normalization':'L2 per sequence','PCA_dimensions':50,'tsne_iterations':2000,'primary_seed':42,'primary_perplexity':30,'caution':'Cluster rarefaction intervals are sampling ranges; ecology intervals are site bootstrap confidence intervals. Neither t-SNE separation nor predicted structure establishes function.'},indent=2))
