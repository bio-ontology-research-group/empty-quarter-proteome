from common import *
import json,random,collections,math,subprocess,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE,trustworthiness
R=N/'results/expanded_analysis';R.mkdir(parents=True,exist_ok=True);C=N/'results/published_clusters';P=N/'results/published_helio'
meta=tab(N/'results/f3_finer/taxonomy_unique.tsv');by={r['id']:r for r in meta};pub={r.id:str(r.seq)for r in SeqIO.parse(P/'site_conserving_references.faa','fasta')};published_hashes={hashlib.sha256(s.encode()).hexdigest()for s in pub.values()};reads={r['id']:r for r in tab(N/'results/read_support/loci.tsv')};allmeta=tab(N/'results/f3_comparison/taxonomy.tsv');passed_hashes={r['sequence_sha256']for r in allmeta if r['source']=='EQ'and reads[r['original_id']]['coverage_support_gate']=='1'}
out=[]
for kind in ['full','core']:
    for lev in ['50','80']:
        cm={m:r for r,m in(l.split()for l in open(C/(kind+lev+'_cluster.tsv')))};refclus={cm[i]for i in pub}|{cm[r['id']]for r in meta if r['source']=='curated_reference'};groups={s:[r['id']for r in meta if r['source']==s and r['class']=='Saccharimonadia'and(s!='EQ'or r['sequence_sha256']in passed_hashes)]for s in ['EQ','GEM']};n=min(map(len,groups.values()));assert n>=5
        for src,ids in groups.items():
            vals=[len({cm[x]for x in random.Random(20261015+j).sample(ids,n)}-refclus)for j in range(1000)];out.append({'kind':kind,'identity':lev,'source':src,'read_supported_EQ_unique':len(groups['EQ']),'equal_n':n,'mean_external_clusters':float(np.mean(vals)),'q025':float(np.quantile(vals,.025)),'q975':float(np.quantile(vals,.975))})
write(R/'read_supported_rarefaction.tsv',out)
# Recompute all-source embeddings in one common GPU run; exact-sequence dedup for projection only.
z=np.load(N/'expanded_embeddings/embeddings.npz',allow_pickle=True);ids=list(map(str,z['ids']));emb=z['embeddings'];hashes=list(map(str,z['clean_sequence_sha256']))if'clean_sequence_sha256'in z.files else None;z.close()
seq={r.id:str(r.seq)for r in SeqIO.parse(C/'full.faa','fasta')};assert set(ids)==set(seq)
sources={i:(by[i]['source']if i in by else('curated_reference'if i.startswith('Helio__curated_reference__')else('EQ'if i.startswith('Helio__EQ__')else('SMAG'if i.startswith('Helio__SMAG__')else('GEM'if i.startswith('Helio__GEM__')else'Published')))))for i in ids}
ii={x:i for i,x in enumerate(ids)};refs=[i for i in ids if sources[i]in {'Published','curated_reference'}];query=[r for r in allmeta if r['source']in {'EQ','SMAG','GEM'}];dd=cdist(emb[[ii[r['id']]for r in query]],emb[[ii[i]for i in refs]],'cosine');distrows=[]
for k,r in enumerate(query):
    j=int(dd[k].argmin());distrows.append({'id':r['id'],'source':r['source'],'phylum':r['phylum'],'cosine_nearest_expanded_reference':float(dd[k,j]),'nearest_reference':refs[j],'exact_sequence_in_published_panel':int(r['sequence_sha256']in published_hashes)})
write(R/'native_distances.tsv',distrows)
# Stable source-priority dedup: retain EQ point if identical to a published reference, flag overlap explicitly.
priority={'EQ':0,'SMAG':1,'GEM':2,'curated_reference':3,'Published':4};seen=set();plotids=[]
for pid in sorted(ids,key=lambda pid:(priority[sources[pid]],pid)):
    h=hashlib.sha256(seq[pid].encode()).hexdigest()
    if h not in seen:seen.add(h);plotids.append(pid)
x=emb[[ii[i]for i in plotids]];x=x/np.linalg.norm(x,axis=1,keepdims=True);pc=PCA(n_components=50,random_state=42).fit_transform(x);coords=[];quality=[]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
col={'EQ':'#b54232','SMAG':'#007f7b','GEM':'#456daf','curated_reference':'#555555','Published':'#b8b8b8'}
def scatter(ax,y):
    for src in ['Published','curated_reference','GEM','SMAG','EQ']:
        jj=[j for j,pid in enumerate(plotids)if sources[pid]==src];ax.scatter(y[jj,0],y[jj,1],s=7 if src=='Published'else 16,c=col[src],alpha=.45 if src=='Published'else.9,label=src.replace('curated_reference','Reference anchors'),linewidths=0)
    ax.set(xlabel='t-SNE 1',ylabel='t-SNE 2');ax.legend(frameon=False,fontsize=6,markerscale=1.3)
for seed in [42,123,20261015]:
    for perp in [5,30,50]:
        initialization='pca'if seed==42 else'random'
        y=TSNE(n_components=2,perplexity=perp,random_state=seed,init=initialization,learning_rate='auto',max_iter=2000).fit_transform(pc)
        coords += [{'id':pid,'source':sources[pid],'seed':seed,'perplexity':perp,'initialization':initialization,'x':float(v[0]),'y':float(v[1])}for pid,v in zip(plotids,y)];quality.append({'seed':seed,'perplexity':perp,'initialization':initialization,'trustworthiness_cosine_k5':float(trustworthiness(x,y,n_neighbors=5,metric='cosine'))})
        if seed==42 and perp==30:primary=y
        fig,ax=plt.subplots(figsize=(5,4),layout='constrained');scatter(ax,y);fig.savefig(R/f'expanded_tsne_seed{seed}_p{perp}.pdf');plt.close(fig)
write(R/'coordinates.tsv',coords);write(R/'projection_quality.tsv',quality)
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
fig.savefig(R/'fig_helio_expanded.pdf');fig.savefig(R/'fig_helio_expanded.png',dpi=200);plt.close(fig)
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'input_records':len(ids),'projection_unique_sequences':len(plotids),'native_source_medians':{s:float(np.median([r['cosine_nearest_expanded_reference']for r in distrows if r['source']==s]))for s in ['EQ','SMAG','GEM']},'exact_published_overlap':{s:sum(r['exact_sequence_in_published_panel']for r in distrows if r['source']==s)for s in ['EQ','SMAG','GEM']},'EQ_read_supported_unique_sequences':len(passed_hashes),'limitations':['Reference sequences overlap public resources; comparator proximity is descriptive','t-SNE preserves local relationships imperfectly and is not a novelty or training-data test','Read quality sensitivity applies only to EQ; external raw reads were not evaluated']},indent=2))
