from common import *
import json,subprocess,random,collections,math,numpy as np
from Bio import SeqIO,AlignIO
R=N/'results/published_clusters';R.mkdir(parents=True,exist_ok=True);P=N/'results/published_helio';C=N/'results/f3_comparison'
pub={r.id:str(r.seq)for r in SeqIO.parse(P/'site_conserving_references.faa','fasta')};al=AlignIO.read(P/'profile.sto','stockholm');rf=al.column_annotations['reference_annotation'];match=[i for i,c in enumerate(rf)if c not in '.- '];core={}
for r in al:
    if r.id not in pub:continue
    a=str(r.seq);jj=[i for i in match if a[i]not in '.-'];core[r.id]=''.join(c for c in a[min(jj):max(jj)+1]if c not in '.-').upper();assert core[r.id]in pub[r.id]
meta=tab(N/'results/f3_finer/taxonomy_unique.tsv');out=[]
for kind in ['full','core']:
    original={r.id:str(r.seq)for r in SeqIO.parse(C/('full_same_set.faa'if kind=='full'else'core.faa'),'fasta')};seq={**original,**(pub if kind=='full'else core)}
    with open(R/(kind+'.faa'),'w')as f:
        for pid,s in seq.items():f.write('>'+pid+'\n'+s+'\n')
    for lev in ['50','80']:
        prefix=R/(kind+lev)
        with open(R/(kind+lev+'.log'),'w')as log:subprocess.run(['/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs','easy-cluster',str(R/(kind+'.faa')),str(prefix),str(R/('tmp'+kind+lev)),'--min-seq-id',str(int(lev)/100),'-c','0.8','--cov-mode','0','--threads','4'],check=True,stdout=log,stderr=log)
        cm={m:r for r,m in(l.split()for l in open(str(prefix)+'_cluster.tsv'))};refclus={cm[i]for i in pub}|{cm[r['id']]for r in meta if r['source']=='curated_reference'}
        for rank,val in [('phylum','Patescibacteriota'),('class','Saccharimonadia'),('class','Acidimicrobiia')]:
            groups={s:[r['id']for r in meta if r['source']==s and r[rank]==val]for s in ['EQ','SMAG','GEM']}
            for comp in ['SMAG','GEM']:
                n=min(len(groups['EQ']),len(groups[comp]))
                if n<5:continue
                for src in ['EQ',comp]:
                    ids=groups[src];vals=[len({cm[x]for x in random.Random(20261015+j).sample(ids,n)}-refclus)for j in range(1000)];counts=collections.Counter(cm[i]for i in ids);den=math.comb(len(ids),n);ex=sum(1-(math.comb(len(ids)-k,n)/den if len(ids)-k>=n else 0)for c,k in counts.items()if c not in refclus)
                    out.append({'rank':rank,'taxon':val,'kind':kind,'identity':lev,'comparison':'EQ+'+comp,'source':src,'equal_n':n,'mean_external_clusters':float(np.mean(vals)),'q025':float(np.quantile(vals,.025)),'q975':float(np.quantile(vals,.975)),'analytical_expected_clusters':ex,'clusters_absent_expanded_reference':len(set(counts)-refclus)})
write(R/'rarefaction.tsv',out)
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'published_site_conserving_references':len(pub),'curated_panel_records':25,'reference_definition':'union of the1566profile/site-filtered literature sequences and25site-conserving curated panel sequences; redundant references retained as IDs in joint cluster partition','limits':['Operational clustering and sampled comparison panels','Conditional rarefaction ranges, not population confidence intervals','Expanded panel is a sensitivity, not an exhaustive global reference']},indent=2))
