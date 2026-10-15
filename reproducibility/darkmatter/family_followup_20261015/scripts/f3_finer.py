from common import *
import collections,json,random,numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner,substitution_matrices
R=N/'results/f3_finer';R.mkdir(parents=True,exist_ok=True);C=N/'results/f3_comparison'
meta=tab(C/'taxonomy.tsv');seen=set();rr=[]
for r in sorted(meta,key=lambda r:r['id']):
    k=(r['source'],r['sequence_sha256'])
    if k in seen:continue
    seen.add(k);r['class']=next((v[3:]for v in r['published_taxonomy'].split(';')if v.startswith('c__')),'Unresolved');r['order']=next((v[3:]for v in r['published_taxonomy'].split(';')if v.startswith('o__')),'Unresolved');rr.append(r)
write(R/'taxonomy_unique.tsv',rr);out=[];counts=[]
for rank in ['class','order']:
    for value in sorted({r[rank]for r in rr if r['source']=='EQ'}):
        groups={s:[r['id']for r in rr if r['source']==s and r[rank]==value]for s in ['EQ','SMAG','GEM']}
        counts.append({'rank':rank,'taxon':value,**{s:len(ids)for s,ids in groups.items()}})
        for kind in ['full','core']:
            for lev in ['50','80']:
                cm={mem:rep for rep,mem in (l.split()for l in open(C/(kind+lev+'_cluster.tsv')))};refs={cm[r['id']]for r in rr if r['source']=='curated_reference'}
                for comp in ['SMAG','GEM']:
                    n=min(len(groups['EQ']),len(groups[comp]))
                    if n<5 or value=='Unresolved':continue
                    for src in ['EQ',comp]:
                        vals=[len({cm[x]for x in random.Random(20261015+j).sample(groups[src],n)}-refs)for j in range(1000)]
                        other=comp if src=='EQ'else'EQ';otherclus={cm[x]for x in groups[other]}
                        out.append({'rank':rank,'taxon':value,'kind':kind,'identity':lev,'comparison':'EQ+'+comp,'source':src,'equal_n':n,'mean_external_clusters':float(np.mean(vals)),'q025':float(np.quantile(vals,.025)),'q975':float(np.quantile(vals,.975)),'full_source_clusters_absent_other_panel':len({cm[x]for x in groups[src]}-otherclus)})
write(R/'counts.tsv',counts);write(R/'rarefaction.tsv',out)
# All-pairs global core identities in the focal matched class, independently of embeddings.
seq={r.id:str(r.seq)for r in SeqIO.parse(C/'core.faa','fasta')};pa=PairwiseAligner();pa.mode='global';pa.substitution_matrix=substitution_matrices.load('BLOSUM62');pa.open_gap_score=-10;pa.extend_gap_score=-.5
focal=[r for r in rr if r['class']=='Saccharimonadia'and r['source']in {'EQ','SMAG','GEM'}];pairs=[];nearest=[]
for i,r in enumerate(focal):
    for t in focal[i+1:]:
        a=pa.align(seq[r['id']],seq[t['id']])[0];pid=100*a.counts().identities/a.shape[1];pairs.append({'query':r['id'],'target':t['id'],'identity':pid})
for r in focal:
    for src in ['EQ','SMAG','GEM']:
        targets={t['id']for t in focal if t['source']==src and t['genome']!=r['genome']};hits=[(x['identity'],x['target']if x['query']==r['id']else x['query'])for x in pairs if(x['query']==r['id']and x['target']in targets)or(x['target']==r['id']and x['query']in targets)]
        if hits:
            val,pid=max(hits);nearest.append({'id':r['id'],'source':r['source'],'target_source':src,'n_comparators':len(targets),'max_core_global_identity':val,'target':pid})
write(R/'Saccharimonadia_pairwise.tsv',pairs);write(R/'Saccharimonadia_nearest.tsv',nearest)
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'focal_class':'Saccharimonadia','counts':dict(collections.Counter(r['source']for r in focal)),'limits':['Finer taxonomy from source versions, no genome reclassification','Global core identity uses aligned length including gaps','Sequence diversity among sampled panels, not endemism or functional novelty']},indent=2))
