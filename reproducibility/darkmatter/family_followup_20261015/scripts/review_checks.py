"""Root checks prompted by Argon's independent results review."""
from common import *
import json,collections,subprocess,random,numpy as np
from Bio import SeqIO
R=N/'results/review_checks';R.mkdir(parents=True,exist_ok=True)
meta=tab(N/'results/f3_finer/taxonomy_unique.tsv');sac=[r for r in meta if r['class']=='Saccharimonadia'and r['source']in {'EQ','GEM','SMAG'}]
rows=[]
for kind in ['full','core']:
    seq={r.id:str(r.seq)for r in SeqIO.parse(N/'results/published_clusters'/(kind+'.faa'),'fasta')}
    with open(R/(kind+'.faa'),'w')as f:
        for r in sac:f.write('>'+r['id']+'\n'+seq[r['id']]+'\n')
    for level in [50,80]:
        stem=R/(kind+str(level))
        with open(R/(kind+str(level)+'.log'),'w')as log:subprocess.run(['/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs','easy-cluster',str(R/(kind+'.faa')),str(stem),str(R/('tmp'+kind+str(level))),'--min-seq-id',str(level/100),'-c','0.8','--cov-mode','0','--threads','4'],stdout=log,stderr=log,check=True)
        cm={m:r for r,m in(l.split()for l in open(str(stem)+'_cluster.tsv'))}
        for source in ['EQ','GEM']:
            ids=[r['id']for r in sac if r['source']==source];vals=[len({cm[i]for i in random.Random(20261015+j).sample(ids,10)})for j in range(1000)]
            rows.append({'kind':kind,'identity':level,'source':source,'unique_sequences':len(ids),'clusters':len({cm[i]for i in ids}),'equal_n':10,'mean_total_clusters':float(np.mean(vals)),'q025':float(np.quantile(vals,.025)),'q975':float(np.quantile(vals,.975))})
write(R/'reference_free_rarefaction.tsv',rows)
f4=tab(N/'results/f4_integrity/loci.tsv');dual=[r for r in f4 if r['strict_ABC']=='1'and r['complete_A']==r['complete_B']];domains=[]
for r in dual:
    hits=collections.defaultdict(list)
    for line in open(D/'f345/results/full_v1/f4_ectc'/('f4_domtbl_'+r['source']+'.tsv')):
        if line.startswith('#'):continue
        x=line.split()
        if x[0]==r['complete_A']and x[4].split('.')[0]in {'TIGR02406','TIGR02407'}:
            hits[x[4].split('.')[0]].append({'protein':x[0],'length':int(x[2]),'profile':x[4].split('.')[0],'score':float(x[13]),'ievalue':float(x[12]),'profile_coverage':(int(x[16])-int(x[15])+1)/int(x[5]),'start':int(x[17]),'end':int(x[18])})
    chosen=[max(hits[k],key=lambda x:x['score'])for k in ['TIGR02406','TIGR02407']]
    overlap=max(0,min(x['end']for x in chosen)-max(x['start']for x in chosen)+1)
    assert all(x['profile_coverage']>=.7 and x['ievalue']<=1e-5 for x in chosen)
    for h in chosen:domains.append({'EctC_id':r['id'],'source':r['source'],'domain_overlap_aa':overlap,**h})
write(R/'dual_AB_domains.tsv',domains)
# Correct the misleading A/B order tie when both annotations belong to one ORF.
corrected=[]
for r in f4:
    if r['strict_ABC']!='1':continue
    order=r['gene_order']
    if r['complete_A']==r['complete_B']:order=order.replace('A-B','AB').replace('B-A','AB')
    corrected.append({'id':r['id'],'source':r['source'],'A_B_distinct_ORFs':int(r['complete_A']!=r['complete_B']),'marker_order':order})
write(R/'neighborhood_gene_orders.tsv',corrected)
pub={str(r.seq)for r in SeqIO.parse(N/'results/published_helio/site_conserving_references.faa','fasta')};seq={r.id:str(r.seq)for r in SeqIO.parse(N/'results/published_clusters/full.faa','fasta')}
overlap=[]
for r in meta:
    if r['source']not in {'EQ','GEM','SMAG'}:continue
    s=seq[r['id']];overlap.append({'id':r['id'],'source':r['source'],'class':r['class'],'exact_published':int(s in pub),'substring_published':int(any(s in p or p in s for p in pub))})
write(R/'published_overlap_unique.tsv',overlap)
(R/'summary.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'reference_free_input_sequences':len(sac),'EQ_ABC_distinct_ORFs':sum(r['source']=='EQ_MAG'and r['A_B_distinct_ORFs']for r in corrected),'EQ_ABC_single_AB_ORF':sum(r['source']=='EQ_MAG'and not r['A_B_distinct_ORFs']for r in corrected),'dual_AB_domain_max_overlap':max(r['domain_overlap_aa']for r in domains),'limits':['Reference-free clustering still depends on heuristic partition and sampled sequences','Dual-profile AB ORFs are putative fusions, without targeted junction validation','Corrected marker order merges A/B assignments to the same ORF; old tie-sorted labels are not physical gene order']},indent=2))
