"""Source-symmetric retinal-site sensitivity; preserves primary profile census and raw analyses."""
import csv,json,hashlib,os,random,statistics,subprocess
from pathlib import Path
from Bio import AlignIO
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import binomtest
D=Path('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean');A=D/'root_audit';C=A/'f345_comparison';R=A/'f3_site_conserving';R.mkdir(exist_ok=True)
def tab(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
 with open(p,'w')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
z0=np.load(A/'f345_embeddings/embeddings.npz',allow_pickle=True);z={k:z0[k]for k in z0.files};z0.close();ix={str(x):i for i,x in enumerate(z['ids'])};summary={};kept_eq={}
for fam,alpath,acc,pos in [('Type1',A/'f345_anchor_alignment/Type1.sto','P02945',229),('Helio',A/'f3_reference_sensitivity/helio.sto','AVZ43932.1',241)]:
 al=AlignIO.read(alpath,'stockholm');seqs={r.id:str(r.seq)for r in al};rid=fam+'__curated_reference__'+acc;count=0
 for col,c in enumerate(seqs[rid]):
  if c not in '.-':count+=1
  if count==pos:break
 assert al.column_annotations['reference_annotation'][col]not in '.- ';keep={pid for pid,s in seqs.items()if s[col].upper()=='K'};meta=tab(C/fam/'metadata.tsv');rr=[r for r in meta if r['id']in keep];kept_eq[fam]={r['original_id']for r in rr if r['source']=='EQ'};out=R/fam;out.mkdir(exist_ok=True);write(out/'metadata.tsv',rr);write(out/'excluded.tsv',[r for r in meta if r['id']not in keep])
 sources=['EQ','SMAG','GEM','curated_reference'];record={'source_counts':{s:sum(r['source']==s for r in rr)for s in sources},'closest_reference':{},'cluster_sensitivity':{}}
 byid={r['id']:r for r in rr};best={}
 for l in open(C/fam/'hits.tsv'):
  q,t,pid,qs,qe,ts,te,ev,bit,ql,tl=l.split()
  if q not in keep or t not in keep or q==t or (int(qe)-int(qs)+1)/int(ql)<.7 or (int(te)-int(ts)+1)/int(tl)<.7:continue
  key=(q,byid[t]['source']);rec={'query':q,'target_source':key[1],'target':t,'identity':float(pid),'bitscore':float(bit)}
  if key not in best or rec['bitscore']>best[key]['bitscore']:best[key]=rec
 write(out/'closest_hits.tsv',list(best.values()))
 for src in sources[1:]:
  vals=[v['identity']for k,v in best.items()if k[1]==src];record['closest_reference'][src]={'coverage_passing_hits':len(vals),'median_identity':statistics.median(vals)if vals else None}
 refs=[r for r in rr if r['source']=='curated_reference'];dd=cdist(z['embeddings'][[ix[r['id']]for r in rr]],z['embeddings'][[ix[r['id']]for r in refs]],'cosine');dist=[]
 for i,r in enumerate(rr):
  for j,ref in enumerate(refs):
   if r['id']==ref['id']:dd[i,j]=np.inf
  dist.append({'id':r['id'],'source':r['source'],'cosine_nearest_reference':float(dd[i].min())})
 write(out/'distances.tsv',dist);record['median_cosine']={src:statistics.median(v['cosine_nearest_reference']for v in dist if v['source']==src)for src in sources}
 for level in ['50','80']:
  cm={mem:rep for rep,mem in (l.split()for l in open(C/fam/('cluster'+level+'_cluster.tsv')))};refclus={cm[r['id']]for r in refs};groups={}
  for src in sources[:3]:
   seen=set();gg=[]
   for r in sorted(rr,key=lambda r:r['id']):
    if r['source']==src and r['sequence_sha256']not in seen:seen.add(r['sequence_sha256']);gg.append(r['id'])
   groups[src]=gg
  n=min(100,*map(len,groups.values()));res={'equal_addition_n':n,'fixed_joint_clustering':True,'sources':{}}
  for src,ii in groups.items():res['sources'][src]={'sequences':len(ii),'clusters_outside_reference':len({cm[i]for i in ii}-refclus),'mean_equal_addition_clusters_outside_reference':statistics.mean(len({cm[i]for i in random.Random(20261015+k).sample(ii,n)}-refclus)for k in range(100))if n else None}
  record['cluster_sensitivity'][level]=res
 ii=[ix[r['id']]for r in rr];np.savez_compressed(out/'embeddings.npz',**{k:z[k][ii]for k in z});summary[fam]=record
 subprocess.run(['/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python','/home/hohndor/eq-family-discovery-20261015/cpu/family_embedding_analysis.py','--embeddings',str(out/'embeddings.npz'),'--metadata',str(out/'metadata.tsv'),'--output-dir',str(out/'plots'),'--seeds','42,123,20261015','--perplexities','5,30,50','--tsne-max-iter','1000'],check=True)
# Quality-conditioned context sensitivity; same matching rule, locus selection independent of marker outcome.
context=tab(A/'f3_context_audit/all_loci.tsv');cr={};selected_all=[]
for fam,full in [('Type1','Type-1_Microbial_Rhodopsin'),('Helio','Heliorhodopsin_Type-3')]:
 eligible=[r for r in context if r['family']==full and r['id']in kept_eq[fam]and r['eligible']=='1'];seen=set();selected=[]
 for r in sorted(eligible,key=lambda r:hashlib.sha256(r['id'].encode()).hexdigest()):
  if r['mag']not in seen:seen.add(r['mag']);selected.append(r)
 selected_all+=selected;b=sum(r['rho_any']=='1'and r['ctrl_any']=='0'for r in selected);c=sum(r['rho_any']=='0'and r['ctrl_any']=='1'for r in selected)
 cr[fam]={'n_pairs':len(selected),'rho_positive':sum(r['rho_any']=='1'for r in selected),'control_positive':sum(r['ctrl_any']=='1'for r in selected),'discordant_rho_only':b,'discordant_control_only':c,'p_exact':binomtest(b,b+c,.5).pvalue if b+c else 1}
last=0
for k,(fam,r)in enumerate(sorted(cr.items(),key=lambda kv:kv[1]['p_exact'])):last=max(last,min(1,r['p_exact']*(2-k)));r['p_Holm']=last
write(R/'quality_conditioned_context_loci.tsv',selected_all)
(R/'summary.json').write_text(json.dumps({'job_id':os.environ['SLURM_JOB_ID'],'families':summary,'quality_conditioned_context':cr,'limits':['Post-census reference-site sensitivity, not an independently preregistered validation.','Same retinal-Lys rule applied to every source; exact annotation of all30Type1referencepositions previously calibrated.','Site conservation does not demonstrate retinal binding or sensory function.','Fixed original joint cluster partition retained; no taxonomy-adjusted exceptional-diversity test.','Helio excludes short131aaA0A2C8YND9 reference lacking the reference-aligned site.']},indent=2))
(R/'COMPLETE').write_text(os.environ['SLURM_JOB_ID']+'\n')
