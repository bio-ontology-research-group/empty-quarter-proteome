from common import *
import pandas as pd,numpy as np,gzip
from scipy.stats import false_discovery_control
R=N/'results';O=R/'hs2';d=pd.read_csv(O/'thermo_genes.tsv.gz',sep='\t',keep_default_na=False);d['partial']=d.partial.astype(str).str.zfill(2);d=d.sort_values(['genome','contig','start','protein']);d['pos']=d.groupby(['genome','contig']).cumcount();d['length']=d.end-d.start+1
loci=pd.read_csv(O/'tested_loci.tsv.gz',sep='\t');groups=d.groupby(['genome','contig']);pr=[]
for key,g in loci.groupby(['genome','contig']):
 b=groups.get_group(key);n=len(b);anchor=b.anchor.to_numpy(bool);near=np.zeros(n,bool)
 for j in np.where(anchor)[0]:near[max(0,j-2):min(n,j+3)]=True
 for r in g.itertuples():
  target=b[b.protein==r.protein].iloc[0];eligible=(b.partial=='00')&(b.anchor==0)&(b.pos>=2)&(b.pos<n-2)&b.length.between(target.length/2,target.length*2)
  assert eligible.sum()>0
  pr.append({'family':r.family,'genome':r.genome,'protein':r.protein,'near':r.near_anchor,'p_length_matched':near[eligible].mean(),'genus':';'.join(r.taxonomy.split(';')[:6]),'order':';'.join(r.taxonomy.split(';')[:4])})
p=pd.DataFrame(pr);rows=[]
for mode in ['species_length_matched','one_representative_per_genus_label_stratum']:
 dd=p.sort_values('genome')
 if mode=='one_representative_per_genus_label_stratum':dd=dd.drop_duplicates(['family','genus'])
 for family,g in dd.groupby('family'):
  dist=np.array([1.])
  for prob in g.p_length_matched:dist=np.convolve(dist,[1-prob,prob])
  obs=int(g.near.sum());rows.append({'mode':mode,'family':family,'units':len(g),'near_units':obs,'expected':g.p_length_matched.sum(),'p':min(1.,dist[obs:].sum())})
tests=pd.DataFrame(rows);tests['q']=np.nan
for mode,ix in tests.groupby('mode').groups.items():tests.loc[ix,'q']=false_discovery_control(tests.loc[ix,'p'])
tests.to_csv(O/'neighbourhood_sensitivity.tsv',sep='\t',index=False)
exact=pd.read_csv(R/'audit/HS2_exact_tail.tsv',sep='\t');chosen=set(exact.loc[exact.q_exact<.05,'family']);
if not chosen:
 print('No exact-null discoveries; sensitivity table retained');raise SystemExit(0)
cand=d[d.family.isin(chosen)];cand.to_csv(O/'candidate_all_loci.tsv',sep='\t',index=False);targets=set(cand.protein);hits=[];sources=[]
for genome in cand.genome.unique():
 path=T/'functional/pfam_full/domtbl'/f'{genome}.domtbl.gz';assert path.exists(),path;sources.append({'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size})
 with gzip.open(path,'rt')as f:
  for line in f:
   if line.startswith('#')or not line.strip():continue
   z=line.split();assert '::'in z[0],z[:6]
   if z[0]in targets:hits.append({'protein':z[0],'pfam':z[4],'name':z[3],'domain_evalue':float(z[12]),'hmm_coverage':(int(z[16])-int(z[15])+1)/int(z[5]),'description':' '.join(z[22:])})
pd.DataFrame(hits,columns=['protein','pfam','name','domain_evalue','hmm_coverage','description']).to_csv(O/'candidate_pfam_hits.tsv',sep='\t',index=False);write(O/'candidate_annotation_sources.tsv',sources)
markers=pd.read_csv(R/'marker_loci.tsv.gz',sep='\t');anchor_fam=markers.groupby('protein').family.agg(lambda v:';'.join(sorted(set(v)))).to_dict();info=[]
for r in loci[loci.family.isin(chosen)].itertuples():
 b=groups.get_group((r.genome,r.contig));neigh=b[(b.pos-r.pos).abs()<=2];neigh=neigh[neigh.anchor==1]
 info.append({'family':r.family,'genome':r.genome,'protein':r.protein,'genus':';'.join(r.taxonomy.split(';')[:6]),'order':';'.join(r.taxonomy.split(';')[:4]),'near':r.near_anchor,'anchor_families':';'.join(sorted({anchor_fam[x]for x in neigh.protein}))})
pd.DataFrame(info).to_csv(O/'candidate_context.tsv',sep='\t',index=False)
print(tests[tests.family.isin(chosen)].to_string(index=False));print('DONE HS2 refine')
