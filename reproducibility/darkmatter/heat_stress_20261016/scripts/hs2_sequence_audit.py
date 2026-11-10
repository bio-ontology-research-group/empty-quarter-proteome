from common import *
import pandas as pd,gzip
from concurrent.futures import ThreadPoolExecutor
R=N/'results';O=R/'hs2';d=pd.read_csv(O/'candidate_all_loci.tsv',sep='\t',keep_default_na=False);idx=pd.read_csv(R/'genome_index.tsv',sep='\t').set_index('genome_id');seqs={}
def fasta(path):
 op=gzip.open if str(path).endswith('.gz')else open
 with op(path,'rt')as f:
  h=None;s=[]
  for l in f:
   if l.startswith('>'):
    if h is not None:yield h,''.join(s).rstrip('*')
    h=l[1:].strip();s=[]
   else:s.append(l.strip())
  if h is not None:yield h,''.join(s).rstrip('*')
for genome,g in d.groupby('genome'):
 targets=set(g.protein)
 for h,s in fasta(idx.loc[genome,'proteome']):
  pid=h.split()[0]
  if pid in targets:seqs[pid]=s
jobs=[]
for sample,g in d.groupby(d.catalog_id.str.split('::').str[0]):jobs.append((sample,g))
def worker(job):
 sample,g=job;wanted={r.catalog_id.split('::')[1]:r for r in g.itertuples()};found={};path=W/f'04_protein_catalog/assembly_proteins/{sample}/{sample}.faa.gz'
 for h,s in fasta(path):
  pid=h.split()[0].split('::')[-1]
  if pid in wanted:found[pid]=(h,s)
 rows=[]
 for pid,r in wanted.items():
  h,s=found.get(pid,('', ''));current=seqs[r.protein];rows.append({'protein':r.protein,'family':r.family,'catalog_id':r.catalog_id,'current_aa':len(current),'original_aa':len(s),'exact_sequence_match':current==s,'current_sha256':hashlib.sha256(current.encode()).hexdigest(),'original_sha256':hashlib.sha256(s.encode()).hexdigest(),'original_header':h,'source':str(path)})
 return rows
rows=[]
with ThreadPoolExecutor(max_workers=4)as pool:
 for res in pool.map(worker,jobs):rows.extend(res)
pd.DataFrame(rows).to_csv(O/'sequence_identity_audit.tsv',sep='\t',index=False);print(pd.DataFrame(rows).groupby('family').exact_sequence_match.agg(['sum','count']).to_string());print('DONE sequence audit')
