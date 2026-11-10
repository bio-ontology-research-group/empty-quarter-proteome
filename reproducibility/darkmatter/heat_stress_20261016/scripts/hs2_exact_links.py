from common import *
import pandas as pd,gzip
from concurrent.futures import ProcessPoolExecutor
R=N/'results';O=R/'hs2_exact';O.mkdir(exist_ok=True)
def records(path):
 with gzip.open(path,'rb')as f:
  header=None;seq=[]
  for line in f:
   if line.startswith(b'>'):
    if header is not None:yield header,b''.join(seq).rstrip(b'*')
    header=line[1:].split()[0].decode();seq=[]
   else:seq.append(line.strip())
  if header is not None:yield header,b''.join(seq).rstrip(b'*')
def worker(job):
 sample,proteomes=job;targets={};rows={}
 for genome,path in proteomes:
  for pid,seq in records(path):
   raw=pid.split('::')[-1];contig=raw.rsplit('_',1)[0];digest=hashlib.sha256(seq).hexdigest();targets.setdefault((contig,digest),[]).append(pid);rows[pid]={'protein':pid,'genome':genome,'sample':sample,'name_based_catalog_id':sample+'::'+raw,'sequence_sha256':digest,'aa_length':len(seq)}
 contigs={key[0]for key in targets};matches={};path=W/f'04_protein_catalog/assembly_proteins/{sample}/{sample}.faa.gz';before=path.stat()
 def accept(pid,seq):
  digest=hashlib.sha256(b''.join(seq).rstrip(b'*')).hexdigest();contig=pid.rsplit('_',1)[0]
  for target in targets.get((contig,digest),[]):matches.setdefault(target,[]).append(sample+'::'+pid)
 with gzip.open(path,'rb')as f:
  pid=None;seq=[];wanted=False
  for line in f:
   if line.startswith(b'>'):
    if wanted:accept(pid,seq)
    pid=line[1:].split()[0].decode().split('::')[-1];wanted=pid.rsplit('_',1)[0]in contigs;seq=[]
   elif wanted:seq.append(line.strip())
  if wanted:accept(pid,seq)
 for pid,r in rows.items():
  candidates=sorted(set(matches.get(pid,[])));r['catalog_id']=candidates[0]if len(candidates)==1 else '';r['status']='exact_unique'if len(candidates)==1 else 'ambiguous_exact'if candidates else 'no_exact_sequence';r['identifier_changed']=bool(r['catalog_id'] and r['catalog_id']!=r['name_based_catalog_id'])
 digest=sha(path);after=path.stat();assert(before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
 return list(rows.values()),{'path':str(path),'bytes':after.st_size,'sha256':digest}
if __name__=='__main__':
 idx=pd.read_csv(R/'genome_index.tsv',sep='\t');idx=idx[idx.gtdb_classification.str.contains('p__Thermoplasmatota;')];jobs=[(sample,list(zip(g.genome_id,g.proteome)))for sample,g in idx.groupby('source_sample')];rows=[];sources=[]
 with ProcessPoolExecutor(max_workers=8)as pool:
  for i,(rr,src)in enumerate(pool.map(worker,jobs)):
   rows.extend(rr);sources.append(src);print('sources_done',i+1,'of',len(jobs),'genes',len(rows),flush=True)
 d=pd.DataFrame(rows);assert d.protein.is_unique;d.to_csv(O/'exact_locus_links.tsv.gz',sep='\t',index=False);write(O/'original_proteome_sources.tsv',sources);summary={'job':os.environ['SLURM_JOB_ID'],'source_samples':len(jobs),'loci':len(d),'status':d.status.value_counts().to_dict(),'changed_identifiers':int(d.identifier_changed.sum()),'matching':'identical amino-acid sequence after terminal-stop removal, same original contig, unique original locus; no fuzzy gene-number matching'};(O/'summary.json').write_text(json.dumps(summary,indent=2));print(summary)
