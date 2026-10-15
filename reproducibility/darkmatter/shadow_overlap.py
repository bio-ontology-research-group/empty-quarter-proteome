#!/usr/bin/env python3
"""Flag >=50% opposite-strand overlap of admitted dark members with non-dark calls.

All coordinates are original Pyrodigal protein FASTA headers. A flag is a candidate overlap,
not proof of a spurious ORF. Non-dark status is resolved through cluster90 and
the complete sequence-dark ID set, not inferred from absence in the MAG subset.
"""
import csv,gzip,hashlib,json,sys,os
from pathlib import Path
from collections import defaultdict,Counter
from concurrent.futures import ProcessPoolExecutor
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');R=W/'06_darkmatter/darkfam_basis_clean/results';out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
fams=set((R/'FAMILIES_FINAL.txt').read_text().split());member={}
for line in (R/'fam_cluster.tsv').open():
 f,m=line.rstrip().split('\t')
 if f in fams:member[m]=f
assert len(member)==536437,len(member)
queries=defaultdict(lambda:defaultdict(list));seen=0
for line in (R/'dark_reps.faa').open():
 if not line.startswith('>'):continue
 parts=line[1:].strip().split(' # ');pid=parts[0]
 if pid not in member:continue
 sample=pid.split('::')[0];contig=pid.split('::')[-1].rsplit('_',1)[0];a,b,strand=map(int,parts[1:4]);queries[sample][contig].append((pid,member[pid],a,b,strand));seen+=1
assert seen==len(member),(seen,len(member))
def process_sample(args):
 sample,query=args;p=W/f'04_protein_catalog/assembly_proteins/{sample}/{sample}.faa.gz';genes=defaultdict(list);count=0
 if not p.exists():return [],{'sample':sample,'error':'missing protein FASTA','query_members':sum(map(len,query.values()))}
 try:
  with gzip.open(p,'rt') as f:
   for line in f:
    if not line.startswith('>'):continue
    parts=line[1:].strip().split(' # ');pid=parts[0].split()[0]
    if '::' not in pid:pid=sample+'::'+pid
    contig=pid.split('::')[-1].rsplit('_',1)[0]
    if contig not in query:continue
    a,b,strand=map(int,parts[1:4]);genes[contig].append((pid,a,b,strand));count+=1
 except (EOFError,OSError) as e:
  return [],{'sample':sample,'error':str(e),'query_members':sum(map(len,query.values()))}
 pairs=[];matched=0
 for contig,qs in query.items():
  gs=genes[contig];keys={(pid,a,b,s) for pid,a,b,s in gs}
  for pid,fam,a,b,s in qs:
   if (pid,a,b,s) in keys:matched+=1
   for other,c,d,t in gs:
    if s==t:continue
    overlap=max(0,min(b,d)-max(a,c)+1)
    if overlap/(b-a+1)>=.5:pairs.append([fam,pid,other,overlap,b-a+1,d-c+1])
 h=hashlib.sha256()
 with p.open('rb') as f:
  for blob in iter(lambda:f.read(8*1024*1024),b''):h.update(blob)
 return pairs,{'sample':sample,'path':str(p),'sha256':h.hexdigest(),'query_members':sum(map(len,query.values())),'exact_coordinate_matches':matched,'overlap_pairs':len(pairs),'error':''}
allpairs=[];checks=[]
with ProcessPoolExecutor(max_workers=12) as ex:
 for pairs,check in ex.map(process_sample,[(s,dict(q)) for s,q in queries.items()]):allpairs+=pairs;checks.append(check);print(check['sample'],check.get('exact_coordinate_matches'),len(pairs),flush=True)
(out/'gff_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
assert sum(x.get('exact_coordinate_matches',0) for x in checks)==len(member),'Incomplete coordinate join; inspect gff_checks.json'
WANT={p[2].encode() for p in allpairs};cluster=W/'04_protein_catalog/all_linclust90/all_cluster90.tsv'
def lookup(bounds):
 start,end=bounds;hits={}
 with cluster.open('rb') as f:
  f.seek(start)
  if start:f.readline()
  pos=f.tell()
  while pos<end:
   line=f.readline()
   if not line:break
   pos+=len(line);rep,_,gene=line.rstrip(b'\n').partition(b'\t')
   if gene in WANT:hits[gene.decode()]=rep.decode()
 return hits
mapping={}
if WANT:
 size=cluster.stat().st_size
 with ProcessPoolExecutor(max_workers=12) as ex:
  for found in ex.map(lookup,[(size*i//12,size*(i+1)//12) for i in range(12)]):mapping.update(found)
 assert len(mapping)==len(WANT),(len(mapping),len(WANT))
 targets=set(mapping.values());dark=set();darkfile=W/'04_protein_catalog/darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz'
 with gzip.open(darkfile,'rt') as f:
  for line in f:
   ident=line.strip()
   if ident in targets:dark.add(ident)
rows=[]
for row in allpairs:rows.append(row+[mapping[row[2]],int(mapping[row[2]] in dark)])
with (out/'overlaps.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['family_id','dark_member','opposite_gene','overlap_nt','dark_nt','opposite_nt','opposite_cluster90','opposite_is_dark']);w.writerows(rows)
flagged={r[1] for r in rows if r[-1]==0};flagfams={r[0] for r in rows if r[-1]==0}
with (out/'family_flags.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['family_id','members','opposite_nondark_overlap_members']);counts=Counter(member.values());flags=Counter(member[m] for m in flagged)
 for fam in sorted(fams):w.writerow([fam,counts[fam],flags[fam]])
summary={'admitted_families':len(fams),'admitted_members':len(member),'overlap_pairs_any':len(rows),'members_opposite_nondark_overlap':len(flagged),'families_with_opposite_nondark_overlap':len(flagfams)}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
