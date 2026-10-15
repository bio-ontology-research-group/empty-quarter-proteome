#!/usr/bin/env python3
"""Trace eligible culture-comparison dark calls through cluster90 to admitted families."""
import csv,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');R=W/'06_darkmatter/p4_review_20261013/results';O=R/'culture_family_membership';O.mkdir(exist_ok=True)
calls=[r for r in csv.DictReader((R/'culture_nucleotide_cigar/per_gene.tsv').open(),delimiter='\t') if r['dark']=='1' and r['eligible']=='1'];want={r['protein'].encode() for r in calls};p=W/'04_protein_catalog/all_linclust90/all_cluster90.tsv'
def scan(bounds):
 start,end=bounds;out={}
 with p.open('rb') as f:
  f.seek(start)
  if start:f.readline()
  pos=f.tell()
  while pos<end:
   line=f.readline()
   if not line:break
   pos+=len(line);rep,_,gene=line.rstrip(b'\n').partition(b'\t')
   if gene in want:out[gene.decode()]=rep.decode()
 return out
mapping={};size=p.stat().st_size
with ProcessPoolExecutor(max_workers=4) as ex:
 for found in ex.map(scan,[(size*i//4,size*(i+1)//4) for i in range(4)]):mapping.update(found)
assert len(mapping)==len(want),(len(mapping),len(want))
famroot=W/'06_darkmatter/darkfam_basis_clean/results';admitted=set((famroot/'FAMILIES_FINAL.txt').read_text().split());targets=set(mapping.values());fm={}
for line in (famroot/'fam_cluster.tsv').open():
 fam,member=line.rstrip().split('\t')
 if member in targets and fam in admitted:fm[member]=fam
for r in calls:r.update(cluster90=mapping[r['protein']],admitted_family=fm.get(mapping[r['protein']],''))
with (O/'per_gene.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(calls[0]),delimiter='\t');w.writeheader();w.writerows(calls)
summary={'eligible_dark_calls':len(calls),'admitted_family_members':sum(bool(r['admitted_family']) for r in calls),'distinct_admitted_families':len({r['admitted_family'] for r in calls if r['admitted_family']})};(O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
