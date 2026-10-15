#!/usr/bin/env python3
"""Reaggregate raw Foldseek rows, preserving retrieval failure and alignment gates."""
import csv,json,hashlib,sys
from pathlib import Path
from collections import defaultdict
base=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
best=defaultdict(lambda:[0.,0.,0.,0.,0]);hashes={}
for name in ['pdb_hits.raw.tsv','afdb_hits.raw.tsv']:
 p=base/'screen'/name;h=hashlib.sha256()
 with p.open('rb') as f:
  for line in f:
   h.update(line);c=line.decode().rstrip().split('\t');q=c[0].removesuffix('.pdb');b=best[q];b[4]+=1;a=float(c[2]);b[0]=max(b[0],float(c[3]));b[1]=max(b[1],a)
   if int(c[13])>=50:b[2]=max(b[2],a)
   if int(c[13])>=100 and float(c[9])>=.5:b[3]=max(b[3],a)
 hashes[str(p)]=h.hexdigest()
records=[];summaries={}
for arm,name in [('dark','per_rep_structures.tsv'),('control','per_control_structures.tsv')]:
 p=base/name;hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();rows=[]
 for r in csv.DictReader(p.open(),delimiter='\t'):
  if r['confidence_qualified']!='1':continue
  ident=r['id'];b=best.get(ident,[0,0,0,0,0]);assert abs(b[0]-max(float(r['best_pdb_qtm']),float(r['best_afdb_qtm'])))<1e-7;assert abs(b[1]-max(float(r['best_pdb_alntm']),float(r['best_afdb_alntm'])))<1e-7
  rows.append(dict(id=ident,arm=arm,length=int(r['length']),plddt=float(r['mean_plddt']),qtm=b[0],alntm=b[1],alntm_min50=b[2],alntm_min100_cov50=b[3],raw_rows=b[4]))
 records+=rows;retrieved=[r for r in rows if r['raw_rows']>0];k=[r for r in rows if r['qtm']<.5]
 summaries[arm]={'qualified':len(rows),'no_rows':len(rows)-len(retrieved),'qtm_below':len(k),'alntm_below':sum(r['alntm']<.5 for r in rows),'retrieved_qtm_below':sum(r['qtm']<.5 for r in retrieved),'retrieved_alntm_below':sum(r['alntm']<.5 for r in retrieved),'min50_alntm_below':sum(r['alntm_min50']<.5 for r in rows),'qtm_below_but_min100_cov50_alntm_hit':sum(r['alntm_min100_cov50']>=.5 for r in k)}
with (out/'per_model.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
(out/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n');(out/'inputs.sha256.json').write_text(json.dumps(hashes,indent=2)+'\n');print(json.dumps(summaries,indent=2))
