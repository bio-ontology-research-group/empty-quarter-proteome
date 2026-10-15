#!/usr/bin/env python3
"""Extract family CDS covariates on IBEX; run under Slurm, never on login."""
import csv, hashlib, json, statistics, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
base=Path(sys.argv[1]); out=Path(sys.argv[2]); out.mkdir(parents=True,exist_ok=True)
def extract(row):
 p=Path(row['fna']); blob=p.read_bytes(); seqs=[]; current=[]
 for line in blob.decode().splitlines():
  if line.startswith('>'):
   if current:seqs.append(''.join(current).upper())
   current=[]
  else:current.append(line.strip())
 if current:seqs.append(''.join(current).upper())
 valid=[s for s in seqs if len(s)>=90]
 gc=[sum(s.count(x) for x in 'GC')/sum(s.count(x) for x in 'ACGT') for s in valid if sum(s.count(x) for x in 'ACGT')]
 return [row['family_id'],row['set'],len(seqs),statistics.median([len(s)/3 for s in valid]) if valid else '',statistics.median(gc) if gc else '',str(p),hashlib.sha256(blob).hexdigest()]
rows=[]
for arm in ['dark','ctrl']:
 with (base/f'family_index_{arm}.tsv').open() as f: rows+=list(csv.DictReader(f,delimiter='\t'))
with (out/'family_covariates.tsv').open('w') as f,ProcessPoolExecutor(max_workers=8) as ex:
 w=csv.writer(f,delimiter='\t');w.writerow(['family_id','set','input_sequences','median_cds_aa','median_gc','fna','sha256'])
 for row in ex.map(extract,rows,chunksize=200):w.writerow(row)
print(json.dumps({'families':len(rows),'output':str(out)}))
