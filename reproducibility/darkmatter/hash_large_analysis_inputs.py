#!/usr/bin/env python3
"""Complete SHA256 for large files in the clean-analysis inventory, deduplicated.
Read-only Slurm task; records files that change during hashing as errors.
"""
import csv,hashlib,json,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
R=Path('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013');out=R/'results/large_analysis_checksums';out.mkdir(exist_ok=True)
rows={}
for r in csv.DictReader((R/'results/ibex_files.tsv').open(),delimiter='\t'):
 if r['verification']=='metadata_only_large':rows.setdefault(r['resolved_path'],r)
def one(item):
 path,r=item;p=Path(path);v={'path':path,'inventory_path':r['path'],'bytes':r['bytes'],'sha256':'','error':''}
 try:
  a=p.stat();h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
  z=p.stat();assert (a.st_size,a.st_mtime_ns)==(z.st_size,z.st_mtime_ns),'changed during read'
  assert a.st_size==int(r['bytes']),'changed since inventory'
  v['sha256']=h.hexdigest()
 except Exception as e:v['error']=str(e)
 return v
with (out/'sha256.tsv').open('w') as f,ThreadPoolExecutor(max_workers=4) as ex:
 w=csv.DictWriter(f,fieldnames=['path','inventory_path','bytes','sha256','error'],delimiter='\t');w.writeheader()
 for i,r in enumerate(ex.map(one,rows.items())):w.writerow(r);f.flush()
print('Completed',len(rows),'unique large analysis paths')
