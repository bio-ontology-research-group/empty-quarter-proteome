#!/usr/bin/env python3
"""Recompute structural threshold sensitivity and frame attrition from frozen inputs."""
import csv,gzip,hashlib,io,json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parent
manifest={r['relative_source']:r for r in csv.DictReader((root/'current_checks_inputs.tsv').open(),delimiter='\t')}
def read(path):
 m=manifest[path]; b=gzip.decompress((root/m['local_path']).read_bytes())
 if hashlib.sha256(b).hexdigest()!=m['source_sha256']:raise ValueError(path)
 return list(csv.DictReader(io.StringIO(b.decode()),delimiter='\t'))
def main():
 out={'structural':{},'frame':{}}
 for arm,fn in [('dark','per_rep_structures.tsv'),('control','per_control_structures.tsv')]:
  rows=[r for r in read('04_structure/00_structure_clean/05_validate/'+fn) if r['confidence_qualified']=='1']
  q=np.array([max(float(r['best_pdb_qtm']),float(r['best_afdb_qtm'])) for r in rows]);a=np.array([max(float(r['best_pdb_alntm']),float(r['best_afdb_alntm'])) for r in rows])
  zero=(q==0)&(a==0)
  out['structural'][arm]={'n':len(rows),'zero_reported_scores':int(zero.sum()),'thresholds':{str(t):{'qtm_below':int((q<t).sum()),'alntm_below':int((a<t).sum()),'qtm_below_with_score':int(((q<t)&~zero).sum()),'alntm_below_with_score':int(((a<t)&~zero).sum())} for t in [.4,.45,.5,.55,.6]},'qtm_045_to_05':int(((q>=.45)&(q<.5)).sum())}
 base='06_darkmatter/darkfam_basis_clean/results/selection/'
 n={r['family_id']:r for r in read(base+'dnds_null_frame1.tsv')};r={r['family_id']:r for r in read(base+'dnds_dark_subset.tsv')};both=sorted(n.keys()&r.keys())
 out['frame']={'input_families':10000,'null':len(n),'annotated':len(r),'both':len(both),'null_only':len(n.keys()-r.keys()),'annotated_only':len(r.keys()-n.keys()),'neither':10000-len(n.keys()|r.keys()),'paired_median_difference':float(np.median([float(n[k]['median_dNdS'])-float(r[k]['median_dNdS']) for k in both]))}
 (root/'round2_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
