#!/usr/bin/env python3
"""Curator-only: freeze retained review inputs after provenance inspection.
Do not run during validation: verify_review_inputs.py consumes this fixed ledger.
"""
import csv,hashlib
from pathlib import Path
H=Path(__file__).resolve().parent
paths=list((H/'inputs').rglob('*.gz'))
paths += [H/'selection_covariates/family_covariates.tsv.gz']
for folder in ['environment_inputs','culture_nucleotide_cigar','culture_locus','culture_locus_all','culture_family_membership','structural_raw_audit','shadow_overlap','antifam','environments']:
 for p in (H/folder).rglob('*'):
  if not p.is_file() or p.suffix not in ['.tsv','.json','.txt','.sha256','.tbl','.domtbl','.hmm','.gz','.py']:continue
  if p.name in ['selection_sensitivity.json','posthoc_recovery_subset.json','per_pair_positive_control.tsv']:continue
  paths.append(p)
with (H/'review_inputs.sha256.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256'])
 for p in sorted(set(paths)):w.writerow([str(p.relative_to(H)),p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()])
print('Pinned',len(set(paths)),'review inputs')
