#!/usr/bin/env python3
"""Fail closed if a frozen input or promoted upstream script has changed."""
import csv,hashlib
from pathlib import Path
root=Path(__file__).resolve().parent
n=0
for manifest,column in [('review_inputs.sha256.tsv','path'),('upstream_code_manifest.tsv','local_path')]:
 for r in csv.DictReader((root/manifest).open(),delimiter='\t'):
  name=r[column]
  p=root/name
  if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise SystemExit('Checksum mismatch: '+str(p))
  n+=1
print('Verified',n,'frozen inputs and upstream scripts')
