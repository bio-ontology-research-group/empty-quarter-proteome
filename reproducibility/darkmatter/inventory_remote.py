#!/usr/bin/env python3
"""Inventory P4 inputs without changing them; run under Slurm on IBEX.

Hash files up to 16 MiB. Larger files are explicitly metadata-only pending
claim-level selection for full hashing. Never follow directory symlinks.
"""
import csv
import argparse
import hashlib
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('output_directory')
parser.add_argument('--roots-file')
parser.add_argument('--metadata-only', action='store_true')
parser.add_argument('--name', default='ibex_files.tsv')
args = parser.parse_args()

root = Path('/ibex/scratch/projects/c2014/rob/metagenomics')
roots = [root / p for p in (
    'metadata', 'mimag_rrna_trna_clean', '04_functional/repko_clean',
    '04_protein_catalog/darkmatter_clean', '04_structure/00_structure_clean',
    '04_structure/pfam_candidates_clean', '06_darkmatter',
    '06_reference_bias/refbias_clean', '05_refine/_drep',
    '06_metabolism/uv_survival_clean', '06_metabolism/osmolyte_clean')]
if args.roots_file:
    roots = [Path(x) for x in Path(args.roots_file).read_text().splitlines() if x.strip()]
out = Path(args.output_directory)
out.mkdir(parents=True, exist_ok=True)
with (out / args.name).open('w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['path', 'resolved_path', 'bytes', 'mtime_ns', 'sha256', 'verification', 'error'])
    for base in roots:
        if not base.exists():
            writer.writerow([str(base), '', '', '', '', 'missing_root', ''])
            continue
        for parent, dirs, files in os.walk(base, followlinks=False):
            dirs[:] = [d for d in dirs if d != 'p4_review_20261013']
            dirs.sort()
            for name in sorted(files):
                p = Path(parent) / name
                try:
                    s = p.stat()
                    digest = ''
                    state = 'metadata_only' if args.metadata_only else 'metadata_only_large'
                    if not args.metadata_only and s.st_size <= 16 * 1024**2:
                        h = hashlib.sha256()
                        with p.open('rb') as source:
                            for chunk in iter(lambda: source.read(1024**2), b''):
                                h.update(chunk)
                        digest = h.hexdigest()
                        state = 'sha256'
                    writer.writerow([str(p), str(p.resolve()), s.st_size, s.st_mtime_ns, digest, state, ''])
                except OSError as e:
                    writer.writerow([str(p), '', '', '', '', 'error', str(e)])
            f.flush()
