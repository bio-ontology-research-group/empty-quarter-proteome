#!/usr/bin/env python3
"""Distinguish terminal and internal table-11 stops in the frame-test inputs."""
import csv
import hashlib
import json
from pathlib import Path
import sys
from collections import Counter

index = Path(sys.argv[1])
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
counts = Counter()
stops = {'TAA', 'TAG', 'TGA'}
with index.open() as source, (out / 'frame_stop_inputs.tsv').open('w') as ledger, (out / 'frame_stop_families.tsv').open('w') as results:
    w = csv.writer(ledger, delimiter='\t'); w.writerow(['family_id', 'path', 'sha256'])
    q = csv.writer(results, delimiter='\t'); q.writerow(['family_id', 'sequences', 'terminal_stop', 'internal_stop', 'eligible_original', 'eligible_masked'])
    for row in csv.DictReader(source, delimiter='\t'):
        p = Path(row['fna'])
        blob = p.read_bytes()
        w.writerow([row['family_id'], str(p), hashlib.sha256(blob).hexdigest()])
        seqs = []
        current = []
        for line in blob.decode().splitlines():
            if line.startswith('>'):
                if current: seqs.append(''.join(current).upper())
                current = []
            else: current.append(line.strip())
        if current: seqs.append(''.join(current).upper())
        c = Counter()
        for nt in seqs:
            c['sequences'] += 1
            nt = nt[:len(nt)//3*3]
            codons = [nt[i:i+3] for i in range(0, len(nt), 3)]
            terminal = bool(codons and codons[-1] in stops)
            internal = any(x in stops for x in codons[:-1])
            c['terminal_stop'] += terminal
            c['internal_stop'] += internal
            if len(nt) >= 90:
                c['eligible_masked'] += 1
                if not internal: c['eligible_original'] += 1
        counts.update(c)
        counts['families'] += 1
        counts['families_with_internal_stop'] += c['internal_stop'] > 0
        q.writerow([row['family_id']] + [c[k] for k in ('sequences', 'terminal_stop', 'internal_stop', 'eligible_original', 'eligible_masked')])
(out / 'frame_stop_summary.json').write_text(json.dumps(dict(counts), indent=2) + '\n')
print(json.dumps(dict(counts), indent=2))
