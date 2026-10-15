#!/usr/bin/env python3
"""Join existing private ENA receipts to the frozen WGS run table; no submissions."""
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from collections import Counter

here = Path(__file__).resolve().parent
repo = here.parents[1]
with (repo / 'code/ena_deposit/wgs_runs.tsv').open() as f:
    runs = list(csv.DictReader(f, delimiter='\t'))
receipts = {}
for folder in ('receipts', 'title_receipts'):
    for p in sorted((here / 'ena_evidence' / folder).glob('*.xml')):
        tree = ET.parse(p).getroot()
        run = tree.find('RUN')
        if tree.get('success') != 'true' or run is None:
            continue
        alias = run.get('alias', '').removeprefix('webin-reads-')
        accession = run.get('accession')
        if alias in receipts and receipts[alias][0] != accession:
            raise ValueError('Conflicting run accession for ' + alias)
        receipts[alias] = (accession, str(p.relative_to(repo)), tree.get('receiptDate'))
with (here / 'ena_run_map.tsv').open('w') as f:
    cols = ['sample_id', 'run_alias', 'biosample_accession', 'instrument_model', 'run_accession', 'receipt_path', 'receipt_date', 'status']
    w = csv.DictWriter(f, fieldnames=cols, delimiter='\t'); w.writeheader()
    for r in runs:
        accession, path, date = receipts.get(r['run_alias'], ('', '', ''))
        w.writerow({**{k: r[k] for k in cols[:4]}, 'run_accession': accession, 'receipt_path': path, 'receipt_date': date, 'status': 'private_submission_receipt' if accession else 'receipt_missing'})
summary = {'runs': len(runs), 'samples': len({r['sample_id'] for r in runs}), 'matched_receipts': sum(r['run_alias'] in receipts for r in runs), 'instrument_counts': dict(Counter(r['instrument_model'] for r in runs)), 'public_run_rows_on_2026_10_13': 0, 'public_status_basis': 'ena_evidence/public_filereport.tsv; private receipts do not establish current public availability'}
(here / 'ena_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
