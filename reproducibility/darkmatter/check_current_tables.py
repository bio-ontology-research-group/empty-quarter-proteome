#!/usr/bin/env python3
"""Independent compact-table extraction for the current P4 review.

This checks table arithmetic, not upstream sequence searches or biological validity.
"""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from statistics import median

root = Path(__file__).resolve().parent
with (root / 'current_checks_inputs.tsv').open() as f:
    inputs = {x['relative_source']: x for x in csv.DictReader(f, delimiter='\t')}

def rows(path):
    entry = inputs[path]
    compressed = (root / entry['local_path']).read_bytes()
    if hashlib.sha256(compressed).hexdigest() != entry['compressed_sha256']:
        raise ValueError('Compressed checksum mismatch: ' + path)
    data = gzip.decompress(compressed)
    if hashlib.sha256(data).hexdigest() != entry['source_sha256']:
        raise ValueError('Source checksum mismatch: ' + path)
    return list(csv.DictReader(io.StringIO(data.decode()), delimiter='\t'))

out = {}
selection = '06_darkmatter/darkfam_basis_clean/results/selection/'
all_dnds = rows(selection + 'dnds.tsv')
for arm in ('dark', 'ctrl'):
    a = [float(x['median_dNdS']) for x in all_dnds if x['set'] == arm]
    out['selection_' + arm] = {'families': len(a), 'median': median(a), 'fraction_below_one': sum(x < 1 for x in a) / len(a)}
null = {x['family_id']: float(x['median_dNdS']) for x in rows(selection + 'dnds_null_frame1.tsv')}
real = {x['family_id']: float(x['median_dNdS']) for x in rows(selection + 'dnds_dark_subset.tsv')}
common = set(null) & set(real)
out['frame_contrast'] = {'null_n': len(null), 'real_n': len(real), 'paired_n': len(common), 'null_median': median(null.values()), 'real_median': median(real.values()), 'null_higher_fraction': sum(null[k] > real[k] for k in common) / len(common)}
reps = rows('06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv')
thermo = [x for x in reps if x['phylum'] == 'Thermoplasmatota']
out['representatives'] = {'n': len(reps), 'median_dark': median(float(x['dark_fraction']) for x in reps), 'thermoplasmatota_n': len(thermo), 'thermoplasmatota_median_dark': median(float(x['dark_fraction']) for x in thermo)}
for label, file in [('candidate', 'per_rep_structures.tsv'), ('control', 'per_control_structures.tsv')]:
    records = rows('04_structure/00_structure_clean/05_validate/' + file)
    qualified = [x for x in records if float(x['length']) >= 100 and float(x['mean_plddt']) >= 0.70]
    out['structure_' + label] = {'all_n': len(records), 'qualified_n': len(qualified), 'qtm_nohit': sum(max(float(x['best_pdb_qtm']), float(x['best_afdb_qtm'])) < 0.5 for x in qualified), 'alntm_nohit': sum(max(float(x['best_pdb_alntm']), float(x['best_afdb_alntm'])) < 0.5 for x in qualified)}
out['length_table'] = rows('06_darkmatter/cultured_comparison_clean/results/step3_length_summary.tsv')
paired = rows('06_darkmatter/cultured_comparison_clean/results/step3_paired_by_species.tsv')
out['culture_comparison'] = {'pairs': len(paired), 'culture_lower': sum(float(x['delta']) < 0 for x in paired), 'median_delta': median(float(x['delta']) for x in paired)}
(root / 'current_table_checks.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
