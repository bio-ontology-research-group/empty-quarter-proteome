#!/usr/bin/env python3
"""Focused round-1 reaggregation of frozen tables; no new sequence searches."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from statistics import median
from scipy.stats import wilcoxon

here = Path(__file__).resolve().parent
with (here / 'current_checks_inputs.tsv').open() as f:
    inputs = {r['relative_source']: r for r in csv.DictReader(f, delimiter='\t')}

def read(path):
    r = inputs[path]
    b = gzip.decompress((here / r['local_path']).read_bytes())
    assert hashlib.sha256(b).hexdigest() == r['source_sha256']
    return list(csv.DictReader(io.StringIO(b.decode()), delimiter='\t'))

out = {}
ecoli = 'V42PRr2__single_SemiBin_33'
rows = read('06_darkmatter/cultured_comparison_clean/results/pair_anatomy/pair_anatomy_per_pair.tsv')
for label, data in [('all_pairs', rows), ('exclude_confirmed_ecoli', [r for r in rows if r['mag_rep'] != ecoli])]:
    sums = {k: sum(int(r[k]) for r in data) for k in ('mag_annot', 'mag_annot_shared', 'mag_dark', 'mag_dark_shared', 'mag_dark_full', 'mag_dark_full_shared')}
    out[label] = {'pairs': len(data), **sums, 'annotated_share': sums['mag_annot_shared']/sums['mag_annot'], 'dark_share': sums['mag_dark_shared']/sums['mag_dark'], 'median_shared_identity': median(float(r['shared_median_identity']) for r in data), 'pairs_annotated_more_shared': sum(float(r['mag_annot_with_partner_pct']) > float(r['mag_dark_with_partner_pct']) for r in data)}
rates = [r for r in read('06_darkmatter/cultured_comparison_clean/results/step3_paired_by_species.tsv') if r['mag_rep'] != ecoli]
delta = [float(r['delta']) for r in rates]
out['culture_dark_rate_exclude_ecoli'] = {'pairs': len(rates), 'culture_lower': sum(x < 0 for x in delta), 'median_delta': median(delta), 'wilcoxon_two_sided_p': float(wilcoxon(delta).pvalue)}
selection = read('06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv')
strata = []
for axis, bins in [('mean_dS', [(0,.1),(.1,.25),(.25,.5),(.5,1),(1,2)]), ('mean_pid', [(0,.7),(.7,.8),(.8,.9),(.9,.95),(.95,1)])]:
    for lo, hi in bins:
        for arm in ('dark','ctrl'):
            x = [float(r['median_dNdS']) for r in selection if r['set'] == arm and lo < float(r[axis]) <= hi]
            strata.append({'axis':axis, 'lower_exclusive':lo, 'upper_inclusive':hi, 'arm':arm, 'families':len(x), 'median_dnds':median(x) if x else None})
out['divergence_strata'] = strata
reps = read('06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv')
t = [float(r['dark_fraction']) for r in reps if r['phylum']=='Thermoplasmatota' and r['first_unnamed_rank']=='2']
c = [float(r['dark_fraction']) for r in reps if r['phylum']!='Thermoplasmatota' and r['first_unnamed_rank']=='2']
out['class_placeholder_contrast'] = {'thermo_n':len(t), 'other_n':len(c), 'thermo_median':median(t), 'other_median':median(c), 'ratio':median(t)/median(c)}
(here / 'round1_sensitivity.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='divergence_strata'},indent=2))
