#!/usr/bin/env python3
"""Part 3 read-only audit: exact stats for SNTC-flagged reps and Thermoplasmatota classes SW-10-69-26 vs E2."""
from pathlib import Path
import json
import gzip
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

ROOT = Path('/home/leechuck/Documents/papers/empty-quarter-mags')
ROUND5 = ROOT / 'reviews/darkmatter-2026-10-13/round5'
OUT = ROUND5 / 'argon_checks'
snap = ROOT / 'reproducibility/darkmatter/source_snapshot'

reps = pd.read_csv(gzip.open(ROOT / 'reproducibility/darkmatter/inputs/06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv.gz'), sep='\t')
contam_hits = pd.read_csv(snap / '06_darkmatter/contamination_clean/results/control_contig_hits.tsv', sep='\t')
tr = pd.read_csv(snap / '06_darkmatter/lineage_darkness_clean/results/v2_thermo_relatives_pairs.tsv', sep='\t')

out = {}

# 1. SNTC read-profiled and contig-hit genomes in v2_reps
sntc_read_genomes = ['V42PRr2__single_SemiBin_33', '51Dr3__single_SemiBin_390']
sntc_contig_genomes = list(contam_hits[contam_hits['control'] == 'SNTC']['genome'].unique())
all_12_sntc = sorted(set(sntc_read_genomes + sntc_contig_genomes))

sub_read = reps[reps['genome'].isin(sntc_read_genomes)]
sub_contig = reps[reps['genome'].isin(sntc_contig_genomes)]
sub_all12 = reps[reps['genome'].isin(all_12_sntc)]

out['sntc_contamination_in_reps'] = {
    'read_profiled_2': sub_read[['genome', 'phylum', 'genus', 'n_proteins', 'n_dark', 'dark_fraction']].to_dict(orient='records'),
    'read_profiled_2_total_proteins': int(sub_read['n_proteins'].sum()),
    'read_profiled_2_total_dark': int(sub_read['n_dark'].sum()),
    'contig_hit_10_n': len(sub_contig),
    'contig_hit_10_total_proteins': int(sub_contig['n_proteins'].sum()),
    'contig_hit_10_total_dark': int(sub_contig['n_dark'].sum()),
    'contig_hit_10_records': sub_contig[['genome', 'phylum', 'genus', 'n_proteins', 'n_dark', 'dark_fraction']].to_dict(orient='records'),
    'all_12_total_proteins': int(sub_all12['n_proteins'].sum()),
    'all_12_total_dark': int(sub_all12['n_dark'].sum()),
    'catalogue_total_proteins_in_reps': int(reps['n_proteins'].sum()),
    'catalogue_total_dark_in_reps': int(reps['n_dark'].sum()),
}

# 2. Thermoplasmatota SW-10-69-26 (n=185) vs E2 (n=51) vs Halobacteriota (n=37) vs Thermoproteota (n=70)
tr_m = tr.merge(reps, on='genome', how='inner')
thermo_classes = {}
for cls_name, g in tr_m.groupby('cls'):
    thermo_classes[cls_name] = {
        'n': len(g),
        'median_length_mb': float(g['length'].median() / 1e6),
        'iqr_length_mb': [float(g['length'].quantile(0.25) / 1e6), float(g['length'].quantile(0.75) / 1e6)],
        'median_n_proteins': float(g['n_proteins'].median()),
        'median_n_dark': float(g['n_dark'].median()),
        'median_dark_fraction': float(g['dark_fraction'].median()),
        'iqr_dark_fraction': [float(g['dark_fraction'].quantile(0.25)), float(g['dark_fraction'].quantile(0.75))],
        'median_completeness': float(g['completeness'].median()),
        'median_contamination': float(g['contamination'].median()),
        'median_N50': float(g['N50'].median()),
        'median_best_external_aai': float(g['aai'].median()),
        'median_best_desert_aai': float(g['aai_desert'].median()),
        'median_rbh_frac': float(g['aai_frac_rbh'].median()),
        'spearman_ext_aai_vs_dark': [float(x) for x in spearmanr(g['aai'], g['dark_fraction'])],
        'spearman_desert_aai_vs_dark': [float(x) for x in spearmanr(g['aai_desert'], g['dark_fraction'])],
        'spearman_rbh_frac_vs_dark': [float(x) for x in spearmanr(g['aai_frac_rbh'], g['dark_fraction'])],
    }

sw_g = tr_m[tr_m['cls'] == 'SW-10-69-26']
e2_g = tr_m[tr_m['cls'] == 'E2']
thermo_classes['sw_vs_e2_mannwhitney'] = {
    'dark_fraction_p': float(mannwhitneyu(sw_g['dark_fraction'], e2_g['dark_fraction']).pvalue),
    'length_p': float(mannwhitneyu(sw_g['length'], e2_g['length']).pvalue),
    'n_proteins_p': float(mannwhitneyu(sw_g['n_proteins'], e2_g['n_proteins']).pvalue),
    'aai_p': float(mannwhitneyu(sw_g['aai'], e2_g['aai']).pvalue),
    'rbh_frac_p': float(mannwhitneyu(sw_g['aai_frac_rbh'], e2_g['aai_frac_rbh']).pvalue),
}
out['thermo_classes'] = thermo_classes

(OUT / 'audit_part3.json').write_text(json.dumps(out, indent=2) + '\n')
print("Saved clean audit_part3.json")
