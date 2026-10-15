#!/usr/bin/env python3
"""Read-only audit of baseline_main.tex, baseline_references.bib, compact tables, and round5_extensions."""
from pathlib import Path
import re
import json
import gzip
import csv
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, spearmanr, mannwhitneyu, fisher_exact

ROOT = Path('/home/leechuck/Documents/papers/empty-quarter-mags')
ROUND5 = ROOT / 'reviews/darkmatter-2026-10-13/round5'
OUT = ROUND5 / 'argon_checks'
OUT.mkdir(parents=True, exist_ok=True)

tex_path = ROUND5 / 'baseline_main.tex'
bib_path = ROUND5 / 'baseline_references.bib'
tex = tex_path.read_text()
lines = tex.splitlines()

report = {}

# 1. Check cross-references, tables, figures, supplementary references in baseline_main.tex
table_mentions = []
for idx, line in enumerate(lines, 1):
    if re.search(r'(Table|Fig|Supplementary|Extended Data|Note)', line, re.I):
        table_mentions.append((idx, line.strip()))

# Extract all \cite* keys
cite_keys = set()
for m in re.finditer(r'\\cite[a-zA-Z]*\{([^}]+)\}', tex):
    for k in m.group(1).split(','):
        cite_keys.add(k.strip())

bib_text = bib_path.read_text()
bib_entries = {}
for m in re.finditer(r'@\w+\{([^,]+),([\s\S]*?)\n\}', bib_text):
    key = m.group(1).strip()
    body = m.group(2)
    bib_entries[key] = body

missing_bib = sorted(cite_keys - set(bib_entries.keys()))
cited_bib_details = {}
for k in sorted(cite_keys):
    body = bib_entries.get(k, '')
    has_vol = bool(re.search(r'\bvolume\s*=', body, re.I))
    has_pages = bool(re.search(r'\bpages\s*=', body, re.I))
    has_note = bool(re.search(r'\bnote\s*=', body, re.I))
    has_others = 'and others' in body
    if not has_vol or not has_pages or has_note or has_others:
        cited_bib_details[k] = {
            'has_volume': has_vol,
            'has_pages': has_pages,
            'has_note': has_note,
            'has_others': has_others,
            'snippet': ' '.join(body.split())[:250]
        }

report['tex_line_count'] = len(lines)
report['cited_keys_count'] = len(cite_keys)
report['missing_bib_keys'] = missing_bib
report['cited_bib_incomplete'] = cited_bib_details

# 2. Audit H1 (culture_caller)
ext_dir = ROOT / 'reproducibility/darkmatter/round5_extensions/results'
cc_dir = ext_dir / 'culture_caller'
cc_df = pd.read_csv(cc_dir / 'per_gene.tsv.gz', sep='\t')
cc_sum = pd.read_csv(cc_dir / 'summary.tsv', sep='\t')
cc_pair = pd.read_csv(cc_dir / 'per_pair.tsv', sep='\t')

# Also join with culture_nucleotide_cigar/per_gene.tsv to see protein_shared status!
cigar_df = pd.read_csv(ROOT / 'reproducibility/darkmatter/culture_nucleotide_cigar/per_gene.tsv', sep='\t')
cc_merged = cc_df.merge(
    cigar_df[['rep', 'protein', 'protein_shared']],
    on=['rep', 'protein'],
    how='left',
    validate='many_to_one'
)

h1_audit = {}
for thresh in [0.5, 0.9]:
    sub = cc_merged[cc_merged['overlap_threshold'] == thresh]
    h1_audit[str(thresh)] = {}
    for dark_val, dark_label in [(1, 'dark'), (0, 'annotated')]:
        for sh_val, sh_label in [('all', [0, 1]), ('unshared_only', [0]), ('shared_only', [1])]:
            g = sub[(sub['dark'] == dark_val) & (sub['protein_shared'].isin(sh_label))]
            h1_audit[str(thresh)][f'{dark_label}_{sh_label}'] = {
                'n': len(g),
                'bakta_no_cds': int((g['bakta'] == 'no_major_CDS_overlap').sum()),
                'pyrodigal_no_cds': int((g['pyrodigal'] == 'no_major_CDS_overlap').sum()),
                'bakta_same_phase': int((g['bakta'] == 'same_orientation_start_phase').sum()),
                'pyrodigal_same_phase': int((g['pyrodigal'] == 'same_orientation_start_phase').sum()),
                'bakta_other_frame': int((g['bakta'] == 'other_frame_or_orientation').sum()),
                'pyrodigal_other_frame': int((g['pyrodigal'] == 'other_frame_or_orientation').sum()),
                'gained_cds': int(((g['bakta'] == 'no_major_CDS_overlap') & (g['pyrodigal'] != 'no_major_CDS_overlap')).sum()),
                'lost_cds': int(((g['bakta'] != 'no_major_CDS_overlap') & (g['pyrodigal'] == 'no_major_CDS_overlap')).sum()),
            }

# Check per-pair distribution of H1 changes (at 0.5 threshold)
sub05_dark = cc_merged[(cc_merged['overlap_threshold'] == 0.5) & (cc_merged['dark'] == 1)]
pair_gains = sub05_dark.groupby('rep').apply(
    lambda g: pd.Series({
        'n_dark_dna': len(g),
        'n_unshared_dna': int((g['protein_shared'] == 0).sum()),
        'bakta_no_cds': int((g['bakta'] == 'no_major_CDS_overlap').sum()),
        'pyro_no_cds': int((g['pyrodigal'] == 'no_major_CDS_overlap').sum()),
        'gained': int(((g['bakta'] == 'no_major_CDS_overlap') & (g['pyrodigal'] != 'no_major_CDS_overlap')).sum()),
        'lost': int(((g['bakta'] != 'no_major_CDS_overlap') & (g['pyrodigal'] == 'no_major_CDS_overlap')).sum()),
        'net_change': int((g['bakta'] == 'no_major_CDS_overlap').sum() - (g['pyrodigal'] == 'no_major_CDS_overlap').sum())
    }),
    include_groups=False
)
h1_audit['per_pair_dark_05'] = pair_gains.reset_index().to_dict(orient='records')
# Concentration of gains across pairs:
total_gained = pair_gains['gained'].sum()
top1_pair = pair_gains.sort_values('gained', ascending=False).iloc[0]
top3_gained = pair_gains.sort_values('gained', ascending=False).iloc[:3]['gained'].sum()
h1_audit['concentration'] = {
    'total_gained': int(total_gained),
    'pairs_with_any_gain': int((pair_gains['gained'] > 0).sum()),
    'pairs_with_net_improvement': int((pair_gains['net_change'] > 0).sum()),
    'pairs_with_zero_change': int((pair_gains['net_change'] == 0).sum()),
    'pairs_with_net_worsening': int((pair_gains['net_change'] < 0).sum()),
    'top1_rep': str(pair_gains.sort_values('gained', ascending=False).index[0]),
    'top1_gained': int(top1_pair['gained']),
    'top3_gained': int(top3_gained),
    'wilcoxon_pair_rate_p': float(wilcoxon(pair_gains['bakta_no_cds']/pair_gains['n_dark_dna'], pair_gains['pyro_no_cds']/pair_gains['n_dark_dna']).pvalue)
}

# Check length distribution of rescued vs still-missing dark calls in H1
sub05_dark_bakta_no = sub05_dark[sub05_dark['bakta'] == 'no_major_CDS_overlap']
rescued_dark = sub05_dark_bakta_no[sub05_dark_bakta_no['pyrodigal'] != 'no_major_CDS_overlap']
still_no_cds_dark = sub05_dark_bakta_no[sub05_dark_bakta_no['pyrodigal'] == 'no_major_CDS_overlap']
h1_audit['length_of_rescued_vs_unrescued'] = {
    'rescued_n': len(rescued_dark),
    'rescued_median_aa': float(rescued_dark['aa'].median()),
    'rescued_complete_pct': float((rescued_dark['partial'].astype(str) == '00').mean() * 100),
    'still_no_cds_n': len(still_no_cds_dark),
    'still_no_cds_median_aa': float(still_no_cds_dark['aa'].median()),
    'still_no_cds_complete_pct': float((still_no_cds_dark['partial'].astype(str) == '00').mean() * 100),
}
report['h1_audit'] = h1_audit

# 3. Audit H2 (temporal permutation null)
temp_sum = pd.read_csv(ext_dir / 'temporal_summary.tsv', sep='\t')
temp_fam = pd.read_csv(ext_dir / 'temporal_per_family.tsv.gz', sep='\t')
meta = pd.read_csv(gzip.open(ROOT / 'reproducibility/darkmatter/inputs/metadata/sample_metadata.tsv.gz'), sep='\t')
meta['trip_num'] = meta['trip_num'].astype(int)

# Check why within-(site x compartment) permutation increases all-5-trip detection!
# Let's inspect the site x compartment blocks and trip distribution
blocks_df = meta.groupby(['site', 'compartment']).agg(
    n_samples=('sample', 'count'),
    n_trips=('trip_num', 'nunique'),
    has_trip2=('trip_num', lambda x: int((x == 2).any())),
    trip2_samples=('trip_num', lambda x: int((x == 2).sum()))
).reset_index()

# In how many familias is null_all_five_count/499 > observed_all_five?
temp_fam['null_prob'] = temp_fam['null_all_five_count'] / 499.0
temp_fam['diff'] = temp_fam['observed_all_five'] - temp_fam['null_prob']
# Families that are ALWAYS all-5 in null (null_all_five_count == 499) or NEVER all-5 in null (== 0):
h2_audit = {
    'total_families': len(temp_fam),
    'always_all5_under_null': int((temp_fam['null_all_five_count'] == 499).sum()),
    'never_all5_under_null': int((temp_fam['null_all_five_count'] == 0).sum()),
    'variable_under_null': int(((temp_fam['null_all_five_count'] > 0) & (temp_fam['null_all_five_count'] < 499)).sum()),
    'obs_1_but_null_lt_005': int(((temp_fam['observed_all_five'] == 1) & (temp_fam['null_prob'] < 0.05)).sum()),
    'obs_0_but_null_gt_095': int(((temp_fam['observed_all_five'] == 0) & (temp_fam['null_prob'] > 0.95)).sum()),
    'blocks_total': len(blocks_df),
    'blocks_with_trip2': int(blocks_df['has_trip2'].sum()),
    'samples_in_blocks_with_trip2': int(blocks_df.loc[blocks_df['has_trip2'] == 1, 'n_samples'].sum()),
    'trip2_total_samples': int((meta['trip_num'] == 2).sum()),
    'trip2_sites': int(meta.loc[meta['trip_num'] == 2, 'site'].nunique()),
}
report['h2_audit'] = h2_audit

# 4. Audit H3 (paired functional-profile null and threshold sensitivity)
prof_sum = json.loads((ext_dir / 'profile_summary.json').read_text())
prof_thresh_sample = pd.read_csv(ext_dir / 'profile_thresholds_per_sample.tsv', sep='\t')
prof_corr = pd.read_csv(ext_dir / 'profile_pairwise_rank_correlations.tsv.gz', sep='\t', index_col=0)

# Check NaN count in recall_of_mag at each threshold
thresh_nan_audit = {}
for cut, g in prof_thresh_sample.groupby('threshold'):
    thresh_nan_audit[str(cut)] = {
        'n_samples': len(g),
        'mag_kos_zero_samples': int((g['mag_kos'] == 0).sum()),
        'recall_nan_samples': int(g['recall_of_mag'].isna().sum()),
        'recall_non_nan_samples': int(g['recall_of_mag'].notna().sum()),
        'median_mag_kos': float(g['mag_kos'].median()),
        'mean_mag_kos': float(g['mag_kos'].mean()),
        'max_mag_kos': int(g['mag_kos'].max()),
        'median_picrust_kos': float(g['picrust_kos'].median()),
        'median_jaccard': float(g['jaccard'].median()),
        'median_recall_ignoring_nan': float(g['recall_of_mag'].median()),
    }

# Also check in H3: for each sample i, how often is its own paired PICRUSt profile its #1 best match among all 448 samples, or within its block?
S = prof_corr.values
samples_h3 = list(prof_corr.index)
diag = np.diag(S)
# Off-diagonal overall mean:
off_diag_mean = float((S.sum() - diag.sum()) / (len(samples_h3) * (len(samples_h3) - 1)))
# Rank of true partner among all 448 PICRUSt profiles for each MAG sample:
true_rank_mag_to_pic = [int((S[i, :] > S[i, i]).sum()) + 1 for i in range(len(samples_h3))]
# Also check why the site-conditioned null advantage is only 0.00232:
# In trip_compartment_site, what are the 146 exchangeable samples?
# They are biological replicates at the SAME site, SAME compartment, and SAME expedition (trip)!
md_h3 = meta.set_index('sample').loc[samples_h3].reset_index()
tcs_blocks = [v for v in md_h3.groupby(['trip', 'compartment', 'site']).indices.values() if len(v) > 1]
site_comp_blocks = [v for v in md_h3.groupby(['compartment', 'site']).indices.values() if len(v) > 1]
# Wait! What if we permute across trips within (site, compartment), or across sites within (trip, compartment)?
# Notice that (trip, compartment, site) has multiple samples ONLY when co-located biological replicates (r1, r2, r3) from the SAME visit were profiled!
h3_audit = {
    'thresh_nan_audit': thresh_nan_audit,
    'overall_paired_mean': float(diag.mean()),
    'overall_paired_median': float(np.median(diag)),
    'overall_unpaired_mean': off_diag_mean,
    'top1_match_count_of_448': int(sum(r == 1 for r in true_rank_mag_to_pic)),
    'median_rank_of_true_partner_out_of_448': float(np.median(true_rank_mag_to_pic)),
    'tcs_blocks_count': len(tcs_blocks),
    'tcs_samples_count': sum(len(v) for v in tcs_blocks),
    'sc_blocks_count': len(site_comp_blocks),
    'sc_samples_count': sum(len(v) for v in site_comp_blocks),
}
report['h3_audit'] = h3_audit

(OUT / 'audit_part1.json').write_text(json.dumps(report, indent=2) + '\n')
print("Saved audit_part1.json")
