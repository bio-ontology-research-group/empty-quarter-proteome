#!/usr/bin/env python3
"""Part 2 read-only audit: manuscript cross-checks, H1/H2/H3 deep-dive, and candidate extension data verification."""
from pathlib import Path
import re
import json
import gzip
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, spearmanr, mannwhitneyu

ROOT = Path('/home/leechuck/Documents/papers/empty-quarter-mags')
ROUND5 = ROOT / 'reviews/darkmatter-2026-10-13/round5'
OUT = ROUND5 / 'argon_checks'

tex = (ROUND5 / 'baseline_main.tex').read_text()
lines = tex.splitlines()

out = {}

# 1. Check all occurrences of Table / Supplementary / Extended Data / Figure in baseline_main.tex
refs_in_tex = []
for i, l in enumerate(lines, 1):
    if any(k in l for k in ['Table', 'Supplementary', 'Extended Data', 'Fig.', 'Figure', '\\ref{', '\\label{', 'todo', 'unreg']):
        refs_in_tex.append({'line': i, 'text': l.strip()})
out['refs_in_tex'] = refs_in_tex

# 2. H1 deep-dive with dtype={'partial': str}
ext_dir = ROOT / 'reproducibility/darkmatter/round5_extensions/results'
cc_dir = ext_dir / 'culture_caller'
cc_df = pd.read_csv(cc_dir / 'per_gene.tsv.gz', sep='\t', dtype={'partial': str})
cigar_df = pd.read_csv(ROOT / 'reproducibility/darkmatter/culture_nucleotide_cigar/per_gene.tsv', sep='\t', dtype={'partial': str})
cc_merged = cc_df.merge(
    cigar_df[['rep', 'protein', 'protein_shared']],
    on=['rep', 'protein'],
    how='left',
    validate='many_to_one'
)

sub05 = cc_merged[cc_merged['overlap_threshold'] == 0.5]
sub05_dark = sub05[sub05['dark'] == 1]
sub05_dark_bakta_no = sub05_dark[sub05_dark['bakta'] == 'no_major_CDS_overlap']
rescued_dark = sub05_dark_bakta_no[sub05_dark_bakta_no['pyrodigal'] != 'no_major_CDS_overlap']
still_no_cds_dark = sub05_dark_bakta_no[sub05_dark_bakta_no['pyrodigal'] == 'no_major_CDS_overlap']

# Breakdown by partiality ('00' = complete, '01'/'10'/'11' = partial/contig-edge)
partial_breakdown_dark = {}
for p_flag, g in sub05_dark.groupby(sub05_dark['partial'] == '00'):
    label = 'complete_00' if p_flag else 'partial_non00'
    partial_breakdown_dark[label] = {
        'n': len(g),
        'median_aa': float(g['aa'].median()),
        'n_ge100aa': int((g['aa'] >= 100).sum()),
        'bakta_no_cds': int((g['bakta'] == 'no_major_CDS_overlap').sum()),
        'bakta_no_cds_pct': float((g['bakta'] == 'no_major_CDS_overlap').mean() * 100),
        'pyro_no_cds': int((g['pyrodigal'] == 'no_major_CDS_overlap').sum()),
        'pyro_no_cds_pct': float((g['pyrodigal'] == 'no_major_CDS_overlap').mean() * 100),
        'gained': int(((g['bakta'] == 'no_major_CDS_overlap') & (g['pyrodigal'] != 'no_major_CDS_overlap')).sum()),
        'lost': int(((g['bakta'] != 'no_major_CDS_overlap') & (g['pyrodigal'] == 'no_major_CDS_overlap')).sum()),
    }

# Per-pair summary statistics for H1
cc_pair = pd.read_csv(cc_dir / 'per_pair.tsv', sep='\t')
dark_pair = cc_pair[cc_pair['dark'] == 1]
annot_pair = cc_pair[cc_pair['dark'] == 0]
out['h1_deep'] = {
    'partial_breakdown_dark_05': partial_breakdown_dark,
    'rescued_n': len(rescued_dark),
    'rescued_median_aa': float(rescued_dark['aa'].median()),
    'rescued_complete_00_n': int((rescued_dark['partial'] == '00').sum()),
    'rescued_complete_00_pct': float((rescued_dark['partial'] == '00').mean() * 100),
    'still_no_cds_n': len(still_no_cds_dark),
    'still_no_cds_median_aa': float(still_no_cds_dark['aa'].median()),
    'still_no_cds_complete_00_n': int((still_no_cds_dark['partial'] == '00').sum()),
    'still_no_cds_complete_00_pct': float((still_no_cds_dark['partial'] == '00').mean() * 100),
    'per_pair_dark_bakta_median': float(dark_pair['bakta_no_cds'].median()),
    'per_pair_dark_bakta_iqr': [float(dark_pair['bakta_no_cds'].quantile(0.25)), float(dark_pair['bakta_no_cds'].quantile(0.75))],
    'per_pair_dark_pyro_median': float(dark_pair['pyrodigal_no_cds'].median()),
    'per_pair_dark_pyro_iqr': [float(dark_pair['pyrodigal_no_cds'].quantile(0.25)), float(dark_pair['pyrodigal_no_cds'].quantile(0.75))],
    'per_pair_dark_delta_median': float((dark_pair['pyrodigal_no_cds'] - dark_pair['bakta_no_cds']).median()),
}

# 3. H3 additional block structures from profile_pairwise_rank_correlations.tsv.gz
prof_corr = pd.read_csv(ext_dir / 'profile_pairwise_rank_correlations.tsv.gz', sep='\t', index_col=0)
meta = pd.read_csv(gzip.open(ROOT / 'reproducibility/darkmatter/inputs/metadata/sample_metadata.tsv.gz'), sep='\t')
samples_h3 = list(prof_corr.index)
S = prof_corr.values
md_h3 = meta.set_index('sample').loc[samples_h3].reset_index()

rng = np.random.default_rng(20261016)
h3_extra_nulls = {}
for label, keys in [
    ('site_compartment_across_trips', ['site', 'compartment']),
    ('trip_site_across_compartments', ['trip', 'site']),
    ('site_only', ['site']),
]:
    blocks = [np.asarray(v) for v in md_h3.groupby(keys).indices.values() if len(v) > 1]
    ix = np.concatenate(blocks)
    obs_mean = float(np.diag(S)[ix].mean())
    null_draws = []
    for _ in range(999):
        perm = np.arange(len(samples_h3))
        for v in blocks:
            perm[v] = rng.permutation(v)
        null_draws.append(float(S[ix, perm[ix]].mean()))
    null_arr = np.array(null_draws)
    h3_extra_nulls[label] = {
        'exchangeable_samples': int(len(ix)),
        'blocks': int(len(blocks)),
        'observed_mean': obs_mean,
        'null_mean': float(null_arr.mean()),
        'advantage': float(obs_mean - null_arr.mean()),
        'p_upper': float((1 + int((null_arr >= obs_mean).sum())) / 1000)
    }

# Also check within `trip_compartment_site`: what are the off-diagonal correlations inside those 66 replicate blocks?
# Notice that inside a `trip_compartment_site` block, all samples are biological replicates (r1, r2, r3) from the same site, compartment, and trip!
out['h3_extra_nulls'] = h3_extra_nulls

# 4. Check existing data in source_snapshot for candidate falsifiable extensions!
snap = ROOT / 'reproducibility/darkmatter/source_snapshot'

# Extension Candidate A: Thermoplasmatota known relatives AAI vs dark_fraction and class SW-10-69-26 vs E2 reductive evolution
thermo_rel_path = snap / '06_darkmatter/lineage_darkness_clean/results/v2_thermo_relatives_pairs.tsv'
reps_path = ROOT / 'reproducibility/darkmatter/inputs/06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv.gz'
if thermo_rel_path.exists():
    tr = pd.read_csv(thermo_rel_path, sep='\t')
    reps = pd.read_csv(gzip.open(reps_path), sep='\t')
    tr_merged = tr.merge(reps, on='genome', how='inner')
    rho_ext, p_ext = spearmanr(tr_merged['aai'], tr_merged['dark_fraction'])
    rho_desert, p_desert = spearmanr(tr_merged['aai_desert'], tr_merged['dark_fraction'])
    rho_rbh, p_rbh = spearmanr(tr_merged['aai_frac_rbh'], tr_merged['dark_fraction'])
    sw = tr_merged[tr_merged['cls'] == 'SW-10-69-26']
    e2 = tr_merged[tr_merged['cls'] == 'E2']
    out['candidate_thermo_aai'] = {
        'n_merged': len(tr_merged),
        'columns': list(tr.columns),
        'spearman_best_external_aai_vs_dark': {'rho': float(rho_ext), 'p': float(p_ext)},
        'spearman_best_desert_aai_vs_dark': {'rho': float(rho_desert), 'p': float(p_desert)},
        'spearman_rbh_frac_vs_dark': {'rho': float(rho_rbh), 'p': float(p_rbh)},
        'aai_ge65_n': int((tr_merged['aai'] >= 65).sum()),
        'aai_lt65_n': int((tr_merged['aai'] < 65).sum()),
        'aai_ge65_median_dark': float(tr_merged.loc[tr_merged['aai'] >= 65, 'dark_fraction'].median()),
        'aai_lt65_median_dark': float(tr_merged.loc[tr_merged['aai'] < 65, 'dark_fraction'].median()),
        'sw_n': len(sw),
        'sw_median_dark': float(sw['dark_fraction'].median()),
        'sw_median_genome_mb': float(sw['genome_size_bp'].median() / 1e6) if 'genome_size_bp' in sw.columns else None,
        'sw_median_cds': float(sw['total_genes'].median()) if 'total_genes' in sw.columns else None,
        'sw_spearman_aai_vs_dark': {'rho': float(spearmanr(sw['aai'], sw['dark_fraction'])[0]), 'p': float(spearmanr(sw['aai'], sw['dark_fraction'])[1])},
        'e2_n': len(e2),
        'e2_median_dark': float(e2['dark_fraction'].median()),
        'e2_median_genome_mb': float(e2['genome_size_bp'].median() / 1e6) if 'genome_size_bp' in e2.columns else None,
        'e2_median_cds': float(e2['total_genes'].median()) if 'total_genes' in e2.columns else None,
        'e2_spearman_aai_vs_dark': {'rho': float(spearmanr(e2['aai'], e2['dark_fraction'])[0]), 'p': float(spearmanr(e2['aai'], e2['dark_fraction'])[1])},
    }

# Extension Candidate B: SMAG (40,039 soil MAGs, 129.55M proteins) vs GEM (52,515 MAGs, 122.98M proteins) reference-biome test + B6 tables
smag_cvd_path = snap / '06_darkmatter/specificity_clean/results/smag_content_vs_depth.tsv'
smag_union_path = snap / '06_darkmatter/specificity_clean/results/union_specificity_with_smag_summary.tsv'
b6_len_path = snap / '06_darkmatter/specificity_clean/results/b6_length_confound_with_smag.tsv'
b6_occ_path = snap / '06_darkmatter/specificity_clean/results/b6_absence_occupancy_by_length_with_smag.tsv'
if smag_cvd_path.exists():
    out['candidate_smag'] = {
        'smag_content_vs_depth': pd.read_csv(smag_cvd_path, sep='\t').to_dict(orient='records'),
        'union_with_smag': pd.read_csv(smag_union_path, sep='\t').to_dict(orient='records') if smag_union_path.exists() else None,
        'b6_length_confound': pd.read_csv(b6_len_path, sep='\t').to_dict(orient='records') if b6_len_path.exists() else None,
        'b6_absence_occupancy_by_length': pd.read_csv(b6_occ_path, sep='\t').to_dict(orient='records') if b6_occ_path.exists() else None,
    }

# Extension Candidate C: Pfam domain architecture on K_members (4,880) vs QnotK_reps (1,807) vs controls (3,013) + stratified pLDDT floor
pfam_sum_path = snap / '04_structure/pfam_candidates_clean/results/set_summary.tsv'
pfam_duf_path = snap / '04_structure/pfam_candidates_clean/results/duf_share.tsv'
strat_floor_path = snap / '04_structure/00_structure_clean/05_validate/stratified_floor.tsv'
if pfam_sum_path.exists():
    out['candidate_pfam_and_structure'] = {
        'pfam_set_summary': pd.read_csv(pfam_sum_path, sep='\t').to_dict(orient='records'),
        'pfam_duf_share': pd.read_csv(pfam_duf_path, sep='\t').to_dict(orient='records') if pfam_duf_path.exists() else None,
        'stratified_floor': pd.read_csv(strat_floor_path, sep='\t').to_dict(orient='records') if strat_floor_path.exists() else None
    }

# Extension Candidate D: Contamination check (Cutibacterium acnes 51Dr3__single_SemiBin_390 + 10 plasmid-matching reps in contamination_clean)
contam_prof_path = snap / '06_darkmatter/contamination_clean/results/control_profiles.tsv'
contam_hits_path = snap / '06_darkmatter/contamination_clean/results/control_contig_hits.tsv'
if contam_prof_path.exists():
    out['candidate_contamination'] = {
        'control_profiles': pd.read_csv(contam_prof_path, sep='\t').to_dict(orient='records'),
        'control_contig_hits_n': len(pd.read_csv(contam_hits_path, sep='\t')) if contam_hits_path.exists() else None,
        'control_contig_hits_reps': list(pd.read_csv(contam_hits_path, sep='\t')['rep'].unique()) if contam_hits_path.exists() and 'rep' in pd.read_csv(contam_hits_path, sep='\t').columns else None
    }

(OUT / 'audit_part2.json').write_text(json.dumps(out, indent=2) + '\n')
print("Saved audit_part2.json")
