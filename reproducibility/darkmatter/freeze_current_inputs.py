#!/usr/bin/env python3
"""Freeze compact current claim tables; no upstream recomputation is implied."""
import csv
import gzip
import hashlib
from pathlib import Path

here = Path(__file__).resolve().parent
paths = [
    "06_darkmatter/lineage_darkness_clean/results/v2_dark_composition_by_phylum.tsv",
    "06_darkmatter/lineage_darkness_clean/results/v2_thermo_ecology_summary.tsv",
    "mimag_rrna_trna_clean/results/mimag_tiers.tsv",
    "06_reference_bias/refbias_clean/results/per_sample_comparison.tsv",
    "06_darkmatter/specificity_clean/results/content_vs_depth.tsv",
    "06_darkmatter/specificity_clean/results/b6_absence_occupancy_by_length.tsv",
    "metadata/sample_metadata.tsv",
    "06_darkmatter/lineage_darkness_clean/results/v2_matched_null.tsv",
    "06_darkmatter/lineage_darkness_clean/results/v2_phylum_dark.tsv",
    '06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv',
    '06_darkmatter/darkfam_basis_clean/results/selection/dnds_null_frame1.tsv',
    '06_darkmatter/darkfam_basis_clean/results/selection/dnds_dark_subset.tsv',
    '06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv',
    '04_structure/00_structure_clean/05_validate/per_rep_structures.tsv',
    '04_structure/00_structure_clean/05_validate/per_control_structures.tsv',
    '06_darkmatter/cultured_comparison_clean/results/step3_length_summary.tsv',
    '06_darkmatter/cultured_comparison_clean/results/step3_paired_by_species.tsv',
    '06_darkmatter/cultured_comparison_clean/results/pair_anatomy/pair_anatomy_per_pair.tsv',
]
with (here / 'current_checks_inputs.tsv').open('w') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['relative_source', 'ibex_path', 'local_path', 'source_sha256', 'compressed_sha256', 'source_bytes'])
    for path in paths:
        source = here / 'source_snapshot' / path
        blob = source.read_bytes()
        compressed = gzip.compress(blob, mtime=0)
        output = here / 'inputs' / (path + '.gz')
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(compressed)
        w.writerow([path, '/ibex/scratch/projects/c2014/rob/metagenomics/' + path, str(output.relative_to(here)), hashlib.sha256(blob).hexdigest(), hashlib.sha256(compressed).hexdigest(), len(blob)])
