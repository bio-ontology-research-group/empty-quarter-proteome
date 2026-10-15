#!/usr/bin/env python3
"""Part 2b of Round 4 computational audit."""
import os
import gzip
import json
import pandas as pd
import numpy as np
from scipy import stats

ROOT = "/home/leechuck/Documents/papers/empty-quarter-mags"
DM = os.path.join(ROOT, "reproducibility/darkmatter")
INP = os.path.join(DM, "inputs")
SNAP = os.path.join(DM, "source_snapshot")
OUT_DIR = os.path.join(ROOT, "reviews/darkmatter-2026-10-13/round4/reviewer_checks")

out_log = open(os.path.join(OUT_DIR, "r4_audit_part2b.txt"), "w")


def pr(*args, **kwargs):
    print(*args, **kwargs)
    print(*args, **kwargs, file=out_log)


def section(title):
    pr("\n" + "=" * 78)
    pr(title)
    pr("=" * 78)


# Find all files matching key ecology / structure / metabolism names
section("1. FIND ALL ECOLOGY, STRUCTURE & METABOLISM FILES IN REPO")
targets = [
    "partial_summary.tsv",
    "PARTIAL.md",
    "indicators_summary.tsv",
    "INDICATORS.md",
    "permanova_blocks.tsv",
    "varpart.tsv",
    "community_summary.txt",
    "compartment_families_summary.tsv",
    "compartment_families_by_phylum.tsv",
    "set_summary.tsv",
    "k_family_summary.json",
    "duf_share.tsv",
    "top_pfam_K_members.tsv",
    "hit_rate_by_length.tsv",
    "v2_thermo_ecology_summary.tsv",
    "v2_thermo_metabolism_summary.tsv",
    "v2_thermo_metabolism_report.txt",
    "v2_model_report.txt",
    "v2_matched_null_summary.txt",
    "step3_dark_with_cultured_homolog.tsv",
    "step3_cultured_dark_rate_summary.tsv",
]

found = []
for dp, _, fns in os.walk(DM):
    for fn in fns:
        base = fn[:-3] if fn.endswith(".gz") else fn
        if base in targets:
            found.append(os.path.join(dp, fn))

for p in sorted(found):
    pr("FOUND:", os.path.relpath(p, DM))

for p in sorted(found):
    rel = os.path.relpath(p, DM)
    if any(k in rel for k in ["partial_summary", "PARTIAL.md", "indicators_summary", "permanova_blocks", "varpart", "community_summary", "compartment_families_summary", "set_summary", "k_family_summary", "duf_share", "top_pfam_K_members", "hit_rate_by_length", "step3_dark_with_cultured", "step3_cultured_dark_rate"]):
        pr(f"\n--- {rel} ---")
        op = gzip.open if p.endswith(".gz") else open
        with op(p, "rt") as fh:
            pr(fh.read()[:4000])

out_log.close()
