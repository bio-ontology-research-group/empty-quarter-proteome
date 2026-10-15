#!/usr/bin/env python3
"""
scripts/04_exact_length_matching.py (Candidate for Additional Harmonized Draws)

Exact Integer-Length (1 aa) Common-Support Quota Matching:
  q_s(L) = min(N_SMAG_uniq(L), N_GEM_uniq_s(L))
Constructs matched reference database pairs with IDENTICAL unique protein count,
IDENTICAL total amino-acid residues, and IDENTICAL integer-length distributions.

Configurations:
1. Primary Harmonized: SMAG (Pyrodigal v3.7.1) vs GEM (Harmonized Pyrodigal v3.7.1) [seed 20261014]
2. Overlap Excluded: Primary pair with cross-release exact shared sequences removed [seed 20261014_no_overlap]
3. Harmonized Draw Sensitivity 1: SMAG vs GEM (Harmonized Pyrodigal v3.7.1) [seed 20261015]
4. Harmonized Draw Sensitivity 2: SMAG vs GEM (Harmonized Pyrodigal v3.7.1) [seed 20261016]
5. Caller Sensitivity: SMAG (Pyrodigal v3.7.1) vs GEM (Published Prodigal v2.6.3) [seed 20261014]
6. Published Draw Sensitivity 1: SMAG vs GEM (Published Prodigal v2.6.3) [seed 20261015]
7. Published Draw Sensitivity 2: SMAG vs GEM (Published Prodigal v2.6.3) [seed 20261016]

Emits:
- matched FASTA pairs (results/matched_smag_*.faa, results/matched_gem_*.faa)
- results/cross_release_shared_sequences.tsv
- results/matching_closure_audit.tsv
"""

import os
import sys
import random
import hashlib
import numpy as np
from collections import defaultdict

RESULTS_DIR = os.environ.get("RESULTS_DIR", "results")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")

LENGTH_BINS = [
    ("<100", lambda l: l < 100),
    ("100-149", lambda l: 100 <= l <= 149),
    ("150-249", lambda l: 150 <= l <= 249),
    ("250-399", lambda l: 250 <= l <= 399),
    (">=400", lambda l: l >= 400),
]


def bin_index(length_aa):
    if length_aa < 100: return 0
    elif length_aa <= 149: return 1
    elif length_aa <= 249: return 2
    elif length_aa <= 399: return 3
    else: return 4


def load_canonical_unique_fasta(path):
    """Load a deduplicated unique FASTA, canonicalize and return dict of sha256 -> (id, seq)."""
    seen = {}
    audit = {
        "raw_count": 0,
        "raw_residues": 0,
        "trailing_stars": 0,
        "internal_stars": 0,
        "lowercase": 0,
        "invalid_chars": 0,
    }
    cur_id = None
    cur_lines = []

    def flush_prot(sid, lines):
        raw_seq = "".join(lines)
        audit["raw_count"] += 1
        audit["raw_residues"] += len(raw_seq)
        if raw_seq.endswith("*"):
            audit["trailing_stars"] += 1
        stripped = raw_seq.rstrip("*")
        if "*" in stripped:
            audit["internal_stars"] += 1
            return
        audit["lowercase"] += sum(1 for c in stripped if c.islower())
        canon = stripped.upper()
        bad = sum(1 for c in canon if c not in VALID_AA)
        audit["invalid_chars"] += bad
        if not canon or bad > 0:
            return
        h = hashlib.sha256(canon.encode("ascii")).hexdigest()
        if h not in seen:
            seen[h] = (sid, canon)

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    flush_prot(cur_id, cur_lines)
                cur_id = line[1:].strip().split()[0]
                cur_lines = []
            else:
                cur_lines.append(line.strip())
        if cur_id is not None:
            flush_prot(cur_id, cur_lines)

    return seen, audit


def match_pair(smag_dict, gem_dict, seed, exclude_shared, label):
    """Perform exact integer-length matching between SMAG and GEM unique pools."""
    rng = random.Random(seed)

    # Bin by exact integer length
    s_by_len = defaultdict(list)
    g_by_len = defaultdict(list)

    shared_hashes = set(smag_dict.keys()).intersection(set(gem_dict.keys()))
    print(f"  Pre-matching cross-release shared identical sequences: {len(shared_hashes):,}")

    for h, (sid, seq) in smag_dict.items():
        if exclude_shared and h in shared_hashes:
            continue
        s_by_len[len(seq)].append((sid, seq, h))

    for h, (sid, seq) in gem_dict.items():
        if exclude_shared and h in shared_hashes:
            continue
        g_by_len[len(seq)].append((sid, seq, h))

    common_lengths = sorted(set(s_by_len.keys()).intersection(set(g_by_len.keys())))
    print(f"  Exact length bins with common support: {len(common_lengths):,} distinct lengths")

    s_selected = []
    g_selected = []

    for l in common_lengths:
        s_list = s_by_len[l]
        g_list = g_by_len[l]
        quota = min(len(s_list), len(g_list))
        if quota == 0:
            continue

        # Sort for deterministic shuffle
        s_list.sort(key=lambda x: x[0])
        g_list.sort(key=lambda x: x[0])

        # Sample quota without replacement
        if len(s_list) > quota:
            s_chosen = rng.sample(s_list, quota)
        else:
            s_chosen = list(s_list)

        if len(g_list) > quota:
            g_chosen = rng.sample(g_list, quota)
        else:
            g_chosen = list(g_list)

        s_selected.extend(s_chosen)
        g_selected.extend(g_chosen)

    # Sort deterministically
    s_selected.sort(key=lambda x: x[0])
    g_selected.sort(key=lambda x: x[0])

    # Write matched FASTAs
    s_faa = os.path.join(RESULTS_DIR, f"matched_smag_{label}.faa")
    g_faa = os.path.join(RESULTS_DIR, f"matched_gem_{label}.faa")

    tmp_s_faa = s_faa + f".tmp_{os.getpid()}"
    tmp_g_faa = g_faa + f".tmp_{os.getpid()}"

    with open(tmp_s_faa, "w", encoding="utf-8") as out_s:
        for sid, seq, h in s_selected:
            out_s.write(f">{sid}\n{seq}\n")
    os.replace(tmp_s_faa, s_faa)

    with open(tmp_g_faa, "w", encoding="utf-8") as out_g:
        for sid, seq, h in g_selected:
            out_g.write(f">{sid}\n{seq}\n")
    os.replace(tmp_g_faa, g_faa)

    # Programmatic verification
    s_lens = [len(x[1]) for x in s_selected]
    g_lens = [len(x[1]) for x in g_selected]
    s_uniq_h = {x[2] for x in s_selected}
    g_uniq_h = {x[2] for x in g_selected}
    post_shared = len(s_uniq_h.intersection(g_uniq_h))

    assert len(s_selected) == len(g_selected), f"Count mismatch: {len(s_selected)} vs {len(g_selected)}"
    assert len(s_uniq_h) == len(s_selected), f"SMAG duplicate sequences remain: {len(s_selected) - len(s_uniq_h)}"
    assert len(g_uniq_h) == len(g_selected), f"GEM duplicate sequences remain: {len(g_selected) - len(g_uniq_h)}"
    assert sum(s_lens) == sum(g_lens), f"Total residues mismatch: {sum(s_lens)} vs {sum(g_lens)}"
    assert sorted(s_lens) == sorted(g_lens), "Integer length distributions do not match!"
    if exclude_shared:
        assert post_shared == 0, f"Expected 0 shared sequences, got {post_shared}"

    arr = np.array(s_lens, dtype=np.float64)
    bin_counts = [0, 0, 0, 0, 0]
    for l in s_lens:
        bin_counts[bin_index(l)] += 1

    return {
        "label": label,
        "seed": seed,
        "exclude_shared": exclude_shared,
        "matched_proteins": len(s_selected),
        "matched_residues": int(sum(s_lens)),
        "pre_matching_shared_seqs": len(shared_hashes),
        "post_matching_shared_seqs": post_shared,
        "mean_length_aa": round(float(np.mean(arr)), 2),
        "sd_length_aa": round(float(np.std(arr, ddof=1)), 2),
        "median_length_aa": round(float(np.median(arr)), 2),
        "bin_lt100": bin_counts[0],
        "bin_100_149": bin_counts[1],
        "bin_150_249": bin_counts[2],
        "bin_250_399": bin_counts[3],
        "bin_gte400": bin_counts[4],
        "s_faa_path": s_faa,
        "g_faa_path": g_faa,
    }


def main():
    print("=== H3 Step 4: Exact Integer-Length (1 aa) Common-Support Quota Matching ===")
    smag_path = os.path.join(RESULTS_DIR, "smag_pool_unique.faa")
    print(f"Loading unique SMAG pool from {smag_path}...")
    if not os.path.exists(smag_path):
        raise FileNotFoundError(f"SMAG unique pool missing: {smag_path}")
    s_dict, s_aud = load_canonical_unique_fasta(smag_path)
    print(f"  SMAG unique proteins: {len(s_dict):,}, total residues: {sum(len(v[1]) for v in s_dict.values()):,}")

    configurations = [
        # (label, gem_unique_path, seed, exclude_shared)
        ("harmonized_seed20261014", os.path.join(RESULTS_DIR, "gem_harmonized_pool_unique_seed20261014.faa"), 20261014, False),
        ("harmonized_seed20261014_no_overlap", os.path.join(RESULTS_DIR, "gem_harmonized_pool_unique_seed20261014.faa"), 20261014, True),
        ("harmonized_seed20261015", os.path.join(RESULTS_DIR, "gem_harmonized_pool_unique_seed20261015.faa"), 20261015, False),
        ("harmonized_seed20261016", os.path.join(RESULTS_DIR, "gem_harmonized_pool_unique_seed20261016.faa"), 20261016, False),
        ("published_seed20261014", os.path.join(RESULTS_DIR, "gem_published_pool_unique_seed20261014.faa"), 20261014, False),
        ("published_seed20261015", os.path.join(RESULTS_DIR, "gem_published_pool_unique_seed20261015.faa"), 20261015, False),
        ("published_seed20261016", os.path.join(RESULTS_DIR, "gem_published_pool_unique_seed20261016.faa"), 20261016, False),
    ]

    closure_rows = []
    shared_tsv_path = os.path.join(RESULTS_DIR, "cross_release_shared_sequences.tsv")
    with open(shared_tsv_path, "w", encoding="utf-8") as out_sh:
        out_sh.write("comparison_label\tseed\tshared_seq_sha256\tlength_aa\tsmag_seq_id\tgem_seq_id\n")

    for label, gem_path, seed, exclude_shared in configurations:
        print(f"\n--- Matching Pair: {label} ---")
        if not os.path.exists(gem_path):
            raise FileNotFoundError(f"Required comparison {label} input missing: {gem_path}")
        g_dict, g_aud = load_canonical_unique_fasta(gem_path)
        print(f"  GEM unique proteins: {len(g_dict):,}, total residues: {sum(len(v[1]) for v in g_dict.values()):,}")

        # Record shared sequences
        shared_hashes = set(s_dict.keys()).intersection(set(g_dict.keys()))
        with open(shared_tsv_path, "a", encoding="utf-8") as out_sh:
            for sh in sorted(shared_hashes):
                sid_s, seq_s = s_dict[sh]
                sid_g, seq_g = g_dict[sh]
                out_sh.write(f"{label}\t{seed}\t{sh}\t{len(seq_s)}\t{sid_s}\t{sid_g}\n")

        info = match_pair(s_dict, g_dict, seed, exclude_shared, label)
        closure_rows.append(info)
        print(f"  SUCCESS: Matched {info['matched_proteins']:,} proteins ({info['matched_residues']:,} aa)")
        print(f"  Mean length: {info['mean_length_aa']} aa, 5-bin counts: {info['bin_lt100']}, {info['bin_100_149']}, {info['bin_150_249']}, {info['bin_250_399']}, {info['bin_gte400']}")

    audit_path = os.path.join(RESULTS_DIR, "matching_closure_audit.tsv")
    tmp_audit_path = audit_path + f".tmp_{os.getpid()}"
    with open(tmp_audit_path, "w", encoding="utf-8") as out:
        out.write("label\tseed\texclude_shared\tmatched_proteins\tmatched_residues\tpre_matching_shared_seqs\t"
                  "post_matching_shared_seqs\tmean_length_aa\tsd_length_aa\tmedian_length_aa\t"
                  "bin_lt100\tbin_100_149\tbin_150_249\tbin_250_399\tbin_gte400\ts_faa_path\tg_faa_path\n")
        for r in closure_rows:
            out.write("\t".join(str(r[k]) for k in [
                "label", "seed", "exclude_shared", "matched_proteins", "matched_residues",
                "pre_matching_shared_seqs", "post_matching_shared_seqs", "mean_length_aa",
                "sd_length_aa", "median_length_aa", "bin_lt100", "bin_100_149", "bin_150_249",
                "bin_250_399", "bin_gte400", "s_faa_path", "g_faa_path"
            ]) + "\n")
    os.replace(tmp_audit_path, audit_path)

    print(f"\nWrote matching closure audit to {audit_path}")


if __name__ == "__main__":
    main()
