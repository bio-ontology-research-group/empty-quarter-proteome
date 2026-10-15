#!/usr/bin/env python3
"""
scripts/05_diamond_search_full_census.py
DIAMOND Searches & Statistical Census Across ALL 212,084 Queries:
- 106,042 dark families (FAMILIES_FINAL.txt)
- 106,042 annotated controls (CONTROL_FAMILIES.txt)
- Query FASTA: query_all.faa

Executes paired searches:
1. Primary Harmonized: SMAG (Pyrodigal v3.7.1) vs GEM (Harmonized Pyrodigal v3.7.1) [seed 20261014]
2. Overlap Excluded: Primary pair with cross-release shared sequences removed [seed 20261014_no_overlap]
3. Identity Sensitivity: Primary pair searched with --id 30 [seed 20261014_id30]
4. Caller Sensitivity: SMAG (Pyrodigal v3.7.1) vs GEM (Published Prodigal v2.6.3) [seed 20261014]
5. Draw Sensitivity 1: SMAG vs GEM (Published) [seed 20261015]
6. Draw Sensitivity 2: SMAG vs GEM (Published) [seed 20261016]

Computes exact 2x2 contingency tables (Both, SMAG-only, GEM-only, Neither),
paired differences (b-c)/N, exact binomial McNemar p-values, discordant ratios b/c,
and per-length-stratum breakdowns. Separate dark and control arms (never pooled).

Usage:
    python3 scripts/05_diamond_search_full_census.py [threads]
"""

import os
import sys
import subprocess
import time
from collections import defaultdict
from scipy import stats

RESULTS_DIR = "results"
QUERY_FAA = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/results/query_all.faa"
DARK_IDS_PATH = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results/FAMILIES_FINAL.txt"
CTRL_IDS_PATH = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results/CONTROL_FAMILIES.txt"
DIAMOND_BIN = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"

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


def load_query_lengths(path):
    lens = {}
    cur_id = None
    cur_len = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    lens[cur_id] = cur_len
                cur_id = line[1:].strip().split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip().rstrip("*"))
        if cur_id is not None:
            lens[cur_id] = cur_len
    return lens


def run_diamond_search(faa_path, db_path, m8_path, threads, id_cutoff=0):
    """Build diamond db if needed and run search."""
    if not os.path.exists(db_path + ".dmnd"):
        print(f"  Building DIAMOND DB: {db_path}...")
        t0 = time.time()
        cmd_mk = [DIAMOND_BIN, "makedb", "--in", faa_path, "-d", db_path, "-p", str(threads)]
        subprocess.run(cmd_mk, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        print(f"  DB built in {time.time() - t0:.1f}s.")

    if not os.path.exists(m8_path) or os.path.getsize(m8_path) == 0:
        print(f"  Running DIAMOND search: {m8_path} (id_cutoff={id_cutoff})...")
        t0 = time.time()
        cmd_srch = [
            DIAMOND_BIN, "blastp",
            "-q", QUERY_FAA,
            "-d", db_path,
            "--sensitive",
            "--evalue", "1e-3",
            "-k", "1",
            "-p", str(threads),
            "-f", "6", "qseqid", "sseqid", "pident", "evalue",
            "-o", m8_path
        ]
        if id_cutoff > 0:
            cmd_srch.extend(["--id", str(id_cutoff)])
        subprocess.run(cmd_srch, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        print(f"  Search completed in {time.time() - t0:.1f}s ({os.path.getsize(m8_path) / 1e6:.1f} MB).")


def load_hit_ids(m8_path):
    hits = set()
    with open(m8_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                hits.add(parts[0])
    return hits


def score_cohort(q_ids, s_hits, g_hits, q_lens):
    total = len(q_ids)
    a = 0  # both
    b = 0  # smag only
    c = 0  # gem only
    d = 0  # neither

    bin_counts = [defaultdict(int) for _ in range(5)]

    for qid in q_ids:
        in_s = qid in s_hits
        in_g = qid in g_hits
        l_idx = bin_index(q_lens[qid])

        if in_s and in_g:
            a += 1
            bin_counts[l_idx]["a"] += 1
        elif in_s and not in_g:
            b += 1
            bin_counts[l_idx]["b"] += 1
        elif not in_s and in_g:
            c += 1
            bin_counts[l_idx]["c"] += 1
        else:
            d += 1
            bin_counts[l_idx]["d"] += 1

    p_s = 100.0 * (a + b) / total
    p_g = 100.0 * (a + c) / total
    delta = 100.0 * (b - c) / total
    disc_total = b + c
    pval = stats.binomtest(b, disc_total, 0.5).pvalue if disc_total > 0 else 1.0
    disc_ratio = round(b / c, 2) if c > 0 else "undefined"

    stratum_results = []
    for idx, (label, _) in enumerate(LENGTH_BINS):
        sb = bin_counts[idx]
        tot_h = sb["a"] + sb["b"] + sb["c"] + sb["d"]
        sh_s = 100.0 * (sb["a"] + sb["b"]) / tot_h if tot_h > 0 else 0.0
        sh_g = 100.0 * (sb["a"] + sb["c"]) / tot_h if tot_h > 0 else 0.0
        d_h = 100.0 * (sb["b"] - sb["c"]) / tot_h if tot_h > 0 else 0.0
        disc_h = sb["b"] + sb["c"]
        p_h = stats.binomtest(sb["b"], disc_h, 0.5).pvalue if disc_h > 0 else 1.0
        ratio_h = round(sb["b"] / sb["c"], 2) if sb["c"] > 0 else "undefined"
        stratum_results.append({
            "stratum": label,
            "total": tot_h,
            "a": sb["a"],
            "b": sb["b"],
            "c": sb["c"],
            "d": sb["d"],
            "smag_hit_pct": round(sh_s, 2),
            "gem_hit_pct": round(sh_g, 2),
            "delta_pct": round(d_h, 2),
            "mcnemar_pval": p_h,
            "disc_ratio": ratio_h,
        })

    return {
        "total": total,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "smag_hit_pct": round(p_s, 2),
        "gem_hit_pct": round(p_g, 2),
        "delta_pct": round(delta, 2),
        "mcnemar_pval": pval,
        "disc_ratio": disc_ratio,
        "strata": stratum_results,
    }


def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    print(f"=== H3 Step 5: Paired DIAMOND Searches & Statistical Census Across ALL 212,084 Queries ({threads} threads) ===")

    print("Loading query identifiers and lengths...")
    with open(DARK_IDS_PATH, "r") as f:
        dark_ids = sorted(line.strip() for line in f if line.strip())
    with open(CTRL_IDS_PATH, "r") as f:
        ctrl_ids = sorted(line.strip() for line in f if line.strip())
    print(f"  Dark queries: {len(dark_ids):,}")
    print(f"  Control queries: {len(ctrl_ids):,}")

    q_lens = load_query_lengths(QUERY_FAA)
    print(f"  Loaded lengths for {len(q_lens):,} queries.")

    searches = [
        # (run_label, s_faa, g_faa, id_cutoff, description)
        (
            "harmonized_seed20261014",
            os.path.join(RESULTS_DIR, "matched_smag_harmonized_seed20261014.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_harmonized_seed20261014.faa"),
            0,
            "Primary Harmonized (Pyrodigal v3.7.1 both arms, exact length quota matched)"
        ),
        (
            "harmonized_seed20261014_no_overlap",
            os.path.join(RESULTS_DIR, "matched_smag_harmonized_seed20261014_no_overlap.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_harmonized_seed20261014_no_overlap.faa"),
            0,
            "Primary Harmonized Overlap-Excluded Sensitivity (0 shared sequences)"
        ),
        (
            "harmonized_seed20261014_id30",
            os.path.join(RESULTS_DIR, "matched_smag_harmonized_seed20261014.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_harmonized_seed20261014.faa"),
            30,
            "Primary Harmonized Identity-Threshold Sensitivity (--id 30 during search)"
        ),
        (
            "published_seed20261014",
            os.path.join(RESULTS_DIR, "matched_smag_published_seed20261014.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_published_seed20261014.faa"),
            0,
            "Caller Sensitivity: SMAG (Pyrodigal v3.7.1) vs GEM (Published Prodigal v2.6.3)"
        ),
        (
            "published_seed20261015",
            os.path.join(RESULTS_DIR, "matched_smag_published_seed20261015.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_published_seed20261015.faa"),
            0,
            "Reference Draw Sensitivity 1: Seed 20261015"
        ),
        (
            "published_seed20261016",
            os.path.join(RESULTS_DIR, "matched_smag_published_seed20261016.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_published_seed20261016.faa"),
            0,
            "Reference Draw Sensitivity 2: Seed 20261016"
        ),
    ]

    summary_rows = []
    stratum_rows = []

    for run_lbl, s_faa, g_faa, id_cut, desc in searches:
        print(f"\n==================================================================")
        print(f"Executing Search Pair: {run_lbl} ({desc})")
        print(f"==================================================================")
        if not os.path.exists(s_faa) or not os.path.exists(g_faa):
            print(f"  Missing input FASTA: {s_faa} or {g_faa} — skipping.")
            continue

        s_db = s_faa.replace(".faa", "_db")
        g_db = g_faa.replace(".faa", "_db")
        suffix = f"_id{id_cut}" if id_cut > 0 else ""
        s_m8 = os.path.join(RESULTS_DIR, f"hits_smag_{run_lbl}{suffix}.m8")
        g_m8 = os.path.join(RESULTS_DIR, f"hits_gem_{run_lbl}{suffix}.m8")

        run_diamond_search(s_faa, s_db, s_m8, threads, id_cutoff=id_cut)
        run_diamond_search(g_faa, g_db, g_m8, threads, id_cutoff=id_cut)

        print("  Scoring census hits against 106,042 dark and 106,042 control queries...")
        s_hits = load_hit_ids(s_m8)
        g_hits = load_hit_ids(g_m8)
        print(f"    Raw hit queries: SMAG={len(s_hits):,}, GEM={len(g_hits):,}")

        # Dark cohort scoring
        dark_res = score_cohort(dark_ids, s_hits, g_hits, q_lens)
        ctrl_res = score_cohort(ctrl_ids, s_hits, g_hits, q_lens)

        for cohort_name, res in [("Dark", dark_res), ("Control", ctrl_res)]:
            summary_rows.append((
                run_lbl, desc, cohort_name, id_cut, res["total"],
                res["a"], res["b"], res["c"], res["d"],
                res["smag_hit_pct"], res["gem_hit_pct"], res["delta_pct"],
                res["disc_ratio"], f"{res['mcnemar_pval']:.3e}"
            ))
            print(f"    [{cohort_name}] N={res['total']:,} | Both={res['a']:,} | SMAG-only(b)={res['b']:,} | GEM-only(c)={res['c']:,} | Neither={res['d']:,}")
            print(f"      SMAG hit rate: {res['smag_hit_pct']}%, GEM hit rate: {res['gem_hit_pct']}%, Delta: {res['delta_pct']:+0.2f}%")
            print(f"      Discordant ratio: {res['disc_ratio']}, McNemar p-val: {res['mcnemar_pval']:.3e}")

            for st in res["strata"]:
                stratum_rows.append((
                    run_lbl, cohort_name, st["stratum"], id_cut, st["total"],
                    st["a"], st["b"], st["c"], st["d"],
                    st["smag_hit_pct"], st["gem_hit_pct"], st["delta_pct"],
                    st["disc_ratio"], f"{st['mcnemar_pval']:.3e}"
                ))

    # Write summary table
    summary_path = os.path.join(RESULTS_DIR, "matched_full_census_summary.tsv")
    with open(summary_path, "w") as out:
        out.write("run_label\tdescription\tcohort\tid_cutoff\ttotal_queries\tboth_hits_a\tsmag_only_b\t"
                  "gem_only_c\tneither_d\tsmag_hit_pct\tgem_hit_pct\tpaired_delta_pct\tdiscordant_ratio\tmcnemar_pvalue\n")
        for r in summary_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    # Write stratum table
    stratum_path = os.path.join(RESULTS_DIR, "matched_length_stratified_results.tsv")
    with open(stratum_path, "w") as out:
        out.write("run_label\tcohort\tlength_stratum\tid_cutoff\tstratum_queries\tboth_hits_a\tsmag_only_b\t"
                  "gem_only_c\tneither_d\tsmag_hit_pct\tgem_hit_pct\tpaired_delta_pct\tdiscordant_ratio\tmcnemar_pvalue\n")
        for r in stratum_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    print(f"\nWrote full census summary to {summary_path}")
    print(f"Wrote length-stratified results to {stratum_path}")


if __name__ == "__main__":
    main()
