#!/usr/bin/env python3
"""
scripts/09_v3_match_gem_and_search.py
v3 Repair for Hypothesis 3 (Critical Findings 1, 2, & 4):
1. Canonicalize all sequences (strip trailing '*', uppercase, audit internal '*' and invalid chars).
2. Sample GEM genomes across 3 deterministic seeds (20261013, 20261014, 20261015) until each seed's
   deduplicated canonical GEM pool reaches >= the deduplicated SMAG pool size (692,194 unique proteins).
   Write full file-level GEM genome manifests (with .faa.gz byte size and SHA-256) for every seed.
3. Audit within-pool redundancy and exact cross-release sequence overlap (saving shared sequence IDs).
4. Enforce REAL exact integer-length (1 aa) common-support quota matching:
     q_s(L) = min(N_SMAG_uniq(L), N_GEM_uniq_s(L))
   so that matched SMAG and GEM databases have 0% internal duplicates, IDENTICAL protein count,
   IDENTICAL total amino acid residues, and IDENTICAL exact length distributions.
   Verify closure programmatically via fail-closed assertions.
5. Physically build DIAMOND databases and execute paired searches for:
   - Primary rule (E <= 1e-3, --sensitive -e 1e-3 -k 1) across all 3 seeds (20261013, 20261014, 20261015)
   - Overlap-excluded sensitivity (seed 20261013 with cross-release exact shared sequences excluded)
   - Secondary hit-rule sensitivity (seed 20261013 with --id 30 applied during DIAMOND search)
6. Compute separate dark and control per-stratum and Horvitz-Thompson / stratified population-weighted
   marginal estimates (with finite-population-corrected SEs), never pooling dark+ctrl, and using
   'undefined' (never 'inf') when denominator is 0.
"""

import os
import sys
import gzip
import math
import random
import hashlib
import subprocess
import time
from collections import defaultdict
import numpy as np
from scipy import stats

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")
LENGTH_BINS = [
    ("<100", lambda l: l < 100),
    ("100-149", lambda l: 100 <= l <= 149),
    ("150-249", lambda l: 150 <= l <= 249),
    ("250-399", lambda l: 250 <= l <= 399),
    (">=400", lambda l: l >= 400),
]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def bin_index(length_aa):
    if length_aa < 100:
        return 0
    elif length_aa <= 149:
        return 1
    elif length_aa <= 249:
        return 2
    elif length_aa <= 399:
        return 3
    else:
        return 4


def load_canonical_fasta(path):
    """Load a FASTA file, canonicalizing sequences and auditing stars/case."""
    records = []  # list of (seq_id, canon_seq, seq_sha256)
    raw_count = 0
    raw_residues = 0
    trailing_stars = 0
    internal_stars = 0
    lowercase_chars = 0
    invalid_chars = 0

    cur_id = None
    cur_lines = []

    def flush(sid, lines):
        nonlocal raw_count, raw_residues, trailing_stars, internal_stars, lowercase_chars, invalid_chars
        raw_seq = "".join(lines)
        raw_count += 1
        raw_residues += len(raw_seq)
        if raw_seq.endswith("*"):
            trailing_stars += 1
        stripped = raw_seq.rstrip("*")
        if "*" in stripped:
            internal_stars += 1
            return
        lowercase_chars += sum(1 for c in stripped if c.islower())
        canon = stripped.upper()
        bad = sum(1 for c in canon if c not in VALID_AA)
        invalid_chars += bad
        if not canon or bad > 0:
            return
        sh = hashlib.sha256(canon.encode("ascii")).hexdigest()
        records.append((sid, canon, sh))

    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    flush(cur_id, cur_lines)
                cur_id = line[1:].strip().split()[0]
                cur_lines = []
            else:
                cur_lines.append(line.strip())
        if cur_id is not None:
            flush(cur_id, cur_lines)

    audit = {
        "raw_count": raw_count,
        "raw_residues": raw_residues,
        "trailing_stars": trailing_stars,
        "internal_stars": internal_stars,
        "lowercase_chars": lowercase_chars,
        "invalid_chars": invalid_chars,
        "canonical_count": len(records),
        "canonical_residues": sum(len(r[1]) for r in records),
    }
    return records, audit


def sample_gem_pool_for_seed(gem_faa_dir, all_gem_files, target_unique_count, seed, out_dir):
    """
    Deterministically shuffle GEM .faa.gz files with `seed` and read whole genomes until
    unique canonical sequences >= target_unique_count.
    Write `v3_gem_genome_manifest_seed{seed}.tsv`.
    """
    rng = random.Random(seed)
    shuffled_files = list(all_gem_files)
    rng.shuffle(shuffled_files)

    seen_sha256 = {}  # seq_sha256 -> (seq_id, canon_seq)
    raw_prot_total = 0
    raw_res_total = 0
    trailing_star_total = 0
    internal_star_total = 0
    lowercase_total = 0
    invalid_total = 0
    canon_prot_total = 0
    canon_res_total = 0
    manifest_rows = []

    for fname in shuffled_files:
        fpath = os.path.join(gem_faa_dir, fname)
        f_bytes = os.path.getsize(fpath)
        f_sha256 = sha256_file(fpath)
        genome_id = fname[:-7] if fname.endswith(".faa.gz") else os.path.splitext(fname)[0]

        g_raw_prots = 0
        g_raw_res = 0
        g_trail_stars = 0
        g_intern_stars = 0
        g_lower = 0
        g_invalid = 0
        g_canon_prots = 0
        g_canon_res = 0
        g_within_uniq = set()
        g_new_global_uniq = 0

        cur_id = None
        cur_lines = []

        def flush_gem(sid, lines):
            nonlocal g_raw_prots, g_raw_res, g_trail_stars, g_intern_stars
            nonlocal g_lower, g_invalid, g_canon_prots, g_canon_res, g_new_global_uniq
            raw_seq = "".join(lines)
            g_raw_prots += 1
            g_raw_res += len(raw_seq)
            if raw_seq.endswith("*"):
                g_trail_stars += 1
            stripped = raw_seq.rstrip("*")
            if "*" in stripped:
                g_intern_stars += 1
                return
            g_lower += sum(1 for c in stripped if c.islower())
            canon = stripped.upper()
            bad = sum(1 for c in canon if c not in VALID_AA)
            g_invalid += bad
            if not canon or bad > 0:
                return
            g_canon_prots += 1
            g_canon_res += len(canon)
            sh = hashlib.sha256(canon.encode("ascii")).hexdigest()
            g_within_uniq.add(sh)
            if sh not in seen_sha256:
                seen_sha256[sh] = (sid, canon)
                g_new_global_uniq += 1

        with gzip.open(fpath, "rt") as gz:
            for line in gz:
                if line.startswith(">"):
                    if cur_id is not None:
                        flush_gem(cur_id, cur_lines)
                    cur_id = line[1:].strip().split()[0]
                    cur_lines = []
                else:
                    cur_lines.append(line.strip())
            if cur_id is not None:
                flush_gem(cur_id, cur_lines)

        raw_prot_total += g_raw_prots
        raw_res_total += g_raw_res
        trailing_star_total += g_trail_stars
        internal_star_total += g_intern_stars
        lowercase_total += g_lower
        invalid_total += g_invalid
        canon_prot_total += g_canon_prots
        canon_res_total += g_canon_res

        manifest_rows.append((
            len(manifest_rows) + 1,
            fname,
            genome_id,
            f_bytes,
            f_sha256,
            g_raw_prots,
            g_raw_res,
            g_trail_stars,
            g_intern_stars,
            g_lower,
            g_invalid,
            g_canon_prots,
            g_canon_res,
            len(g_within_uniq),
            g_new_global_uniq,
            len(seen_sha256),
        ))

        if len(seen_sha256) >= target_unique_count:
            break

    man_path = os.path.join(out_dir, f"v3_gem_genome_manifest_seed{seed}.tsv")
    with open(man_path, "w") as out:
        out.write(
            "sample_order\tfilename\tgenome_id\tfile_bytes\tfile_sha256\traw_protein_count\t"
            "raw_residue_count\ttrailing_star_count\tinternal_star_count\tlowercase_char_count\t"
            "invalid_aa_char_count\tcanonical_protein_count\tcanonical_residue_count\t"
            "within_genome_unique_proteins\tnew_pool_unique_proteins\tcumulative_pool_unique_proteins\n"
        )
        for r in manifest_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    pool_stats = {
        "seed": seed,
        "genomes_sampled": len(manifest_rows),
        "raw_proteins": raw_prot_total,
        "raw_residues": raw_res_total,
        "trailing_stars": trailing_star_total,
        "internal_stars": internal_star_total,
        "lowercase_chars": lowercase_total,
        "invalid_chars": invalid_total,
        "canonical_proteins": canon_prot_total,
        "canonical_residues": canon_res_total,
        "unique_canonical_proteins": len(seen_sha256),
        "unique_canonical_residues": sum(len(v[1]) for v in seen_sha256.values()),
        "internal_duplicates": canon_prot_total - len(seen_sha256),
        "internal_redundancy_pct": 100.0 * (canon_prot_total - len(seen_sha256)) / canon_prot_total,
        "manifest_path": man_path,
    }
    return seen_sha256, pool_stats


def build_exact_length_matched_pair(smag_uniq_dict, gem_uniq_dict, seed, exclude_shared, out_dir, label):
    """
    Construct exact integer-length (1 aa) common-support matched FASTA databases for SMAG and GEM.
    Verify programmatically that protein count, residue count, length histogram, and 5-bin counts
    are 100% identical and internal duplicates are 0.
    """
    shared_hashes = set(smag_uniq_dict.keys()).intersection(set(gem_uniq_dict.keys()))

    s_by_len = defaultdict(list)
    for sh, (sid, seq) in smag_uniq_dict.items():
        if exclude_shared and sh in shared_hashes:
            continue
        s_by_len[len(seq)].append((sid, seq, sh))

    g_by_len = defaultdict(list)
    for sh, (sid, seq) in gem_uniq_dict.items():
        if exclude_shared and sh in shared_hashes:
            continue
        g_by_len[len(seq)].append((sid, seq, sh))

    all_lengths = sorted(set(s_by_len.keys()).intersection(set(g_by_len.keys())))
    s_selected = []
    g_selected = []

    for length_aa in all_lengths:
        s_list = sorted(s_by_len[length_aa], key=lambda x: x[0])
        g_list = sorted(g_by_len[length_aa], key=lambda x: x[0])
        q_len = min(len(s_list), len(g_list))
        if q_len == 0:
            continue
        rng_s = random.Random(seed * 100000 + length_aa * 2 + 1)
        rng_g = random.Random(seed * 100000 + length_aa * 2 + 2)
        s_chosen = rng_s.sample(s_list, q_len) if len(s_list) > q_len else s_list
        g_chosen = rng_g.sample(g_list, q_len) if len(g_list) > q_len else g_list
        s_selected.extend(s_chosen)
        g_selected.extend(g_chosen)

    # Sort deterministically by seq_id before writing FASTA
    s_selected.sort(key=lambda x: x[0])
    g_selected.sort(key=lambda x: x[0])

    s_faa = os.path.join(out_dir, f"v3_smag_matched_{label}.faa")
    g_faa = os.path.join(out_dir, f"v3_gem_matched_{label}.faa")

    with open(s_faa, "w") as fs:
        for sid, seq, _ in s_selected:
            fs.write(f">{sid}\n{seq}\n")

    with open(g_faa, "w") as fg:
        for sid, seq, _ in g_selected:
            fg.write(f">{sid}\n{seq}\n")

    # Programmatic closure verification by re-reading the written FASTA files
    s_recs, s_aud = load_canonical_fasta(s_faa)
    g_recs, g_aud = load_canonical_fasta(g_faa)

    s_lens = sorted(len(r[1]) for r in s_recs)
    g_lens = sorted(len(r[1]) for r in g_recs)
    s_uniq_h = {r[2] for r in s_recs}
    g_uniq_h = {r[2] for r in g_recs}
    post_shared = len(s_uniq_h.intersection(g_uniq_h))

    # Fail-closed closure assertions
    assert len(s_recs) == len(g_recs), f"Count mismatch: {len(s_recs)} vs {len(g_recs)}"
    assert len(s_uniq_h) == len(s_recs), f"SMAG duplicates remain: {len(s_recs) - len(s_uniq_h)}"
    assert len(g_uniq_h) == len(g_recs), f"GEM duplicates remain: {len(g_recs) - len(g_uniq_h)}"
    assert s_aud["trailing_stars"] == 0 and g_aud["trailing_stars"] == 0, "Trailing stars remain!"
    assert s_aud["internal_stars"] == 0 and g_aud["internal_stars"] == 0, "Internal stars remain!"
    assert s_aud["canonical_residues"] == g_aud["canonical_residues"], "Residue total mismatch!"
    assert s_lens == g_lens, "Exact integer-length distributions do not match!"
    if exclude_shared:
        assert post_shared == 0, f"Expected 0 shared sequences in overlap-excluded mode, got {post_shared}"

    bins_counts = [0, 0, 0, 0, 0]
    for l in s_lens:
        bins_counts[bin_index(l)] += 1

    arr = np.array(s_lens, dtype=np.float64)
    closure_info = {
        "label": label,
        "seed": seed,
        "exclude_shared": exclude_shared,
        "matched_proteins": len(s_recs),
        "matched_residues": s_aud["canonical_residues"],
        "mean_length_aa": float(np.mean(arr)),
        "sd_length_aa": float(np.std(arr, ddof=1)),
        "median_length_aa": float(np.median(arr)),
        "q10_aa": float(np.percentile(arr, 10)),
        "q25_aa": float(np.percentile(arr, 25)),
        "q75_aa": float(np.percentile(arr, 75)),
        "q90_aa": float(np.percentile(arr, 90)),
        "bin_lt100": bins_counts[0],
        "bin_100_149": bins_counts[1],
        "bin_150_249": bins_counts[2],
        "bin_250_399": bins_counts[3],
        "bin_ge400": bins_counts[4],
        "pre_match_shared_exact_seqs": len(shared_hashes),
        "post_match_shared_exact_seqs": post_shared,
        "smag_faa_path": s_faa,
        "gem_faa_path": g_faa,
        "smag_faa_bytes": os.path.getsize(s_faa),
        "gem_faa_bytes": os.path.getsize(g_faa),
        "smag_faa_sha256": sha256_file(s_faa),
        "gem_faa_sha256": sha256_file(g_faa),
    }
    return closure_info


def run_diamond_paired_search(diamond_bin, s_faa, g_faa, query_faa, tmp_dir, out_dir, run_label, min_id, threads):
    """Build DIAMOND DBs and run paired blastp searches under identical parameters."""
    s_db = os.path.join(tmp_dir, f"v3_smag_{run_label}")
    g_db = os.path.join(tmp_dir, f"v3_gem_{run_label}")

    subprocess.run([diamond_bin, "makedb", "--in", s_faa, "-d", s_db, "-p", str(threads)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([diamond_bin, "makedb", "--in", g_faa, "-d", g_db, "-p", str(threads)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    s_m8 = os.path.join(out_dir, f"v3_smag_hits_{run_label}.m8")
    g_m8 = os.path.join(out_dir, f"v3_gem_hits_{run_label}.m8")

    outfmt = ["6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    cmd_s = [diamond_bin, "blastp", "-q", query_faa, "-d", s_db, "--sensitive", "-e", "1e-3", "-k", "1", "-p", str(threads)]
    cmd_g = [diamond_bin, "blastp", "-q", query_faa, "-d", g_db, "--sensitive", "-e", "1e-3", "-k", "1", "-p", str(threads)]
    if min_id is not None:
        cmd_s.extend(["--id", str(min_id)])
        cmd_g.extend(["--id", str(min_id)])
    cmd_s.extend(["-f"] + outfmt + ["-o", s_m8])
    cmd_g.extend(["-f"] + outfmt + ["-o", g_m8])

    subprocess.run(cmd_s, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run(cmd_g, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    for dm in [s_db + ".dmnd", g_db + ".dmnd"]:
        if os.path.exists(dm):
            os.remove(dm)

    return s_m8, g_m8


def evaluate_paired_search(s_m8, g_m8, queries, pop_stratum_counts, run_label, seed, hit_rule_desc, overlap_policy):
    """
    Evaluate paired outcomes separately for dark and ctrl cohorts.
    Compute per-stratum statistics and Horvitz-Thompson / stratified population-weighted
    marginal statistics with finite-population-corrected SEs.
    Never pool dark+ctrl; format 0/0 ratios as 'undefined'.
    """
    s_hits = {}
    with open(s_m8) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            s_hits[parts[0]] = (parts[1], float(parts[2]), float(parts[10]), float(parts[11]))

    g_hits = {}
    with open(g_m8) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            g_hits[parts[0]] = (parts[1], float(parts[2]), float(parts[10]), float(parts[11]))

    stratum_rows = []
    cohort_rows = []

    for cohort in ["dark", "ctrl"]:
        N_pop_total = sum(pop_stratum_counts[cohort].values())
        wt_p_smag = 0.0
        wt_p_gem = 0.0
        wt_p_both = 0.0
        wt_p_smag_only = 0.0
        wt_p_gem_only = 0.0
        wt_p_neither = 0.0
        wt_delta = 0.0
        wt_var_delta = 0.0
        wt_var_smag = 0.0
        wt_var_gem = 0.0

        panel_both = 0
        panel_s_only = 0
        panel_g_only = 0
        panel_neither = 0
        panel_n = 0

        for bin_name, _ in LENGTH_BINS:
            q_bin = [q for q in queries if q["cohort"] == cohort and q["lbin"] == bin_name]
            n_h = len(q_bin)
            N_h = pop_stratum_counts[cohort][bin_name]
            pi_h = n_h / N_h
            W_h = N_h / N_pop_total
            fpc_h = 1.0 - pi_h

            a_h = sum(1 for q in q_bin if q["qid"] in s_hits and q["qid"] in g_hits)
            b_h = sum(1 for q in q_bin if q["qid"] in s_hits and q["qid"] not in g_hits)
            c_h = sum(1 for q in q_bin if q["qid"] not in s_hits and q["qid"] in g_hits)
            d_h = sum(1 for q in q_bin if q["qid"] not in s_hits and q["qid"] not in g_hits)

            p_s_h = (a_h + b_h) / n_h
            p_g_h = (a_h + c_h) / n_h
            delta_h = (b_h - c_h) / n_h

            # Sample variance of paired differences d_hi in {-1, 0, +1}
            # sum(d_hi^2) = b_h + c_h
            s2_d_h = ((b_h + c_h) - n_h * (delta_h ** 2)) / (n_h - 1) if n_h > 1 else 0.0
            s2_d_h = max(0.0, s2_d_h)
            var_delta_h = fpc_h * s2_d_h / n_h
            se_delta_h = math.sqrt(var_delta_h)

            s2_s_h = n_h * p_s_h * (1.0 - p_s_h) / (n_h - 1) if n_h > 1 else 0.0
            s2_g_h = n_h * p_g_h * (1.0 - p_g_h) / (n_h - 1) if n_h > 1 else 0.0
            var_s_h = fpc_h * s2_s_h / n_h
            var_g_h = fpc_h * s2_g_h / n_h

            # Exact binomial McNemar two-sided p-value and continuity-corrected chi2
            if (b_h + c_h) > 0:
                binom_p = float(stats.binomtest(b_h, b_h + c_h, 0.5, alternative="two-sided").pvalue)
                chi2_cc = ((abs(b_h - c_h) - 1.0) ** 2) / (b_h + c_h)
            else:
                binom_p = 1.0
                chi2_cc = 0.0

            if c_h > 0:
                disc_ratio_str = f"{b_h / c_h:.2f}"
            else:
                disc_ratio_str = "undefined"

            stratum_rows.append({
                "run_label": run_label,
                "seed": seed,
                "hit_rule": hit_rule_desc,
                "overlap_policy": overlap_policy,
                "cohort": cohort,
                "length_bin": bin_name,
                "n_h": n_h,
                "N_h": N_h,
                "inclusion_prob_pi_h": f"{pi_h:.6f}",
                "pop_weight_W_h": f"{W_h:.6f}",
                "both": a_h,
                "smag_only": b_h,
                "gem_only": c_h,
                "neither": d_h,
                "smag_hits": a_h + b_h,
                "gem_hits": a_h + c_h,
                "smag_hit_rate_pct": f"{100.0 * p_s_h:.2f}",
                "gem_hit_rate_pct": f"{100.0 * p_g_h:.2f}",
                "delta_pct": f"{100.0 * delta_h:.2f}",
                "se_delta_fpc_pct": f"{100.0 * se_delta_h:.2f}",
                "ci95_low_pct": f"{100.0 * (delta_h - 1.96 * se_delta_h):.2f}",
                "ci95_high_pct": f"{100.0 * (delta_h + 1.96 * se_delta_h):.2f}",
                "mcnemar_chi2_cc": f"{chi2_cc:.2f}",
                "mcnemar_exact_pval": f"{binom_p:.4e}",
                "discordant_ratio_b_over_c": disc_ratio_str,
            })

            # Accumulate population-weighted marginal estimators
            wt_p_smag += W_h * p_s_h
            wt_p_gem += W_h * p_g_h
            wt_p_both += W_h * (a_h / n_h)
            wt_p_smag_only += W_h * (b_h / n_h)
            wt_p_gem_only += W_h * (c_h / n_h)
            wt_p_neither += W_h * (d_h / n_h)
            wt_delta += W_h * delta_h
            wt_var_delta += (W_h ** 2) * var_delta_h
            wt_var_smag += (W_h ** 2) * var_s_h
            wt_var_gem += (W_h ** 2) * var_g_h

            panel_both += a_h
            panel_s_only += b_h
            panel_g_only += c_h
            panel_neither += d_h
            panel_n += n_h

        wt_se_delta = math.sqrt(wt_var_delta)
        wt_se_smag = math.sqrt(wt_var_smag)
        wt_se_gem = math.sqrt(wt_var_gem)
        wt_z = (wt_delta / wt_se_delta) if wt_se_delta > 0 else 0.0
        wt_pval = float(2.0 * stats.norm.sf(abs(wt_z))) if wt_se_delta > 0 else 1.0
        wt_disc_ratio = f"{wt_p_smag_only / wt_p_gem_only:.2f}" if wt_p_gem_only > 0 else "undefined"

        panel_s_hits = panel_both + panel_s_only
        panel_g_hits = panel_both + panel_g_only
        panel_s_rate = 100.0 * panel_s_hits / panel_n
        panel_g_rate = 100.0 * panel_g_hits / panel_n
        panel_delta = 100.0 * (panel_s_only - panel_g_only) / panel_n
        if (panel_s_only + panel_g_only) > 0:
            panel_chi2 = ((abs(panel_s_only - panel_g_only) - 1.0) ** 2) / (panel_s_only + panel_g_only)
            panel_exact_p = float(stats.binomtest(panel_s_only, panel_s_only + panel_g_only, 0.5).pvalue)
        else:
            panel_chi2 = 0.0
            panel_exact_p = 1.0
        panel_disc_ratio = f"{panel_s_only / panel_g_only:.2f}" if panel_g_only > 0 else "undefined"

        cohort_rows.append({
            "run_label": run_label,
            "seed": seed,
            "hit_rule": hit_rule_desc,
            "overlap_policy": overlap_policy,
            "cohort": cohort,
            "N_pop": N_pop_total,
            "n_panel": panel_n,
            "pop_wt_smag_hit_rate_pct": f"{100.0 * wt_p_smag:.2f}",
            "pop_wt_smag_se_pct": f"{100.0 * wt_se_smag:.2f}",
            "pop_wt_gem_hit_rate_pct": f"{100.0 * wt_p_gem:.2f}",
            "pop_wt_gem_se_pct": f"{100.0 * wt_se_gem:.2f}",
            "pop_wt_both_pct": f"{100.0 * wt_p_both:.2f}",
            "pop_wt_smag_only_pct": f"{100.0 * wt_p_smag_only:.2f}",
            "pop_wt_gem_only_pct": f"{100.0 * wt_p_gem_only:.2f}",
            "pop_wt_neither_pct": f"{100.0 * wt_p_neither:.2f}",
            "pop_wt_delta_pct": f"{100.0 * wt_delta:.2f}",
            "pop_wt_se_delta_fpc_pct": f"{100.0 * wt_se_delta:.2f}",
            "pop_wt_ci95_low_pct": f"{100.0 * (wt_delta - 1.96 * wt_se_delta):.2f}",
            "pop_wt_ci95_high_pct": f"{100.0 * (wt_delta + 1.96 * wt_se_delta):.2f}",
            "pop_wt_wald_z": f"{wt_z:.2f}",
            "pop_wt_wald_pval": f"{wt_pval:.4e}",
            "pop_wt_discordant_ratio": wt_disc_ratio,
            "panel_unwt_both": panel_both,
            "panel_unwt_smag_only": panel_s_only,
            "panel_unwt_gem_only": panel_g_only,
            "panel_unwt_neither": panel_neither,
            "panel_unwt_smag_hits": panel_s_hits,
            "panel_unwt_gem_hits": panel_g_hits,
            "panel_unwt_smag_hit_rate_pct": f"{panel_s_rate:.2f}",
            "panel_unwt_gem_hit_rate_pct": f"{panel_g_rate:.2f}",
            "panel_unwt_delta_pct": f"{panel_delta:.2f}",
            "panel_mcnemar_chi2_cc": f"{panel_chi2:.2f}",
            "panel_mcnemar_exact_pval": f"{panel_exact_p:.4e}",
            "panel_unwt_discordant_ratio": panel_disc_ratio,
        })

    return stratum_rows, cohort_rows


def decompose_top1_vs_dbsize_bias(out_dir):
    """
    Quantify how much of the 03_pilot_subsample_search.py discrepancy (110 dark, 895 ctrl)
    is genuine top-1 masking at full-database E-value threshold vs database-size E-value rescaling.
    Full GEM letters: ~39.2B (122,978,191 proteins); sub2500 letters: 1,835,451,488 (5,815,177 proteins),
    protein ratio = 122978191 / 5815177 = 21.1478x.
    """
    comp_path = os.path.join(out_dir, "pilot_query_comparison_details.tsv")
    if not os.path.exists(comp_path):
        return
    scale_factor = 122978191.0 / 5815177.0  # 21.1478x
    strict_e_cutoff = 1e-3 / scale_factor    # 4.7286e-5

    rows_by_cohort = defaultdict(list)
    with open(comp_path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            d = dict(zip(header, parts))
            rows_by_cohort[d["cohort"]].append(d)

    out_path = os.path.join(out_dir, "v3_top1_vs_dbsize_decomposition.tsv")
    with open(out_path, "w") as out:
        out.write(
            "cohort\tn_queries\tempirical_sub2500_hits\tposthoc_filtered_hits\ttotal_discrepancy_hits\t"
            "discrepancy_with_other_full_gem_top1\tdiscrepancy_with_no_full_gem_hit\t"
            "strict_top1_masked_e_le_4p73e5\tborderline_dbsize_rescaled_e_gt_4p73e5\t"
            "strict_top1_masked_share_pct\tborderline_dbsize_rescaled_share_pct\n"
        )
        for cohort in ["dark", "ctrl"]:
            c_rows = rows_by_cohort[cohort]
            n_tot = len(c_rows)
            emp = [r for r in c_rows if r["has_empirical"] == "True"]
            post = [r for r in c_rows if r["has_posthoc"] == "True"]
            disc = [r for r in c_rows if r["category"] == "masked_alternative_hit"]
            has_other_full = [r for r in disc if r["full_gem_subject"] != "NONE"]
            no_full = [r for r in disc if r["full_gem_subject"] == "NONE"]
            strict_masked = [r for r in disc if float(r["empirical_evalue"]) <= strict_e_cutoff]
            borderline = [r for r in disc if float(r["empirical_evalue"]) > strict_e_cutoff]
            n_disc = len(disc)
            s_pct = 100.0 * len(strict_masked) / n_disc if n_disc > 0 else 0.0
            b_pct = 100.0 * len(borderline) / n_disc if n_disc > 0 else 0.0
            out.write(
                f"{cohort}\t{n_tot}\t{len(emp)}\t{len(post)}\t{n_disc}\t"
                f"{len(has_other_full)}\t{len(no_full)}\t{len(strict_masked)}\t{len(borderline)}\t"
                f"{s_pct:.2f}\t{b_pct:.2f}\n"
            )


def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    base_dir = "/ibex/scratch/projects/c2014/rob/metagenomics"
    gem_faa_dir = os.path.join(base_dir, "databases/gem/faa/faa")
    spec_dir = os.path.join(base_dir, "06_darkmatter/specificity_clean/results")
    basis_dir = os.path.join(base_dir, "06_darkmatter/darkfam_basis_clean/results")
    h3_dir = os.path.join(base_dir, "06_darkmatter/p6_hypotheses_clean/h3_references")
    out_dir = os.path.join(h3_dir, "results")
    tmp_dir = os.path.join(h3_dir, "tmp")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    diamond_bin = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"
    if not os.path.exists(diamond_bin):
        diamond_bin = "diamond"
    dmnd_ver = subprocess.check_output([diamond_bin, "version"], text=True).strip()
    print(f"=== H3 v3 Step 2: Canonicalization, Exact Length Matching, and Paired Searches ({dmnd_ver}) ===")

    # 1. Compute exact population stratum counts N_h from FAMILIES_FINAL.txt, CONTROL_FAMILIES.txt, and query_all.faa
    print("Computing exact population length-stratum counts N_h from upstream query_all.faa...")
    q_all_recs, q_all_aud = load_canonical_fasta(os.path.join(spec_dir, "query_all.faa"))
    q_all_len = {sid: len(seq) for sid, seq, _ in q_all_recs}
    print(f"  query_all.faa audit: {q_all_aud}")

    with open(os.path.join(basis_dir, "FAMILIES_FINAL.txt")) as f:
        dark_pop_ids = [l.strip() for l in f if l.strip()]
    with open(os.path.join(basis_dir, "CONTROL_FAMILIES.txt")) as f:
        ctrl_pop_ids = [l.strip() for l in f if l.strip()]

    pop_stratum_counts = {"dark": {}, "ctrl": {}}
    for bname, bfunc in LENGTH_BINS:
        pop_stratum_counts["dark"][bname] = sum(1 for qid in dark_pop_ids if bfunc(q_all_len[qid]))
        pop_stratum_counts["ctrl"][bname] = sum(1 for qid in ctrl_pop_ids if bfunc(q_all_len[qid]))
    print("  Dark population N_h:", pop_stratum_counts["dark"], "Total:", sum(pop_stratum_counts["dark"].values()))
    print("  Ctrl population N_h:", pop_stratum_counts["ctrl"], "Total:", sum(pop_stratum_counts["ctrl"].values()))

    # Verify and canonicalize pilot_query_2000.faa
    orig_query_faa = os.path.join(out_dir, "pilot_query_2000.faa")
    q2000_recs, q2000_aud = load_canonical_fasta(orig_query_faa)
    print(f"  pilot_query_2000.faa audit: {q2000_aud}")
    v3_query_faa = os.path.join(out_dir, "v3_pilot_query_2000_canonical.faa")
    with open(v3_query_faa, "w") as out_q:
        for sid, seq, _ in q2000_recs:
            out_q.write(f">{sid}\n{seq}\n")

    # Load pilot query manifest
    queries = []
    with open(os.path.join(out_dir, "pilot_query_manifest.tsv")) as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            queries.append({"qid": parts[0], "cohort": parts[1], "lbin": parts[2], "laa": int(parts[3])})

    # 2. Load SMAG canonical and unique pools
    smag_canon_path = os.path.join(out_dir, "v3_smag_pool_canonical.faa")
    smag_uniq_path = os.path.join(out_dir, "v3_smag_pool_unique.faa")
    smag_canon_recs, smag_canon_aud = load_canonical_fasta(smag_canon_path)
    smag_uniq_recs, smag_uniq_aud = load_canonical_fasta(smag_uniq_path)
    smag_uniq_dict = {sh: (sid, seq) for sid, seq, sh in smag_uniq_recs}
    target_uniq = len(smag_uniq_dict)
    print(f"Loaded SMAG unique pool: {target_uniq} unique canonical proteins "
          f"(from {smag_canon_aud['canonical_count']} canonical proteins, "
          f"{smag_canon_aud['canonical_count'] - target_uniq} duplicates).")

    # 3. Sample GEM candidate pools for 3 deterministic seeds and match exact length distributions
    all_gem_files = sorted(f for f in os.listdir(gem_faa_dir) if f.endswith(".faa.gz"))
    assert len(all_gem_files) == 52515, f"Expected 52515 GEM .faa.gz files, found {len(all_gem_files)}"

    seeds = [20261013, 20261014, 20261015]
    shared_seq_rows = []
    closure_rows = []
    all_stratum_rows = []
    all_cohort_rows = []
    gem_pool_stats_by_seed = {}

    for seed in seeds:
        print(f"--- Sampling GEM candidate pool for seed={seed} (target >= {target_uniq} unique proteins) ---")
        t_s = time.time()
        gem_uniq_dict, g_pool_st = sample_gem_pool_for_seed(gem_faa_dir, all_gem_files, target_uniq, seed, out_dir)
        gem_pool_stats_by_seed[seed] = g_pool_st
        print(f"  [Seed {seed}] Sampled {g_pool_st['genomes_sampled']} GEM MAGs in {time.time() - t_s:.1f}s: "
              f"raw={g_pool_st['raw_proteins']}, trailing_stars={g_pool_st['trailing_stars']}, "
              f"internal_stars={g_pool_st['internal_stars']}, unique_canonical={g_pool_st['unique_canonical_proteins']}, "
              f"internal_redundancy={g_pool_st['internal_redundancy_pct']:.2f}%")

        # Audit exact cross-release canonical sequence overlap
        shared_hashes = sorted(set(smag_uniq_dict.keys()).intersection(set(gem_uniq_dict.keys())))
        print(f"  [Seed {seed}] Exact canonical shared sequences with SMAG pool: {len(shared_hashes)}")
        for sh in shared_hashes:
            s_id, s_seq = smag_uniq_dict[sh]
            g_id, _ = gem_uniq_dict[sh]
            shared_seq_rows.append((seed, sh, len(s_seq), s_id, g_id))

        # Construct primary exact-length-matched pair (overlap retained)
        run_label = f"seed{seed}"
        c_info = build_exact_length_matched_pair(
            smag_uniq_dict, gem_uniq_dict, seed=seed, exclude_shared=False, out_dir=out_dir, label=run_label
        )
        c_info["gem_genomes_sampled"] = g_pool_st["genomes_sampled"]
        c_info["gem_pool_raw_proteins"] = g_pool_st["raw_proteins"]
        c_info["gem_pool_trailing_stars"] = g_pool_st["trailing_stars"]
        c_info["gem_pool_internal_stars"] = g_pool_st["internal_stars"]
        c_info["gem_pool_unique_proteins"] = g_pool_st["unique_canonical_proteins"]
        c_info["gem_pool_internal_redundancy_pct"] = g_pool_st["internal_redundancy_pct"]
        closure_rows.append(c_info)
        print(f"  [Seed {seed} Exact-Length Match] Matched proteins={c_info['matched_proteins']} (SMAG==GEM), "
              f"Matched residues={c_info['matched_residues']} (SMAG==GEM), "
              f"Mean={c_info['mean_length_aa']:.2f} aa, Median={c_info['median_length_aa']:.1f} aa.")

        # Run Primary Paired DIAMOND Search (E <= 1e-3)
        s_m8, g_m8 = run_diamond_paired_search(
            diamond_bin, c_info["smag_faa_path"], c_info["gem_faa_path"],
            v3_query_faa, tmp_dir, out_dir, run_label=run_label, min_id=None, threads=threads
        )
        st_rows, co_rows = evaluate_paired_search(
            s_m8, g_m8, queries, pop_stratum_counts,
            run_label=run_label, seed=seed, hit_rule_desc="E<=1e-3", overlap_policy="overlap_retained"
        )
        all_stratum_rows.extend(st_rows)
        all_cohort_rows.extend(co_rows)
        for cr in co_rows:
            print(f"    [{run_label} {cr['cohort']}] Pop-Weighted: SMAG={cr['pop_wt_smag_hit_rate_pct']}%, "
                  f"GEM={cr['pop_wt_gem_hit_rate_pct']}%, Delta={cr['pop_wt_delta_pct']}% "
                  f"(SE={cr['pop_wt_se_delta_fpc_pct']}%, 95%CI=[{cr['pop_wt_ci95_low_pct']}, {cr['pop_wt_ci95_high_pct']}%]) | "
                  f"Panel(n=1000): b={cr['panel_unwt_smag_only']}, c={cr['panel_unwt_gem_only']}, p={cr['panel_mcnemar_exact_pval']}")

        # For primary seed 20261013, also run (a) overlap-excluded sensitivity and (b) --id 30 search sensitivity
        if seed == 20261013:
            # (a) Overlap-excluded sensitivity
            no_ov_label = "seed20261013_no_overlap"
            c_info_no_ov = build_exact_length_matched_pair(
                smag_uniq_dict, gem_uniq_dict, seed=seed, exclude_shared=True, out_dir=out_dir, label=no_ov_label
            )
            c_info_no_ov["gem_genomes_sampled"] = g_pool_st["genomes_sampled"]
            c_info_no_ov["gem_pool_raw_proteins"] = g_pool_st["raw_proteins"]
            c_info_no_ov["gem_pool_trailing_stars"] = g_pool_st["trailing_stars"]
            c_info_no_ov["gem_pool_internal_stars"] = g_pool_st["internal_stars"]
            c_info_no_ov["gem_pool_unique_proteins"] = g_pool_st["unique_canonical_proteins"]
            c_info_no_ov["gem_pool_internal_redundancy_pct"] = g_pool_st["internal_redundancy_pct"]
            closure_rows.append(c_info_no_ov)

            s_m8_no, g_m8_no = run_diamond_paired_search(
                diamond_bin, c_info_no_ov["smag_faa_path"], c_info_no_ov["gem_faa_path"],
                v3_query_faa, tmp_dir, out_dir, run_label=no_ov_label, min_id=None, threads=threads
            )
            st_no, co_no = evaluate_paired_search(
                s_m8_no, g_m8_no, queries, pop_stratum_counts,
                run_label=no_ov_label, seed=seed, hit_rule_desc="E<=1e-3", overlap_policy="overlap_excluded"
            )
            all_stratum_rows.extend(st_no)
            all_cohort_rows.extend(co_no)

            # (b) Secondary Hit Rule Sensitivity: E <= 1e-3 AND pident >= 30% applied during DIAMOND search (--id 30)
            id30_label = "seed20261013_id30"
            s_m8_id30, g_m8_id30 = run_diamond_paired_search(
                diamond_bin, c_info["smag_faa_path"], c_info["gem_faa_path"],
                v3_query_faa, tmp_dir, out_dir, run_label=id30_label, min_id=30, threads=threads
            )
            st_id30, co_id30 = evaluate_paired_search(
                s_m8_id30, g_m8_id30, queries, pop_stratum_counts,
                run_label=id30_label, seed=seed, hit_rule_desc="E<=1e-3_and_pident>=30", overlap_policy="overlap_retained"
            )
            all_stratum_rows.extend(st_id30)
            all_cohort_rows.extend(co_id30)

    # 4. Write output tables
    shared_path = os.path.join(out_dir, "v3_cross_release_shared_sequences.tsv")
    with open(shared_path, "w") as out:
        out.write("seed\tseq_sha256\tlength_aa\tsmag_representative_id\tgem_representative_id\n")
        for r in shared_seq_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    closure_path = os.path.join(out_dir, "v3_matching_closure_audit.tsv")
    closure_cols = [
        "label", "seed", "exclude_shared", "gem_genomes_sampled",
        "gem_pool_raw_proteins", "gem_pool_trailing_stars", "gem_pool_internal_stars",
        "gem_pool_unique_proteins", "gem_pool_internal_redundancy_pct",
        "pre_match_shared_exact_seqs", "post_match_shared_exact_seqs",
        "matched_proteins", "matched_residues", "mean_length_aa", "sd_length_aa",
        "median_length_aa", "q10_aa", "q25_aa", "q75_aa", "q90_aa",
        "bin_lt100", "bin_100_149", "bin_150_249", "bin_250_399", "bin_ge400",
        "smag_faa_bytes", "smag_faa_sha256", "gem_faa_bytes", "gem_faa_sha256",
    ]
    with open(closure_path, "w") as out:
        out.write("\t".join(closure_cols) + "\n")
        for r in closure_rows:
            vals = []
            for c in closure_cols:
                v = r[c]
                if isinstance(v, float):
                    vals.append(f"{v:.2f}")
                else:
                    vals.append(str(v))
            out.write("\t".join(vals) + "\n")

    strat_path = os.path.join(out_dir, "v3_matched_length_stratified_results.tsv")
    strat_cols = list(all_stratum_rows[0].keys())
    with open(strat_path, "w") as out:
        out.write("\t".join(strat_cols) + "\n")
        for r in all_stratum_rows:
            out.write("\t".join(str(r[c]) for c in strat_cols) + "\n")

    cohort_path = os.path.join(out_dir, "v3_matched_cohort_summary.tsv")
    cohort_cols = list(all_cohort_rows[0].keys())
    with open(cohort_path, "w") as out:
        out.write("\t".join(cohort_cols) + "\n")
        for r in all_cohort_rows:
            out.write("\t".join(str(r[c]) for c in cohort_cols) + "\n")

    # 5. Decompose top-1 vs database-size rescaling from 03_pilot_subsample_search.py
    decompose_top1_vs_dbsize_bias(out_dir)

    # 6. Compute SHA-256 checksums of all v3 and preserved v2 files in results/
    checksums_path = os.path.join(out_dir, "v3_artifact_checksums.tsv")
    all_res_files = sorted(f for f in os.listdir(out_dir) if os.path.isfile(os.path.join(out_dir, f)) and f != "v3_artifact_checksums.tsv")
    with open(checksums_path, "w") as out:
        out.write("filename\tbytes\tsha256\n")
        for fname in all_res_files:
            fp = os.path.join(out_dir, fname)
            out.write(f"{fname}\t{os.path.getsize(fp)}\t{sha256_file(fp)}\n")

    print("DONE_09_V3_MATCH_GEM_AND_SEARCH")


if __name__ == "__main__":
    main()
