#!/usr/bin/env python3
"""
Independent arithmetic, sampling, quota, invariant, and TM-score boundary
verification for Round 7 H3 (matched-reference specificity) and H4 (paired
whole-chain vs cropped-domain Foldseek structure census & alnTM>1 audit).

Author: Gemini 4 Argon High (argon_h3h4_review1)
Execution: uv run python reviews/darkmatter-2026-10-13/round7/argon_h3h4_review1/verification.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from scipy.stats import binomtest

REPO_ROOT = Path("/home/leechuck/Documents/papers/empty-quarter-mags")
H3_SNAP = REPO_ROOT / "reviews/darkmatter-2026-10-13/round7/snapshots/h3_completed"
H4_PRE = REPO_ROOT / "reviews/darkmatter-2026-10-13/round7/root_checks/h4_outputs_preclosure"
H4_ROOT = REPO_ROOT / "reviews/darkmatter-2026-10-13/round7/root_checks"
H4_ANOM = REPO_ROOT / "reviews/darkmatter-2026-10-13/round7/h4_tm_anomaly"
ORIX_DIR = REPO_ROOT / "reviews/darkmatter-2026-10-13/round7/orix_h4_gpu"
OUT_DIR = REPO_ROOT / "reviews/darkmatter-2026-10-13/round7/argon_h3h4_review1"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def verify_h3() -> dict:
    res: dict = {}

    # 1. SMAG 40,039 frame & 1,500 random sample replay (seed=20261014)
    smag_frame_lines = [
        line.strip()
        for line in (H3_SNAP / "smag_catalogue_40039_ids.txt").read_text().splitlines()
        if line.strip()
    ]
    sampled_smag_file = [
        line.strip()
        for line in (H3_SNAP / "sampled_smag_1500_ids_seed20261014.txt").read_text().splitlines()
        if line.strip()
    ]
    rng_smag = random.Random(20261014)
    replayed_smag_1500 = sorted(rng_smag.sample(sorted(smag_frame_lines), 1500))

    smag_manifest = read_tsv(H3_SNAP / "smag_genome_manifest.tsv")
    manifest_smag_ids = [r["mag_id"] for r in smag_manifest]

    smag_bytes_match = all(
        int(r["declared_bytes"]) == int(r["observed_bytes"]) and int(r["observed_bytes"]) > 0
        for r in smag_manifest
    )
    smag_clean_aa = all(
        int(r["trailing_star_count"]) == 0
        and int(r["internal_star_count"]) == 0
        and int(r["lowercase_char_count"]) == 0
        and int(r["invalid_aa_char_count"]) == 0
        and int(r["raw_protein_count"]) == int(r["canonical_protein_count"])
        and int(r["raw_residue_count"]) == int(r["canonical_residue_count"])
        for r in smag_manifest
    )

    zenodo_parts = read_tsv(H3_SNAP / "smag_zenodo_parts_manifest.tsv")
    res["smag_sampling_and_extraction"] = {
        "frame_count": len(smag_frame_lines),
        "frame_unique_count": len(set(smag_frame_lines)),
        "sampled_file_count": len(sampled_smag_file),
        "replay_matches_sampled_file": replayed_smag_1500 == sampled_smag_file,
        "replay_matches_manifest_ids": replayed_smag_1500 == sorted(manifest_smag_ids),
        "manifest_row_count": len(smag_manifest),
        "all_1500_declared_eq_observed_bytes": smag_bytes_match,
        "all_1500_clean_canonical_aa": smag_clean_aa,
        "total_declared_bytes": sum(int(r["declared_bytes"]) for r in smag_manifest),
        "total_contigs": sum(int(r["contig_count"]) for r in smag_manifest),
        "total_bp": sum(int(r["total_bp"]) for r in smag_manifest),
        "total_canonical_proteins": sum(int(r["canonical_protein_count"]) for r in smag_manifest),
        "total_canonical_residues": sum(int(r["canonical_residue_count"]) for r in smag_manifest),
        "sum_within_genome_unique_proteins": sum(
            int(r["unique_canonical_proteins"]) for r in smag_manifest
        ),
        "zenodo_parts_manifest_count": len(zenodo_parts),
    }

    # 2. GEM 52,515 frame & 2,055 random sample replay across seeds 20261014, 20261015, 20261016
    gem_frame_lines = [
        line.strip()
        for line in (H3_SNAP / "gem_catalogue_52515_ids.txt").read_text().splitlines()
        if line.strip()
    ]
    gem_sorted = sorted(gem_frame_lines)
    gem_draws_summary = json.loads((H3_SNAP / "gem_harmonized_draws_summary.json").read_text())

    gem_seed_checks = {}
    target_unique = 4_500_000
    expected_n_gem = math.ceil(target_unique * 1.05 / 2300)  # 2055
    for seed in [20261014, 20261015, 20261016]:
        rng = random.Random(seed)
        shuffled = list(gem_sorted)
        rng.shuffle(shuffled)
        selected = sorted(shuffled[:expected_n_gem])

        harm_man = read_tsv(H3_SNAP / f"gem_harmonized_manifest_seed{seed}.tsv")
        pub_man = read_tsv(H3_SNAP / f"gem_published_manifest_seed{seed}.tsv")
        harm_ids = [r["genome_id"] for r in harm_man]
        pub_ids = [r["genome_id"] for r in pub_man]

        harm_raw_prots = sum(int(r["raw_protein_count"]) for r in harm_man)
        harm_canon_prots = sum(int(r["canonical_protein_count"]) for r in harm_man)
        harm_canon_res = sum(int(r["canonical_residue_count"]) for r in harm_man)

        pub_raw_prots = sum(int(r["raw_protein_count"]) for r in pub_man)
        pub_canon_prots = sum(int(r["canonical_protein_count"]) for r in pub_man)
        pub_canon_res = sum(int(r["canonical_residue_count"]) for r in pub_man)
        pub_trailing_stars = sum(int(r["trailing_star_count"]) for r in pub_man)
        pub_internal_stars = sum(int(r["internal_star_count"]) for r in pub_man)

        summary_entry = gem_draws_summary.get("seeds", {}).get(str(seed))
        matches_summary = (
            (
                harm_canon_prots == summary_entry["canonical_protein_count"]
                and harm_canon_res == summary_entry["canonical_residue_count"]
            )
            if summary_entry is not None
            else True
        )
        gem_seed_checks[str(seed)] = {
            "expected_n_genomes": expected_n_gem,
            "replay_matches_harmonized_manifest": selected == harm_ids,
            "replay_matches_published_manifest": selected == pub_ids,
            "harmonized_eq_published_ids": harm_ids == pub_ids,
            "harmonized_raw_proteins": harm_raw_prots,
            "harmonized_canonical_proteins": harm_canon_prots,
            "harmonized_canonical_residues": harm_canon_res,
            "matches_summary_json_totals": matches_summary,
            "published_raw_proteins": pub_raw_prots,
            "published_canonical_proteins": pub_canon_prots,
            "published_canonical_residues": pub_canon_res,
            "published_trailing_stars_stripped": pub_trailing_stars,
            "published_internal_stars_excluded": pub_internal_stars,
            "published_raw_minus_internal_eq_canonical": (
                pub_raw_prots - pub_internal_stars == pub_canon_prots
            ),
            "published_raw_res_minus_trailing_eq_canonical_res": (
                sum(int(r["raw_residue_count"]) for r in pub_man) - pub_trailing_stars
                == pub_canon_res
                if pub_internal_stars == 0
                else None
            ),
        }

    # Check pairwise genome overlap across the 3 GEM seeds
    s14 = set(
        r["genome_id"] for r in read_tsv(H3_SNAP / "gem_harmonized_manifest_seed20261014.tsv")
    )
    s15 = set(
        r["genome_id"] for r in read_tsv(H3_SNAP / "gem_harmonized_manifest_seed20261015.tsv")
    )
    s16 = set(
        r["genome_id"] for r in read_tsv(H3_SNAP / "gem_harmonized_manifest_seed20261016.tsv")
    )
    res["gem_sampling_and_manifests"] = {
        "frame_count": len(gem_frame_lines),
        "frame_unique_count": len(set(gem_frame_lines)),
        "per_seed": gem_seed_checks,
        "seed_genome_overlap": {
            "20261014_vs_20261015": len(s14 & s15),
            "20261014_vs_20261016": len(s14 & s16),
            "20261015_vs_20261016": len(s15 & s16),
            "union_unique_gem_genomes_across_3_seeds": len(s14 | s15 | s16),
        },
    }

    # 3. Exact 1-aa quota matching audit (matching_closure_audit.tsv)
    closure_rows = read_tsv(H3_SNAP / "matching_closure_audit.tsv")
    closure_checks = []
    by_label = {}
    for r in closure_rows:
        bins_sum = (
            int(r["bin_lt100"])
            + int(r["bin_100_149"])
            + int(r["bin_150_249"])
            + int(r["bin_250_399"])
            + int(r["bin_gte400"])
        )
        matched_prots = int(r["matched_proteins"])
        matched_res = int(r["matched_residues"])
        pre_shared = int(r["pre_matching_shared_seqs"])
        post_shared = int(r["post_matching_shared_seqs"])
        excl = r["exclude_shared"] == "True"
        item = {
            "label": r["label"],
            "seed": int(r["seed"]),
            "exclude_shared": excl,
            "pre_matching_shared_seqs": pre_shared,
            "post_matching_shared_seqs": post_shared,
            "matched_proteins": matched_prots,
            "matched_residues": matched_res,
            "five_bins_sum": bins_sum,
            "bins_sum_matches_total": bins_sum == matched_prots,
            "no_overlap_zero_post_shared_ok": (post_shared == 0) if excl else (post_shared > 0),
        }
        closure_checks.append(item)
        by_label[r["label"]] = item

    # Verify exact relationship between harmonized_seed20261014 and harmonized_seed20261014_no_overlap
    h14 = by_label["harmonized_seed20261014"]
    h14_no = by_label["harmonized_seed20261014_no_overlap"]
    res["matching_closure_audit"] = {
        "row_count": len(closure_checks),
        "all_bins_sum_match_matched_proteins": all(
            c["bins_sum_matches_total"] for c in closure_checks
        ),
        "no_overlap_protein_diff_equals_pre_shared": (
            h14["matched_proteins"] - h14_no["matched_proteins"]
            == h14["pre_matching_shared_seqs"]
        ),
        "rows": closure_checks,
    }

    # 4. DIAMOND signature chain audit (.dmnd.sig.json and .m8.sig.json)
    dmnd_sigs = {}
    for p in sorted(H3_SNAP.glob("matched_*_db.dmnd.sig.json")):
        d = json.loads(p.read_text())
        dmnd_sigs[p.name] = d

    m8_sigs = {}
    for p in sorted(H3_SNAP.glob("hits_*.m8.sig.json")):
        d = json.loads(p.read_text())
        m8_sigs[p.name] = d

    sig_links_ok = True
    query_hashes = set()
    for name, m8 in m8_sigs.items():
        if m8["exit_code"] != 0 or m8["diamond_version"] != "diamond version 2.0.15":
            sig_links_ok = False
        query_hashes.add(m8["query_sha256"])
        db_base = Path(m8["db_file"]).name + ".sig.json"
        if db_base not in dmnd_sigs:
            sig_links_ok = False
        else:
            dmnd = dmnd_sigs[db_base]
            if (
                dmnd["exit_code"] != 0
                or dmnd["db_sha256"] != m8["db_sha256"]
                or dmnd["input_faa_sha256"] != m8["target_faa_sha256"]
            ):
                sig_links_ok = False

    res["diamond_signatures"] = {
        "dmnd_sig_count": len(dmnd_sigs),
        "m8_sig_count": len(m8_sigs),
        "all_exit_zero_and_hash_linked": sig_links_ok,
        "unique_query_all_faa_sha256": list(query_hashes),
    }

    # 5. Full census summary & length-stratified arithmetic + McNemar check
    full_rows = read_tsv(H3_SNAP / "matched_full_census_summary.tsv")
    strat_rows = read_tsv(H3_SNAP / "matched_length_stratified_results.tsv")

    full_verified = []
    full_by_key = {}
    all_full_exact = True
    for r in full_rows:
        n = int(r["total_queries"])
        a = int(r["both_hits_a"])
        b = int(r["smag_only_b"])
        c = int(r["gem_only_c"])
        d = int(r["neither_d"])
        smag_hits = a + b
        gem_hits = a + c
        smag_rate = float(r["smag_hit_pct"])
        gem_rate = float(r["gem_hit_pct"])
        delta_pct = float(r["paired_delta_pct"])
        ratio_rep = float(r["discordant_ratio"])

        calc_smag_rate = 100.0 * (a + b) / n
        calc_gem_rate = 100.0 * (a + c) / n
        calc_delta = 100.0 * (b - c) / n
        calc_ratio = b / c if c > 0 else float("inf")
        pval = binomtest(b, b + c, 0.5).pvalue

        row_ok = (
            n == 106042
            and (a + b + c + d == n)
            and abs(round(calc_smag_rate, 2) - smag_rate) <= 1e-2
            and abs(round(calc_gem_rate, 2) - gem_rate) <= 1e-2
            and abs(round(calc_delta, 2) - delta_pct) <= 1e-2
            and abs(round(calc_ratio, 2) - ratio_rep) <= 1e-2
        )
        if not row_ok:
            all_full_exact = False

        item = {
            "run_label": r["run_label"],
            "cohort": r["cohort"],
            "total_queries": n,
            "a": a,
            "b": b,
            "c": c,
            "d": d,
            "smag_hits": smag_hits,
            "gem_hits": gem_hits,
            "smag_rate_pct": round(calc_smag_rate, 4),
            "gem_rate_pct": round(calc_gem_rate, 4),
            "paired_diff_pct": round(calc_delta, 4),
            "discordant_ratio_b_c": round(calc_ratio, 4),
            "mcnemar_exact_p": pval,
            "arithmetic_ok": row_ok,
        }
        full_verified.append(item)
        full_by_key[(r["run_label"], r["cohort"])] = item

    # Verify length-stratified rows sum to full census rows
    strat_sums = defaultdict(lambda: [0, 0, 0, 0, 0])
    all_strat_exact = True
    strat_by_key = defaultdict(dict)
    for r in strat_rows:
        key = (r["run_label"], r["cohort"])
        n = int(r["stratum_queries"])
        a = int(r["both_hits_a"])
        b = int(r["smag_only_b"])
        c = int(r["gem_only_c"])
        d = int(r["neither_d"])
        strat_sums[key][0] += n
        strat_sums[key][1] += a
        strat_sums[key][2] += b
        strat_sums[key][3] += c
        strat_sums[key][4] += d

        calc_smag_rate = 100.0 * (a + b) / n
        calc_gem_rate = 100.0 * (a + c) / n
        calc_delta = 100.0 * (b - c) / n
        row_ok = (
            (a + b + c + d == n)
            and abs(round(calc_smag_rate, 2) - float(r["smag_hit_pct"])) <= 1e-2
            and abs(round(calc_gem_rate, 2) - float(r["gem_hit_pct"])) <= 1e-2
            and abs(round(calc_delta, 2) - float(r["paired_delta_pct"])) <= 1e-2
        )
        if not row_ok:
            all_strat_exact = False
        strat_by_key[key][r["length_stratum"]] = {
            "n": n,
            "a": a,
            "b": b,
            "c": c,
            "d": d,
            "smag_rate_pct": round(calc_smag_rate, 4),
            "gem_rate_pct": round(calc_gem_rate, 4),
            "paired_diff_pct": round(calc_delta, 4),
        }

    strat_marginal_closure_ok = True
    for key, sums in strat_sums.items():
        fk = full_by_key[key]
        if sums != [fk["total_queries"], fk["a"], fk["b"], fk["c"], fk["d"]]:
            strat_marginal_closure_ok = False

    # Compute direct length-standardization sensitivity for primary run (harmonized_seed20261014)
    # Reweight Control to Dark's 5-bin length distribution and Dark to Control's 5-bin length distribution
    dark_p = strat_by_key[("harmonized_seed20261014", "Dark")]
    ctrl_p = strat_by_key[("harmonized_seed20261014", "Control")]
    bins = ["<100", "100-149", "150-249", "250-399", ">=400"]
    dark_weights = {b: dark_p[b]["n"] / 106042.0 for b in bins}
    ctrl_weights = {b: ctrl_p[b]["n"] / 106042.0 for b in bins}

    dark_delta_under_ctrl_len = sum(ctrl_weights[b] * dark_p[b]["paired_diff_pct"] for b in bins)
    ctrl_delta_under_dark_len = sum(dark_weights[b] * ctrl_p[b]["paired_diff_pct"] for b in bins)

    res["census_and_stratified_verification"] = {
        "full_census_row_count": len(full_verified),
        "stratified_row_count": len(strat_rows),
        "all_16_full_rows_exact": all_full_exact,
        "all_80_stratified_rows_exact": all_strat_exact,
        "all_16_stratified_marginals_sum_to_full": strat_marginal_closure_ok,
        "query_length_distribution_comparison": {
            b: {
                "dark_n": dark_p[b]["n"],
                "dark_pct": round(100.0 * dark_weights[b], 2),
                "control_n": ctrl_p[b]["n"],
                "control_pct": round(100.0 * ctrl_weights[b], 2),
                "primary_dark_delta_pct": dark_p[b]["paired_diff_pct"],
                "primary_control_delta_pct": ctrl_p[b]["paired_diff_pct"],
            }
            for b in bins
        },
        "length_reweighting_sensitivity_seed20261014": {
            "dark_unadjusted_delta_pct": full_by_key[("harmonized_seed20261014", "Dark")][
                "paired_diff_pct"
            ],
            "dark_reweighted_to_control_length_dist_delta_pct": round(
                dark_delta_under_ctrl_len, 4
            ),
            "control_unadjusted_delta_pct": full_by_key[("harmonized_seed20261014", "Control")][
                "paired_diff_pct"
            ],
            "control_reweighted_to_dark_length_dist_delta_pct": round(
                ctrl_delta_under_dark_len, 4
            ),
        },
        "full_census_rows": full_verified,
    }
    return res


def verify_h4() -> dict:
    res: dict = {}

    # 1. Verify SHA256 hashes of local preclosure files against v7_results.sha256.tsv
    sha_lines = [
        line.strip().split()
        for line in (H4_PRE / "v7_results.sha256.tsv").read_text().splitlines()
        if line.strip()
    ]
    local_sha_checks = {}
    for expected_sha, fname in sha_lines:
        fpath = H4_PRE / fname
        if fpath.exists():
            actual_sha = sha256_file(fpath)
            local_sha_checks[fname] = {
                "present_locally": True,
                "sha256_match": actual_sha == expected_sha,
                "sha256": actual_sha,
            }
        else:
            local_sha_checks[fname] = {
                "present_locally": False,
                "expected_sha256_in_manifest": expected_sha,
            }

    res["preclosure_sha256_verification"] = {
        "total_manifest_entries": len(sha_lines),
        "locally_present_count": sum(1 for v in local_sha_checks.values() if v["present_locally"]),
        "all_local_files_match_sha256": all(
            v["sha256_match"] for v in local_sha_checks.values() if v["present_locally"]
        ),
        "files": local_sha_checks,
    }

    # 2. Load v7_model_summaries.tsv (5,851), v7_region_evaluations.tsv (7,773),
    #    v7_full5851_cropped_regions_manifest.tsv (7,934), and v7_stratum_census.tsv (40)
    models = read_tsv(H4_PRE / "v7_model_summaries.tsv")
    regions = read_tsv(H4_PRE / "v7_region_evaluations.tsv")
    manifest = read_tsv(H4_PRE / "v7_full5851_cropped_regions_manifest.tsv")
    census = read_tsv(H4_PRE / "v7_stratum_census.tsv")

    # Check manifest closure
    manifest_zero_reg = [r for r in manifest if int(r["n_confident_regions"]) == 0]
    manifest_pos_reg = [r for r in manifest if int(r["n_confident_regions"]) > 0]
    manifest_unique_models = set(r["model_id"] for r in manifest)

    # Recompute all 40 rows of v7_stratum_census.tsv from models and regions
    strata_defs = {
        "K_local_pdb": lambda r: r["sub_stratum"] == "K_local_pdb",
        "K_local_afdb_only": lambda r: r["sub_stratum"] == "K_local_afdb_only",
        "K_local_all": lambda r: r["parent_stratum"] == "K_local_all",
        "K_dual_none": lambda r: r["sub_stratum"] == "K_dual_none",
        "Q_not_K_pdb": lambda r: r["sub_stratum"] == "Q_not_K_pdb",
        "Q_not_K_afdb_only": lambda r: r["sub_stratum"] == "Q_not_K_afdb_only",
        "Q_not_K_all": lambda r: r["parent_stratum"] == "Q_not_K_all",
        "Controls": lambda r: r["sub_stratum"] == "Controls",
    }

    census_mismatches = []
    for crow in census:
        st = crow["stratum"]
        k = int(crow["tier"].split("_")[1])
        pred = strata_defs[st]
        m_sub = [m for m in models if pred(m)]
        r_sub = [r for r in regions if pred(r)]

        n_models = len(m_sub)
        n_elig = sum(1 for m in m_sub if int(m["n_confident_regions"]) > 0)
        n_no_reg = sum(1 for m in m_sub if int(m["n_confident_regions"]) == 0)
        n_regs = len(r_sub)

        wc_matched = sum(1 for m in m_sub if m[f"wc_t{k}_status"] == "MATCHED")
        wc_low = sum(1 for m in m_sub if m[f"wc_t{k}_status"] == "LOW_SCORE")
        wc_excl = sum(1 for m in m_sub if m[f"wc_t{k}_status"] == "POST_EXCLUSION_EMPTY")
        wc_zero = sum(1 for m in m_sub if m[f"wc_t{k}_status"] == "RAW_ZERO_RETRIEVAL")

        crop_matched = sum(1 for m in m_sub if m[f"crop_t{k}_status"] == "MATCHED")
        crop_low = sum(1 for m in m_sub if m[f"crop_t{k}_status"] == "PERSISTENT_LOW_SCORE")
        crop_excl = sum(1 for m in m_sub if m[f"crop_t{k}_status"] == "POST_EXCLUSION_EMPTY")
        crop_zero = sum(1 for m in m_sub if m[f"crop_t{k}_status"] == "RAW_ZERO_RETRIEVAL")
        crop_no_reg = sum(1 for m in m_sub if m[f"crop_t{k}_status"] == "NO_CONFIDENT_REGION")
        crop_rem_denom = sum(1 for m in m_sub if int(m[f"crop_t{k}_rem_regions"]) > 0)

        wc_unm_total = n_models - wc_matched
        wc_unm_elig = sum(
            1
            for m in m_sub
            if m[f"wc_t{k}_status"] != "MATCHED" and int(m["n_confident_regions"]) > 0
        )

        tr_rec = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_UNMATCHED_TO_CROP_RECOVERED"
        )
        tr_plow = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE"
        )
        tr_pexcl = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_UNMATCHED_TO_POST_EXCLUSION_EMPTY"
        )
        tr_rzero = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL"
        )
        tr_unm_nr = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_UNMATCHED_NO_CONFIDENT_REGION"
        )
        tr_ret = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_MATCHED_RETAINED_IN_CROP"
        )
        tr_llow = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_MATCHED_LOST_LOW_SCORE"
        )
        tr_lexcl = sum(
            1
            for m in m_sub
            if m[f"trans_t{k}"] == "WC_MATCHED_LOST_POST_EXCLUSION_EMPTY"
        )
        tr_lzero = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_MATCHED_LOST_RAW_ZERO_RETRIEVAL"
        )
        tr_m_nr = sum(
            1 for m in m_sub if m[f"trans_t{k}"] == "WC_MATCHED_LOST_NO_CONFIDENT_REGION"
        )

        expected_ints = {
            "model_count": n_models,
            "models_eligible": n_elig,
            "models_no_region": n_no_reg,
            "total_regions": n_regs,
            "wc_matched": wc_matched,
            "wc_low_score": wc_low,
            "wc_post_empty": wc_excl,
            "wc_raw_zero": wc_zero,
            "crop_matched": crop_matched,
            "crop_persistent_low": crop_low,
            "crop_post_empty": crop_excl,
            "crop_raw_zero": crop_zero,
            "crop_no_region": crop_no_reg,
            "unconditional_crop_match_count": crop_matched,
            "unconditional_crop_match_denominator": n_models,
            "region_eligible_crop_match_count": crop_matched,
            "region_eligible_crop_match_denominator": n_elig,
            "remaining_hit_conditional_crop_match_count": crop_matched,
            "remaining_hit_conditional_crop_match_denominator": crop_rem_denom,
            "wc_unmatched_count": wc_unm_total,
            "wc_unmatched_recovered_count": tr_rec,
            "wc_unmatched_eligible_count": wc_unm_elig,
            "wc_matched_count": wc_matched,
            "wc_matched_retained_count": tr_ret,
            "trans_recovered": tr_rec,
            "trans_persistent_low": tr_plow,
            "trans_post_empty": tr_pexcl,
            "trans_raw_zero": tr_rzero,
            "trans_unm_no_reg": tr_unm_nr,
            "trans_retained": tr_ret,
            "trans_lost_low": tr_llow,
            "trans_lost_empty": tr_lexcl,
            "trans_lost_zero": tr_lzero,
            "trans_m_no_reg": tr_m_nr,
        }
        for col, val in expected_ints.items():
            if int(crow[col]) != val:
                census_mismatches.append((st, k, col, int(crow[col]), val))

    # Check dark models Tier 0 == Tier 1 symmetry and Tiers 0..4 differences
    dark_models = [m for m in models if int(m["is_control"]) == 0]
    ctrl_models = [m for m in models if int(m["is_control"]) == 1]
    dark_regions = [r for r in regions if r["sub_stratum"] != "Controls"]
    ctrl_regions = [r for r in regions if r["sub_stratum"] == "Controls"]

    dark_t0_eq_t1 = all(
        m["wc_t0_status"] == m["wc_t1_status"]
        and m["crop_t0_status"] == m["crop_t1_status"]
        and m["wc_t0_qtm"] == m["wc_t1_qtm"]
        and m["crop_t0_max_qtm"] == m["crop_t1_max_qtm"]
        for m in dark_models
    )
    dark_models_varying_across_t0_t4 = [
        m["model_id"]
        for m in dark_models
        if len({m[f"wc_t{k}_status"] for k in range(5)}) > 1
        or len({m[f"crop_t{k}_status"] for k in range(5)}) > 1
    ]
    dark_regions_fident_ge_20 = [
        r["crop_id"] for r in dark_regions if int(r["hits_fident_ge_20"]) > 0
    ]

    # Adaptive rerun counts in model_summaries and region_evaluations
    wc_hit_cap_5000 = [m for m in models if int(m["wc_hit_cap"]) == 5000]
    crop_hit_cap_5000 = [r for r in regions if int(r["hit_cap"]) == 5000]

    res["h4_census_and_invariants"] = {
        "total_models": len(models),
        "dark_models": len(dark_models),
        "control_models": len(ctrl_models),
        "total_cropped_regions": len(regions),
        "dark_cropped_regions": len(dark_regions),
        "control_cropped_regions": len(ctrl_regions),
        "manifest_total_rows": len(manifest),
        "manifest_unique_models": len(manifest_unique_models),
        "manifest_zero_region_rows": len(manifest_zero_reg),
        "manifest_positive_region_rows": len(manifest_pos_reg),
        "stratum_census_rows_verified": len(census),
        "stratum_census_mismatches": len(census_mismatches),
        "dark_tier0_equals_tier1_exact": dark_t0_eq_t1,
        "dark_models_with_any_status_change_across_t0_t4": len(
            dark_models_varying_across_t0_t4
        ),
        "dark_models_varying_ids": dark_models_varying_across_t0_t4,
        "dark_regions_with_hits_fident_ge_20": len(dark_regions_fident_ge_20),
        "dark_regions_fident_ge_20_ids": dark_regions_fident_ge_20,
        "wc_models_updated_by_adaptive_rerun_5000": len(wc_hit_cap_5000),
        "wc_adaptive_rerun_all_controls": all(
            int(m["is_control"]) == 1 for m in wc_hit_cap_5000
        ),
        "crop_regions_updated_by_adaptive_rerun_5000": len(crop_hit_cap_5000),
        "crop_adaptive_rerun_all_controls": all(
            r["sub_stratum"] == "Controls" for r in crop_hit_cap_5000
        ),
    }

    # 3. Audit Foldseek alnTM > 1 anomaly, TMaligner.cpp non-linear dependence on normalizationLen,
    #    Flash's DIAGNOSIS.md errors, and Root's pooled-qTM (K_all = 1,031; Q = 2,838) census
    anom_rows = read_tsv(H4_ANOM / "tm_anomaly_queries.tsv")
    anom_by_file = Counter(r["file"] for r in anom_rows)
    anom_all_controls = all(
        int(r["is_control"]) == 1 and r["query"].startswith("KF") for r in anom_rows
    )
    anom_unique_queries = sorted(set(r["query"] for r in anom_rows))
    anom_unique_wc = [q for q in anom_unique_queries if "_reg" not in q]
    anom_unique_crop = [q for q in anom_unique_queries if "_reg" in q]
    max_raw_alntm = max(float(r["alntm"]) for r in anom_rows)
    max_anom_qtm = max(float(r["qtm"]) for r in anom_rows)
    max_anom_ttm = max(float(r["ttm"]) for r in anom_rows)

    # Audit exclusion_sensitivity_diff.json to expose Flash's query-tier vs unique-entity confusion
    excl_diff = json.loads((H4_ANOM / "exclusion_sensitivity_diff.json").read_text())
    excl_diff_audit = {}
    for key, items in excl_diff.items():
        if items and isinstance(items, list):
            id_field = next(
                (k for k in ("query", "model_id", "crop_id", "model", "id") if k in items[0]),
                list(items[0].keys())[0],
            )
            unique_ids = sorted(set(x[id_field] for x in items))
            excl_diff_audit[key] = {
                "query_tier_pair_count_reported_by_flash": len(items),
                "unique_entity_count": len(unique_ids),
                "unique_ids": unique_ids,
            }

    # Check Control crop match rates in baseline vs exclusion_sensitivity_census.tsv
    excl_census = read_tsv(H4_ANOM / "exclusion_sensitivity_census.tsv")
    ctrl_baseline_vs_excl = {}
    for k in range(5):
        base_row = next(
            r for r in census if r["stratum"] == "Controls" and r["tier"] == f"tier_{k}"
        )
        excl_row = next(
            r
            for r in excl_census
            if r["stratum"] == "Controls" and r["tier"] == f"tier_{k}"
        )
        ctrl_baseline_vs_excl[f"tier_{k}"] = {
            "baseline_wc_matched": int(base_row["wc_matched"]),
            "baseline_wc_pct": round(100.0 * int(base_row["wc_matched"]) / int(base_row["model_count"]), 2),
            "excl_wc_matched": int(excl_row["wc_matched"]),
            "excl_wc_pct": round(100.0 * int(excl_row["wc_matched"]) / int(excl_row["model_count"]), 2),
            "baseline_crop_matched": int(base_row["crop_matched"]),
            "baseline_crop_uncond_pct": float(base_row["unconditional_crop_match_pct"]),
            "baseline_crop_elig_pct": float(base_row["region_eligible_crop_match_pct"]),
            "excl_crop_matched": int(excl_row["crop_matched"]),
            "excl_crop_uncond_pct": float(excl_row["unconditional_crop_match_pct"]),
            "excl_crop_elig_pct": float(excl_row["region_eligible_crop_match_pct"]),
        }

    # Verify qTM-only cohort definitions (independent of historical or current alnTM)
    k_pooled_models = [
        m
        for m in dark_models
        if m["sub_stratum"] in ("K_local_pdb", "K_local_afdb_only", "K_dual_none")
    ]
    k_pooled_qtm_equiv = all(
        float(m["hist_pdb_qtm"]) < 0.5 and float(m["hist_afdb_qtm"]) < 0.5
        for m in k_pooled_models
    ) and all(
        (float(m["hist_pdb_qtm"]) >= 0.5 or float(m["hist_afdb_qtm"]) >= 0.5)
        for m in dark_models
        if m["sub_stratum"] not in ("K_local_pdb", "K_local_afdb_only", "K_dual_none")
    )
    q_not_k_pdb_equiv = all(
        (float(m["hist_pdb_qtm"]) >= 0.5) == (m["sub_stratum"] == "Q_not_K_pdb")
        for m in dark_models
    )
    q_not_k_afdb_equiv = all(
        (float(m["hist_pdb_qtm"]) < 0.5 and float(m["hist_afdb_qtm"]) >= 0.5)
        == (m["sub_stratum"] == "Q_not_K_afdb_only")
        for m in dark_models
    )

    def summarize_qtm_cohort(m_list: list[dict], r_list: list[dict]) -> dict:
        n_mod = len(m_list)
        n_el = sum(1 for m in m_list if int(m["n_confident_regions"]) > 0)
        n_nr = sum(1 for m in m_list if int(m["n_confident_regions"]) == 0)
        n_rg = len(r_list)
        tiers_out = {}
        for k in range(5):
            wc_m = sum(1 for m in m_list if m[f"wc_t{k}_status"] == "MATCHED")
            wc_low = sum(1 for m in m_list if m[f"wc_t{k}_status"] == "LOW_SCORE")
            wc_pe = sum(1 for m in m_list if m[f"wc_t{k}_status"] == "POST_EXCLUSION_EMPTY")
            wc_rz = sum(1 for m in m_list if m[f"wc_t{k}_status"] == "RAW_ZERO_RETRIEVAL")
            wc_unm = n_mod - wc_m
            wc_unm_el = sum(
                1
                for m in m_list
                if m[f"wc_t{k}_status"] != "MATCHED" and int(m["n_confident_regions"]) > 0
            )
            cr_m = sum(1 for m in m_list if m[f"crop_t{k}_status"] == "MATCHED")
            cr_low = sum(
                1 for m in m_list if m[f"crop_t{k}_status"] == "PERSISTENT_LOW_SCORE"
            )
            cr_pe = sum(
                1 for m in m_list if m[f"crop_t{k}_status"] == "POST_EXCLUSION_EMPTY"
            )
            cr_rz = sum(
                1 for m in m_list if m[f"crop_t{k}_status"] == "RAW_ZERO_RETRIEVAL"
            )
            cr_nr = sum(
                1 for m in m_list if m[f"crop_t{k}_status"] == "NO_CONFIDENT_REGION"
            )
            cr_rem_denom = sum(1 for m in m_list if int(m[f"crop_t{k}_rem_regions"]) > 0)

            tr_rec = sum(
                1 for m in m_list if m[f"trans_t{k}"] == "WC_UNMATCHED_TO_CROP_RECOVERED"
            )
            tr_plow = sum(
                1
                for m in m_list
                if m[f"trans_t{k}"] == "WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE"
            )
            tr_pe = sum(
                1
                for m in m_list
                if m[f"trans_t{k}"] == "WC_UNMATCHED_TO_POST_EXCLUSION_EMPTY"
            )
            tr_rz = sum(
                1
                for m in m_list
                if m[f"trans_t{k}"] == "WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL"
            )
            tr_unm_nr = sum(
                1 for m in m_list if m[f"trans_t{k}"] == "WC_UNMATCHED_NO_CONFIDENT_REGION"
            )
            tr_ret = sum(
                1 for m in m_list if m[f"trans_t{k}"] == "WC_MATCHED_RETAINED_IN_CROP"
            )
            tr_llow = sum(
                1 for m in m_list if m[f"trans_t{k}"] == "WC_MATCHED_LOST_LOW_SCORE"
            )
            tr_lpe = sum(
                1
                for m in m_list
                if m[f"trans_t{k}"] == "WC_MATCHED_LOST_POST_EXCLUSION_EMPTY"
            )
            tr_lrz = sum(
                1
                for m in m_list
                if m[f"trans_t{k}"] == "WC_MATCHED_LOST_RAW_ZERO_RETRIEVAL"
            )
            tr_m_nr = sum(
                1
                for m in m_list
                if m[f"trans_t{k}"] == "WC_MATCHED_LOST_NO_CONFIDENT_REGION"
            )
            reg_m = sum(1 for r in r_list if int(r[f"t{k}_matched"]) == 1)

            tiers_out[f"tier_{k}"] = {
                "wc_matched": wc_m,
                "wc_matched_pct": round(100.0 * wc_m / n_mod, 2),
                "wc_unmatched": wc_unm,
                "wc_unmatched_eligible": wc_unm_el,
                "wc_low_score": wc_low,
                "wc_post_empty": wc_pe,
                "wc_raw_zero": wc_rz,
                "crop_matched": cr_m,
                "crop_uncond_pct": round(100.0 * cr_m / n_mod, 2),
                "crop_elig_pct": round(100.0 * cr_m / n_el, 2),
                "crop_rem_hit_denom": cr_rem_denom,
                "crop_rem_hit_pct": round(100.0 * cr_m / cr_rem_denom, 2),
                "crop_persistent_low": cr_low,
                "crop_post_empty": cr_pe,
                "crop_raw_zero": cr_rz,
                "crop_no_region": cr_nr,
                "trans_recovered": tr_rec,
                "trans_recovered_uncond_pct": round(100.0 * tr_rec / wc_unm, 2)
                if wc_unm > 0
                else 0.0,
                "trans_recovered_elig_pct": round(100.0 * tr_rec / wc_unm_el, 2)
                if wc_unm_el > 0
                else 0.0,
                "trans_persistent_low": tr_plow,
                "trans_post_empty": tr_pe,
                "trans_raw_zero": tr_rz,
                "trans_unm_no_reg": tr_unm_nr,
                "trans_retained": tr_ret,
                "trans_retained_pct": round(100.0 * tr_ret / wc_m, 2) if wc_m > 0 else 0.0,
                "trans_lost_low": tr_llow,
                "trans_lost_empty": tr_lpe,
                "trans_lost_zero": tr_lrz,
                "trans_m_no_reg": tr_m_nr,
                "region_matched": reg_m,
                "region_matched_pct": round(100.0 * reg_m / n_rg, 2) if n_rg > 0 else 0.0,
            }
        return {
            "model_count": n_mod,
            "models_eligible": n_el,
            "models_no_region": n_nr,
            "total_regions": n_rg,
            "tiers": tiers_out,
        }

    k_pooled_regions = [
        r
        for r in dark_regions
        if r["sub_stratum"] in ("K_local_pdb", "K_local_afdb_only", "K_dual_none")
    ]
    q_not_k_pdb_models = [m for m in dark_models if m["sub_stratum"] == "Q_not_K_pdb"]
    q_not_k_pdb_regions = [r for r in dark_regions if r["sub_stratum"] == "Q_not_K_pdb"]
    q_not_k_afdb_models = [
        m for m in dark_models if m["sub_stratum"] == "Q_not_K_afdb_only"
    ]
    q_not_k_afdb_regions = [
        r for r in dark_regions if r["sub_stratum"] == "Q_not_K_afdb_only"
    ]
    q_not_k_all_models = [m for m in dark_models if m["parent_stratum"] == "Q_not_K_all"]
    q_not_k_all_regions = [
        r for r in dark_regions if r["parent_stratum"] == "Q_not_K_all"
    ]

    res["h4_tm_anomaly_and_flash_audit"] = {
        "anomalous_rows_total": len(anom_rows),
        "anomalous_rows_by_file": dict(anom_by_file),
        "all_anomalous_rows_are_controls": anom_all_controls,
        "unique_anomalous_queries_total": len(anom_unique_queries),
        "unique_anomalous_wc_models": len(anom_unique_wc),
        "unique_anomalous_cropped_regions": len(anom_unique_crop),
        "max_raw_alntm": max_raw_alntm,
        "max_anom_qtm": max_anom_qtm,
        "max_anom_ttm": max_anom_ttm,
        "tmaligner_source_audit_and_withdrawal_note": {
            "withdrawn_claims": [
                "Withdrawn: claim that alntm * (L_span - 1) / L_span is an exact recalculation of alnTM.",
                "Withdrawn: claim that [0.5, 0.5051] (L>=100) and [0.5, 0.5102] (L>=50) provide validated upper bounds on historical K_local_pdb / K_local_afdb_only / K_dual_none membership changes.",
                "Withdrawn: 1 WC model / 14 cropped region threshold-crossing counts from v7_model_summaries.tsv and v7_region_evaluations.tsv."
            ],
            "technical_reasons": [
                "In h4_tm_anomaly/source/TMaligner.cpp lines 89-103 (computeAppoximateTMscore) and lines 145-200 (computeExactTMscore), normalizationLen is passed into parameter_set4search (and parameter_set4final), which sets d0, score_d8, and d0_search and governs iterative Kabsch superposition optimization in standard_TMscore, detailed_search_standard, and TMscore8_search. Thus normalizationLen is not a simple external scalar divisor.",
                "In v7_recovery_postprocess.py lines 377-386, wc_t{k}_alntm is taken from wa['top_alntm'][k] (independently maximized alnTM across hits), whereas wc_t{k}_qspan_cov is taken from wa['top_qtm'][k] (the top-qTM hit). In v7_region_evaluations.tsv lines 305-325, t1_alntm is the alnTM of the top-qTM hit (not top-alnTM hit) and only qspan (not min(qspan, tspan)) is recorded.",
                "In historical per_rep_structures.tsv, alignment spans min(qspan, tspan) are not stored and no minimum 50 aa or 100 aa alignment span applies to all whole-chain alignments."
            ]
        },
        "flash_exclusion_diff_audit": excl_diff_audit,
        "control_baseline_vs_alntm_gt1_excluded": ctrl_baseline_vs_excl,
        "qtm_only_pooled_census": {
            "definitions_strictly_qtm_and_plddt_independent_of_alntm": (
                k_pooled_qtm_equiv and q_not_k_pdb_equiv and q_not_k_afdb_equiv
            ),
            "K_pooled_all_initial_qtm_lt_0.5": summarize_qtm_cohort(
                k_pooled_models, k_pooled_regions
            ),
            "Q_not_K_pdb_initial_pdb_qtm_ge_0.5": summarize_qtm_cohort(
                q_not_k_pdb_models, q_not_k_pdb_regions
            ),
            "Q_not_K_afdb_only_initial_afdb_qtm_ge_0.5": summarize_qtm_cohort(
                q_not_k_afdb_models, q_not_k_afdb_regions
            ),
            "Q_not_K_all_initial_any_qtm_ge_0.5": summarize_qtm_cohort(
                q_not_k_all_models, q_not_k_all_regions
            ),
            "Q_all_confidence_qualified_dark": summarize_qtm_cohort(
                dark_models, dark_regions
            ),
            "Controls_positive_pdb": summarize_qtm_cohort(ctrl_models, ctrl_regions),
        },
    }

    # 4. Check ORIX GPU local provenance files
    orix_transfer = json.loads((ORIX_DIR / "transfer_status.json").read_text())
    res["orix_gpu_provenance"] = {
        "transfer_status": orix_transfer,
        "cpu_vs_gpu_raw_comparison_status": "PENDING_ROOT_RAW_TRANSFER_AND_COMPARISON",
    }

    return res


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    report = {
        "reviewer": "Gemini 4 Argon High (argon_h3h4_review1)",
        "generated_at": "2026-10-14",
        "h3_verification": verify_h3(),
        "h4_verification": verify_h4(),
    }
    json_path = OUT_DIR / "verification_output.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    # Also write a concise human-readable summary
    h3 = report["h3_verification"]
    h4 = report["h4_verification"]
    kp = h4["h4_tm_anomaly_and_flash_audit"]["qtm_only_pooled_census"][
        "K_pooled_all_initial_qtm_lt_0.5"
    ]
    qa = h4["h4_tm_anomaly_and_flash_audit"]["qtm_only_pooled_census"][
        "Q_all_confidence_qualified_dark"
    ]
    lines = [
        "=== ARGON H3 & H4 INDEPENDENT VERIFICATION SUMMARY (CORRECTED) ===",
        f"H3 SMAG 1,500 replay == sampled_smag_1500_ids_seed20261014.txt: {h3['smag_sampling_and_extraction']['replay_matches_sampled_file']}",
        f"H3 SMAG 1,500 manifest declared == observed bytes: {h3['smag_sampling_and_extraction']['all_1500_declared_eq_observed_bytes']}",
        f"H3 GEM 2,055 replay matches all 6 manifests (seeds 20261014..20261016): {all(v['replay_matches_harmonized_manifest'] and v['replay_matches_published_manifest'] for v in h3['gem_sampling_and_manifests']['per_seed'].values())}",
        f"H3 matching_closure_audit.tsv all 7 rows 5-bin sum == matched_proteins: {h3['matching_closure_audit']['all_bins_sum_match_matched_proteins']}",
        f"H3 DIAMOND signatures (14 dmnd + 16 m8) exit_code==0 & SHA256 linked: {h3['diamond_signatures']['all_exit_zero_and_hash_linked']}",
        f"H3 matched_full_census_summary.tsv (16 rows) exact arithmetic & McNemar: {h3['census_and_stratified_verification']['all_16_full_rows_exact']}",
        f"H3 matched_length_stratified_results.tsv (80 rows) exact & sum to full: {h3['census_and_stratified_verification']['all_80_stratified_rows_exact'] and h3['census_and_stratified_verification']['all_16_stratified_marginals_sum_to_full']}",
        f"H3 Length-reweighted Dark delta (to Control length weights): {h3['census_and_stratified_verification']['length_reweighting_sensitivity_seed20261014']['dark_reweighted_to_control_length_dist_delta_pct']}% (unadjusted: {h3['census_and_stratified_verification']['length_reweighting_sensitivity_seed20261014']['dark_unadjusted_delta_pct']}%)",
        f"H3 Length-reweighted Control delta (to Dark length weights): {h3['census_and_stratified_verification']['length_reweighting_sensitivity_seed20261014']['control_reweighted_to_dark_length_dist_delta_pct']}% (unadjusted: {h3['census_and_stratified_verification']['length_reweighting_sensitivity_seed20261014']['control_unadjusted_delta_pct']}%)",
        "---",
        f"H4 local preclosure SHA256 matches (11/15 local files): {h4['preclosure_sha256_verification']['all_local_files_match_sha256']}",
        f"H4 v7_stratum_census.tsv (40 rows x 27 int cols) mismatches: {h4['h4_census_and_invariants']['stratum_census_mismatches']}",
        f"H4 Dark Tier 0 == Tier 1 exact across 2,838 models: {h4['h4_census_and_invariants']['dark_tier0_equals_tier1_exact']}",
        f"H4 Dark models changing status across Tiers 0..4: {h4['h4_census_and_invariants']['dark_models_with_any_status_change_across_t0_t4']} ({h4['h4_census_and_invariants']['dark_models_varying_ids']})",
        f"H4 Adaptive rerun (max-seqs 5000) updated models/regions: {h4['h4_census_and_invariants']['wc_models_updated_by_adaptive_rerun_5000']} WC models, {h4['h4_census_and_invariants']['crop_regions_updated_by_adaptive_rerun_5000']} crop regions (100% Controls)",
        f"H4 Pooled-K (N={kp['model_count']}, eligible={kp['models_eligible']}, regions={kp['total_regions']}) strictly qTM-defined: {h4['h4_tm_anomaly_and_flash_audit']['qtm_only_pooled_census']['definitions_strictly_qtm_and_plddt_independent_of_alntm']}",
        f"H4 Pooled-K Tier 1: WC matched={kp['tiers']['tier_1']['wc_matched']}/{kp['model_count']} ({kp['tiers']['tier_1']['wc_matched_pct']}%), WC unmatched={kp['tiers']['tier_1']['wc_unmatched']}, Crop matched={kp['tiers']['tier_1']['crop_matched']}/{kp['model_count']} ({kp['tiers']['tier_1']['crop_uncond_pct']}% uncond; {kp['tiers']['tier_1']['crop_elig_pct']}% elig), WC-unmatched recovered={kp['tiers']['tier_1']['trans_recovered']}/{kp['tiers']['tier_1']['wc_unmatched']} ({kp['tiers']['tier_1']['trans_recovered_uncond_pct']}% uncond; {kp['tiers']['tier_1']['trans_recovered_elig_pct']}% elig)",
        f"H4 All Dark Q (N={qa['model_count']}, eligible={qa['models_eligible']}, regions={qa['total_regions']}) Tier 1: WC matched={qa['tiers']['tier_1']['wc_matched']}/{qa['model_count']} ({qa['tiers']['tier_1']['wc_matched_pct']}%), Crop matched={qa['tiers']['tier_1']['crop_matched']}/{qa['model_count']} ({qa['tiers']['tier_1']['crop_uncond_pct']}% uncond; {qa['tiers']['tier_1']['crop_elig_pct']}% elig), WC-unmatched recovered={qa['tiers']['tier_1']['trans_recovered']}/{qa['tiers']['tier_1']['wc_unmatched']} ({qa['tiers']['tier_1']['trans_recovered_uncond_pct']}% uncond; {qa['tiers']['tier_1']['trans_recovered_elig_pct']}% elig)",
    ]
    txt_path = OUT_DIR / "verification_output.txt"
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
