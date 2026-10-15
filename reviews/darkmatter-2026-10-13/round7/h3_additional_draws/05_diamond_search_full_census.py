#!/usr/bin/env python3
"""
scripts/05_diamond_search_full_census.py (Candidate for Additional Harmonized Draws)

DIAMOND Searches & Statistical Census Across ALL 212,084 Queries:
- 106,042 dark families (FAMILIES_FINAL.txt)
- 106,042 annotated controls (CONTROL_FAMILIES.txt)
- Query FASTA: query_all.faa

Executes paired searches:
1. Primary Harmonized: SMAG (Pyrodigal v3.7.1) vs GEM (Harmonized Pyrodigal v3.7.1) [seed 20261014]
2. Overlap Excluded: Primary pair with cross-release shared sequences removed [seed 20261014_no_overlap]
3. Identity Sensitivity: Primary pair searched with --id 30 [seed 20261014_id30]
4. Harmonized Draw Sensitivity 1: SMAG vs GEM (Harmonized Pyrodigal v3.7.1) [seed 20261015]
5. Harmonized Draw Sensitivity 2: SMAG vs GEM (Harmonized Pyrodigal v3.7.1) [seed 20261016]
6. Caller Sensitivity: SMAG (Pyrodigal v3.7.1) vs GEM (Published Prodigal v2.6.3) [seed 20261014]
7. Published Draw Sensitivity 1: SMAG vs GEM (Published Prodigal v2.6.3) [seed 20261015]
8. Published Draw Sensitivity 2: SMAG vs GEM (Published Prodigal v2.6.3) [seed 20261016]

Repair Specifications:
- Requires completion signatures pinning:
    * input FAA SHA256
    * query FASTA SHA256
    * DIAMOND binary version
    * all search parameters (mode, sensitivity, evalue, max_target_seqs, outfmt, id_cutoff)
    * completed output SHA256
- Distinguishes database signature (.dmnd.sig.json) from search signature (.m8.sig.json).
  Legitimate reuse of the primary database by the id30 search is supported via verified DB signature.
- Atomic builds & search executions: writes to temporary files and atomically renames ONLY after exit code 0.
- Fail-closed cache semantics: existing output files lacking valid matching signatures fail closed
  with actionable repair messages; files are never silently overwritten or blindly unlinked.
- Valid zero-hit results are recognized as complete when accompanied by a valid exit-0 signature.
- Diagnostic stdout/stderr is preserved and emitted on subprocess failure.
- Census validation: validates exactly 212,084 unique positive-length IDs, disjoint dark and control arms
  (106,042 each), union exactly matching query FASTA IDs, no duplicate FASTA or arm IDs.
- Parsed hit IDs are verified to be strict subsets of query census, and all hit rows are validated
  against declared search thresholds (evalue <= 1e-3, pident >= id_cutoff).
- Scientific DIAMOND settings and contingency/McNemar estimands are strictly preserved.
"""

import os
import sys
import json
import time
import shutil
import hashlib
import tempfile
import subprocess
from collections import defaultdict
from scipy import stats

RESULTS_DIR = os.environ.get("RESULTS_DIR", "results")
QUERY_FAA = os.environ.get("QUERY_FAA", "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/results/query_all.faa")
DARK_IDS_PATH = os.environ.get("DARK_IDS_PATH", "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results/FAMILIES_FINAL.txt")
CTRL_IDS_PATH = os.environ.get("CTRL_IDS_PATH", "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results/CONTROL_FAMILIES.txt")
DIAMOND_BIN = os.environ.get("DIAMOND_BIN", "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")

LENGTH_BINS = [
    ("<100", lambda l: l < 100),
    ("100-149", lambda l: 100 <= l <= 149),
    ("150-249", lambda l: 150 <= l <= 249),
    ("250-399", lambda l: 250 <= l <= 399),
    (">=400", lambda l: l >= 400),
]


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


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def atomic_write_json(target_path, data):
    dir_name = os.path.dirname(os.path.abspath(target_path))
    os.makedirs(dir_name, exist_ok=True)
    tmp_path = target_path + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"
    with open(tmp_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp_path, target_path)


def get_diamond_version(diamond_bin):
    if not os.path.exists(diamond_bin):
        raise FileNotFoundError(f"DIAMOND binary not found at '{diamond_bin}'.")
    try:
        res = subprocess.run([diamond_bin, "version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        ver = res.stdout.strip()
        if not ver:
            ver = res.stderr.strip()
        return ver
    except Exception as e:
        raise RuntimeError(f"Failed to query DIAMOND version from '{diamond_bin}': {e}")


def validate_and_load_query_census(query_faa_path, dark_ids_path, ctrl_ids_path):
    for p in (dark_ids_path, ctrl_ids_path, query_faa_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"Required query census file missing: '{p}'.")

    # 1. Dark query IDs
    dark_ids = []
    with open(dark_ids_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            val = line.strip()
            if val:
                dark_ids.append(val)
    if len(dark_ids) != 106042:
        raise ValueError(f"Expected exactly 106,042 dark IDs in '{dark_ids_path}', found {len(dark_ids):,}.")
    dark_set = set(dark_ids)
    if len(dark_set) != len(dark_ids):
        raise ValueError(f"Duplicate dark IDs detected in '{dark_ids_path}': {len(dark_ids) - len(dark_set)} duplicate entries.")
    for did in dark_ids:
        if len(did) == 0:
            raise ValueError(f"Encountered empty dark ID in '{dark_ids_path}'.")

    # 2. Control query IDs
    ctrl_ids = []
    with open(ctrl_ids_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            val = line.strip()
            if val:
                ctrl_ids.append(val)
    if len(ctrl_ids) != 106042:
        raise ValueError(f"Expected exactly 106,042 control IDs in '{ctrl_ids_path}', found {len(ctrl_ids):,}.")
    ctrl_set = set(ctrl_ids)
    if len(ctrl_set) != len(ctrl_ids):
        raise ValueError(f"Duplicate control IDs detected in '{ctrl_ids_path}': {len(ctrl_ids) - len(ctrl_set)} duplicate entries.")
    for cid in ctrl_ids:
        if len(cid) == 0:
            raise ValueError(f"Encountered empty control ID in '{ctrl_ids_path}'.")

    # 3. Disjointness & arm union check
    overlap = dark_set.intersection(ctrl_set)
    if overlap:
        raise ValueError(f"Dark and control query ID sets are NOT disjoint! Overlap count: {len(overlap):,}; example: {list(overlap)[:3]}")
    expected_census_set = dark_set.union(ctrl_set)
    if len(expected_census_set) != 212084:
        raise ValueError(f"Expected combined census to have exactly 212,084 unique IDs, got {len(expected_census_set):,}.")

    # 4. Load query FASTA lengths and headers
    query_lens = {}
    cur_id = None
    cur_len = 0
    with open(query_faa_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            if line.startswith(">"):
                if cur_id is not None:
                    if cur_len == 0:
                        raise ValueError(f"Query sequence for '{cur_id}' in '{query_faa_path}' has length 0.")
                    query_lens[cur_id] = cur_len
                cur_id = line[1:].strip().split()[0]
                if not cur_id:
                    raise ValueError(f"Line {line_no} in '{query_faa_path}': Empty FASTA ID.")
                if cur_id in query_lens:
                    raise ValueError(f"Duplicate FASTA ID '{cur_id}' detected in '{query_faa_path}'.")
                cur_len = 0
            else:
                cur_len += len(line.strip().rstrip("*"))
        if cur_id is not None:
            if cur_len == 0:
                raise ValueError(f"Query sequence for '{cur_id}' in '{query_faa_path}' has length 0.")
            query_lens[cur_id] = cur_len

    # 5. Union equality check
    fasta_ids = set(query_lens.keys())
    if len(fasta_ids) != 212084:
        raise ValueError(f"Query FASTA contains {len(fasta_ids):,} sequences; expected exactly 212,084 unique sequences.")
    if fasta_ids != expected_census_set:
        diff_dark = dark_set - fasta_ids
        diff_ctrl = ctrl_set - fasta_ids
        diff_extra = fasta_ids - expected_census_set
        raise ValueError(
            f"Query FASTA IDs do not match census arm IDs. Missing dark: {len(diff_dark)}, "
            f"missing ctrl: {len(diff_ctrl)}, extra in FASTA: {len(diff_extra)}."
        )

    return dark_ids, ctrl_ids, query_lens, expected_census_set


def ensure_diamond_db(faa_path, db_path, threads, diamond_bin):
    db_file = db_path + ".dmnd"
    db_sig_file = db_path + ".dmnd.sig.json"

    if not os.path.exists(faa_path):
        raise FileNotFoundError(f"Input FAA file missing for database creation: '{faa_path}'.")

    faa_sha256 = sha256_file(faa_path)
    diamond_ver = get_diamond_version(diamond_bin)

    if os.path.exists(db_file):
        if not os.path.exists(db_sig_file):
            raise RuntimeError(
                f"[SEARCH CACHE INTEGRITY ERROR] DIAMOND database '{db_file}' exists on disk, "
                f"but completion signature file '{db_sig_file}' is missing! "
                f"Cannot verify database integrity or input provenance. Failing closed without overwriting or deleting. "
                f"Repair action: inspect or remove '{db_file}' before retrying."
            )
        try:
            with open(db_sig_file, "r", encoding="utf-8") as f:
                sig = json.load(f)
        except Exception as e:
            raise RuntimeError(
                f"[SEARCH CACHE INTEGRITY ERROR] Corrupted signature file '{db_sig_file}': {e}. "
                f"Repair action: inspect or remove '{db_file}' and '{db_sig_file}'."
            )

        db_sha256 = sha256_file(db_file)
        if (
            sig.get("input_faa_sha256") != faa_sha256
            or sig.get("diamond_version") != diamond_ver
            or sig.get("db_sha256") != db_sha256
            or sig.get("exit_code") != 0
        ):
            raise RuntimeError(
                f"[SEARCH CACHE INTEGRITY ERROR] Database '{db_file}' signature mismatch:\n"
                f"  Recorded FAA SHA256: {sig.get('input_faa_sha256')} vs Current: {faa_sha256}\n"
                f"  Recorded DIAMOND Ver: {sig.get('diamond_version')} vs Current: {diamond_ver}\n"
                f"  Recorded DB SHA256:  {sig.get('db_sha256')} vs Current: {db_sha256}\n"
                f"  Exit code:           {sig.get('exit_code')}\n"
                f"Failing closed to prevent using stale or corrupt database. "
                f"Repair action: verify provenance or remove '{db_file}' and '{db_sig_file}'."
            )
        print(f"  [CACHE HIT] Verified DIAMOND DB: '{db_file}' (SHA256: {db_sha256[:12]}...).")
        return db_sha256

    print(f"  Building DIAMOND DB: {db_path}...")
    t0 = time.time()
    tmp_prefix = f"{db_path}.tmp_{os.getpid()}_{int(time.time() * 1000)}"
    tmp_db_file = tmp_prefix + ".dmnd"

    cmd_mk = [diamond_bin, "makedb", "--in", faa_path, "-d", tmp_prefix, "-p", str(threads)]
    res = subprocess.run(cmd_mk, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        if os.path.exists(tmp_db_file):
            try:
                os.remove(tmp_db_file)
            except OSError:
                pass
        print(f"DIAMOND makedb failed with code {res.returncode}:\nSTDERR:\n{res.stderr}\nSTDOUT:\n{res.stdout}", file=sys.stderr)
        raise RuntimeError(f"DIAMOND makedb failed for '{db_path}' (code {res.returncode}): {res.stderr.strip()}")

    if not os.path.exists(tmp_db_file):
        raise FileNotFoundError(f"DIAMOND makedb exited 0 but expected output '{tmp_db_file}' was not created!")

    db_sha256 = sha256_file(tmp_db_file)
    os.replace(tmp_db_file, db_file)

    sig_data = {
        "signature_type": "diamond_makedb_v1",
        "db_file": db_file,
        "input_faa_path": os.path.abspath(faa_path),
        "input_faa_sha256": faa_sha256,
        "diamond_version": diamond_ver,
        "makedb_args": {
            "subcommand": "makedb",
            "in": faa_path,
            "db": db_path,
            "threads": threads,
        },
        "db_sha256": db_sha256,
        "exit_code": 0,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    atomic_write_json(db_sig_file, sig_data)
    print(f"  DB built and verified in {time.time() - t0:.1f}s (DB SHA256: {db_sha256[:12]}...).")
    return db_sha256


def run_diamond_search(faa_path, db_path, m8_path, query_faa, query_sha256, threads, diamond_bin, id_cutoff=0):
    db_sha256 = ensure_diamond_db(faa_path, db_path, threads, diamond_bin)
    m8_sig_file = m8_path + ".sig.json"
    faa_sha256 = sha256_file(faa_path)
    diamond_ver = get_diamond_version(diamond_bin)

    declared_search_params = {
        "subcommand": "blastp",
        "query_faa": os.path.abspath(query_faa),
        "query_sha256": query_sha256,
        "db_file": os.path.abspath(db_path + ".dmnd"),
        "db_sha256": db_sha256,
        "target_faa": os.path.abspath(faa_path),
        "target_faa_sha256": faa_sha256,
        "sensitive": True,
        "evalue": "1e-3",
        "k": 1,
        "outfmt": ["6", "qseqid", "sseqid", "pident", "evalue"],
        "id_cutoff": id_cutoff,
    }

    if os.path.exists(m8_path):
        if not os.path.exists(m8_sig_file):
            raise RuntimeError(
                f"[SEARCH CACHE INTEGRITY ERROR] DIAMOND search hits file '{m8_path}' exists on disk, "
                f"but completion signature file '{m8_sig_file}' is missing! "
                f"Refusing to reuse unverified hits or blindly overwrite. "
                f"Repair action: inspect or remove '{m8_path}' before retrying."
            )
        try:
            with open(m8_sig_file, "r", encoding="utf-8") as f:
                sig = json.load(f)
        except Exception as e:
            raise RuntimeError(
                f"[SEARCH CACHE INTEGRITY ERROR] Corrupted signature file '{m8_sig_file}': {e}. "
                f"Repair action: inspect or remove '{m8_path}' and '{m8_sig_file}'."
            )

        m8_sha256 = sha256_file(m8_path)
        mismatches = []
        if sig.get("query_sha256") != query_sha256:
            mismatches.append(f"query_sha256 (sig: {sig.get('query_sha256')} vs cur: {query_sha256})")
        if sig.get("db_sha256") != db_sha256:
            mismatches.append(f"db_sha256 (sig: {sig.get('db_sha256')} vs cur: {db_sha256})")
        if sig.get("target_faa_sha256") != faa_sha256:
            mismatches.append(f"target_faa_sha256 (sig: {sig.get('target_faa_sha256')} vs cur: {faa_sha256})")
        if sig.get("diamond_version") != diamond_ver:
            mismatches.append(f"diamond_version (sig: {sig.get('diamond_version')} vs cur: {diamond_ver})")
        if sig.get("id_cutoff") != id_cutoff:
            mismatches.append(f"id_cutoff (sig: {sig.get('id_cutoff')} vs cur: {id_cutoff})")
        if sig.get("output_sha256") != m8_sha256:
            mismatches.append(f"output_sha256 (sig: {sig.get('output_sha256')} vs cur: {m8_sha256})")
        if sig.get("exit_code") != 0:
            mismatches.append(f"exit_code (sig: {sig.get('exit_code')})")

        if mismatches:
            raise RuntimeError(
                f"[SEARCH CACHE INTEGRITY ERROR] DIAMOND search hits file '{m8_path}' signature mismatch:\n  "
                + "\n  ".join(mismatches)
                + f"\nFailing closed. Repair action: verify or remove '{m8_path}' and '{m8_sig_file}'."
            )

        hit_count = sum(1 for line in open(m8_path, "r", encoding="utf-8") if line.strip())
        print(f"  [CACHE HIT] Verified search results: '{m8_path}' ({hit_count:,} hits, SHA256: {m8_sha256[:12]}...).")
        return

    print(f"  Running DIAMOND search: {m8_path} (id_cutoff={id_cutoff})...")
    t0 = time.time()
    tmp_m8 = f"{m8_path}.tmp_{os.getpid()}_{int(time.time() * 1000)}"

    cmd_srch = [
        diamond_bin, "blastp",
        "-q", query_faa,
        "-d", db_path,
        "--sensitive",
        "--evalue", "1e-3",
        "-k", "1",
        "-p", str(threads),
        "-f", "6", "qseqid", "sseqid", "pident", "evalue",
        "-o", tmp_m8,
    ]
    if id_cutoff > 0:
        cmd_srch.extend(["--id", str(id_cutoff)])

    res = subprocess.run(cmd_srch, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        if os.path.exists(tmp_m8):
            try:
                os.remove(tmp_m8)
            except OSError:
                pass
        print(f"DIAMOND search failed with code {res.returncode}:\nSTDERR:\n{res.stderr}\nSTDOUT:\n{res.stdout}", file=sys.stderr)
        raise RuntimeError(f"DIAMOND search failed for '{m8_path}' (code {res.returncode}): {res.stderr.strip()}")

    if not os.path.exists(tmp_m8):
        raise FileNotFoundError(f"DIAMOND search exited 0 but expected output '{tmp_m8}' was not found!")

    output_sha256 = sha256_file(tmp_m8)
    os.replace(tmp_m8, m8_path)

    sig_data = {
        "signature_type": "diamond_search_v1",
        "m8_path": os.path.abspath(m8_path),
        "query_faa": os.path.abspath(query_faa),
        "query_sha256": query_sha256,
        "db_file": os.path.abspath(db_path + ".dmnd"),
        "db_sha256": db_sha256,
        "target_faa": os.path.abspath(faa_path),
        "target_faa_sha256": faa_sha256,
        "diamond_version": diamond_ver,
        "search_params": declared_search_params,
        "id_cutoff": id_cutoff,
        "output_sha256": output_sha256,
        "exit_code": 0,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    atomic_write_json(m8_sig_file, sig_data)
    print(f"  Search completed and verified in {time.time() - t0:.1f}s ({os.path.getsize(m8_path) / 1e6:.1f} MB).")


def load_and_validate_hit_ids(m8_path, valid_query_ids, id_cutoff=0):
    hits = set()
    with open(m8_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            line_str = line.strip()
            if not line_str:
                continue
            parts = line_str.split("\t")
            if len(parts) != 4:
                raise ValueError(
                    f"Line {line_no} in '{m8_path}': expected 4 columns (qseqid, sseqid, pident, evalue), got {len(parts)}: '{line_str}'."
                )
            qid, sid, pident_str, evalue_str = parts
            if qid not in valid_query_ids:
                raise ValueError(
                    f"Line {line_no} in '{m8_path}': hit query ID '{qid}' is NOT in the declared query census!"
                )
            try:
                pident = float(pident_str)
            except ValueError:
                raise ValueError(f"Line {line_no} in '{m8_path}': invalid float for pident '{pident_str}'.")
            try:
                evalue = float(evalue_str)
            except ValueError:
                raise ValueError(f"Line {line_no} in '{m8_path}': invalid float for evalue '{evalue_str}'.")

            if not (0.0 <= pident <= 100.0):
                raise ValueError(f"Line {line_no} in '{m8_path}': pident {pident} out of range [0, 100].")
            if evalue < 0.0 or evalue > (1e-3 + 1e-12):
                raise ValueError(f"Line {line_no} in '{m8_path}': evalue {evalue} violates declared threshold <= 1e-3.")
            if id_cutoff > 0 and pident < (float(id_cutoff) - 1e-6):
                raise ValueError(f"Line {line_no} in '{m8_path}': pident {pident} violates declared threshold >= {id_cutoff}.")

            hits.add(qid)
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
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("Validating full query census across dark, control, and FASTA inputs...")
    dark_ids, ctrl_ids, q_lens, census_id_set = validate_and_load_query_census(QUERY_FAA, DARK_IDS_PATH, CTRL_IDS_PATH)
    query_sha256 = sha256_file(QUERY_FAA)
    print(f"  Dark queries:    {len(dark_ids):,} (validated disjoint, non-empty, unique)")
    print(f"  Control queries: {len(ctrl_ids):,} (validated disjoint, non-empty, unique)")
    print(f"  Total census:    {len(census_id_set):,} queries (FASTA SHA256: {query_sha256[:12]}...)")

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
            "harmonized_seed20261015",
            os.path.join(RESULTS_DIR, "matched_smag_harmonized_seed20261015.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_harmonized_seed20261015.faa"),
            0,
            "Reference Draw Sensitivity 1 (Harmonized Pyrodigal v3.7.1): Seed 20261015"
        ),
        (
            "harmonized_seed20261016",
            os.path.join(RESULTS_DIR, "matched_smag_harmonized_seed20261016.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_harmonized_seed20261016.faa"),
            0,
            "Reference Draw Sensitivity 2 (Harmonized Pyrodigal v3.7.1): Seed 20261016"
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
            "Reference Draw Sensitivity 1 (Published Prodigal v2.6.3): Seed 20261015"
        ),
        (
            "published_seed20261016",
            os.path.join(RESULTS_DIR, "matched_smag_published_seed20261016.faa"),
            os.path.join(RESULTS_DIR, "matched_gem_published_seed20261016.faa"),
            0,
            "Reference Draw Sensitivity 2 (Published Prodigal v2.6.3): Seed 20261016"
        ),
    ]

    summary_rows = []
    stratum_rows = []

    for run_lbl, s_faa, g_faa, id_cut, desc in searches:
        print(f"\n==================================================================")
        print(f"Executing Search Pair: {run_lbl} ({desc})")
        print(f"==================================================================")
        if not os.path.exists(s_faa) or not os.path.exists(g_faa):
            raise FileNotFoundError(f"Required comparison {run_lbl} input missing: '{s_faa}' or '{g_faa}'.")

        s_db = s_faa.replace(".faa", "_db")
        g_db = g_faa.replace(".faa", "_db")
        suffix = f"_id{id_cut}" if id_cut > 0 else ""
        s_m8 = os.path.join(RESULTS_DIR, f"hits_smag_{run_lbl}{suffix}.m8")
        g_m8 = os.path.join(RESULTS_DIR, f"hits_gem_{run_lbl}{suffix}.m8")

        run_diamond_search(s_faa, s_db, s_m8, QUERY_FAA, query_sha256, threads, DIAMOND_BIN, id_cutoff=id_cut)
        run_diamond_search(g_faa, g_db, g_m8, QUERY_FAA, query_sha256, threads, DIAMOND_BIN, id_cutoff=id_cut)

        print("  Scoring census hits against 106,042 dark and 106,042 control queries...")
        s_hits = load_and_validate_hit_ids(s_m8, census_id_set, id_cutoff=id_cut)
        g_hits = load_and_validate_hit_ids(g_m8, census_id_set, id_cutoff=id_cut)
        print(f"    Raw hit queries: SMAG={len(s_hits):,}, GEM={len(g_hits):,}")

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
    tmp_summary_path = summary_path + f".tmp_{os.getpid()}"
    with open(tmp_summary_path, "w", encoding="utf-8") as out:
        out.write("run_label\tdescription\tcohort\tid_cutoff\ttotal_queries\tboth_hits_a\tsmag_only_b\t"
                  "gem_only_c\tneither_d\tsmag_hit_pct\tgem_hit_pct\tpaired_delta_pct\tdiscordant_ratio\tmcnemar_pvalue\n")
        for r in summary_rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp_summary_path, summary_path)

    # Write stratum table
    stratum_path = os.path.join(RESULTS_DIR, "matched_length_stratified_results.tsv")
    tmp_stratum_path = stratum_path + f".tmp_{os.getpid()}"
    with open(tmp_stratum_path, "w", encoding="utf-8") as out:
        out.write("run_label\tcohort\tlength_stratum\tid_cutoff\tstratum_queries\tboth_hits_a\tsmag_only_b\t"
                  "gem_only_c\tneither_d\tsmag_hit_pct\tgem_hit_pct\tpaired_delta_pct\tdiscordant_ratio\tmcnemar_pvalue\n")
        for r in stratum_rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp_stratum_path, stratum_path)

    print(f"\nWrote full census summary to {summary_path}")
    print(f"Wrote length-stratified results to {stratum_path}")


if __name__ == "__main__":
    main()
