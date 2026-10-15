#!/usr/bin/env python3
"""
scripts/07_prepare_full_cohorts.py
Prepares and asserts complete, immutable cohorts for full family discovery across F3, F4, and F5:
1. EQ Representatives (all 5,438 MAGs):
   - Reuses Root's verified full cohort:
     /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f1_full/cohort.faa
     and cohort.tsv (20,451,986 proteins across all 5,438 genomes).
   - Asserts Wdb.csv contains exactly 5,438 genomes and 5,438 distinct secondary clusters.
   - Symlinks or verifies cohorts/eq_reps_5438.faa and cohorts/eq_reps_5438_manifest.tsv.
2. Environmental Same-Caller Cohorts (Pyrodigal v3.7.1 harmonized):
   - Asserts SMAG canonical FASTA:
     /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/results/smag_pool_canonical.faa
     and smag_protein_coordinate_mapping.tsv exist and are > 1 GB.
   - Asserts GEM canonical FASTA:
     /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/results/gem_harmonized_pool_canonical_seed20261014.faa
     and gem_harmonized_protein_coordinate_mapping_seed20261014.tsv exist and are > 1 GB.
3. Asserts GTDB-Tk R232 summaries (bac120 and ar53), CheckM genomeInformation.csv,
   and verified profile collections.
4. Emits cohorts/cohorts_preflight.json.
"""

import csv
import json
import os
import shutil
import sys
import time
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
COHORTS_DIR = BASE_DIR / "cohorts"
COHORTS_DIR.mkdir(parents=True, exist_ok=True)
REFS_DIR = BASE_DIR / "references"
PROFILES_DIR = REFS_DIR / "profiles"
META_DIR = BASE_DIR / "metadata"

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
WDB_CSV = W / "05_refine/_drep/out/data_tables/Wdb.csv"
ROOT_F1_FULL_FAA = W / "06_darkmatter/family_discovery_20261015_clean/root_audit/f1_full/cohort.faa"
ROOT_F1_FULL_TSV = W / "06_darkmatter/family_discovery_20261015_clean/root_audit/f1_full/cohort.tsv"

H3_DIR = W / "06_darkmatter/p7_analyses_clean/h3_references/results"
SMAG_CANONICAL = H3_DIR / "smag_pool_canonical.faa"
SMAG_COORD_MAP = H3_DIR / "smag_protein_coordinate_mapping.tsv"
GEM_CANONICAL = H3_DIR / "gem_harmonized_pool_canonical_seed20261014.faa"
GEM_COORD_MAP = H3_DIR / "gem_harmonized_protein_coordinate_mapping_seed20261014.tsv"

GTDB_BAC = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
GTDB_AR = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
CHECKM_CSV = META_DIR / "genomeInformation.csv"

MMSEQS_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs")
DIAMOND_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")


def main():
    t0 = time.time()
    print("=" * 80)
    print("ASSERTING AND PREPARING FULL COHORTS FOR F3-F5 DISCOVERY")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print(f"Cohorts directory: {COHORTS_DIR}")
    print("=" * 80)

    # 1. Assert Wdb.csv exists and has exactly 5,438 genomes
    assert WDB_CSV.exists(), f"CRITICAL: Missing Wdb.csv at {WDB_CSV}"
    with open(WDB_CSV, "r", encoding="utf-8") as f:
        wdb_rows = list(csv.DictReader(f))
    assert len(wdb_rows) == 5438, f"Assertion failed: expected 5438 genomes in Wdb.csv, got {len(wdb_rows)}"
    unique_clusters = {r["cluster"] for r in wdb_rows}
    assert len(unique_clusters) == 5438, f"Assertion failed: expected 5438 unique clusters, got {len(unique_clusters)}"
    print(f"[OK] Wdb.csv: exactly {len(wdb_rows)} genomes and {len(unique_clusters)} unique secondary clusters.")

    # 2. Assert Root's verified full EQ cohort
    assert ROOT_F1_FULL_FAA.exists(), f"CRITICAL: Missing Root cohort.faa at {ROOT_F1_FULL_FAA}"
    assert ROOT_F1_FULL_FAA.stat().st_size > 8_000_000_000, f"Root cohort.faa too small: {ROOT_F1_FULL_FAA.stat().st_size} bytes"
    assert ROOT_F1_FULL_TSV.exists(), f"CRITICAL: Missing Root cohort.tsv at {ROOT_F1_FULL_TSV}"

    with open(ROOT_F1_FULL_TSV, "r", encoding="utf-8") as f:
        tsv_rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(tsv_rows) == 5438, f"Assertion failed: expected 5438 manifest rows in cohort.tsv, got {len(tsv_rows)}"
    total_eq_proteins = sum(int(r["proteins"]) for r in tsv_rows)
    assert total_eq_proteins == 20451986, f"Assertion failed: expected 20,451,986 EQ proteins, got {total_eq_proteins}"
    print(f"[OK] Root EQ Cohort: exactly 5,438 genomes, 20,451,986 proteins ({ROOT_F1_FULL_FAA.stat().st_size / 1e9:.2f} GB).")

    # Symlink to local cohorts/ for transparent access
    local_eq_faa = COHORTS_DIR / "eq_reps_5438.faa"
    local_eq_tsv = COHORTS_DIR / "eq_reps_5438_manifest.tsv"
    if not local_eq_faa.exists():
        local_eq_faa.unlink(missing_ok=True)
        local_eq_faa.symlink_to(ROOT_F1_FULL_FAA)
    if not local_eq_tsv.exists():
        local_eq_tsv.unlink(missing_ok=True)
        local_eq_tsv.symlink_to(ROOT_F1_FULL_TSV)

    local_smag_faa = COHORTS_DIR / "smag_pool_canonical.faa"
    local_gem_faa = COHORTS_DIR / "gem_harmonized_pool_canonical_seed20261014.faa"
    if not local_smag_faa.exists():
        local_smag_faa.unlink(missing_ok=True)
        local_smag_faa.symlink_to(SMAG_CANONICAL)
    if not local_gem_faa.exists():
        local_gem_faa.unlink(missing_ok=True)
        local_gem_faa.symlink_to(GEM_CANONICAL)

    smag_short = COHORTS_DIR / "smag_canonical.faa"
    gem_short = COHORTS_DIR / "gem_canonical.faa"
    if not smag_short.exists():
        smag_short.unlink(missing_ok=True)
        smag_short.symlink_to(SMAG_CANONICAL)
    if not gem_short.exists():
        gem_short.unlink(missing_ok=True)
        gem_short.symlink_to(GEM_CANONICAL)

    # 3. Assert Environmental Cohorts
    assert SMAG_CANONICAL.exists(), f"CRITICAL: Missing SMAG canonical FASTA: {SMAG_CANONICAL}"
    assert SMAG_CANONICAL.stat().st_size > 1_000_000_000, f"SMAG FASTA too small: {SMAG_CANONICAL.stat().st_size}"
    assert SMAG_COORD_MAP.exists(), f"CRITICAL: Missing SMAG coordinate mapping: {SMAG_COORD_MAP}"
    assert SMAG_COORD_MAP.stat().st_size > 1_000_000_000, f"SMAG coord map too small: {SMAG_COORD_MAP.stat().st_size}"
    print(f"[OK] SMAG Environmental Cohort: {SMAG_CANONICAL} ({SMAG_CANONICAL.stat().st_size / 1e9:.2f} GB)")

    assert GEM_CANONICAL.exists(), f"CRITICAL: Missing GEM canonical FASTA: {GEM_CANONICAL}"
    assert GEM_CANONICAL.stat().st_size > 1_000_000_000, f"GEM FASTA too small: {GEM_CANONICAL.stat().st_size}"
    assert GEM_COORD_MAP.exists(), f"CRITICAL: Missing GEM coordinate mapping: {GEM_COORD_MAP}"
    assert GEM_COORD_MAP.stat().st_size > 1_000_000_000, f"GEM coord map too small: {GEM_COORD_MAP.stat().st_size}"
    print(f"[OK] GEM Environmental Cohort: {GEM_CANONICAL} ({GEM_CANONICAL.stat().st_size / 1e9:.2f} GB)")

    # 4. Assert GTDB-Tk R232 summaries
    assert GTDB_BAC.exists(), f"CRITICAL: Missing GTDB-Tk bac120: {GTDB_BAC}"
    assert GTDB_AR.exists(), f"CRITICAL: Missing GTDB-Tk ar53: {GTDB_AR}"
    print(f"[OK] GTDB-Tk R232 summaries verified (bac120 & ar53).")

    # 5. Assert CheckM metadata
    assert CHECKM_CSV.exists(), f"CRITICAL: Missing CheckM metadata: {CHECKM_CSV}"
    print(f"[OK] CheckM metadata verified ({CHECKM_CSV}).")

    # 6. Assert Profile collections and DIAMOND references
    f3_hmm = PROFILES_DIR / "f3_rhodopsins.hmm"
    f4_hmm = PROFILES_DIR / "f4_ectc.hmm"
    f5_hmm = PROFILES_DIR / "f5_carbonic_anhydrase.hmm"
    for hmm in [f3_hmm, f4_hmm, f5_hmm]:
        assert hmm.exists(), f"CRITICAL: Missing profile {hmm}"
        assert hmm.stat().st_size > 100_000, f"Profile too small: {hmm}"
    print(f"[OK] Verified profile collections: f3_rhodopsins.hmm, f4_ectc.hmm, f5_carbonic_anhydrase.hmm.")

    all_refs_dmnd = REFS_DIR / "all_f345_references.dmnd"
    all_refs_faa = REFS_DIR / "all_f345_references.faa"
    assert all_refs_dmnd.exists(), f"CRITICAL: Missing DIAMOND reference db: {all_refs_dmnd}"
    assert all_refs_faa.exists(), f"CRITICAL: Missing reference FASTA: {all_refs_faa}"
    print(f"[OK] Verified reference databases: all_f345_references.dmnd ({all_refs_dmnd.stat().st_size} bytes).")

    # 7. Assert Binaries
    assert MMSEQS_BIN.exists(), f"CRITICAL: Missing MMseqs2 binary at {MMSEQS_BIN}"
    assert DIAMOND_BIN.exists(), f"CRITICAL: Missing DIAMOND binary at {DIAMOND_BIN}"
    print(f"[OK] Verified binaries: MMseqs2 ({MMSEQS_BIN}) and DIAMOND ({DIAMOND_BIN}).")

    # Freeze executed code BEFORE recording preflight
    frozen_dir = BASE_DIR / "results" / "full_v1" / "frozen_code"
    frozen_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(__file__, frozen_dir / "executed_07_prepare_full_cohorts.py")

    preflight = {
        "preflight_verified": True,
        "missing_inputs": 0,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "eq_reps": {
            "genomes": len(wdb_rows),
            "unique_clusters": len(unique_clusters),
            "proteins": total_eq_proteins,
            "cohort_faa": str(ROOT_F1_FULL_FAA),
            "cohort_tsv": str(ROOT_F1_FULL_TSV),
            "faa_size_bytes": ROOT_F1_FULL_FAA.stat().st_size,
        },
        "environmental_cohorts": {
            "smag_canonical_faa": str(SMAG_CANONICAL),
            "smag_faa_size_bytes": SMAG_CANONICAL.stat().st_size,
            "smag_coordinate_map": str(SMAG_COORD_MAP),
            "gem_canonical_faa": str(GEM_CANONICAL),
            "gem_faa_size_bytes": GEM_CANONICAL.stat().st_size,
            "gem_coordinate_map": str(GEM_COORD_MAP),
        },
        "taxonomy": {
            "gtdb_bac120": str(GTDB_BAC),
            "gtdb_ar53": str(GTDB_AR),
        },
        "metadata": {
            "checkm": str(CHECKM_CSV),
            "wdb": str(WDB_CSV),
        },
        "references": {
            "f3_hmm": str(f3_hmm),
            "f4_hmm": str(f4_hmm),
            "f5_hmm": str(f5_hmm),
            "all_refs_dmnd": str(all_refs_dmnd),
            "all_refs_faa": str(all_refs_faa),
        },
        "binaries": {
            "mmseqs2": str(MMSEQS_BIN),
            "diamond": str(DIAMOND_BIN),
        },
        "duration_seconds": round(time.time() - t0, 2),
    }

    out_json = COHORTS_DIR / "cohorts_preflight.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(preflight, f, indent=2)
    print("\nCohort Preflight check complete. Saved to " + str(out_json))
    print("All preflight cohort assertions verified.")


if __name__ == "__main__":
    main()
