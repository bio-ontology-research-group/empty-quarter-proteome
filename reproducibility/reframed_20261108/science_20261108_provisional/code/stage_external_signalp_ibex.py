#!/usr/bin/env python3
"""Prepare bounded external comparator SignalP 6.0h FASTA chunks, manifest, and IBEX Slurm scripts.

Implements `coordinated-revision-20261108/external-signalp-design.txt`:
- Uses the exact 301 external Thermoplasmatota comparators and deduplicates the 54 MGnify duplicates
  of GTDB r232 species representatives -> 247 deduplicated external comparators.
- Extracts their exact protein FASTAs from `/project/k10018/work/eq-annot-20261030/results/c7_thermo_evol/thermo_proteomes.tar.gz`,
  recording SHA-256 and protein counts in `external_signalp_manifest.tsv`.
- Packages all proteins of the 247 deduplicated external comparators into bounded chunks of ~8,000 proteins
  (`chunks/signalp/chunk_0001.faa` .. `chunk_NNNN.faa`) with `>{genome_id}|{raw_header}` identifiers matching
  the EQ 237 SignalP 6.0h run on IBEX (`gpu_annot_20261030`).
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import tarfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_science import identify_mgnify_gtdb_duplicates, require_file, sha256_file

P_ANNOT = Path("/project/k10018/work/eq-annot-20261030")
P_OUT_ROOT = Path(os.environ.get("EQ_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261108"))
P_IN = P_OUT_ROOT / "inputs"
MAX_PROTEINS_PER_CHUNK = 8000


def main() -> None:
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else (P_OUT_ROOT / "external_signalp_ibex_pkg")
    chunks_dir = out_dir / "chunks" / "signalp"
    logs_dir = out_dir / "logs"
    out_sig_dir = out_dir / "out" / "signalp"
    chunks_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    out_sig_dir.mkdir(parents=True, exist_ok=True)

    comp_df = pd.read_csv(require_file(P_IN / "v2_thermo_relatives_comparators.tsv", "v2_thermo_relatives_comparators.tsv"), sep="\t")
    skani_df = pd.read_csv(require_file(P_IN / "v2_thermo_relatives_skani.tsv.gz", "v2_thermo_relatives_skani.tsv.gz"), sep="\t")
    dup_ids = identify_mgnify_gtdb_duplicates(comp_df, skani_df)
    comp_df["dedup_retained"] = ~comp_df["genome_id"].astype(str).isin(dup_ids)
    assert int(comp_df["dedup_retained"].sum()) == 247, f"Expected 247 deduplicated external comparators, got {comp_df['dedup_retained'].sum()}"

    ext_all_set = set(comp_df["genome_id"].astype(str))
    ext_retained_set = set(comp_df.loc[comp_df["dedup_retained"], "genome_id"].astype(str))

    tar_path = require_file(P_ANNOT / "results/c7_thermo_evol/thermo_proteomes.tar.gz", "thermo_proteomes.tar.gz")
    tar_sha = sha256_file(tar_path)

    genome_faa_info: dict[str, dict[str, object]] = {}
    cur_chunk_idx = 1
    cur_chunk_proteins = 0
    cur_chunk_genomes: set[str] = set()
    cur_chunk_lines: list[str] = []
    chunk_summaries: list[dict[str, object]] = []

    def _flush_chunk(c_idx: int, lines: list[str], n_prot: int, g_set: set[str]) -> None:
        if not lines:
            return
        c_name = f"chunk_{c_idx:04d}"
        c_path = chunks_dir / f"{c_name}.faa"
        content = "".join(lines).encode("utf-8")
        c_path.write_bytes(content)
        chunk_summaries.append(
            {
                "chunk_id": c_name,
                "chunk_index": c_idx,
                "n_proteins": n_prot,
                "n_genomes": len(g_set),
                "sha256": hashlib.sha256(content).hexdigest(),
            }
        )

    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf:
            if not member.isfile() or not member.name.endswith(".faa"):
                continue
            gid = Path(member.name).stem
            if gid not in ext_all_set:
                continue
            f = tf.extractfile(member)
            if f is None:
                raise ValueError(f"Could not extract {member.name}")
            raw_bytes = f.read()
            faa_sha = hashlib.sha256(raw_bytes).hexdigest()
            text = raw_bytes.decode("utf-8", errors="replace")

            n_prot_genome = 0
            g_chunks: list[str] = []
            retained = gid in ext_retained_set

            for raw_line in text.splitlines(keepends=True):
                if raw_line.startswith(">"):
                    n_prot_genome += 1
                    if retained:
                        if cur_chunk_proteins >= MAX_PROTEINS_PER_CHUNK:
                            _flush_chunk(cur_chunk_idx, cur_chunk_lines, cur_chunk_proteins, cur_chunk_genomes)
                            cur_chunk_idx += 1
                            cur_chunk_proteins = 0
                            cur_chunk_genomes = set()
                            cur_chunk_lines = []
                        c_label = f"chunk_{cur_chunk_idx:04d}"
                        if not g_chunks or g_chunks[-1] != c_label:
                            g_chunks.append(c_label)
                        cur_chunk_proteins += 1
                        cur_chunk_genomes.add(gid)
                        cur_chunk_lines.append(f">{gid}|{raw_line[1:]}")
                else:
                    if retained:
                        cur_chunk_lines.append(raw_line)

            genome_faa_info[gid] = {
                "n_proteins": n_prot_genome,
                "input_faa_tar_member": member.name,
                "input_sha256": faa_sha,
                "chunk_ids": ",".join(g_chunks) if retained else "excluded_duplicate_of_gtdb_rep",
            }

    if cur_chunk_lines:
        _flush_chunk(cur_chunk_idx, cur_chunk_lines, cur_chunk_proteins, cur_chunk_genomes)

    if len(genome_faa_info) != 301:
        raise ValueError(f"Expected 301 external genomes in thermo_proteomes.tar.gz, found {len(genome_faa_info)}")

    manifest_rows = []
    for _, row in comp_df.iterrows():
        gid = str(row["genome_id"])
        info = genome_faa_info[gid]
        manifest_rows.append(
            {
                "genome_id": gid,
                "source_cohort": str(row["source"]),
                "dedup_retained": bool(row["dedup_retained"]),
                "same_assembly_as_gtdb_rep": str(row["same_assembly_as_gtdb_rep"]),
                "lineage": str(row["lineage"]),
                "completeness": float(row["completeness"]),
                "contamination": float(row["contamination"]),
                "n_proteins": int(info["n_proteins"]),
                "input_tar_path": str(tar_path),
                "input_tar_sha256": tar_sha,
                "input_faa_tar_member": str(info["input_faa_tar_member"]),
                "input_sha256": str(info["input_sha256"]),
                "chunk_ids": str(info["chunk_ids"]),
                "prior_ibex_signalp_status": "not_in_gpu_annot_20261030_eq237_only",
            }
        )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(out_dir / "external_signalp_manifest.tsv", sep="\t", index=False)

    chunks_df = pd.DataFrame(chunk_summaries)
    chunks_df.to_csv(out_dir / "external_signalp_chunks_summary.tsv", sep="\t", index=False)

    n_chunks = len(chunks_df)
    total_retained_proteins = int(chunks_df["n_proteins"].sum())

    sbatch_content = f"""#!/bin/bash
# Bounded External Comparator SignalP 6.0h Array on IBEX CPU nodes
# Matches exact settings of /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/gpu_annot_20261030/signalp6_ibex.sbatch
# Total chunks: {n_chunks} ({total_retained_proteins} proteins across 247 deduplicated external Thermoplasmatota comparators)
# Smoke test chunk 1 first:
#   sbatch --array=1 signalp6_external_ibex.sbatch
# Then submit remaining chunks with %32 cap:
#   sbatch --array=2-{n_chunks}%32 signalp6_external_ibex.sbatch
#SBATCH --job-name=eq_ext_sigp6
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp/logs/%x_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp/logs/%x_%A_%a.err

set -euo pipefail
module load signalp/6.0h
B=/ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp
C=$(printf 'chunk_%04d' "$SLURM_ARRAY_TASK_ID")
IN="$B/chunks/signalp/$C.faa"
OUT="$B/out/signalp/$C"
if [ -f "$OUT.done" ] && [ -s "$OUT.prediction_results.txt.gz" ]; then
  echo "[IBEX SignalP $C] Already completed. Skipping."
  exit 0
fi
mkdir -p "$B/out/signalp" "$B/logs"
T=$(mktemp -d "$OUT.tmp.XXXX")
trap 'rm -rf "$T"' EXIT
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
signalp6 --fastafile "$IN" --organism other --output_dir "$T" --format none --mode fast \\
         --write_procs "$SLURM_CPUS_PER_TASK"
gzip -c "$T/prediction_results.txt" > "$OUT.prediction_results.txt.gz.tmp"
mv "$OUT.prediction_results.txt.gz.tmp" "$OUT.prediction_results.txt.gz"
if [ -f "$T/output.gff3" ]; then
  gzip -c "$T/output.gff3" > "$OUT.gff3.gz.tmp"
  mv "$OUT.gff3.gz.tmp" "$OUT.gff3.gz"
fi
rm -rf "$T"
date -Is > "$OUT.done"
echo "[IBEX SignalP $C] Completed at $(date -Is)."
"""
    (out_dir / "signalp6_external_ibex.sbatch").write_text(sbatch_content)

    collect_py = """#!/usr/bin/env python3
import csv, gzip
from pathlib import Path

B = Path("/ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp")
with open(B / "external_signalp_manifest.tsv", newline="") as fh:
    manifest_rows = list(csv.DictReader(fh, delimiter="\\t"))
with open(B / "external_signalp_chunks_summary.tsv", newline="") as fh:
    chunks_rows = list(csv.DictReader(fh, delimiter="\\t"))

out_gz = B / "out" / "external247_signalp6_predictions.tsv.gz"
counts_tot = {}
counts_sp = {}
n_total_scored = 0

with gzip.open(out_gz, "wt", newline="") as out_fh:
    w = csv.writer(out_fh, delimiter="\\t")
    w.writerow(["genome_id", "gene_id", "signalp_prediction", "has_signal_peptide", "signalp_non_other_prob", "chunk_id"])
    for r in chunks_rows:
        c_id = r["chunk_id"]
        p_gz = B / "out" / "signalp" / f"{c_id}.prediction_results.txt.gz"
        done_f = B / "out" / "signalp" / f"{c_id}.done"
        if not (p_gz.exists() and done_f.exists()):
            raise FileNotFoundError(f"Missing completed chunk output for {c_id}")
        with gzip.open(p_gz, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\\n").split("\\t")
                if len(parts) < 8:
                    continue
                hdr = parts[0]
                pred = parts[1]
                p_other = float(parts[2])
                tok = hdr.split()[0]
                gid, gene_id = tok.split("|", 1)
                has_sp = (pred != "OTHER")
                w.writerow([gid, gene_id, pred, "True" if has_sp else "False", f"{1.0 - p_other:.6f}", c_id])
                counts_tot[gid] = counts_tot.get(gid, 0) + 1
                if has_sp:
                    counts_sp[gid] = counts_sp.get(gid, 0) + 1
                n_total_scored += 1

with open(B / "out" / "external247_signalp6_per_genome_summary.tsv", "w", newline="") as fh:
    fieldnames = list(manifest_rows[0].keys()) + ["n_proteins_scored", "n_signalp_positive", "pct_signalp_positive"]
    w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\\t")
    w.writeheader()
    for m in manifest_rows:
        if m["dedup_retained"] != "True":
            continue
        gid = m["genome_id"]
        n_sc = counts_tot.get(gid, 0)
        n_pos = counts_sp.get(gid, 0)
        if int(m["n_proteins"]) != n_sc:
            raise ValueError(f"Protein count mismatch for {gid}: manifest={m['n_proteins']}, scored={n_sc}")
        row_out = dict(m)
        row_out["n_proteins_scored"] = n_sc
        row_out["n_signalp_positive"] = n_pos
        row_out["pct_signalp_positive"] = f"{100.0 * n_pos / n_sc:.4f}" if n_sc else "0.0"
        w.writerow(row_out)

print(f"Successfully bundled {n_total_scored} predictions across {len(counts_tot)} deduplicated external genomes.")
"""
    (out_dir / "collect_external_signalp.py").write_text(collect_py)

    collect_sbatch = f"""#!/bin/bash
#SBATCH --job-name=eq_ext_sigp6_col
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp/logs/%x_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp/logs/%x_%j.err
set -euo pipefail
python3 -s /ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp/collect_external_signalp.py
"""
    (out_dir / "collect_external_signalp.sbatch").write_text(collect_sbatch)

    summary = {
        "n_external_comparators_total": 301,
        "n_external_comparators_dedup_retained": 247,
        "n_external_comparators_excluded_duplicates": 54,
        "n_proteins_dedup_retained": total_retained_proteins,
        "n_chunks": n_chunks,
        "max_proteins_per_chunk": MAX_PROTEINS_PER_CHUNK,
    }
    (out_dir / "package_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
