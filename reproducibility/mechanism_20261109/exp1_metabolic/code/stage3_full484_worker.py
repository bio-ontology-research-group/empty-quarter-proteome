#!/usr/bin/env python3
"""Stage 3 Array Worker of `exp1_metabolic`: Evaluates an assigned shard (`task_id` in `0..num_tasks-1`)
of the full 484 deduplicated genomes (237 EQ MAGs + 247 external comparators) using `evaluate_single_genome_full_protocol`.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from stage2_pilot24_ensemble import evaluate_single_genome_full_protocol


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage 3 Array Worker: Full 484-genome MILP ensemble shard.")
    ap.add_argument("--inventory-484-tsv", type=Path, required=True)
    ap.add_argument("--pool-pkl", type=Path, required=True)
    ap.add_argument("--gapseq-out-dir", type=Path, default=Path("/project/k10018/work/eq-annot-20261030/results/c11_gapseq/out"))
    ap.add_argument("--supp-models-dir", type=Path, default=Path("/scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models"))
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--task-id", type=int, required=True)
    ap.add_argument("--num-tasks", type=int, default=8)
    ap.add_argument("--n-target-solutions", type=int, default=20)
    ap.add_argument("--milp-time-limit-s", type=float, default=6.0)
    ap.add_argument("--workers", type=int, default=32)
    args = ap.parse_args()

    t0 = time.time()
    shard_dir = args.out_dir / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)

    inv_df = pd.read_csv(args.inventory_484_tsv, sep="\t")
    eligible_df = inv_df[
        (inv_df["has_draft_sbml"] == True)
        & (inv_df["qc_error"].fillna("") == "")
    ].sort_values("genome").reset_index(drop=True)

    shard_df = eligible_df.iloc[args.task_id :: args.num_tasks].copy().reset_index(drop=True)
    print(
        f"[Stage 3 Shard {args.task_id}/{args.num_tasks}] Assigned {len(shard_df)} of {len(eligible_df)} genomes",
        flush=True,
    )

    tasks: list[dict[str, Any]] = []
    for _, row in shard_df.iterrows():
        tasks.append(
            {
                "genome": str(row["genome"]),
                "cohort": str(row["cohort"]),
                "canonical_order": str(row["canonical_order"]),
                "canonical_family": str(row["canonical_family"]),
                "canonical_genus": str(row["canonical_genus"]),
                "completeness": row["completeness"],
                "contamination": row["contamination"],
                "gapseq_out_dir": str(args.gapseq_out_dir),
                "supp_models_dir": str(args.supp_models_dir),
                "pool_pkl_path": str(args.pool_pkl),
                "n_target_solutions": args.n_target_solutions,
                "milp_time_limit_s": args.milp_time_limit_s,
                "run_masking_control": False,
            }
        )

    all_summary: list[dict[str, Any]] = []
    all_solution_ledger: list[dict[str, Any]] = []
    all_rxn_freq: list[dict[str, Any]] = []
    all_wkf: list[dict[str, Any]] = []
    all_supp: list[dict[str, Any]] = []

    with ProcessPoolExecutor(max_workers=min(args.workers, max(1, len(tasks)), os.cpu_count() or 32)) as ex:
        futs = {ex.submit(evaluate_single_genome_full_protocol, t): t["genome"] for t in tasks}
        for idx, fut in enumerate(as_completed(futs), 1):
            g_name = futs[fut]
            res = fut.result()
            all_summary.extend(res["summary_rows"])
            all_solution_ledger.extend(res.get("solution_ledger_rows", []))
            all_rxn_freq.extend(res["rxn_freq_rows"])
            all_wkf.extend(res["wkf_rows"])
            all_supp.extend(res["supplement_rows"])
            print(f"[Stage 3 Shard {args.task_id}: {idx}/{len(tasks)}] Completed {g_name}", flush=True)

    tag = f"shard_{args.task_id:02d}"
    pd.DataFrame(all_summary).to_csv(shard_dir / f"exp1_{tag}_summary.tsv", sep="\t", index=False)
    pd.DataFrame(all_solution_ledger).to_csv(shard_dir / f"exp1_{tag}_solution_ledger.tsv", sep="\t", index=False)
    pd.DataFrame(all_rxn_freq).to_csv(shard_dir / f"exp1_{tag}_rxn_freq.tsv", sep="\t", index=False)
    pd.DataFrame(all_wkf).to_csv(shard_dir / f"exp1_{tag}_wkf.tsv", sep="\t", index=False)
    pd.DataFrame(all_supp).to_csv(shard_dir / f"exp1_{tag}_supplement.tsv", sep="\t", index=False)

    dt = time.time() - t0
    print(f"[Stage 3 Shard {args.task_id} Complete in {dt:.1f}s] Wrote {len(shard_df)} genomes to {shard_dir}", flush=True)


if __name__ == "__main__":
    main()
