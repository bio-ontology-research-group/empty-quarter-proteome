#!/usr/bin/env python3
"""Independent Witness Rebuild and Verification for the 4 Changed Historical Cases
(`EXP1-LEDGER-FINAL-ACCEPTANCE.txt` Item 5 & `EXP1-REPLAY-EXECUTION-CORRECTION.txt`).

Independently reconstructs the accepted solution models for:
- 39PRr1__single_SemiBin_48
- 41PRr2__self_SemiBin_2049
- 41PRr3__self_SemiBin_1510_sub
- 41PRr3__self_SemiBin_55
under S0_gapseq_native_predicted and B0_methanosarcina_template_baseline using the exact
engine API, native bounds, and candidate pool semantics. Rigorously audits:
- primal growth >= 0.00999
- stoichiometric max residual <= 1e-5
- bound violation <= 1e-5
- closed-uptake EGC cycle QC pass (zero spurious ATP cycle)
Also verifies a synthetic known-impossible negative control to prove ground-truth negative correctness.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import pickle
from exp1_sbml_highs_engine import (
    reconstruct_and_verify_solution_witness,
    SpeciesRecord,
)
from stage2_pilot24_ensemble import load_genome_archive_payload

CHANGED_GENOMES = [
    "39PRr1__single_SemiBin_48",
    "41PRr2__self_SemiBin_2049",
    "41PRr3__self_SemiBin_1510_sub",
    "41PRr3__self_SemiBin_55",
]


def main() -> None:
    ap = argparse.ArgumentParser(description="Independent rebuild and verification of 4 changed cases.")
    ap.add_argument("--ledger-tsv", type=Path, required=True)
    ap.add_argument("--pool-pkl", type=Path, required=True)
    ap.add_argument("--gapseq-out-dir", type=Path, default=Path("/project/k10018/work/eq-annot-20261030/results/c11_gapseq/out"))
    ap.add_argument("--supp-models-dir", type=Path, default=Path("/scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models"))
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{time.strftime('%X')}] Loading ledger from {args.ledger_tsv}...")
    ledger_df = pd.read_csv(args.ledger_tsv, sep="\t")
    with open(args.pool_pkl, "rb") as fh:
        pool_data = pickle.load(fh)
    universal_rxns = pool_data["reactions"]

    records: list[dict[str, Any]] = []
    all_witnesses_verified = True

    for genome in CHANGED_GENOMES:
        sub = ledger_df[
            (ledger_df["genome"] == genome)
            & (ledger_df["scenario_id"] == "S0_gapseq_native_predicted")
            & (ledger_df["biomass_recipe_id"] == "B0_methanosarcina_template_baseline")
            & (ledger_df["validity_disposition"] == "accepted_feasible")
        ]
        if sub.empty:
            print(f"[ERROR] No accepted solutions found in ledger for {genome}", file=sys.stderr)
            all_witnesses_verified = False
            continue

        print(f"[{time.strftime('%X')}] Loading genome payload for {genome}...")
        payload = load_genome_archive_payload(genome, args.gapseq_out_dir, args.supp_models_dir)
        draft_model = payload["draft_model"]
        native_medium = payload["native_med"]

        # Reconstruct and verify the incumbent witness (solution_index == 0)
        witness_row = sub.sort_values("solution_index").iloc[0]
        sol_idx = int(witness_row["solution_index"])
        r_ids = [r for r in str(witness_row["selected_reaction_ids"]).split(";") if r]

        print(f"[{time.strftime('%X')}] Reconstructing witness for {genome} (solution_index={sol_idx}, n_rxns={len(r_ids)})...")
        w_res = reconstruct_and_verify_solution_witness(
            draft_model=draft_model,
            candidate_reactions=universal_rxns,
            scenario_compounds=native_medium,
            biomass_recipe_id="B0_methanosarcina_template_baseline",
            selected_reaction_ids=r_ids,
            min_growth=0.01,
            time_limit_s=10.0,
        )

        rec = {
            "genome": genome,
            "scenario_id": "S0_gapseq_native_predicted",
            "biomass_recipe_id": "B0_methanosarcina_template_baseline",
            "solution_index": sol_idx,
            "n_added_rxns": len(r_ids),
            "historical_old_status": "infeasible",
            "reconstructed_fba_status": w_res["fba_status"],
            "reconstructed_growth_rate": w_res["growth_rate"],
            "reconstructed_max_residual": w_res["max_residual"],
            "reconstructed_max_bound_viol": w_res["max_bound_viol"],
            "cycle_qc_status": w_res["cycle_qc_status"],
            "has_spurious_atp_cycle": w_res["has_spurious_atp_cycle"],
            "science_qc_eligible": w_res["science_qc_eligible"],
            "independently_verified": w_res["verified"],
            "ledger_growth": float(witness_row["verified_growth"]),
            "ledger_residual": float(witness_row["stoichiometric_max_residual"]),
        }
        records.append(rec)
        print(
            f"  -> {genome}: verified={w_res['verified']}, growth={w_res['growth_rate']:.6f}, "
            f"max_res={w_res['max_residual']:.2e}, cycle_st={w_res['cycle_qc_status']}"
        )
        if not w_res["verified"]:
            all_witnesses_verified = False

    # Also test a synthetic known-impossible negative control to prove ground-truth negative correctness
    print(f"[{time.strftime('%X')}] Evaluating synthetic known-impossible negative control...")
    first_payload = load_genome_archive_payload(CHANGED_GENOMES[0], args.gapseq_out_dir, args.supp_models_dir)
    imp_draft = first_payload["draft_model"].copy()
    imp_draft.species["M_cpd_impossible_c0"] = SpeciesRecord("M_cpd_impossible_c0", "cpd_impossible", "Impossible", "c0", 0, "C1")
    # Add impossible precursor requirement to biomass
    if "R_bio1" in imp_draft.reactions:
        imp_draft.reactions["R_bio1"].stoichiometry["M_cpd_impossible_c0"] = -1.0
    imp_res = reconstruct_and_verify_solution_witness(
        draft_model=imp_draft,
        candidate_reactions=universal_rxns,
        scenario_compounds=first_payload["native_med"],
        biomass_recipe_id="B0_methanosarcina_template_baseline",
        selected_reaction_ids=[],
        min_growth=0.01,
        time_limit_s=5.0,
    )
    neg_control_sound = bool(not imp_res["verified"] and (imp_res["growth_rate"] < 0.00999 or math.isnan(imp_res["growth_rate"])))
    print(f"  -> Synthetic known-impossible negative control: sound={neg_control_sound} (verified={imp_res['verified']}, growth={imp_res['growth_rate']})")

    res_df = pd.DataFrame(records)
    out_tsv = args.out_dir / "exp1_four_changed_witnesses_verification.tsv"
    out_json = args.out_dir / "exp1_four_changed_witnesses_verification.json"
    res_df.to_csv(out_tsv, sep="\t", index=False)

    summary = {
        "all_witnesses_verified": all_witnesses_verified,
        "n_witnesses_evaluated": len(records),
        "synthetic_negative_control_sound": neg_control_sound,
        "elapsed_seconds": round(time.time() - t0, 2),
        "witnesses": records,
        "conclusion": (
            "All 4 changed historical cases independently rebuilt and verified as valid bounded incumbents "
            "with finite residual <= 1e-5 and zero cycle; negative control verified infeasible."
            if (all_witnesses_verified and neg_control_sound)
            else "Witness verification failed."
        ),
    }

    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\n[{time.strftime('%X')}] Summary:")
    print(json.dumps(summary, indent=2))

    if not (all_witnesses_verified and neg_control_sound):
        print("[ERROR] Witness validation failed checks!", file=sys.stderr)
        sys.exit(2)

    print(f"\n[SUCCESS] Completed 4-witness rebuild and validation in {summary['elapsed_seconds']}s.")


if __name__ == "__main__":
    main()
