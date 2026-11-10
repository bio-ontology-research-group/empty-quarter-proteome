#!/usr/bin/env python3
"""Stage 1 of `exp1_metabolic`: Full 538/484-genome Gapseq SBML inventory, Methanosarcina biomass template audit,
mass/charge & closed-uptake EGC QC, deterministic cohort-union candidate reaction pool assembly with conflict/balance
auditing (`EXP1-INVENTORY-REVIEW.txt`), and Pilot24 + Reference8 selection.
"""
from __future__ import annotations

import argparse
import csv
import os
import pickle
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from exp1_sbml_highs_engine import (
    NUTRIENT_SCENARIOS,
    OBSERVED_COHORT_538_POOL_SCOPE,
    SCEN_S3_ANOXIC_DESERT,
    SCEN_S4_OXIC_DESERT,
    ReactionRecord,
    SpeciesRecord,
    aggregate_candidate_pool_deterministically,
    apply_nutrient_scenario_bounds,
    audit_model_mass_charge_balance,
    detect_closed_uptake_egc,
    parse_gapseq_sbml_bytes,
    parse_genome_reaction_and_transporter_evidence,
    solve_fba_highs,
    validate_archaeal_biomass_and_audit_template,
)


def load_frozen_external_247_manifest(
    manifest_tsv: Path,
    comparators_df: pd.DataFrame,
    expected_total_external: int | None = 301,
    expected_retained_external: int | None = 247,
) -> tuple[set[str], set[str]]:
    """Load the frozen external comparator membership from `external_signalp_manifest.tsv` (`EXP1-INVENTORY-REVIEW.txt`).
    Verifies exact set equality of all external comparators and exact `dedup_retained == True` membership.
    """
    man_df = pd.read_csv(manifest_tsv, sep="\t")
    all_man_ids = set(man_df["genome_id"].astype(str))
    comp_col = "genome_id" if "genome_id" in comparators_df.columns else "genome"
    comp_ids = set(comparators_df[comp_col].astype(str))
    if all_man_ids != comp_ids:
        diff1 = sorted(all_man_ids - comp_ids)[:5]
        diff2 = sorted(comp_ids - all_man_ids)[:5]
        raise RuntimeError(
            f"External comparator set mismatch between external_signalp_manifest.tsv ({len(all_man_ids)}) "
            f"and v2_thermo_relatives_comparators.tsv ({len(comp_ids)}): man-comp={diff1}, comp-man={diff2}"
        )
    if expected_total_external is not None and len(all_man_ids) != expected_total_external:
        raise RuntimeError(
            f"Expected {expected_total_external} total external genomes in {manifest_tsv}, got {len(all_man_ids)}"
        )
    retained_mask = man_df["dedup_retained"].astype(str).str.lower().isin({"true", "1", "yes"})
    retained_247 = set(man_df.loc[retained_mask, "genome_id"].astype(str))
    if expected_retained_external is not None and len(retained_247) != expected_retained_external:
        raise RuntimeError(
            f"Expected exactly {expected_retained_external} dedup_retained genomes in {manifest_tsv}, got {len(retained_247)}"
        )
    return all_man_ids, retained_247


def build_cohort_metadata_538_and_484(
    ref_rep_df: pd.DataFrame,
    legacy_eq_df: pd.DataFrame,
    comp_df: pd.DataFrame,
    external_manifest_tsv: Path,
    expected_eq_count: int | None = 237,
    expected_total_external: int | None = 301,
    expected_retained_external: int | None = 247,
) -> pd.DataFrame:
    """Build the unified cohort metadata table across canonical EQ representatives (`refined_representatives.tsv`),
    legacy EQ taxonomy (`thermoplasmatota.tsv`), and frozen external comparators (`v2_thermo_relatives_comparators.tsv`
    + `external_signalp_manifest.tsv`).
    """
    _all_ext, frozen_retained_ext = load_frozen_external_247_manifest(
        external_manifest_tsv,
        comp_df,
        expected_total_external=expected_total_external,
        expected_retained_external=expected_retained_external,
    )
    eq_g_col = "genome" if "genome" in legacy_eq_df.columns else "Name"
    eq_tax_col = "gtdb_classification" if "gtdb_classification" in legacy_eq_df.columns else "classification"
    legacy_tax_map = dict(zip(legacy_eq_df[eq_g_col].astype(str), legacy_eq_df[eq_tax_col].astype(str)))
    eq_ids = set(legacy_eq_df[eq_g_col].astype(str))
    if expected_eq_count is not None and len(eq_ids) != expected_eq_count:
        raise RuntimeError(f"Expected {expected_eq_count} EQ genomes in thermoplasmatota.tsv, got {len(eq_ids)}")

    ref_g_col = "mag" if "mag" in ref_rep_df.columns else "Name"
    ref_tax_col = "gtdb_r232_classification" if "gtdb_r232_classification" in ref_rep_df.columns else "classification"
    ref_eq_df = ref_rep_df[ref_rep_df[ref_g_col].astype(str).isin(eq_ids)].copy()
    if len(ref_eq_df) != len(eq_ids):
        raise RuntimeError(
            f"Expected all {len(eq_ids)} EQ genomes in refined_representatives.tsv, found {len(ref_eq_df)}"
        )

    cohort_rows: list[dict[str, Any]] = []
    for _, r in ref_eq_df.iterrows():
        g = str(r[ref_g_col])
        canon_cls = str(r[ref_tax_col])
        if "p__Thermoplasmatota" not in canon_cls:
            raise RuntimeError(f"EQ genome {g} not in p__Thermoplasmatota: {canon_cls}")
        leg_cls = legacy_tax_map.get(g, canon_cls)
        parts = {p.split("__")[0]: p for p in canon_cls.split(";") if "__" in p}
        cohort_rows.append(
            {
                "genome": g,
                "cohort": "EQ_desert",
                "source": "EQ_MAG",
                "is_deduplicated_484": True,
                "canonical_classification": canon_cls,
                "legacy_classification": leg_cls,
                "taxonomy_changed_canonical_vs_legacy": canon_cls != leg_cls,
                "canonical_order": parts.get("o", "o__"),
                "canonical_family": parts.get("f", "f__"),
                "canonical_genus": parts.get("g", "g__"),
                "completeness": float(r["Completeness"]),
                "contamination": float(r["Contamination"]),
                "quality_tier": str(r.get("quality_tier", "UNKNOWN")),
                "total_contigs": int(r.get("Total_Contigs", 0)) if pd.notna(r.get("Total_Contigs")) else 0,
                "contig_n50": float(r.get("Contig_N50", 0.0)) if pd.notna(r.get("Contig_N50")) else 0.0,
            }
        )

    comp_g_col = "genome_id" if "genome_id" in comp_df.columns else "genome"
    comp_tax_col = "lineage" if "lineage" in comp_df.columns else "gtdb_taxonomy"
    comp_src_col = "source_cohort" if "source_cohort" in comp_df.columns else "source"
    for _, r in comp_df.iterrows():
        g = str(r[comp_g_col])
        cls_str = str(r[comp_tax_col])
        parts = {p.split("__")[0]: p for p in cls_str.split(";") if "__" in p}
        comp_val = float(r["completeness"]) if pd.notna(r.get("completeness")) else float("nan")
        cont_val = float(r["contamination"]) if pd.notna(r.get("contamination")) else float("nan")
        cohort_rows.append(
            {
                "genome": g,
                "cohort": "external_comparator",
                "source": str(r[comp_src_col]),
                "is_deduplicated_484": g in frozen_retained_ext,
                "canonical_classification": cls_str,
                "legacy_classification": cls_str,
                "taxonomy_changed_canonical_vs_legacy": False,
                "canonical_order": parts.get("o", "o__"),
                "canonical_family": parts.get("f", "f__"),
                "canonical_genus": parts.get("g", "g__"),
                "completeness": comp_val,
                "contamination": cont_val,
                "quality_tier": "HQ" if (pd.notna(comp_val) and comp_val >= 90.0 and pd.notna(cont_val) and cont_val <= 5.0) else "MQ_or_unknown",
                "total_contigs": 0,
                "contig_n50": 0.0,
            }
        )

    cohort_meta_df = pd.DataFrame(cohort_rows)
    expected_dedup_set = eq_ids | frozen_retained_ext
    actual_dedup_set = set(cohort_meta_df.loc[cohort_meta_df["is_deduplicated_484"], "genome"].astype(str))
    if actual_dedup_set != expected_dedup_set:
        raise RuntimeError(f"Deduplicated set equality check failed: len={len(actual_dedup_set)} vs {len(expected_dedup_set)}")
    return cohort_meta_df


def load_excluded_reactions(excluded_tsv: Path) -> set[str]:
    excl: set[str] = set()
    if not excluded_tsv.exists():
        return excl
    with open(excluded_tsv, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("id"):
                continue
            rid = line.split("\t")[0].strip()
            if rid.startswith("rxn"):
                excl.add(rid)
    return excl


def load_gapseq_seed_reference_tables(
    seed_reactions_tsv: Path,
    seed_metabolites_tsv: Path,
) -> tuple[dict[str, str], dict[str, tuple[str | None, int | None]]]:
    """Load canonical reaction reversibility (`seed_reactions_corrected.tsv`) and metabolite formula/charge
    (`seed_metabolites_edited.tsv`) from the frozen gapseq reference database.
    """
    rxn_dirs: dict[str, str] = {}
    if seed_reactions_tsv.exists():
        with open(seed_reactions_tsv, "r", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rid = (row.get("id") or "").strip()
                rev = (row.get("reversibility") or "").strip()
                if rid.startswith("rxn") and rev:
                    rxn_dirs[rid] = rev

    met_meta: dict[str, tuple[str | None, int | None]] = {}
    if seed_metabolites_tsv.exists():
        with open(seed_metabolites_tsv, "r", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                cid = (row.get("id") or "").strip()
                if not cid.startswith("cpd"):
                    continue
                f_raw = (row.get("formula") or "").strip()
                f_val = f_raw if (f_raw and f_raw.lower() not in {"null", "na", "none", ""}) else None
                q_raw = (row.get("charge") or "").strip()
                try:
                    q_val: int | None = int(float(q_raw)) if (q_raw and q_raw.lower() not in {"null", "na", "none", ""}) else None
                except ValueError:
                    q_val = None
                met_meta[cid] = (f_val, q_val)

    return rxn_dirs, met_meta


def inspect_single_genome(task: dict[str, Any]) -> dict[str, Any]:
    """Worker function to parse and QC a single genome's gapseq archive or supplementary directory."""
    g = task["genome"]
    tar_path = Path(task["tar_path"]) if task.get("tar_path") else None
    supp_dir = Path(task["supp_dir"]) if task.get("supp_dir") else None

    draft_bytes: bytes | None = None
    desert_bytes: bytes | None = None
    native_bytes: bytes | None = None
    rxn_tbl_text = ""
    trans_tbl_text = ""
    archive_source = "missing"

    if tar_path and tar_path.exists():
        archive_source = "c11_gapseq_tar_gz"
        with tarfile.open(tar_path, "r:gz") as tf:
            names = set(tf.getnames())
            d_key = f"{g}/{g}-draft.xml"
            f_des_key = f"{g}/fill_desert/{g}.xml"
            f_nat_key = f"{g}/fill_native/{g}.xml"
            r_tbl_key = f"{g}/{g}-all-Reactions.tbl"
            t_tbl_key = f"{g}/{g}-Transporter.tbl"
            if d_key in names:
                draft_bytes = tf.extractfile(d_key).read()  # type: ignore[union-attr]
            if f_des_key in names:
                desert_bytes = tf.extractfile(f_des_key).read()  # type: ignore[union-attr]
            if f_nat_key in names:
                native_bytes = tf.extractfile(f_nat_key).read()  # type: ignore[union-attr]
            if r_tbl_key in names:
                rxn_tbl_text = tf.extractfile(r_tbl_key).read().decode("utf-8", errors="replace")  # type: ignore[union-attr]
            if t_tbl_key in names:
                trans_tbl_text = tf.extractfile(t_tbl_key).read().decode("utf-8", errors="replace")  # type: ignore[union-attr]
    elif supp_dir and supp_dir.exists():
        archive_source = "supplementary_recovered_draft"
        d_file = supp_dir / f"{g}-draft.xml"
        r_file = supp_dir / f"{g}-all-Reactions.tbl"
        t_file = supp_dir / f"{g}-Transporter.tbl"
        if d_file.exists():
            draft_bytes = d_file.read_bytes()
        if r_file.exists():
            rxn_tbl_text = r_file.read_text(encoding="utf-8", errors="replace")
        if t_file.exists():
            trans_tbl_text = t_file.read_text(encoding="utf-8", errors="replace")

    if draft_bytes is None:
        return {
            "genome": g,
            "archive_source": archive_source,
            "has_draft_sbml": False,
            "has_fill_desert_sbml": False,
            "has_fill_native_sbml": False,
            "qc_error": "missing_draft_sbml",
            "extracted_reactions": {},
            "extracted_species": {},
            "draft_reaction_ids": set(),
            "fill_desert_reaction_ids": set(),
            "fill_native_reaction_ids": set(),
        }

    rxn_evi, trans_evi = parse_genome_reaction_and_transporter_evidence(rxn_tbl_text, trans_tbl_text)
    draft_mod = parse_gapseq_sbml_bytes(draft_bytes, provenance_default="draft_sequence_supported")
    bio_audit = validate_archaeal_biomass_and_audit_template(draft_mod, rxn_evidence_map=rxn_evi)
    draft_bal = audit_model_mass_charge_balance(draft_mod)
    draft_egc = detect_closed_uptake_egc(draft_mod, time_limit_s=5.0)

    s3_cpds = NUTRIENT_SCENARIOS[SCEN_S3_ANOXIC_DESERT]["compounds"]
    s4_cpds = NUTRIENT_SCENARIOS[SCEN_S4_OXIC_DESERT]["compounds"]
    draft_s3 = solve_fba_highs(apply_nutrient_scenario_bounds(draft_mod, s3_cpds), time_limit_s=5.0)
    draft_s4 = solve_fba_highs(apply_nutrient_scenario_bounds(draft_mod, s4_cpds), time_limit_s=5.0)

    extracted_rxns: dict[str, ReactionRecord] = {}
    extracted_species: dict[str, SpeciesRecord] = dict(draft_mod.species)
    draft_rids: set[str] = set()
    des_rids: set[str] = set()
    nat_rids: set[str] = set()

    for rid, rxn in draft_mod.reactions.items():
        if not rxn.is_biomass:
            draft_rids.add(rid)
            rc = ReactionRecord(**vars(rxn))
            rc.gpr_ast = None
            rc.genes = set()
            rc.provenance = "observed_in_cohort_draft_models"
            extracted_rxns[rid] = rc

    des_n_rxns = 0
    des_egc_cycle = None
    des_s3_growth = float("nan")
    des_s4_growth = float("nan")
    if desert_bytes is not None:
        des_mod = parse_gapseq_sbml_bytes(desert_bytes, provenance_default="fill_desert")
        des_n_rxns = len(des_mod.reactions)
        des_egc = detect_closed_uptake_egc(des_mod, time_limit_s=5.0)
        des_egc_cycle = des_egc["has_spurious_atp_cycle"]
        des_s3 = solve_fba_highs(apply_nutrient_scenario_bounds(des_mod, s3_cpds), time_limit_s=5.0)
        des_s4 = solve_fba_highs(apply_nutrient_scenario_bounds(des_mod, s4_cpds), time_limit_s=5.0)
        des_s3_growth = float(des_s3["growth_rate"])
        des_s4_growth = float(des_s4["growth_rate"])
        extracted_species.update(des_mod.species)
        for rid, rxn in des_mod.reactions.items():
            if not rxn.is_biomass:
                des_rids.add(rid)
                if rid not in extracted_rxns:
                    rc = ReactionRecord(**vars(rxn))
                    rc.gpr_ast = None
                    rc.genes = set()
                    rc.provenance = "observed_only_in_cohort_gapfilled_models"
                    extracted_rxns[rid] = rc

    nat_n_rxns = 0
    nat_egc_cycle = None
    if native_bytes is not None:
        nat_mod = parse_gapseq_sbml_bytes(native_bytes, provenance_default="fill_native")
        nat_n_rxns = len(nat_mod.reactions)
        nat_egc = detect_closed_uptake_egc(nat_mod, time_limit_s=5.0)
        nat_egc_cycle = nat_egc["has_spurious_atp_cycle"]
        extracted_species.update(nat_mod.species)
        for rid, rxn in nat_mod.reactions.items():
            if not rxn.is_biomass:
                nat_rids.add(rid)
                if rid not in extracted_rxns:
                    rc = ReactionRecord(**vars(rxn))
                    rc.gpr_ast = None
                    rc.genes = set()
                    rc.provenance = "observed_only_in_cohort_gapfilled_models"
                    extracted_rxns[rid] = rc

    return {
        "genome": g,
        "archive_source": archive_source,
        "has_draft_sbml": True,
        "has_fill_desert_sbml": desert_bytes is not None,
        "has_fill_native_sbml": native_bytes is not None,
        "qc_error": "",
        "n_draft_species": len(draft_mod.species),
        "n_draft_reactions": len(draft_mod.reactions),
        "n_draft_genes": len(draft_mod.gene_products),
        "n_fill_desert_reactions": des_n_rxns,
        "n_fill_native_reactions": nat_n_rxns,
        "n_blast_rxn_evidence_hits": len(rxn_evi),
        "n_blast_core_candidates": sum(1 for v in rxn_evi.values() if v["is_core_candidate"]),
        "n_transporter_evidence_hits": len(trans_evi),
        **draft_bal,
        "draft_egc_status": draft_egc["egc_check_status"],
        "draft_has_spurious_atp_cycle": draft_egc["has_spurious_atp_cycle"],
        "draft_closed_uptake_atp_flux": draft_egc["closed_uptake_atp_dissipation_flux"],
        "fill_desert_has_spurious_atp_cycle": des_egc_cycle,
        "fill_native_has_spurious_atp_cycle": nat_egc_cycle,
        "draft_s3_anoxic_desert_growth": float(draft_s3["growth_rate"]),
        "draft_s4_oxic_desert_growth": float(draft_s4["growth_rate"]),
        "fill_desert_s3_anoxic_desert_growth": des_s3_growth,
        "fill_desert_s4_oxic_desert_growth": des_s4_growth,
        **bio_audit,
        "extracted_reactions": extracted_rxns,
        "extracted_species": extracted_species,
        "draft_reaction_ids": draft_rids,
        "fill_desert_reaction_ids": des_rids,
        "fill_native_reaction_ids": nat_rids,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage 1: Full 538/484-genome Gapseq SBML inventory, QC, and Pilot24 selection.")
    ap.add_argument("--science-inputs", type=Path, default=Path("/scratch/project/k10018/eq-science-20261108/inputs"))
    ap.add_argument(
        "--external-manifest-tsv",
        type=Path,
        default=Path("/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/inputs/external_signalp_manifest.tsv"),
    )
    ap.add_argument("--gapseq-out-dir", type=Path, default=Path("/project/k10018/work/eq-annot-20261030/results/c11_gapseq/out"))
    ap.add_argument("--supp-models-dir", type=Path, default=Path("/scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models"))
    ap.add_argument("--excluded-rxns-tsv", type=Path, default=Path("/project/k10018/sw/envs/gapseq/share/gapseq/dat/biomass/excluded_reactions.tsv"))
    ap.add_argument("--seed-reactions-tsv", type=Path, default=Path("/project/k10018/sw/envs/gapseq/share/gapseq/dat/seed_reactions_corrected.tsv"))
    ap.add_argument("--seed-metabolites-tsv", type=Path, default=Path("/project/k10018/sw/envs/gapseq/share/gapseq/dat/seed_metabolites_edited.tsv"))
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=64)
    args = ap.parse_args()

    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load canonical EQ 237 metadata + legacy EQ taxonomy + frozen 247 external comparator manifest
    ref_rep_df = pd.read_csv(args.science_inputs / "refined_representatives.tsv", sep="\t", low_memory=False)
    legacy_eq_df = pd.read_csv(args.science_inputs / "thermoplasmatota.tsv", sep="\t")
    comp_df = pd.read_csv(args.science_inputs / "v2_thermo_relatives_comparators.tsv", sep="\t")

    cohort_meta_df = build_cohort_metadata_538_and_484(
        ref_rep_df=ref_rep_df,
        legacy_eq_df=legacy_eq_df,
        comp_df=comp_df,
        external_manifest_tsv=args.external_manifest_tsv,
        expected_eq_count=237,
        expected_total_external=301,
        expected_retained_external=247,
    )
    assert len(cohort_meta_df) == 538, f"Expected 538 total cohort genomes, got {len(cohort_meta_df)}"
    assert int(cohort_meta_df["is_deduplicated_484"].sum()) == 484

    # 2. Dispatch parallel inspection across all 538 genomes
    tasks: list[dict[str, Any]] = []
    for g in sorted(cohort_meta_df["genome"].astype(str)):
        tar_p = args.gapseq_out_dir / f"{g}.tar.gz"
        supp_p = args.supp_models_dir / g
        tasks.append(
            {
                "genome": g,
                "tar_path": str(tar_p) if tar_p.exists() else "",
                "supp_dir": str(supp_p) if supp_p.exists() else "",
            }
        )

    excluded_base_ids = load_excluded_reactions(args.excluded_rxns_tsv)
    seed_rxn_dirs, seed_met_meta = load_gapseq_seed_reference_tables(args.seed_reactions_tsv, args.seed_metabolites_tsv)

    per_genome_pool_extractions: list[dict[str, Any]] = []
    qc_records: list[dict[str, Any]] = []

    with ProcessPoolExecutor(max_workers=min(args.workers, os.cpu_count() or 32)) as ex:
        futs = {ex.submit(inspect_single_genome, t): t["genome"] for t in tasks}
        for fut in as_completed(futs):
            res = fut.result()
            per_genome_pool_extractions.append(
                {
                    "genome": res["genome"],
                    "extracted_reactions": res.pop("extracted_reactions", {}),
                    "extracted_species": res.pop("extracted_species", {}),
                    "draft_reaction_ids": res.pop("draft_reaction_ids", set()),
                    "fill_desert_reaction_ids": res.pop("fill_desert_reaction_ids", set()),
                    "fill_native_reaction_ids": res.pop("fill_native_reaction_ids", set()),
                }
            )
            qc_records.append(res)

    # 3. Deterministic candidate reaction pool aggregation, conflict check, and balance audit (`EXP1-INVENTORY-REVIEW.txt`)
    pool_agg = aggregate_candidate_pool_deterministically(
        per_genome_extractions=per_genome_pool_extractions,
        seed_rxn_directions=seed_rxn_dirs,
        seed_metabolite_metadata=seed_met_meta,
        excluded_base_ids=excluded_base_ids,
    )
    eligible_rxns: dict[str, ReactionRecord] = pool_agg["reactions"]
    canonical_species: dict[str, SpeciesRecord] = pool_agg["species"]

    pd.DataFrame(pool_agg["pool_audit_rows"]).to_csv(
        args.out_dir / "exp1_candidate_reaction_pool_audit.tsv", sep="\t", index=False
    )
    rxn_conf_df = (
        pd.DataFrame(pool_agg["reaction_conflicts"])
        if pool_agg["reaction_conflicts"]
        else pd.DataFrame(columns=["rxn_id", "donor_a", "donor_b", "conflict_reason"])
    )
    rxn_conf_df.to_csv(args.out_dir / "exp1_candidate_reaction_pool_conflicts.tsv", sep="\t", index=False)
    sp_conf_df = (
        pd.DataFrame(pool_agg["species_conflicts"])
        if pool_agg["species_conflicts"]
        else pd.DataFrame(columns=["species_id", "donor_a", "donor_b", "conflict_reason"])
    )
    sp_conf_df.to_csv(args.out_dir / "exp1_candidate_species_conflicts.tsv", sep="\t", index=False)

    qc_df = pd.DataFrame(qc_records)
    full_538_df = cohort_meta_df.merge(qc_df, on="genome", how="left")
    full_538_df = full_538_df.sort_values(["cohort", "genome"]).reset_index(drop=True)
    dedup_484_df = full_538_df[full_538_df["is_deduplicated_484"]].copy().reset_index(drop=True)

    full_538_df.to_csv(args.out_dir / "exp1_model_inventory_and_qc_538.tsv", sep="\t", index=False)
    dedup_484_df.to_csv(args.out_dir / "exp1_model_inventory_and_qc_484.tsv", sep="\t", index=False)

    # 4. Write dedicated Methanosarcina biomass template audit table for the 484 deduplicated genomes
    bio_cols = [
        "genome",
        "cohort",
        "source",
        "canonical_order",
        "canonical_family",
        "canonical_genus",
        "completeness",
        "contamination",
        "archive_source",
        "biomass_rxn_id",
        "biomass_rxn_name",
        "biomass_provenance",
        "n_biomass_reactants",
        "gam_atp_stoichiometry",
        "f420_demanded_in_template",
        "f420_sequence_supported_synthesis_in_draft",
        "f420_sequence_supported_use_in_draft",
        "n_non_f420_methanosarcina_cofactors_demanded",
        "n_all8_template_trace_cofactors_demanded",
        "template_trace_cofactors_demanded_csv",
        "methanofuran_demanded_in_template",
        "n_template_trace_cofactors_supported_in_draft",
        "n_template_trace_cofactors_unsupported_in_draft",
        "template_trace_cofactors_unsupported_in_draft_csv",
        "n_archaeal_ether_lipids_and_polyamines_demanded",
    ]
    dedup_484_df.reindex(columns=bio_cols).to_csv(args.out_dir / "exp1_biomass_template_audit_484.tsv", sep="\t", index=False)

    # 5. Save deterministic candidate reaction pool pickle
    pool_payload = {
        "candidate_pool_scope": OBSERVED_COHORT_538_POOL_SCOPE,
        "n_candidate_reactions": len(eligible_rxns),
        "n_total_audited_reactions": len(pool_agg["pool_audit_rows"]),
        "n_reaction_conflicts_excluded": len(pool_agg["reaction_conflicts"]),
        "n_species": len(canonical_species),
        "excluded_base_ids": sorted(excluded_base_ids),
        "reactions": eligible_rxns,
        "species": canonical_species,
    }
    with open(args.out_dir / "exp1_candidate_reaction_pool.pkl", "wb") as fh:
        pickle.dump(pool_payload, fh, protocol=pickle.HIGHEST_PROTOCOL)

    # 6. Select Pilot24 EQ MAGs (highest-quality, lineage-stratified across all major EQ orders/families)
    #    plus 8 HQ external comparator reference genomes for the artificial genome-incompleteness masking control
    eq_eligible = dedup_484_df[
        (dedup_484_df["cohort"] == "EQ_desert")
        & (dedup_484_df["has_draft_sbml"] == True)
        & (dedup_484_df["qc_error"] == "")
        & (dedup_484_df["draft_has_spurious_atp_cycle"] == False)
    ].copy()
    eq_eligible["quality_score"] = eq_eligible["completeness"] - 5.0 * eq_eligible["contamination"]
    eq_eligible = eq_eligible.sort_values(["quality_score", "contig_n50", "genome"], ascending=[False, False, True]).reset_index(drop=True)

    selected_pilot24: list[str] = []
    seen_p24: set[str] = set()
    families = eq_eligible["canonical_family"].unique().tolist()
    while len(selected_pilot24) < 24:
        added_in_round = 0
        for fam in families:
            fam_rows = eq_eligible[(eq_eligible["canonical_family"] == fam) & (~eq_eligible["genome"].isin(seen_p24))]
            if not fam_rows.empty:
                g_pick = str(fam_rows.iloc[0]["genome"])
                selected_pilot24.append(g_pick)
                seen_p24.add(g_pick)
                added_in_round += 1
                if len(selected_pilot24) >= 24:
                    break
        if added_in_round == 0:
            break

    pilot24_df = eq_eligible[eq_eligible["genome"].isin(seen_p24)].copy()
    pilot24_df["pilot_role"] = "EQ_pilot24_primary"

    ext_eligible = dedup_484_df[
        (dedup_484_df["cohort"] == "external_comparator")
        & (dedup_484_df["has_draft_sbml"] == True)
        & (dedup_484_df["qc_error"] == "")
        & (dedup_484_df["draft_has_spurious_atp_cycle"] == False)
        & (dedup_484_df["completeness"].notna())
        & (dedup_484_df["completeness"] >= 90.0)
        & (dedup_484_df["contamination"] <= 5.0)
    ].copy()
    ext_eligible["quality_score"] = ext_eligible["completeness"] - 5.0 * ext_eligible["contamination"]
    ext_eligible = ext_eligible.sort_values(["quality_score", "genome"], ascending=[False, True]).head(8).copy()
    ext_eligible["pilot_role"] = "external_HQ_masking_reference"

    manifest_df = pd.concat([pilot24_df, ext_eligible], ignore_index=True)
    manifest_df.to_csv(args.out_dir / "exp1_pilot24_manifest.tsv", sep="\t", index=False)

    dt = time.time() - t0
    print(
        f"[Stage 1 Complete in {dt:.1f}s] 538 archives inspected, 484 deduplicated genomes verified against frozen manifest, "
        f"{len(eligible_rxns)} eligible candidate reactions ({len(pool_agg['reaction_conflicts'])} stoich conflicts excluded), "
        f"{len(pilot24_df)} Pilot24 EQ MAGs + {len(ext_eligible)} HQ masking reference genomes selected."
    )


if __name__ == "__main__":
    main()
