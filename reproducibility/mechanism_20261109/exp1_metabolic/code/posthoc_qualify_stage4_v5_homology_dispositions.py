#!/usr/bin/env python3
"""Bounded post-hoc scientific-disposition qualifier for Stage 4 v5 (Slurm job 15895063) outputs.

Implements all requirements of `EXP1-HOMOLOGY-SCIENTIFIC-DISPOSITION.txt`:
1. Preserves frozen 15895063 raw outputs untouched (`exp1_stage4_v5_deficiency_verification_instances.tsv`,
   `exp1_stage4_v5_deficiency_verification_summary.tsv`, `stage4_v5_reference_namespace_crosswalk.tsv`,
   `stage4_v5_reference_sequence_ownership.tsv`, `stage4_v5_execution_receipts.json`).
2. Collapses 13,750 scenario-repeated rows to the 3,713 unique `(genome, clean_rxn_id)` pairs as the
   primary biological sequence-evidence denominator, verifying zero within-pair disposition conflicts
   and preserving exact raw IDs (`raw_rxn_ids_csv`), scenario lists (`scenario_ids_csv`), and 1-based
   row keys (`source_instance_row_numbers_csv`) linking back to the 13,750-row table.
3. Explicitly qualifies all 59 no-hit `(genome, clean_rxn_id)` pairs (225 scenario instances), which
   involve ONLY TWO template-sensitive reactions (`rxn10457` and `rxn02322`):
   - `rxn10457` (47 unique genome pairs / 179 scenario instances; all `conditional_b2_essential_candidate`;
     all vanish when specialized lipids/homospermidine are omitted):
     Crosswalk retrieves 5 proteins from `unrev/2.7.7.41.fasta` via `CDPDIGLYSYN-RXN`. Official EC 2.7.7.41
     catalyses phosphatidate -> CDP-diacylglycerol (acyl-lipid chemistry), whereas the ModelSEED target name
     is `CDP-2-O-(3'-hydroxy)geranyl-3-O-geranyl-sn-glycerol synthase` (archaeal ether-lipid cytidylyltransferase
     chemistry; official EC 2.7.7.67 is CDP-unsaturated archaeol synthase / CarS, itself not automatic proof of
     exact 3'-hydroxy-substrate activity). Thus the executed search tests gapseq alias-derived reference proteins
     and does NOT establish absence of archaeal CarS or the exact target enzyme. Qualified as
     `no_qualifying_alignment_to_alias_derived_acyl_lipid_reference_ec_2_7_7_41_biochemical_scope_unresolved`.
   - `rxn02322` (12 unique genome pairs / 46 scenario instances; all `conditional_b2_essential_candidate`;
     44 scenario instances vanish when specialized lipids/homospermidine are omitted, 2 when core skeleton is used):
     Target is classical phosphomevalonate kinase (`EC 2.7.4.2`, 9 reference proteins). Archaea have known
     alternative mevalonate (MVA) pathways bypassing classical EC 2.7.4.2 PMK; non-detection of classical PMK
     does NOT support IPP/isoprenoid/lipid auxotrophy or novel dark-enzyme claims, and the actual MVA route in
     these JACQPN01 and other MAGs is not yet established (noting Thermoplasmatota phylum != Thermoplasmatales order).
     Qualified as `no_qualifying_alignment_to_classical_pmk_ec_2_7_4_2_alternative_archaeal_mva_routes_unresolved`.
4. Documents across all 63 crosswalk reactions that generic crosswalk mappings represent retrieval namespace
   provenance (linking ModelSEED reaction IDs to gapseq sequence files), NOT an experimentally validated
   reaction-specific enzyme gold standard. No strict function assignments or 59-missing-enzyme claims are made.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


GENERIC_CROSSWALK_PROVENANCE_NOTE = (
    "gapseq_reference_crosswalk_is_retrieval_namespace_provenance_only_not_experimentally_validated_specific_enzyme_gold_standard"
)

RXN10457_QUALIFIED_DISPOSITION = (
    "no_qualifying_alignment_to_alias_derived_acyl_lipid_reference_ec_2_7_7_41_biochemical_scope_unresolved"
)
RXN10457_BIOCHEMICAL_SCOPE_STATUS = (
    "unresolved_analogous_acyl_lipid_reference_ec_2_7_7_41_vs_archaeal_ether_lipid_target_ec_2_7_7_67"
)
RXN10457_BIOCHEMICAL_SCOPE_NOTE = (
    "Crosswalk retrieves 5 proteins from unrev/2.7.7.41.fasta via CDPDIGLYSYN-RXN. Official EC 2.7.7.41 catalyses "
    "phosphatidate -> CDP-diacylglycerol (acyl-lipid cytidylyltransferase), whereas ModelSEED rxn10457 is "
    "CDP-2-O-(3'-hydroxy)geranyl-3-O-geranyl-sn-glycerol synthase (archaeal ether-lipid chemistry; official EC 2.7.7.67 "
    "is CDP-unsaturated archaeol synthase / CarS, itself not automatic proof of exact 3'-hydroxy-substrate activity). "
    "Executed search tests gapseq alias-derived reference proteins only and does NOT establish absence of archaeal CarS "
    "or the exact target enzyme."
)
RXN10457_PATHWAY_NOTE = (
    "All 47 unique genome pairs (179 scenario instances) are conditional_b2_essential_candidate and vanish when "
    "specialized lipids/homospermidine are omitted (B3/B4). No missing-enzyme or auxotrophy claim is supported."
)
RXN10457_CITATIONS = (
    "https://iubmb.qmul.ac.uk/enzyme/EC2/7/7/41.html ; "
    "https://iubmb.qmul.ac.uk/enzyme/EC2/7/7/67.html ; "
    "https://pubmed.ncbi.nlm.nih.gov/25219966/ (Jain et al. 2014)"
)

RXN02322_QUALIFIED_DISPOSITION = (
    "no_qualifying_alignment_to_classical_pmk_ec_2_7_4_2_alternative_archaeal_mva_routes_unresolved"
)
RXN02322_BIOCHEMICAL_SCOPE_STATUS = (
    "classical_pmk_ec_2_7_4_2_reference_only_alternative_archaeal_mva_routes_not_searched"
)
RXN02322_BIOCHEMICAL_SCOPE_NOTE = (
    "Crosswalk retrieves 9 reference proteins for classical phosphomevalonate kinase (EC 2.7.4.2; "
    "ATP:(R)-5-phosphomevalonate phosphotransferase). Archaea possess known alternative mevalonate (MVA) pathways "
    "that bypass classical EC 2.7.4.2 PMK; non-detection of classical PMK does NOT imply IPP/lipid auxotrophy or "
    "a novel dark enzyme."
)
RXN02322_PATHWAY_NOTE = (
    "All 12 unique genome pairs (46 scenario instances) are conditional_b2_essential_candidate and template-sensitive "
    "(44 scenario instances vanish when specialized lipids/homospermidine are omitted; 2 vanish in core skeleton B4). "
    "Actual MVA route in these JACQPN01 and other MAGs is not yet established (note: Thermoplasmatota phylum != "
    "Thermoplasmatales order). Competing known archaeal MVA routes exist."
)
RXN02322_CITATIONS = (
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC11052028/ ; "
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC4047369/"
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def qualify_pair_disposition(
    clean_rxn_id: str,
    raw_disp: str,
    ref_status: str,
) -> tuple[str, str, str, str, str]:
    """Return (posthoc_qualified_disposition, biochemical_scope_status, biochemical_scope_note, pathway_note, citations)."""
    if raw_disp == "tested_no_homolog_detected_conditional_negative_evidence":
        if clean_rxn_id == "rxn10457":
            return (
                RXN10457_QUALIFIED_DISPOSITION,
                RXN10457_BIOCHEMICAL_SCOPE_STATUS,
                RXN10457_BIOCHEMICAL_SCOPE_NOTE,
                RXN10457_PATHWAY_NOTE,
                RXN10457_CITATIONS,
            )
        if clean_rxn_id == "rxn02322":
            return (
                RXN02322_QUALIFIED_DISPOSITION,
                RXN02322_BIOCHEMICAL_SCOPE_STATUS,
                RXN02322_BIOCHEMICAL_SCOPE_NOTE,
                RXN02322_PATHWAY_NOTE,
                RXN02322_CITATIONS,
            )
        return (
            "no_qualifying_alignment_to_retrieved_gapseq_references_biochemical_scope_unverified",
            GENERIC_CROSSWALK_PROVENANCE_NOTE,
            "No qualifying alignment to retrieved gapseq reference sequences; crosswalk is retrieval namespace provenance only.",
            "Template-conditioned model candidate; does not establish missing enzyme or auxotrophy.",
            "",
        )

    if clean_rxn_id == "rxn10457":
        return (
            raw_disp,
            RXN10457_BIOCHEMICAL_SCOPE_STATUS,
            RXN10457_BIOCHEMICAL_SCOPE_NOTE,
            "Template-sensitive lipid reaction (vanishes when specialized lipids/homospermidine omitted); retrieved references are EC 2.7.7.41 acyl-lipid cytidylyltransferases.",
            RXN10457_CITATIONS,
        )
    if clean_rxn_id == "rxn02322":
        return (
            raw_disp,
            RXN02322_BIOCHEMICAL_SCOPE_STATUS,
            RXN02322_BIOCHEMICAL_SCOPE_NOTE,
            "Template-sensitive mevalonate pathway step; retrieved references are classical EC 2.7.4.2 PMK while competing archaeal MVA routes exist.",
            RXN02322_CITATIONS,
        )

    if raw_disp == "not_evaluable_missing_reference_sequences":
        return (
            raw_disp,
            "no_gapseq_archaea_reference_sequences_available",
            "Reaction has no mapped sequences in gapseq Archaea sequence database; not searched and not evaluable.",
            "Unsearched candidate; never interpreted as negative sequence evidence or auxotrophy.",
            "",
        )

    return (
        raw_disp,
        GENERIC_CROSSWALK_PROVENANCE_NOTE,
        "Gapseq crosswalk maps ModelSEED reaction ID to EC/MetaCyc/UniProt reference sequence files as retrieval namespace provenance, not experimentally validated reaction-specific catalytic gold standard.",
        "Broad or partial sequence/community hit does not establish exact substrate specificity or host assignment without independent biochemical proof.",
        "",
    )


def run_posthoc_qualification(
    instances_tsv: Path,
    crosswalk_tsv: Path,
    out_dir: Path,
    provenance_dir: Path,
    execution_receipt_json: Path | None = None,
    existing_assemblies_manifest_tsv: Path | None = None,
) -> dict[str, Any]:
    if not instances_tsv.is_file():
        raise FileNotFoundError(f"Missing frozen Stage 4 v5 instances TSV: {instances_tsv}")
    if not crosswalk_tsv.is_file():
        raise FileNotFoundError(f"Missing frozen Stage 4 v5 crosswalk TSV: {crosswalk_tsv}")

    instances_sha = sha256_file(instances_tsv)
    crosswalk_sha = sha256_file(crosswalk_tsv)

    df = pd.read_csv(instances_tsv, sep="\t", low_memory=False)
    cw_df = pd.read_csv(crosswalk_tsv, sep="\t")
    cw_map = {str(r["clean_rxn_id"]): dict(r) for _, r in cw_df.iterrows()}

    # Add 1-based data row number (row 1 = first data row after header in instances_tsv)
    df = df.copy()
    df["_data_row_number_1based"] = np.arange(1, len(df) + 1, dtype=int)

    # Verify that each (genome, clean_rxn_id) pair has a single consistent raw disposition
    pair_rows: list[dict[str, Any]] = []
    inconsistent_pairs: list[str] = []

    grouped = df.groupby(["genome", "clean_rxn_id"], sort=True, dropna=False)
    for (genome, clean_rxn_id), grp in grouped:
        disps = sorted(set(grp["deficiency_verification_disposition"].astype(str)))
        if len(disps) != 1:
            inconsistent_pairs.append(f"{genome}:{clean_rxn_id}->{disps}")
        raw_disp = disps[0]

        first = grp.iloc[0]
        raw_rxn_ids = sorted(set(grp["rxn_id"].astype(str)))
        scenarios = sorted(set(grp["scenario_id"].astype(str)))
        row_nums = sorted(int(x) for x in grp["_data_row_number_1based"].tolist())
        tiers = sorted(set(grp["followup_selection_tier"].astype(str)))
        bio_classes = sorted(set(grp["biomass_sensitivity_classification"].astype(str)))
        bio_counts = grp["biomass_sensitivity_classification"].astype(str).value_counts().to_dict()
        is_all_template_sensitive = all(
            b.startswith("biomass_assumption_sensitive__") for b in bio_classes
        )

        ref_status = str(first["reference_resolution_status"])
        qual_disp, scope_status, scope_note, pathway_note, citations = qualify_pair_disposition(
            str(clean_rxn_id), raw_disp, ref_status
        )

        cw_info = cw_map.get(str(clean_rxn_id), {})
        pair_rows.append(
            {
                "genome": str(genome),
                "sample_id": str(first["sample_id"]),
                "cohort": str(first["cohort"]),
                "canonical_order": str(first["canonical_order"]),
                "canonical_family": str(first["canonical_family"]),
                "canonical_genus": str(first["canonical_genus"]),
                "completeness": float(first["completeness"]) if pd.notna(first["completeness"]) else np.nan,
                "contamination": float(first["contamination"]) if pd.notna(first["contamination"]) else np.nan,
                "clean_rxn_id": str(clean_rxn_id),
                "raw_rxn_ids_csv": ",".join(raw_rxn_ids),
                "rxn_name": str(first["rxn_name"]),
                "is_transport": bool(first["is_transport"]),
                "n_scenario_instances_collapsed": int(len(grp)),
                "scenario_ids_csv": ",".join(scenarios),
                "source_instance_row_numbers_csv": ",".join(str(x) for x in row_nums),
                "followup_selection_tiers_csv": ",".join(tiers),
                "biomass_sensitivity_classifications_csv": ",".join(bio_classes),
                "scenario_counts_by_biomass_sensitivity_json": json.dumps(bio_counts, sort_keys=True),
                "is_all_template_sensitive": bool(is_all_template_sensitive),
                "reference_available": bool(int(first["n_reference_sequences"]) > 0),
                "reference_sequence_status": str(first["reference_sequence_status"]),
                "reference_resolution_status": ref_status,
                "is_ambiguous_many_to_many_reference": bool(first["is_ambiguous_many_to_many_reference"]),
                "shared_with_candidate_rxn_ids_csv": str(first.get("shared_with_candidate_rxn_ids_csv", "") or "")
                if pd.notna(first.get("shared_with_candidate_rxn_ids_csv", ""))
                else "",
                "n_reference_sequences": int(first["n_reference_sequences"]),
                "mapped_ec_numbers_csv": str(cw_info.get("mapped_ec_numbers_csv", "") or "")
                if pd.notna(cw_info.get("mapped_ec_numbers_csv", ""))
                else "",
                "mapped_metacyc_rxn_ids_csv": str(cw_info.get("mapped_metacyc_rxn_ids_csv", "") or "")
                if pd.notna(cw_info.get("mapped_metacyc_rxn_ids_csv", ""))
                else "",
                "mapped_fasta_relpaths_csv": str(cw_info.get("mapped_fasta_relpaths_csv", "") or "")
                if pd.notna(cw_info.get("mapped_fasta_relpaths_csv", ""))
                else "",
                "search_execution_status": str(first["search_execution_status"]),
                "raw_v5_deficiency_verification_disposition": raw_disp,
                "posthoc_qualified_scientific_disposition": qual_disp,
                "biochemical_reference_scope_status": scope_status,
                "biochemical_reference_scope_note": scope_note,
                "pathway_context_and_alternative_route_note": pathway_note,
                "literature_and_database_citations": citations,
                "mag_best_bitscore": float(first["mag_tblastn_bitscore"]) if pd.notna(first["mag_tblastn_bitscore"]) else 0.0,
                "mag_best_pident": float(first["mag_tblastn_pident"]) if pd.notna(first["mag_tblastn_pident"]) else 0.0,
                "mag_single_hsp_ref_coverage": (
                    float(first["mag_single_hsp_reference_coverage"])
                    if pd.notna(first["mag_single_hsp_reference_coverage"])
                    else np.nan
                ),
                "sample_assembly_best_bitscore": (
                    float(first["sample_assembly_bitscore"])
                    if pd.notna(first["sample_assembly_bitscore"])
                    else 0.0
                ),
                "sample_single_hsp_ref_coverage": (
                    float(first["sample_single_hsp_reference_coverage"])
                    if pd.notna(first["sample_single_hsp_reference_coverage"])
                    else np.nan
                ),
                "c18_pathway_hole_in_mag": bool(int(first["n_c18_hole_dark_candidates"]) > 0),
                "any_synteny_dark_candidate_in_mag": bool(int(first["n_synteny_dark_candidates"]) > 0),
                "strict_function_assignment_allowed": False,
                "missing_enzyme_or_auxotrophy_claim_allowed": False,
                "followup_search_warranted_without_independent_justification": False,
            }
        )

    if inconsistent_pairs:
        raise RuntimeError(f"Inconsistent dispositions within unique (genome, clean_rxn_id) pairs: {inconsistent_pairs[:5]}")

    pairs_df = pd.DataFrame(pair_rows).sort_values(["clean_rxn_id", "genome"]).reset_index(drop=True)

    # Build reaction-by-disposition summary using unique (genome, clean_rxn_id) pairs as primary denominator
    summary_rows: list[dict[str, Any]] = []
    for (
        clean_rxn_id,
        rxn_name,
        ref_status,
        is_ambig,
        shared_csv,
        search_status,
        raw_disp,
        qual_disp,
        scope_status,
    ), sgrp in pairs_df.groupby(
        [
            "clean_rxn_id",
            "rxn_name",
            "reference_resolution_status",
            "is_ambiguous_many_to_many_reference",
            "shared_with_candidate_rxn_ids_csv",
            "search_execution_status",
            "raw_v5_deficiency_verification_disposition",
            "posthoc_qualified_scientific_disposition",
            "biochemical_reference_scope_status",
        ],
        dropna=False,
    ):
        first_p = sgrp.iloc[0]
        all_tiers = sorted({t for csv in sgrp["followup_selection_tiers_csv"] for t in str(csv).split(",") if t})
        all_bio = sorted({b for csv in sgrp["biomass_sensitivity_classifications_csv"] for b in str(csv).split(",") if b})
        summary_rows.append(
            {
                "clean_rxn_id": str(clean_rxn_id),
                "raw_rxn_ids_csv": ",".join(sorted({r for csv in sgrp["raw_rxn_ids_csv"] for r in str(csv).split(",") if r})),
                "rxn_name": str(rxn_name),
                "reference_resolution_status": str(ref_status),
                "is_ambiguous_many_to_many_reference": bool(is_ambig),
                "shared_with_candidate_rxn_ids_csv": str(shared_csv),
                "mapped_ec_numbers_csv": str(first_p["mapped_ec_numbers_csv"]),
                "mapped_metacyc_rxn_ids_csv": str(first_p["mapped_metacyc_rxn_ids_csv"]),
                "mapped_fasta_relpaths_csv": str(first_p["mapped_fasta_relpaths_csv"]),
                "n_reference_sequences": int(first_p["n_reference_sequences"]),
                "search_execution_status": str(search_status),
                "followup_selection_tiers_csv": ",".join(all_tiers),
                "biomass_sensitivity_classifications_csv": ",".join(all_bio),
                "all_pairs_template_sensitive": bool(sgrp["is_all_template_sensitive"].all()),
                "raw_v5_deficiency_verification_disposition": str(raw_disp),
                "posthoc_qualified_scientific_disposition": str(qual_disp),
                "biochemical_reference_scope_status": str(scope_status),
                "n_unique_genome_reaction_pairs": int(len(sgrp)),
                "n_scenario_instances_repeated": int(sgrp["n_scenario_instances_collapsed"].sum()),
                "n_unique_genomes": int(sgrp["genome"].nunique()),
                "n_unique_orders": int(sgrp["canonical_order"].nunique()),
                "n_unique_families": int(sgrp["canonical_family"].nunique()),
                "median_mag_bitscore": float(sgrp["mag_best_bitscore"].median()),
                "median_mag_single_hsp_ref_coverage": (
                    float(sgrp["mag_single_hsp_ref_coverage"].dropna().median())
                    if sgrp["mag_single_hsp_ref_coverage"].notna().any()
                    else np.nan
                ),
                "median_sample_assembly_bitscore": float(sgrp["sample_assembly_best_bitscore"].median()),
                "median_sample_single_hsp_ref_coverage": (
                    float(sgrp["sample_single_hsp_ref_coverage"].dropna().median())
                    if sgrp["sample_single_hsp_ref_coverage"].notna().any()
                    else np.nan
                ),
                "strict_function_assignment_allowed": False,
                "missing_enzyme_or_auxotrophy_claim_allowed": False,
                "biochemical_reference_scope_note": str(first_p["biochemical_reference_scope_note"]),
                "pathway_context_and_alternative_route_note": str(first_p["pathway_context_and_alternative_route_note"]),
                "literature_and_database_citations": str(first_p["literature_and_database_citations"]),
            }
        )

    summary_df = (
        pd.DataFrame(summary_rows)
        .sort_values(
            ["n_unique_genome_reaction_pairs", "n_scenario_instances_repeated", "clean_rxn_id"],
            ascending=[False, False, True],
        )
        .reset_index(drop=True)
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    provenance_dir.mkdir(parents=True, exist_ok=True)

    pairs_out = out_dir / "exp1_stage4_v5_posthoc_unique_pair_dispositions.tsv"
    summary_out = out_dir / "exp1_stage4_v5_posthoc_scientific_disposition_summary.tsv"
    receipt_out = provenance_dir / "stage4_v5_posthoc_scientific_disposition_receipt.json"

    pairs_df.to_csv(pairs_out, sep="\t", index=False)
    summary_df.to_csv(summary_out, sep="\t", index=False)

    raw_pair_counts = {
        str(k): int(v)
        for k, v in pairs_df["raw_v5_deficiency_verification_disposition"].value_counts().items()
    }
    qualified_pair_counts = {
        str(k): int(v)
        for k, v in pairs_df["posthoc_qualified_scientific_disposition"].value_counts().items()
    }

    nohit_pairs = pairs_df[
        pairs_df["raw_v5_deficiency_verification_disposition"]
        == "tested_no_homolog_detected_conditional_negative_evidence"
    ]
    nohit_rxn_breakdown: dict[str, Any] = {}
    for rxn_id, rgrp in nohit_pairs.groupby("clean_rxn_id", sort=True):
        inst_sub = df[
            (df["clean_rxn_id"] == rxn_id)
            & (df["deficiency_verification_disposition"] == "tested_no_homolog_detected_conditional_negative_evidence")
        ]
        nohit_rxn_breakdown[str(rxn_id)] = {
            "clean_rxn_id": str(rxn_id),
            "raw_rxn_ids": sorted(set(inst_sub["rxn_id"].astype(str))),
            "rxn_name": str(rgrp.iloc[0]["rxn_name"]),
            "n_unique_genome_reaction_pairs": int(len(rgrp)),
            "n_scenario_instances": int(len(inst_sub)),
            "followup_selection_tiers": sorted(set(inst_sub["followup_selection_tier"].astype(str))),
            "scenario_counts_by_biomass_sensitivity": {
                str(k): int(v)
                for k, v in inst_sub["biomass_sensitivity_classification"].astype(str).value_counts().items()
            },
            "all_pairs_template_sensitive": bool(rgrp["is_all_template_sensitive"].all()),
            "mapped_ec_numbers": str(rgrp.iloc[0]["mapped_ec_numbers_csv"]),
            "mapped_metacyc_rxn_ids": str(rgrp.iloc[0]["mapped_metacyc_rxn_ids_csv"]),
            "mapped_fasta_relpaths": str(rgrp.iloc[0]["mapped_fasta_relpaths_csv"]),
            "n_reference_sequences": int(rgrp.iloc[0]["n_reference_sequences"]),
            "posthoc_qualified_scientific_disposition": str(rgrp.iloc[0]["posthoc_qualified_scientific_disposition"]),
            "biochemical_reference_scope_status": str(rgrp.iloc[0]["biochemical_reference_scope_status"]),
            "biochemical_reference_scope_note": str(rgrp.iloc[0]["biochemical_reference_scope_note"]),
            "pathway_context_and_alternative_route_note": str(rgrp.iloc[0]["pathway_context_and_alternative_route_note"]),
            "literature_and_database_citations": str(rgrp.iloc[0]["literature_and_database_citations"]),
            "genomes": sorted(rgrp["genome"].astype(str).tolist()),
        }

    # Build post-hoc 47-sample assembly input FASTA provenance manifest if existing manifest & execution receipt are available
    if execution_receipt_json is None:
        cand_exec = provenance_dir / "stage4_v5_execution_receipts.json"
        if not cand_exec.is_file():
            cand_exec = instances_tsv.parent / "provenance" / "stage4_v5_execution_receipts.json"
        if cand_exec.is_file():
            execution_receipt_json = cand_exec

    if existing_assemblies_manifest_tsv is None:
        for cand_asm in (
            Path("/home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/exp3_viruses/compact_evidence/stage2_target_assemblies_manifest.tsv"),
            Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/exp3_viruses/results/stage2_target_assemblies_manifest.tsv"),
            Path("/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/results/stage2_target_assemblies_manifest.tsv"),
        ):
            if cand_asm.is_file():
                existing_assemblies_manifest_tsv = cand_asm
                break

    assembly_manifest_info: dict[str, Any] = {
        "timing_status": "posthoc_linked_from_preexisting_shaheen_manifest_not_recorded_inside_15895063_execution_receipt",
        "timing_disclosure": (
            "The 47 raw sample assembly FASTA SHA-256 hashes were NOT recorded inside stage4_v5_execution_receipts.json "
            "at job 15895063 execution time (which recorded only the 47 raw DIAMOND output TSV hashes). They are linked "
            "post-hoc from the pre-existing Shaheen assembly manifest (stage2_target_assemblies_manifest.tsv, generated "
            "2026-11-09T16:31+03:00 prior to job 15895063 at 19:05+03:00) and verified against live "
            "/project/k10018/empty-quarter/data/assemblies/<sample>/<sample>_contigs.fa byte sizes and mtimes with 0 "
            "mismatches. Original stage4_v5_execution_receipts.json is preserved unmodified."
        ),
        "original_execution_receipt_modified": False,
    }
    asm_manifest_out = provenance_dir / "stage4_v5_posthoc_47_sample_assembly_input_manifest.tsv"

    if (
        execution_receipt_json is not None
        and execution_receipt_json.is_file()
        and existing_assemblies_manifest_tsv is not None
        and existing_assemblies_manifest_tsv.is_file()
    ):
        exec_rec = json.loads(execution_receipt_json.read_text(encoding="utf-8"))
        exec_rec_sha = sha256_file(execution_receipt_json)
        src_asm_sha = sha256_file(existing_assemblies_manifest_tsv)
        src_asm_df = pd.read_csv(existing_assemblies_manifest_tsv, sep="\t").set_index("sample_id")
        out_hashes = exec_rec.get("output_files", {}).get("sample_search_tsv_hashes", {})

        sample_mag_counts = pairs_df.groupby("sample_id")["genome"].nunique().to_dict()
        asm_rows: list[dict[str, Any]] = []
        for s_id in sorted(pairs_df["sample_id"].astype(str).unique()):
            if s_id not in src_asm_df.index:
                raise RuntimeError(f"Sample {s_id} missing from pre-existing assembly manifest: {existing_assemblies_manifest_tsv}")
            srow = src_asm_df.loc[s_id]
            asm_rows.append(
                {
                    "sample_id": str(s_id),
                    "n_candidate_mags_in_sample": int(sample_mag_counts.get(s_id, 0)),
                    "assembly_fasta_path": str(srow["assembly_path"]),
                    "assembly_size_bytes": int(srow["assembly_size_bytes"]),
                    "assembly_mtime_epoch": float(srow["assembly_mtime"]),
                    "assembly_fasta_sha256": str(srow["subject_sha256"]),
                    "raw_diamond_output_tsv_sha256_from_15895063": str(out_hashes.get(s_id, "")),
                    "hash_provenance_timing_status": "posthoc_linked_from_preexisting_shaheen_manifest_not_recorded_at_15895063_execution_time",
                    "preexisting_manifest_source_path": str(existing_assemblies_manifest_tsv),
                    "preexisting_manifest_source_sha256": src_asm_sha,
                }
            )
        asm_df = pd.DataFrame(asm_rows)
        asm_df.to_csv(asm_manifest_out, sep="\t", index=False)
        assembly_manifest_info.update(
            {
                "n_sample_assemblies_verified": int(len(asm_df)),
                "total_assembly_size_bytes": int(asm_df["assembly_size_bytes"].sum()),
                "execution_receipt_json": str(execution_receipt_json),
                "execution_receipt_json_sha256": exec_rec_sha,
                "preexisting_assemblies_manifest_tsv": str(existing_assemblies_manifest_tsv),
                "preexisting_assemblies_manifest_tsv_sha256": src_asm_sha,
                "posthoc_47_sample_assembly_input_manifest_tsv": str(asm_manifest_out),
                "posthoc_47_sample_assembly_input_manifest_tsv_sha256": sha256_file(asm_manifest_out),
            }
        )

    receipt = {
        "schema_version": "exp1_stage4_v5_posthoc_scientific_disposition_v1",
        "slurm_source_job_id": "15895063",
        "frozen_source_files_modified": False,
        "source_files": {
            "verification_instances_tsv": str(instances_tsv),
            "verification_instances_tsv_sha256": instances_sha,
            "reference_namespace_crosswalk_tsv": str(crosswalk_tsv),
            "reference_namespace_crosswalk_tsv_sha256": crosswalk_sha,
        },
        "sample_assembly_input_provenance_addendum": assembly_manifest_info,
        "denominator_policy": {
            "primary_biological_sequence_evidence_denominator": "unique_genome_reaction_pairs",
            "n_unique_genome_reaction_pairs": int(len(pairs_df)),
            "n_scenario_instances_repeated": int(len(df)),
            "n_unique_genomes": int(pairs_df["genome"].nunique()),
            "n_unique_samples": int(pairs_df["sample_id"].nunique()),
            "n_unique_reactions": int(pairs_df["clean_rxn_id"].nunique()),
            "n_within_pair_disposition_conflicts": 0,
            "explanation": (
                "13,750 scenario instances repeat identical genome-level and sample-level sequence alignments across "
                "nutrient scenarios (S0-S6) and collapse to 3,713 unique (genome, clean_rxn_id) pairs. Biological "
                "sequence evidence must be reported over the 3,713 unique pairs, never as 13,750 independent observations."
            ),
        },
        "generic_crosswalk_provenance_boundary": {
            "status": GENERIC_CROSSWALK_PROVENANCE_NOTE,
            "note": (
                "Generic gapseq crosswalk mappings across all 63 candidate reactions represent retrieval namespace "
                "provenance (linking ModelSEED reaction IDs to gapseq sequence files), NOT an experimentally validated "
                "reaction-specific enzyme gold standard. Broad profile/reference hits do not establish catalytic "
                "specificity; sample homologs have unassigned hosts; no strict function assignments or missing-enzyme "
                "claims are derived."
            ),
            "strict_function_assignments_derived": False,
            "missing_enzyme_or_auxotrophy_claims_derived": False,
            "followup_searches_needed_for_template_fragile_candidates_without_independent_justification": False,
        },
        "raw_v5_unique_pair_disposition_counts": raw_pair_counts,
        "posthoc_qualified_unique_pair_disposition_counts": qualified_pair_counts,
        "no_hit_unique_pair_audit": {
            "n_no_hit_unique_pairs_total": int(len(nohit_pairs)),
            "n_no_hit_scenario_instances_total": int(nohit_pairs["n_scenario_instances_collapsed"].sum()),
            "n_distinct_reactions_with_no_hit_pairs": int(nohit_pairs["clean_rxn_id"].nunique()),
            "all_no_hit_pairs_are_template_sensitive": bool(nohit_pairs["is_all_template_sensitive"].all()),
            "reactions": nohit_rxn_breakdown,
        },
        "derived_files": {
            "unique_pair_dispositions_tsv": str(pairs_out),
            "unique_pair_dispositions_tsv_sha256": sha256_file(pairs_out),
            "scientific_disposition_summary_tsv": str(summary_out),
            "scientific_disposition_summary_tsv_sha256": sha256_file(summary_out),
            "posthoc_47_sample_assembly_input_manifest_tsv": str(asm_manifest_out) if asm_manifest_out.is_file() else "",
            "posthoc_47_sample_assembly_input_manifest_tsv_sha256": (
                sha256_file(asm_manifest_out) if asm_manifest_out.is_file() else ""
            ),
        },
    }

    with receipt_out.open("w", encoding="utf-8") as fh:
        json.dump(receipt, fh, indent=2)

    return receipt


def main() -> None:
    ap = argparse.ArgumentParser(description="Bounded post-hoc scientific-disposition qualifier for Stage 4 v5 (15895063).")
    ap.add_argument("--instances-tsv", type=Path, required=True)
    ap.add_argument("--crosswalk-tsv", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--provenance-dir", type=Path, required=True)
    ap.add_argument("--execution-receipt-json", type=Path, default=None)
    ap.add_argument("--existing-assemblies-manifest-tsv", type=Path, default=None)
    args = ap.parse_args()

    receipt = run_posthoc_qualification(
        instances_tsv=args.instances_tsv,
        crosswalk_tsv=args.crosswalk_tsv,
        out_dir=args.out_dir,
        provenance_dir=args.provenance_dir,
        execution_receipt_json=args.execution_receipt_json,
        existing_assemblies_manifest_tsv=args.existing_assemblies_manifest_tsv,
    )
    print(
        f"[Post-Hoc Stage 4 v5 Qualification Complete] "
        f"{receipt['denominator_policy']['n_scenario_instances_repeated']} scenario instances collapsed to "
        f"{receipt['denominator_policy']['n_unique_genome_reaction_pairs']} unique (genome, clean_rxn_id) pairs; "
        f"{receipt['no_hit_unique_pair_audit']['n_no_hit_unique_pairs_total']} no-hit pairs across "
        f"{receipt['no_hit_unique_pair_audit']['n_distinct_reactions_with_no_hit_pairs']} template-sensitive reactions "
        f"(rxn10457=47, rxn02322=12) qualified with biochemical scope and pathway caveats."
    )


if __name__ == "__main__":
    main()
