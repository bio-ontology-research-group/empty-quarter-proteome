#!/usr/bin/env python3
"""Stage 3: Full 591-Archaeal-MAG Host-Link Evidence Ledger and Historical C1 Audit.

Produces:
1. `results/stage3_all591_archaeal_mag_host_evidence.tsv`:
   Complete 591-row ledger (`237` Thermoplasmatota + `354` other archaea) with
   mutually exclusive `mag_host_disposition` separating zero detected spacers,
   spacers with no viral hit, unverified viral identity (`CheckV Not-determined`
   + no hallmark), unexcised provirus candidates, ambiguous multi-host links,
   technical self-CRISPR-array masks, external habitat viral control matches,
   curated non-viral decoys, failed caller execution, and verified Tier 1 / Tier 2 links.
2. `results/stage3_archaeal_virus_host_pairs_all_tiers.tsv`:
   All candidate archaeal host-virus pairs across all tiers with explicit Axis 1
   (CRISPR/junction) and Axis 2 (viral identity & habitat specificity) columns.
3. `results/stage3_historical_c1_thermoplasmatota_audit.tsv`:
   Row-by-row forensic audit of the 8 historical Thermoplasmatota virus-host
   edges from `edges.tsv` (`host_mag`, `votu_rep`, `host_phylum`, `n_spacers`).
4. `results/stage3_summary.json`.
"""
from __future__ import annotations

import collections
import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exp3_common import (
    classify_host_link_evidence_tier,
    reconcile_host_alias_strict,
    summarize_mag_host_disposition,
)

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
IN_DIR = BRANCH_ROOT / "inputs"
RES_DIR = BRANCH_ROOT / "results"

P_MOB = Path("/project/k10018/empty-quarter/data/mobilome")
P_VOTU_MEM = P_MOB / "catalogues/votu_vclust_1.3.1_ani95_af85/votu_membership.tsv"
P_E2_PROVIR = P_MOB / "results_20261104/curated/e2/provirus_host_links.tsv"
P_SKANI_SP = P_MOB / "results_20261104/curated/e2/mag2sp/strict_mag_to_species.tsv"


def main() -> None:
    t0 = time.time()
    RES_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load canonical 591 archaeal MAGs
    tax_df = pd.read_csv(IN_DIR / "taxonomy_reconciliation_all5438.tsv.gz", sep="\t")
    rep_5438_set = set(tax_df["genome"])
    arch_df = tax_df[tax_df["canonical_domain"] == "d__Archaea"].copy().sort_values("genome")
    arch591_list = list(arch_df["genome"])
    arch591_set = set(arch591_list)
    tax_lookup = tax_df.set_index("genome").to_dict("index")

    # 2. Load non-representative refined bin -> species representative crosswalk
    nonrep_to_species: dict[str, str] = {}
    hcw_path = IN_DIR / "host_crosswalk.tsv"
    if hcw_path.exists():
        for r in csv.DictReader(open(hcw_path), delimiter="\t"):
            if r["host_mag"] not in rep_5438_set and r["species_rep"] in arch591_set:
                nonrep_to_species[r["host_mag"]] = r["species_rep"]

    if P_SKANI_SP.exists():
        with open(P_SKANI_SP) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0] not in rep_5438_set and parts[1] in arch591_set:
                    nonrep_to_species.setdefault(parts[0], parts[1])

    # 3. Load vOTU membership mapping (`contig_id -> votu_id`)
    contig_to_votu: dict[str, str] = {}
    if P_VOTU_MEM.exists():
        with open(P_VOTU_MEM) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    contig_to_votu[parts[0]] = parts[1]
                    contig_to_votu[parts[1]] = parts[1]

    # 4. Load Stage 0, Stage 1, and Stage 2 outputs
    ver_df = pd.read_csv(RES_DIR / "stage0_arch591_fasta_verification.tsv", sep="\t")
    ver_lookup = ver_df.set_index("genome").to_dict("index")

    dup_df = pd.read_csv(RES_DIR / "stage0_arch591_sample_duplicate_contig_audit.tsv", sep="\t")
    dup_lookup = dup_df.set_index("genome").to_dict("index")

    st1_df = pd.read_csv(RES_DIR / "stage1_mag_screening_status.tsv", sep="\t")
    st1_lookup = st1_df.set_index("mag_id").to_dict("index")

    scored_path = RES_DIR / "stage2_all_spacer_alignments_scored.tsv"
    scored_df = (
        pd.read_csv(scored_path, sep="\t")
        if (scored_path.exists() and scored_path.stat().st_size > 10)
        else pd.DataFrame()
    )

    target_id_path = RES_DIR / "stage2_target_viral_identity_table.tsv"
    target_id_df = (
        pd.read_csv(target_id_path, sep="\t")
        if (target_id_path.exists() and target_id_path.stat().st_size > 10)
        else pd.DataFrame()
    )
    target_id_lookup = (
        target_id_df.set_index("target_id").to_dict("index")
        if len(target_id_df) > 0
        else {}
    )

    # 5. Group Stage 2 alignments by `(host_mag, vOTU_or_target)` and compute multi-host specificity
    if len(scored_df) > 0:
        scored_df["votu_or_target_id"] = scored_df["target_id"].map(
            lambda tid: contig_to_votu.get(str(tid), str(tid))
        )
    else:
        scored_df["votu_or_target_id"] = []

    virus_to_genera: dict[str, set[str]] = collections.defaultdict(set)
    virus_to_phyla: dict[str, set[str]] = collections.defaultdict(set)

    if len(scored_df) > 0:
        strict_unmasked = scored_df[
            (scored_df["passes_strict_spacer_gate"] == True)
            & (scored_df["is_technical_crispr_array_self_match"] == 0)
        ]
        for _, r in strict_unmasked.iterrows():
            h_mag = str(r["host_mag"])
            _, rep_mag = reconcile_host_alias_strict(h_mag, arch591_set, nonrep_to_species)
            if rep_mag and rep_mag in tax_lookup:
                v_key = str(r["votu_or_target_id"])
                virus_to_genera[v_key].add(str(tax_lookup[rep_mag]["canonical_genus"]))
                virus_to_phyla[v_key].add(str(tax_lookup[rep_mag]["canonical_phylum"]))

    # 6. Evaluate all `(host_mag, target_id)` pairs
    pair_rows: list[dict[str, Any]] = []
    mag_to_pair_tiers: dict[str, list[str]] = collections.defaultdict(list)
    mag_to_verified_viruses: dict[str, set[str]] = collections.defaultdict(set)

    if len(scored_df) > 0:
        for (h_mag, t_id), grp in scored_df.groupby(["host_mag", "target_id"], as_index=False):
            h_status, rep_mag = reconcile_host_alias_strict(str(h_mag), arch591_set, nonrep_to_species)
            v_key = contig_to_votu.get(str(t_id), str(t_id))

            unmasked_grp = grp[grp["is_technical_crispr_array_self_match"] == 0]
            if len(unmasked_grp) > 0:
                ref_row = unmasked_grp.iloc[0]
                n_strict = int(
                    unmasked_grp[unmasked_grp["passes_strict_spacer_gate"] == True]["unique_spacer_id"].nunique()
                )
                n_exact = int(
                    unmasked_grp[unmasked_grp["match_tier"] == "exact_0mm_full"]["unique_spacer_id"].nunique()
                )
                n_explor = int(
                    unmasked_grp[unmasked_grp["passes_exploratory_spacer_gate"] == True]["unique_spacer_id"].nunique()
                )
                min_eff_mm = unmasked_grp["effective_mismatches"].min()
                max_cov = unmasked_grp["full_spacer_coverage"].max()
                best_tier = (
                    "exact_0mm_full"
                    if n_exact > 0
                    else ("strict_1mm_end2" if n_strict > 0 else "near_2mm_exploratory")
                )
                is_tech_self = False
            else:
                ref_row = grp.iloc[0]
                n_strict = 0
                n_exact = 0
                n_explor = 0
                min_eff_mm = grp["effective_mismatches"].min()
                max_cov = grp["full_spacer_coverage"].max()
                best_tier = "technical_self_crispr_array_masked"
                is_tech_self = True

            v_eval_dict = {
                "viral_identity_status": str(ref_row["viral_identity_status"]),
                "habitat_specificity_role": str(ref_row["habitat_specificity_role"]),
                "is_verified_viral": bool(int(ref_row["is_verified_viral"])),
                "is_provisional_viral": bool(int(ref_row["is_provisional_viral"])),
                "is_unknown_viral_identity": bool(int(ref_row["is_unknown_viral_identity"])),
                "is_cellular_like": bool(int(ref_row["is_cellular_like"])),
                "is_technical_crispr_array_self_match": is_tech_self,
                "is_decoy": bool(int(ref_row["is_decoy"])),
                "is_external_viral_specificity_control": bool(int(ref_row["is_external_viral_specificity_control"])),
                "status_detail": str(ref_row["status_detail"]),
            }

            has_cas = bool(int(ref_row["host_has_cas_operon"]))
            # Only treat provirus junction as verified when `provirus_boundary_status == "verified_provirus_with_host_flanks"`
            # AND `is_same_namespaced_contig == 1` (never from filename coordinates alone!)
            has_prov_junc = (
                str(ref_row["provirus_boundary_status"]) == "verified_provirus_with_host_flanks"
                and bool(int(ref_row["is_same_namespaced_contig"]))
            )

            n_gen = max(1, len(virus_to_genera.get(v_key, set())))
            n_phy = max(1, len(virus_to_phyla.get(v_key, set())))

            tier_res = classify_host_link_evidence_tier(
                host_reconciliation_status=h_status,
                viral_eval=v_eval_dict,
                n_strict_spacers=n_strict,
                n_exact_0mm_spacers=n_exact,
                n_exploratory_spacers=n_explor,
                has_cas_operon_in_mag=has_cas,
                has_verified_provirus_junction_in_mag=has_prov_junc,
                has_classifier_corroboration=False,
                n_distinct_host_genera_for_virus=n_gen,
                n_distinct_host_phyla_for_virus=n_phy,
            )

            t_row_tax = tax_lookup.get(rep_mag, {}) if rep_mag else {}
            pair_entry = {
                "host_mag_raw": str(h_mag),
                "canonical_591_host_mag": rep_mag if rep_mag else "unmapped",
                "host_reconciliation_status": h_status,
                "host_phylum": t_row_tax.get("canonical_phylum", "NA"),
                "host_order": t_row_tax.get("canonical_order", "NA"),
                "host_genus": t_row_tax.get("canonical_genus", "NA"),
                "host_species": t_row_tax.get("canonical_species", "NA"),
                "target_id": str(t_id),
                "target_seq_sha256": str(ref_row.get("target_seq_sha256", "NA")),
                "annotation_provenance": str(ref_row.get("annotation_provenance", "NA")),
                "votu_id": v_key,
                "is_votu_cluster_rep": int(ref_row.get("is_votu_cluster_rep", 0)),
                "pool": str(ref_row["pool"]),
                "best_spacer_match_tier": best_tier,
                "n_exact_0mm_spacers": n_exact,
                "n_strict_spacers": n_strict,
                "n_exploratory_spacers": n_explor,
                "min_effective_mismatches": min_eff_mm,
                "max_full_spacer_coverage": max_cov,
                "host_has_cas_operon": int(has_cas),
                "has_verified_provirus_junction": int(has_prov_junc),
                "is_same_namespaced_contig": int(ref_row["is_same_namespaced_contig"]),
                "viral_identity_status": str(ref_row["viral_identity_status"]),
                "habitat_specificity_role": str(ref_row["habitat_specificity_role"]),
                "checkv_quality": str(ref_row["checkv_quality"]),
                "checkv_viral_genes": ref_row["checkv_viral_genes"],
                "checkv_host_genes": ref_row["checkv_host_genes"],
                "virus_hallmarks": ref_row["virus_hallmarks"],
                "archaeal_viral_hallmarks": ref_row["archaeal_viral_hallmarks"],
                "min_independent_viral_orfs_lower_bound": ref_row.get("min_independent_viral_orfs_lower_bound", "NA"),
                "orf_count_provenance": ref_row.get("orf_count_provenance", "NA"),
                "provirus_boundary_status": str(ref_row["provirus_boundary_status"]),
                "read_junction_support_status": str(ref_row.get("read_junction_support_status", "untested_missing")),
                "n_distinct_host_genera_for_virus": n_gen,
                "n_distinct_host_phyla_for_virus": n_phy,
                "evidence_tier": tier_res["evidence_tier"],
                "eligible_stage4_primary": int(tier_res["eligible_stage4_primary"]),
                "eligible_stage4_sensitivity": int(tier_res["eligible_stage4_sensitivity"]),
                "exclusion_reason": tier_res["exclusion_reason"],
            }
            pair_rows.append(pair_entry)

            if rep_mag and h_status == "canonical_591_direct":
                mag_to_pair_tiers[rep_mag].append(tier_res["evidence_tier"])
                if tier_res["eligible_stage4_sensitivity"]:
                    mag_to_verified_viruses[rep_mag].add(v_key)

    pairs_df = pd.DataFrame(pair_rows)
    if len(pairs_df) > 0:
        pairs_df = pairs_df.sort_values(["evidence_tier", "canonical_591_host_mag", "target_id"])
    pairs_df.to_csv(RES_DIR / "stage3_archaeal_virus_host_pairs_all_tiers.tsv", sep="\t", index=False)

    # 7. Build the complete 591-row archaeal MAG host-link evidence ledger
    ledger_rows: list[dict[str, Any]] = []
    for g in arch591_list:
        t_row = tax_lookup[g]
        v_row = ver_lookup.get(g, {})
        d_row = dup_lookup.get(g, {})
        s_row = st1_lookup.get(g, {})

        fasta_ok = (v_row.get("refined_verification_status") == "verified_exact_match_genome_manifest")
        scr_status = str(s_row.get("screening_status", ""))
        if scr_status:
            screened_ok = (scr_status == "screened_both_callers_ok")
            caller_failed = (scr_status == "failed_caller_execution")
        else:
            screened_ok = bool(s_row.get("minced_ok", 0) and s_row.get("cctyper_ok", 0))
            caller_failed = bool(s_row) and not screened_ok
            scr_status = "screened_both_callers_ok" if screened_ok else ("failed_caller_execution" if caller_failed else "unscreened")
        n_arr = int(s_row.get("n_raw_arrays", 0) or 0)
        n_spc = int(s_row.get("n_independent_spacers", s_row.get("n_valid_spacers", 0)) or 0)
        n_cas = int(s_row.get("n_cas_operons", 0) or 0)

        p_tiers = mag_to_pair_tiers.get(g, [])
        disp = summarize_mag_host_disposition(
            mag_id=g,
            fasta_verified=fasta_ok,
            crispr_screened=screened_ok,
            n_raw_arrays=n_arr,
            n_verified_spacers=n_spc,
            pair_tiers=p_tiers,
            caller_failed=caller_failed,
        )
        t_counts = collections.Counter(p_tiers)
        linked_v = sorted(mag_to_verified_viruses.get(g, set()))

        ledger_rows.append(
            {
                "genome": g,
                "source_sample": g.split("__")[0],
                "canonical_domain": t_row["canonical_domain"],
                "canonical_phylum": t_row["canonical_phylum"],
                "canonical_class": t_row["canonical_class"],
                "canonical_order": t_row["canonical_order"],
                "canonical_family": t_row["canonical_family"],
                "canonical_genus": t_row["canonical_genus"],
                "canonical_species": t_row["canonical_species"],
                "checkm2_completeness": t_row["canonical_completeness"],
                "checkm2_contamination": t_row["canonical_contamination"],
                "pass_gunc": v_row.get("pass_gunc", "unscored"),
                "fasta_verification_status": v_row.get("refined_verification_status", "missing"),
                "catalog_derep_gz_disposition": v_row.get("catalog_derep_gz_disposition", "missing"),
                "n_shared_contigs_with_other_arch591": d_row.get("n_exact_shared_contig_seqs_with_other_arch591", 0),
                "crispr_screening_status": scr_status,
                "crispr_screened": int(screened_ok),
                "n_raw_crispr_arrays": n_arr,
                "n_verified_spacers": n_spc,
                "n_cas_operons": n_cas,
                "c8_defense_systems": s_row.get("c8_defense_systems", "NA"),
                "n_tier1_verified_direct_links": t_counts.get("Tier1_verified_direct", 0),
                "n_tier2_corroborated_links": t_counts.get("Tier2_corroborated_single_spacer_or_provirus", 0),
                "n_tier3_unknown_viral_identity_retained": t_counts.get("Tier3_unverified_viral_identity_retained", 0),
                "n_tier3_unexcised_provirus_retained": t_counts.get("Tier3_unresolved_unexcised_provirus_retained", 0),
                "n_tier3_ambiguous_or_exploratory": (
                    t_counts.get("Tier3_ambiguous_multihost_retained", 0)
                    + t_counts.get("Tier3_exploratory_or_species_projected", 0)
                ),
                "n_external_habitat_viral_control_matches": t_counts.get("External_habitat_viral_control_match_retained", 0),
                "n_cellular_or_plasmid_excluded": t_counts.get("Cellular_or_plasmid_excluded", 0),
                "n_curated_decoy_rejected": t_counts.get("Decoy_control_rejected", 0),
                "n_technical_self_array_masked": t_counts.get("Technical_self_crispr_array_masked", 0),
                "verified_linked_votu_ids": ";".join(linked_v) if linked_v else "none",
                "mag_host_disposition": disp,
            }
        )

    ledger_df = pd.DataFrame(ledger_rows)
    ledger_df.to_csv(RES_DIR / "stage3_all591_archaeal_mag_host_evidence.tsv", sep="\t", index=False)

    # 8. Historical C1 8-Thermoplasmatota virus-host edge forensic audit (using actual `edges.tsv` schema!)
    c1_audit_rows: list[dict[str, Any]] = []
    edges_path = IN_DIR / "edges.tsv"
    if edges_path.exists():
        edf = pd.read_csv(edges_path, sep="\t")
        v_col = "votu_rep" if "votu_rep" in edf.columns else "virus_id"
        spc_col = "n_spacers" if "n_spacers" in edf.columns else "spacer_hit_count"
        thermo_edges = edf[edf["host_phylum"].isin(["Thermoplasmatota", "p__Thermoplasmatota"])].copy()

        for _, er in thermo_edges.iterrows():
            h_mag = str(er["host_mag"])
            v_id = str(er[v_col])
            h_status, rep_mag = reconcile_host_alias_strict(h_mag, arch591_set, nonrep_to_species)
            t_meta = target_id_lookup.get(v_id, {})

            matching_pairs = pd.DataFrame()
            if len(pairs_df) > 0:
                matching_pairs = pairs_df[
                    ((pairs_df["host_mag_raw"] == h_mag) | (pairs_df["canonical_591_host_mag"] == rep_mag))
                    & ((pairs_df["target_id"] == v_id) | (pairs_df["votu_id"] == v_id))
                ]

            if len(matching_pairs) > 0:
                mp = matching_pairs.iloc[0]
                re_recovered = 1
                ev_tier = str(mp["evidence_tier"])
                v_status = str(mp["viral_identity_status"])
                el_prim = int(mp["eligible_stage4_primary"])
                el_sens = int(mp["eligible_stage4_sensitivity"])
                excl = str(mp["exclusion_reason"])
            else:
                re_recovered = 0
                v_status = str(t_meta.get("viral_identity_status", "Unknown_unverified_viral_identity"))
                ev_tier = "Historical_edge_not_re_recovered_at_strict_gate"
                el_prim = 0
                el_sens = 0
                excl = "not_re_recovered_or_nonrep_bin"

            c1_audit_rows.append(
                {
                    "historical_host_mag": h_mag,
                    "canonical_591_species_rep": rep_mag if rep_mag else "unmapped",
                    "host_reconciliation_status": h_status,
                    "historical_votu_rep": v_id,
                    "historical_sample_id": str(er.get("sample_id", "NA")),
                    "historical_compartment": str(er.get("compartment", "NA")),
                    "historical_n_spacers": int(er[spc_col]) if spc_col in er else "NA",
                    "re_recovered_in_stage2": re_recovered,
                    "stage2_viral_identity_status": v_status,
                    "stage2_checkv_quality": t_meta.get("checkv_quality", "NA"),
                    "stage2_checkv_viral_genes": t_meta.get("checkv_viral_genes", "NA"),
                    "stage2_checkv_host_genes": t_meta.get("checkv_host_genes", "NA"),
                    "stage2_virus_hallmarks": t_meta.get("virus_hallmarks", "NA"),
                    "stage2_min_independent_viral_orfs_lower_bound": t_meta.get(
                        "min_independent_viral_orfs_lower_bound", "NA"
                    ),
                    "stage2_provirus_boundary_status": t_meta.get("provirus_boundary_status", "NA"),
                    "stage2_read_junction_support_status": t_meta.get(
                        "read_junction_support_status", "untested_missing"
                    ),
                    "stage3_evidence_tier": ev_tier,
                    "eligible_stage4_primary": el_prim,
                    "eligible_stage4_sensitivity": el_sens,
                    "audit_disposition_reason": excl,
                }
            )

    pd.DataFrame(c1_audit_rows).to_csv(
        RES_DIR / "stage3_historical_c1_thermoplasmatota_audit.tsv", sep="\t", index=False
    )

    thermo_ledger = ledger_df[ledger_df["canonical_phylum"] == "p__Thermoplasmatota"]
    summary = {
        "stage": "stage3_host_link_ledger",
        "status": "completed",
        "code_version_search": "stage2_viral_identity_and_spacers_v2.py (observed snapshot remote exp3_viruses/provenance/root_observed_15892214_20261109_1539, sha256: d714a436dbf89bef30023fe736e6a2874dcef3129d2c79694e5b7645e2f1d324)",
        "code_version_validation_and_collector": "stage2_viral_identity_and_spacers_v3.py",
        "stage3_version": "stage3_host_link_ledger.py",
        "root_observed_snapshot": "provenance/root_observed_15892214_20261109_1539",
        "elapsed_seconds": round(time.time() - t0, 2),
        "n_canonical_591_arch_mags": int(len(ledger_df)),
        "n_canonical_237_thermo_mags": int(len(thermo_ledger)),
        "arch591_dispositions": ledger_df["mag_host_disposition"].value_counts().to_dict(),
        "thermo237_dispositions": thermo_ledger["mag_host_disposition"].value_counts().to_dict(),
        "pair_tier_counts": (
            pairs_df["evidence_tier"].value_counts().to_dict() if len(pairs_df) > 0 else {}
        ),
        "n_historical_c1_thermo_edges_audited": len(c1_audit_rows),
    }
    (RES_DIR / "stage3_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("[Stage 3] Summary:", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
