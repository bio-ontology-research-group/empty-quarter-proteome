#!/usr/bin/env python3
"""Stage 3 v2: Full 591-Archaeal-MAG Host-Link Evidence Ledger and Forensic Audit.

Updates from Root review (16:49):
1. Disambiguates unassigned genus (`g__`): does NOT collapse unrelated unknown lineages
   into a single genus; treats unassigned genera distinctly by species/genome.
2. Identifies cellular-dominated contigs without provirus boundaries (e.g. 70 host genes vs 1 viral gene)
   and properly classifies them as `Unresolved_unexcised_provirus_candidate` in Tier 3.
3. Deduplicates exact identical sequence duplicates (e.g. `48Sr2::48Sr2_190170` and
   `coasm_48_S::k127_1210458` sharing SHA-256 `ca78a3b6ce5eca9cb6608e7bd98477ff1c6b486ee21ba2ace1d45437ece707f1`),
   marking primary sequence representatives and aliases.
4. Preserves baseline outputs and writes versioned v2 artifacts.
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


def load_votu_membership_map(membership_tsv: Path) -> dict[str, str]:
    """Map contig IDs to cluster vOTU representatives."""
    if not membership_tsv.exists():
        return {}
    c2v: dict[str, str] = {}
    with open(membership_tsv) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            c_id = r.get("member_id") or r.get("contig_id") or r.get("sequence_id")
            v_id = r.get("cluster_rep") or r.get("votu_rep") or r.get("votu_id")
            if c_id and v_id:
                c2v[str(c_id).split()[0]] = str(v_id).split()[0]
    return c2v


def main() -> None:
    t0 = time.time()
    print("[Stage 3 v2] Starting evidence ledger and forensic audit...", flush=True)

    # 1. Load canonical 591 Archaeal MAG metadata
    ver_df = pd.read_csv(RES_DIR / "stage0_arch591_fasta_verification.tsv", sep="\t")
    arch591_set = set(ver_df["genome"].astype(str))
    tax_lookup: dict[str, dict[str, str]] = {}
    for _, r in ver_df.iterrows():
        g_id = str(r["genome"])
        tax_lookup[g_id] = {
            "canonical_phylum": str(r.get("gtdbtk_phylum", "p__unknown")),
            "canonical_class": str(r.get("gtdbtk_class", "c__unknown")),
            "canonical_order": str(r.get("gtdbtk_order", "o__unknown")),
            "canonical_family": str(r.get("gtdbtk_family", "f__unknown")),
            "canonical_genus": str(r.get("gtdbtk_genus", "g__unknown")),
            "canonical_species": str(r.get("gtdbtk_species", "s__unknown")),
        }

    # 2. Load species representative reconciliation map
    nonrep_to_species: dict[str, str] = {}
    spec_rep_path = IN_DIR / "canonical_archaeal_species_reps.tsv"
    if spec_rep_path.exists():
        with open(spec_rep_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("nonrep_mag_id") and r.get("species_rep_mag_id"):
                    nonrep_to_species[str(r["nonrep_mag_id"])] = str(r["species_rep_mag_id"])

    # 3. Load vOTU cluster membership
    contig_to_votu = load_votu_membership_map(P_VOTU_MEM)

    # 4. Load Stage 2 scored alignments and target viral identity table
    scored_path = RES_DIR / "stage2_all_spacer_alignments_scored.tsv"
    if not scored_path.exists():
        raise FileNotFoundError(f"Missing required Stage 2 output: {scored_path}")
    scored_df = pd.read_csv(scored_path, sep="\t")

    v_identity_path = RES_DIR / "stage2_target_viral_identity_table.tsv"
    v_identity_df = (
        pd.read_csv(v_identity_path, sep="\t") if v_identity_path.exists() else pd.DataFrame()
    )
    v_identity_map: dict[str, dict[str, Any]] = {}
    if len(v_identity_df) > 0:
        for _, vr in v_identity_df.iterrows():
            v_identity_map[str(vr["target_id"])] = vr.to_dict()

    # 5. Build strict genus / phylum maps for multi-host ambiguity audit
    if len(scored_df) > 0:
        scored_df["votu_or_target_id"] = scored_df["target_id"].apply(
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
                # Disambiguate unassigned genus: do NOT treat 'g__' as a shared genus!
                raw_genus = str(tax_lookup[rep_mag]["canonical_genus"])
                effective_genus = (
                    raw_genus
                    if raw_genus not in ["g__", "", "NA", "nan", "g__unknown"]
                    else f"unassigned_genus_{rep_mag}"
                )
                virus_to_genera[v_key].add(effective_genus)
                virus_to_phyla[v_key].add(str(tax_lookup[rep_mag]["canonical_phylum"]))

    # 6. Evaluate all `(host_mag, target_id)` pairs
    pair_rows: list[dict[str, Any]] = []
    mag_to_pair_tiers: dict[str, list[str]] = collections.defaultdict(list)

    if len(scored_df) > 0:
        grouped = scored_df.groupby(["host_mag", "target_id"])
        for (h_mag, t_id), grp in grouped:
            h_status, rep_mag = reconcile_host_alias_strict(h_mag, arch591_set, nonrep_to_species)
            v_key = contig_to_votu.get(str(t_id), str(t_id))
            unmasked_grp = grp[grp["is_technical_crispr_array_self_match"] == 0]

            if len(unmasked_grp) > 0:
                ref_row = unmasked_grp.iloc[0]
                n_strict = int(
                    unmasked_grp[unmasked_grp["passes_strict_spacer_gate"] == True]["unique_spacer_id"].nunique()
                )
                n_exact = int(
                    unmasked_grp[
                        (unmasked_grp["passes_strict_spacer_gate"] == True)
                        & (unmasked_grp["effective_mismatches"] == 0)
                    ]["unique_spacer_id"].nunique()
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
                "is_decoy": bool(int(ref_row["is_decoy"])),
                "is_technical_crispr_array_self_match": is_tech_self,
                "is_external_viral_specificity_control": bool(int(ref_row["is_external_viral_specificity_control"])),
                "status_detail": str(ref_row["status_detail"]),
            }

            has_cas = bool(int(ref_row["host_has_cas_operon"]))
            has_prov_junc = (
                str(ref_row["provirus_boundary_status"]) == "verified_provirus_with_host_flanks"
                and bool(int(ref_row["is_same_namespaced_contig"]))
            )

            cv_host = int(ref_row.get("checkv_host_genes", 0))
            cv_viral = int(ref_row.get("checkv_viral_genes", 0))
            # Individual evidence classification for 41Sr3::41Sr3_56779 is retained
            # in individual pair audit review (stage3_tier1_pairs_individual_evidence_audit.json)
            # without applying arbitrary global numerical ratio thresholds.

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
                "canonical_591_host_mag": str(rep_mag) if rep_mag else "unmapped_nonrep",
                "host_reconciliation_status": h_status,
                "host_phylum": t_row_tax.get("canonical_phylum", "p__unknown"),
                "host_order": t_row_tax.get("canonical_order", "o__unknown"),
                "host_genus": t_row_tax.get("canonical_genus", "g__unknown"),
                "host_species": t_row_tax.get("canonical_species", "s__unknown"),
                "target_id": str(t_id),
                "target_seq_sha256": str(ref_row["target_seq_sha256"]),
                "annotation_provenance": str(ref_row["annotation_provenance"]),
                "votu_id": str(v_key),
                "is_votu_cluster_rep": str(ref_row["is_votu_cluster_rep"]),
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
                "viral_identity_status": v_eval_dict["viral_identity_status"],
                "habitat_specificity_role": str(ref_row["habitat_specificity_role"]),
                "checkv_quality": str(ref_row.get("checkv_quality", "NA")),
                "checkv_viral_genes": cv_viral,
                "checkv_host_genes": cv_host,
                "virus_hallmarks": int(ref_row.get("virus_hallmarks", 0)),
                "archaeal_viral_hallmarks": int(ref_row.get("archaeal_viral_hallmarks", 0)),
                "min_independent_viral_orfs_lower_bound": ref_row.get("min_independent_viral_orfs_lower_bound", "NA"),
                "orf_count_provenance": ref_row.get("orf_count_provenance", "NA"),
                "provirus_boundary_status": str(ref_row.get("provirus_boundary_status", "NA")),
                "read_junction_support_status": str(ref_row.get("read_junction_support_status", "untested_missing")),
                "n_distinct_host_genera_for_virus": n_gen,
                "n_distinct_host_phyla_for_virus": n_phy,
                "evidence_tier": tier_res["evidence_tier"],
                "eligible_stage4_primary": int(tier_res["eligible_stage4_primary"]),
                "eligible_stage4_sensitivity": int(tier_res["eligible_stage4_sensitivity"]),
                "exclusion_reason": tier_res["exclusion_reason"],
            }
            pair_rows.append(pair_entry)
            mag_to_pair_tiers[str(h_mag)].append(tier_res["evidence_tier"])
            if rep_mag and rep_mag != h_mag:
                mag_to_pair_tiers[str(rep_mag)].append(tier_res["evidence_tier"])

    pairs_df = pd.DataFrame(pair_rows)

    # 7. Deduplicate identical sequence duplicates for same host
    if len(pairs_df) > 0:
        pairs_df["is_primary_sequence_representative"] = 1
        pairs_df["sequence_duplicate_aliases"] = ""

        # Group by (canonical_591_host_mag, target_seq_sha256)
        grouped_seqs = pairs_df.groupby(["canonical_591_host_mag", "target_seq_sha256"])
        for (h_can, s_sha), g_idx in grouped_seqs.groups.items():
            if len(g_idx) > 1:
                # Mark first as primary representative, others as duplicates
                first_i = g_idx[0]
                other_i = g_idx[1:]
                other_aliases = ";".join(pairs_df.loc[other_i, "target_id"].astype(str))
                pairs_df.loc[first_i, "sequence_duplicate_aliases"] = other_aliases
                for dup_i in other_i:
                    pairs_df.loc[dup_i, "is_primary_sequence_representative"] = 0
                    pairs_df.loc[dup_i, "sequence_duplicate_aliases"] = str(pairs_df.loc[first_i, "target_id"])

    pairs_df.to_csv(RES_DIR / "stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv", sep="\t", index=False)

    # Load Stage 0 and Stage 1 verification tables
    ver_df = pd.read_csv(RES_DIR / "stage0_arch591_fasta_verification.tsv", sep="\t")
    ver_lookup = ver_df.set_index("genome").to_dict("index")
    st1_df = pd.read_csv(RES_DIR / "stage1_mag_screening_status.tsv", sep="\t")
    st1_lookup = st1_df.set_index("mag_id").to_dict("index")

    # 8. Full 591 Archaeal MAG Host Evidence Ledger
    ledger_rows: list[dict[str, Any]] = []
    for g_id in sorted(arch591_set):
        t_info = tax_lookup[g_id]
        tiers_seen = mag_to_pair_tiers.get(g_id, [])
        v_row = ver_lookup.get(g_id, {})
        s_row = st1_lookup.get(g_id, {})
        fasta_ok = (v_row.get("refined_verification_status") == "verified_exact_match_genome_manifest")
        scr_status = str(s_row.get("screening_status", ""))
        screened_ok = (scr_status == "screened_both_callers_ok")
        caller_failed = (scr_status == "failed_caller_execution")
        n_arr = int(s_row.get("n_raw_arrays", 0) or 0)
        n_spc = int(s_row.get("n_independent_spacers", s_row.get("n_valid_spacers", 0)) or 0)
        disposition = summarize_mag_host_disposition(
            mag_id=g_id,
            fasta_verified=fasta_ok,
            crispr_screened=screened_ok,
            n_raw_arrays=n_arr,
            n_verified_spacers=n_spc,
            pair_tiers=tiers_seen,
            caller_failed=caller_failed,
        )
        ledger_rows.append(
            {
                "genome": g_id,
                "canonical_phylum": t_info["canonical_phylum"],
                "canonical_class": t_info["canonical_class"],
                "canonical_order": t_info["canonical_order"],
                "canonical_family": t_info["canonical_family"],
                "canonical_genus": t_info["canonical_genus"],
                "canonical_species": t_info["canonical_species"],
                "mag_host_disposition": disposition,
                "n_candidate_virus_pairs": len(tiers_seen),
                "highest_evidence_tier": tiers_seen[0] if tiers_seen else "Zero_spacers_or_no_viral_hits",
            }
        )

    ledger_df = pd.DataFrame(ledger_rows)
    ledger_df.to_csv(RES_DIR / "stage3_all591_archaeal_mag_host_evidence_v2.tsv", sep="\t", index=False)

    # 9. Historical C1 Thermoplasmatota forensic audit v2
    c1_edges_path = IN_DIR / "edges.tsv"
    c1_audit_rows: list[dict[str, Any]] = []
    if c1_edges_path.exists():
        edf = pd.read_csv(c1_edges_path, sep="\t")
        h_col = "host_mag"
        v_col = "votu_rep" if "votu_rep" in edf.columns else "virus_id"
        spc_col = "n_spacers"

        for _, er in edf.iterrows():
            h_mag = str(er[h_col])
            v_id = str(er[v_col])
            h_status, rep_mag = reconcile_host_alias_strict(h_mag, arch591_set, nonrep_to_species)

            matching_pairs = (
                pairs_df[
                    (pairs_df["host_mag_raw"] == h_mag)
                    & (
                        (pairs_df["target_id"] == v_id)
                        | (pairs_df["votu_id"] == v_id)
                        | (pairs_df["canonical_591_host_mag"] == rep_mag)
                        & (pairs_df["target_id"] == v_id)
                    )
                ]
                if len(pairs_df) > 0
                else pd.DataFrame()
            )

            if len(matching_pairs) > 0:
                p_top = matching_pairs.iloc[0]
                re_recovered = 1
                v_status = str(p_top["viral_identity_status"])
                ev_tier = str(p_top["evidence_tier"])
                el_prim = int(p_top["eligible_stage4_primary"])
                el_sens = int(p_top["eligible_stage4_sensitivity"])
                excl = str(p_top["exclusion_reason"])
                t_meta = v_identity_map.get(str(p_top["target_id"]), {})
            else:
                t_meta = v_identity_map.get(v_id, {})
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
                    "stage2_provirus_boundary_status": t_meta.get("provirus_boundary_status", "NA"),
                    "stage3_evidence_tier": ev_tier,
                    "eligible_stage4_primary": el_prim,
                    "eligible_stage4_sensitivity": el_sens,
                    "audit_disposition_reason": excl,
                }
            )

    pd.DataFrame(c1_audit_rows).to_csv(
        RES_DIR / "stage3_historical_c1_thermoplasmatota_audit_v2.tsv", sep="\t", index=False
    )

    # 10. Summary v2
    thermo_ledger = ledger_df[ledger_df["canonical_phylum"] == "p__Thermoplasmatota"]
    t1_pairs = pairs_df[pairs_df["evidence_tier"] == "Tier1_verified_direct"] if len(pairs_df) > 0 else pd.DataFrame()
    t1_unique_seq_pairs = (
        t1_pairs[t1_pairs["is_primary_sequence_representative"] == 1]
        if len(t1_pairs) > 0
        else pd.DataFrame()
    )

    summary = {
        "stage": "stage3_host_link_ledger_v2",
        "status": "completed",
        "code_version_search": "stage2_viral_identity_and_spacers_v2.py (observed snapshot remote exp3_viruses/provenance/root_observed_15892214_20261109_1539, sha256: d714a436dbf89bef30023fe736e6a2874dcef3129d2c79694e5b7645e2f1d324)",
        "code_version_validation_and_collector": "stage2_viral_identity_and_spacers_v3.py",
        "stage3_version": "stage3_host_link_ledger_v2.py",
        "root_observed_snapshot": "provenance/root_observed_15892214_20261109_1539",
        "elapsed_seconds": round(time.time() - t0, 2),
        "n_canonical_591_arch_mags": int(len(ledger_df)),
        "n_canonical_237_thermo_mags": int(len(thermo_ledger)),
        "n_tier1_direct_pair_rows": len(t1_pairs),
        "n_tier1_unique_host_sequence_pairs": len(t1_unique_seq_pairs),
        "n_tier1_unique_host_mags": int(t1_unique_seq_pairs["canonical_591_host_mag"].nunique()) if len(t1_unique_seq_pairs) > 0 else 0,
        "n_tier1_thermoplasmatota_pairs": int((t1_unique_seq_pairs["host_phylum"] == "p__Thermoplasmatota").sum()) if len(t1_unique_seq_pairs) > 0 else 0,
        "arch591_dispositions": ledger_df["mag_host_disposition"].value_counts().to_dict(),
        "thermo237_dispositions": thermo_ledger["mag_host_disposition"].value_counts().to_dict(),
        "pair_tier_counts": (
            pairs_df["evidence_tier"].value_counts().to_dict() if len(pairs_df) > 0 else {}
        ),
        "n_historical_c1_thermo_edges_audited": len(c1_audit_rows),
    }
    (RES_DIR / "stage3_summary_v2.json").write_text(json.dumps(summary, indent=2) + "\n")
    (RES_DIR / "STAGE3_V2_DONE").touch()
    print("[Stage 3 v2] Summary:", json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
