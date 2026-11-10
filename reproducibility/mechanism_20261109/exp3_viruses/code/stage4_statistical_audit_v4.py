#!/usr/bin/env python3
"""Stage 4 v4: scoped viral statistical audit after the v3 inference repair.

v3 is preserved as an invalid-inference baseline.  This module does not import it
or rewrite its outputs.  All inferential tables produced here are exploratory:
DETECTEDLINK is a conditional detection outcome, not infection absence/presence.

Root-design invariants implemented here:
* exact host + sequence-SHA deduplication, no global CheckV gene threshold;
* primary pairs derived from Tier 1 ledger rows minus exact root dispositions;
* zero-spacer genomes are untested, not negatives;
* actual missingness/rank/separation diagnostics without completeness imputation;
* fixed 47-family BH denominators with non-estimable tests retained as NA;
* metadata-only trip x compartment selection and within-dune-field site permutations;
* matched alternative MAGs without fabricated genetic distances or zero filling.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import optimize, stats

BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
DEFAULT_METADATA = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/metadata.tsv")
DEFAULT_MAG_ABUND = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/mag_abundance.tsv.gz")
DEFAULT_VIR_TPM = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/viral_tpm.tsv.gz")
DESIGN_LEAD_FAMILIES = ("22PRr1::22PRr1_2531955_3", "V60PRr2::V60PRr2_1314979_4")
TAX_COLUMNS = [
    "canonical_phylum", "canonical_class", "canonical_order", "canonical_family",
    "canonical_genus", "canonical_species", "canonical_completeness", "canonical_contamination",
]
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def require_columns(df: pd.DataFrame, columns: Iterable[str], name: str) -> None:
    missing = sorted(set(columns) - set(df.columns))
    if missing:
        raise ValueError(f"{name} missing required columns: {missing}")


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"1", "true", "t", "yes"})


def known_taxon(value: Any, prefix: str) -> bool:
    if pd.isna(value):
        return False
    text = str(value).strip().lower()
    return text not in {"", prefix, prefix + "unknown", "unknown", "na", "nan", "none"}


def bh_adjust(pvalues: pd.Series | np.ndarray, family_size: int) -> np.ndarray:
    """BH over a fixed hypothesis family; NA tests remain NA but count in m."""
    p = np.asarray(pvalues, dtype=float)
    if family_size < 1 or family_size < int(np.isfinite(p).sum()):
        raise ValueError("family_size must be >= the number of finite p-values")
    q = np.full(p.shape, np.nan, dtype=float)
    valid = np.flatnonzero(np.isfinite(p))
    if not len(valid):
        return q
    if ((p[valid] < 0) | (p[valid] > 1)).any():
        raise ValueError("p-values outside [0,1]")
    order = valid[np.argsort(p[valid], kind="mergesort")]
    running = 1.0
    for rank in range(len(order), 0, -1):
        idx = order[rank - 1]
        running = min(running, p[idx] * family_size / rank)
        q[idx] = min(1.0, running)
    return q


def load_root_unresolved_targets(path: Path) -> set[str]:
    """Read the individual root disposition; never generalize it to a threshold."""
    payload = json.loads(path.read_text())
    target = payload.get("target")
    disposition = str(payload.get("root_review_disposition", ""))
    if not target or not disposition.startswith("unresolved_"):
        raise ValueError("Root disposition must contain target and unresolved root_review_disposition")
    return {str(target)}


def load_provisional_candidate(path: Path) -> tuple[tuple[str, str, str], dict[str, Any]]:
    """Read the completed root review of the boundary-sensitive second Thermo pair."""
    payload = json.loads(path.read_text())
    required = ["host", "target", "target_sha256", "disposition", "canonical_spacer_sequence",
                "spacer_target_interval_1based_inclusive", "catalog_viral_interval", "fresh_checkv_viral_interval"]
    missing = [key for key in required if key not in payload]
    if missing:
        raise ValueError(f"Provisional disposition missing fields: {missing}")
    if payload["disposition"] != "boundary_sensitive_provisional_not_primary_host_link":
        raise ValueError("Provisional candidate lacks exact boundary-sensitive root disposition")
    if not SHA256_RE.fullmatch(str(payload["target_sha256"])) or not payload.get("exact_sequence_witness_verified"):
        raise ValueError("Provisional candidate lacks exact SHA/sequence witness")
    spacer = payload["spacer_target_interval_1based_inclusive"]
    core = payload["catalog_viral_interval"]
    extension = payload["fresh_checkv_viral_interval"]
    if not (spacer[0] > core[1] and core[0] <= spacer[0] <= spacer[1] <= extension[1]):
        raise ValueError("Provisional spacer is not outside catalog core and inside CheckV extension")
    if len(payload["canonical_spacer_sequence"]) != spacer[1] - spacer[0] + 1:
        raise ValueError("Provisional spacer length and witness interval disagree")
    return (str(payload["host"]), str(payload["target"]), str(payload["target_sha256"])), payload


def resolve_graded_candidate_pairs(
    pairs_df: pd.DataFrame,
    unresolved_targets: set[str],
    provisional_candidates: set[tuple[str, str, str]] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Preserve every raw row and derive primary/provisional sets from the ledger.

    Primary = unique canonical-host/exact-SHA Tier 1 rows not individually unresolved
    and not individually boundary-sensitive provisional. Provisional sensitivity adds
    the root-reviewed boundary-sensitive candidate. No host list or CheckV threshold is embedded.
    """
    provisional_candidates = provisional_candidates or set()
    required = [
        "host_mag_raw", "canonical_591_host_mag", "target_id", "target_seq_sha256",
        "votu_id", "evidence_tier", "checkv_host_genes", "checkv_viral_genes",
    ]
    require_columns(pairs_df, required, "stage3 pairs")
    ledger = pairs_df.copy()
    ledger.insert(0, "raw_stage3_row_number", np.arange(1, len(ledger) + 1))
    for col in ["canonical_591_host_mag", "target_id", "target_seq_sha256"]:
        if ledger[col].isna().any() or (ledger[col].astype(str).str.strip() == "").any():
            raise ValueError(f"stage3 pairs has missing/blank {col}")
    bad_sha = ~ledger["target_seq_sha256"].astype(str).map(lambda x: bool(SHA256_RE.fullmatch(x)))
    if bad_sha.any():
        raise ValueError(f"stage3 pairs has invalid target_seq_sha256 at raw rows {list(ledger.loc[bad_sha, 'raw_stage3_row_number'][:5])}")

    for col in ["is_alias_duplicate", "is_unresolved_root_disposition", "is_boundary_sensitive_provisional",
                "is_primary_detected_link", "is_provisional_detected_link"]:
        ledger[col] = 0
    ledger["review_disposition_reason"] = "non_tier1_not_in_frozen_detection_outcome"
    ledger["host_target_sha_key"] = ledger["canonical_591_host_mag"].astype(str) + "|" + ledger["target_seq_sha256"].astype(str)

    t1 = ledger["evidence_tier"].eq("Tier1_verified_direct")
    canonical = ledger["canonical_591_host_mag"].notna() & ~ledger["canonical_591_host_mag"].isin(["", "unmapped_nonrep"])
    eligible_t1 = t1 & canonical
    representatives: dict[tuple[str, str], int] = {}
    for (host, sha), group in ledger.loc[eligible_t1].groupby(["canonical_591_host_mag", "target_seq_sha256"], sort=True):
        preference = pd.to_numeric(group.get("is_primary_sequence_representative", 0), errors="coerce").fillna(0)
        ordered = group.assign(_preference=preference).sort_values(
            ["_preference", "target_id", "raw_stage3_row_number"], ascending=[False, True, True], kind="mergesort"
        )
        representatives[(str(host), str(sha))] = int(ordered.index[0])

    found_provisional: set[tuple[str, str, str]] = set()
    for idx, row in ledger.loc[t1].iterrows():
        host, target, sha = str(row["canonical_591_host_mag"]), str(row["target_id"]), str(row["target_seq_sha256"])
        key = (host, sha)
        if not canonical.loc[idx]:
            reason = "tier1_rejected_noncanonical_host"
        elif representatives.get(key) != idx:
            ledger.loc[idx, "is_alias_duplicate"] = 1
            reason = "duplicate_host_exact_sequence_sha_alias_collapsed"
        elif target in unresolved_targets:
            ledger.loc[idx, "is_unresolved_root_disposition"] = 1
            reason = "individual_root_disposition_unresolved_not_primary_or_negative_control"
        elif (host, target, sha) in provisional_candidates:
            ledger.loc[idx, "is_boundary_sensitive_provisional"] = 1
            ledger.loc[idx, "is_provisional_detected_link"] = 1
            found_provisional.add((host, target, sha))
            reason = "root_reviewed_boundary_sensitive_provisional_exact_spacer_outside_catalog_core_sensitivity_only"
        else:
            ledger.loc[idx, "is_primary_detected_link"] = 1
            ledger.loc[idx, "is_provisional_detected_link"] = 1
            reason = "root_reviewed_primary_unique_host_exact_sha_detected_link"
        ledger.loc[idx, "review_disposition_reason"] = reason

    if found_provisional != provisional_candidates:
        raise ValueError(f"Provisional candidates not found as unique Tier1 host/target/SHA rows: {sorted(provisional_candidates - found_provisional)}")
    missing_unresolved = sorted(unresolved_targets - set(ledger["target_id"].astype(str)))
    if missing_unresolved:
        raise ValueError(f"Root unresolved targets absent from stage3 ledger: {missing_unresolved}")

    audit_cols = list(ledger.columns)
    primary = ledger.loc[ledger["is_primary_detected_link"].eq(1), audit_cols].copy()
    provisional = ledger.loc[ledger["is_provisional_detected_link"].eq(1), audit_cols].copy()
    return ledger, primary, provisional


def join_candidate_taxonomy(candidate_df: pd.DataFrame, taxonomy_df: pd.DataFrame) -> pd.DataFrame:
    require_columns(candidate_df, ["canonical_591_host_mag"], "candidate pairs")
    require_columns(taxonomy_df, ["genome"] + TAX_COLUMNS, "canonical taxonomy")
    if taxonomy_df["genome"].duplicated().any():
        raise ValueError("canonical taxonomy has duplicate genome rows")
    tax = taxonomy_df[["genome"] + TAX_COLUMNS]
    out = candidate_df.merge(tax, left_on="canonical_591_host_mag", right_on="genome", how="left", validate="many_to_one")
    if out["genome"].isna().any():
        raise ValueError("candidate host missing exact canonical taxonomy join")
    out["taxonomy_join_status"] = "exact_canonical_genome_join_stage3_tax_fields_not_used"
    return out.drop(columns=["genome"])


def partition_detection_opportunity(
    screening_df: pd.DataFrame,
    evidence_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame,
    primary_pairs_df: pd.DataFrame,
    provisional_pairs_df: pd.DataFrame,
) -> pd.DataFrame:
    require_columns(evidence_df, ["genome", "highest_evidence_tier"], "stage3 evidence")
    require_columns(screening_df, ["mag_id", "is_arch591_rep", "n_valid_spacers", "n_raw_spacers", "n_cas_operons", "minced_ok", "cctyper_ok", "screening_status"], "stage1 screening")
    require_columns(taxonomy_df, ["genome"] + TAX_COLUMNS, "canonical taxonomy")
    if evidence_df["genome"].duplicated().any() or screening_df["mag_id"].duplicated().any() or taxonomy_df["genome"].duplicated().any():
        raise ValueError("Duplicate genome keys in detection-partition inputs")
    out = evidence_df[["genome", "highest_evidence_tier"]].merge(
        screening_df, left_on="genome", right_on="mag_id", how="left", validate="one_to_one", indicator=True
    )
    if not out["_merge"].eq("both").all():
        raise ValueError("Canonical evidence genomes missing exact stage1 screening rows")
    if not as_bool(out["is_arch591_rep"]).all():
        raise ValueError("Canonical evidence joined a non-arch591 screening row")
    out = out.drop(columns=["mag_id", "_merge"]).merge(
        taxonomy_df[["genome"] + TAX_COLUMNS], on="genome", how="left", validate="one_to_one"
    )
    out["n_valid_spacers"] = pd.to_numeric(out["n_valid_spacers"], errors="raise")
    if out["n_valid_spacers"].isna().any():
        raise ValueError("Missing valid-spacer counts cannot be imputed")
    assay_ok = as_bool(out["minced_ok"]) & as_bool(out["cctyper_ok"])
    positive_primary = set(primary_pairs_df["canonical_591_host_mag"])
    positive_provisional = set(provisional_pairs_df["canonical_591_host_mag"])
    if not positive_provisional <= set(out["genome"]):
        raise ValueError("Candidate host absent from canonical evidence universe")
    out["is_primary_detected_link_host"] = out["genome"].isin(positive_primary).astype(int)
    out["is_provisional_detected_link_host"] = out["genome"].isin(positive_provisional).astype(int)
    out["is_screen_positive_eligible"] = (assay_ok & out["n_valid_spacers"].gt(0)).astype(int)
    out["primary_DETECTEDLINK"] = np.where(out["is_screen_positive_eligible"].eq(1), out["is_primary_detected_link_host"], np.nan)
    out["provisional_DETECTEDLINK"] = np.where(out["is_screen_positive_eligible"].eq(1), out["is_provisional_detected_link_host"], np.nan)
    out["detection_opportunity_category"] = np.select(
        [~assay_ok, out["n_valid_spacers"].eq(0)],
        ["screen_failed_or_missing_untested", "zero_spacers_untested_not_negative"],
        default="screen_positive_conditional_detection_cohort",
    )
    out["taxonomy_assay_status"] = np.where(out["canonical_phylum"].notna(), "exact_taxonomy_join", "taxonomy_missing_not_imputed")
    if (out["is_primary_detected_link_host"].eq(1) & out["is_screen_positive_eligible"].ne(1)).any():
        raise ValueError("Primary detected-link host is outside the screen-positive eligible cohort")
    return out


def load_and_orient_abundance_matrices(mag_path: Path, vir_path: Path, auth_samples: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    def orient(path: Path, name: str) -> pd.DataFrame:
        df = pd.read_csv(path, sep="\t", index_col=0)
        if df.index.duplicated().any():
            raise ValueError(f"{name} has duplicate entity keys")
        idx = set(map(str, df.index))
        cols = set(map(str, df.columns))
        if idx == auth_samples:
            df = df.T
        elif cols != auth_samples:
            idx_overlap, col_overlap = len(idx & auth_samples), len(cols & auth_samples)
            if idx_overlap > col_overlap:
                df = df.loc[[x for x in df.index if str(x) in auth_samples]].T
            else:
                df = df[[x for x in df.columns if str(x) in auth_samples]]
        df.index = df.index.astype(str)
        df.columns = df.columns.astype(str)
        if df.columns.duplicated().any():
            raise ValueError(f"{name} has duplicate sample columns")
        numeric = df.apply(pd.to_numeric, errors="coerce")
        original_nonmissing = df.notna().to_numpy()
        if (~np.isfinite(numeric.to_numpy(dtype=float)) & original_nonmissing).any():
            raise ValueError(f"{name} has nonnumeric/nonfinite observed cells; missing cells remain NA")
        if (numeric.to_numpy(dtype=float) < 0).any():
            raise ValueError(f"{name} has negative abundances")
        return numeric
    return orient(mag_path, "MAG abundance"), orient(vir_path, "viral TPM")


def compute_abundance_metrics(mag_df: pd.DataFrame, genomes: Iterable[str], samples: list[str]) -> pd.DataFrame:
    rows = []
    for genome in sorted(set(genomes)):
        if genome not in mag_df.index:
            rows.append({"genome": genome, "abundance_assay_status": "missing_matrix_row_not_zero", "n_valid_abundance_samples": 0,
                         "sample_prevalence": np.nan, "mean_sample_abundance": np.nan})
            continue
        vec = mag_df.loc[genome, samples].to_numpy(dtype=float)
        finite = np.isfinite(vec)
        rows.append({"genome": genome, "abundance_assay_status": "observed_matrix_row", "n_valid_abundance_samples": int(finite.sum()),
                     "sample_prevalence": int((vec[finite] > 0).sum()),
                     "mean_sample_abundance": float(np.mean(vec[finite])) if finite.any() else np.nan})
    return pd.DataFrame(rows)


def taxonomic_relatedness(candidate: pd.Series, control: pd.Series) -> str:
    cg, hg = candidate["canonical_genus"], control["canonical_genus"]
    if known_taxon(cg, "g__") and known_taxon(hg, "g__"):
        return "known_different_genus" if str(cg) != str(hg) else "same_genus_excluded"
    # Missing genus is unknown relatedness, never evidence of different genera.
    if known_taxon(candidate["canonical_order"], "o__") and known_taxon(control["canonical_order"], "o__") and str(candidate["canonical_order"]) != str(control["canonical_order"]):
        return "known_different_order_genus_unknown"
    if known_taxon(candidate["canonical_family"], "f__") and known_taxon(control["canonical_family"], "f__") and str(candidate["canonical_family"]) != str(control["canonical_family"]):
        return "known_different_family_genus_unknown"
    return "unknown_genus_relatedness_not_established"


def find_matched_alternative_controls(
    candidate_hosts: Iterable[str],
    primary_hosts: set[str],
    provisional_hosts: set[str],
    unresolved_hosts: set[str],
    partition_df: pd.DataFrame,
    abundance_metrics_df: pd.DataFrame,
    completeness_tolerance: float = 15.0,
    prevalence_tolerance: float = 20.0,
) -> pd.DataFrame:
    """Deterministic matched alternatives; all criteria use observed values only."""
    m = partition_df.merge(abundance_metrics_df, on="genome", how="left", validate="one_to_one")
    excluded = set(primary_hosts) | set(provisional_hosts) | set(unresolved_hosts)
    rows = []
    for host in sorted(set(candidate_hosts)):
        hs = m.loc[m["genome"].eq(host)]
        if len(hs) != 1:
            raise ValueError(f"Candidate host does not have one partition row: {host}")
        h = hs.iloc[0]
        base = {
            "candidate_host_mag": host, "candidate_phylum": h["canonical_phylum"], "candidate_genus": h["canonical_genus"],
            "candidate_completeness": h["canonical_completeness"], "candidate_prevalence": h["sample_prevalence"],
            "candidate_primary_set": int(host in primary_hosts), "candidate_provisional_set": int(host in provisional_hosts),
            "matching_rule": f"same_known_phylum;observed_completeness_delta<={completeness_tolerance};observed_prevalence_delta<={prevalence_tolerance};exclude_all_primary_provisional_unresolved_hosts",
        }
        comp, prev = pd.to_numeric(pd.Series([h["canonical_completeness"]]), errors="coerce").iloc[0], pd.to_numeric(pd.Series([h["sample_prevalence"]]), errors="coerce").iloc[0]
        if not known_taxon(h["canonical_phylum"], "p__") or not np.isfinite(comp) or not np.isfinite(prev):
            rows.append(base | {"matched_alternative_mag": pd.NA, "matching_status": "NA_candidate_missing_known_phylum_completeness_or_prevalence"})
            continue
        pool = m.loc[m["canonical_phylum"].eq(h["canonical_phylum"]) & ~m["genome"].isin(excluded)].copy()
        pool["canonical_completeness"] = pd.to_numeric(pool["canonical_completeness"], errors="coerce")
        pool["sample_prevalence"] = pd.to_numeric(pool["sample_prevalence"], errors="coerce")
        pool = pool.loc[np.isfinite(pool["canonical_completeness"]) & np.isfinite(pool["sample_prevalence"])]
        pool = pool.loc[(pool["canonical_completeness"] - comp).abs().le(completeness_tolerance) & (pool["sample_prevalence"] - prev).abs().le(prevalence_tolerance)]
        if len(pool):
            pool["relatedness_status"] = [taxonomic_relatedness(h, r) for _, r in pool.iterrows()]
            pool = pool.loc[pool["relatedness_status"].ne("same_genus_excluded")]
        if not len(pool):
            rows.append(base | {"matched_alternative_mag": pd.NA, "matching_status": "NA_no_observed_value_match_within_declared_criteria"})
            continue
        pool["delta_completeness"] = pool["canonical_completeness"] - comp
        pool["delta_prevalence"] = pool["sample_prevalence"] - prev
        pool["matching_distance_not_genetic"] = np.sqrt((pool["delta_completeness"] / completeness_tolerance) ** 2 + (pool["delta_prevalence"] / prevalence_tolerance) ** 2)
        best = pool.sort_values(["matching_distance_not_genetic", "genome"], kind="mergesort").iloc[0]
        rows.append(base | {
            "matched_alternative_mag": best["genome"], "matched_alternative_phylum": best["canonical_phylum"],
            "matched_alternative_genus": best["canonical_genus"], "matched_alternative_completeness": best["canonical_completeness"],
            "matched_alternative_prevalence": best["sample_prevalence"], "delta_completeness": best["delta_completeness"],
            "delta_prevalence": best["delta_prevalence"], "matching_distance_not_genetic": best["matching_distance_not_genetic"],
            "relatedness_status": best["relatedness_status"], "matching_status": "matched_alternative_MAG_not_unrelated_taxon_claim",
            "control_detection_opportunity_category": best["detection_opportunity_category"],
        })
    return pd.DataFrame(rows)


def validate_assay_provenance(
    full_loci_df: pd.DataFrame, stage1_summary_path: Path, stage1b_summary_path: Path
) -> dict[str, Any]:
    """Validate independently declared complete protein-family assignment scope.

    Stage 1 screens the original protein-family assignment for all Thermoplasmatota
    genes; Stage 1b retains every candidate-family member locus, including failed
    promotion gates. The union is not being asserted as an assay universe on its own.
    """
    stage1 = json.loads(stage1_summary_path.read_text())
    stage1b = json.loads(stage1b_summary_path.read_text())
    required1 = ["status", "n_eq_thermoplasmatota_genomes", "n_eq_thermoplasmatota_genes"]
    required1b = ["status", "n_refined_genome_fastas_verified_against_manifest", "n_candidate_member_loci_audited"]
    if any(key not in stage1 for key in required1) or any(key not in stage1b for key in required1b):
        raise ValueError("Assay completion receipts lack declared genome/gene/locus scope")
    if stage1["status"] != "completed" or stage1b["status"] != "completed":
        raise ValueError("Assay source pipelines are not complete")
    n_genomes = int(full_loci_df["genome"].nunique())
    if len(full_loci_df) != int(stage1b["n_candidate_member_loci_audited"]):
        raise ValueError("Full-locus table row count differs from Stage 1b completion receipt")
    if n_genomes != int(stage1["n_eq_thermoplasmatota_genomes"]) or n_genomes != int(stage1b["n_refined_genome_fastas_verified_against_manifest"]):
        raise ValueError("Full-locus genome scope differs from independent Stage 1/1b receipts")
    require_columns(full_loci_df, ["family_id", "genome", "final_promotion_gate_pass"], "full candidate-family membership")
    gate_false = int((~as_bool(full_loci_df["final_promotion_gate_pass"])).sum())
    return {
        "status": "validated_independent_complete_original_family_assignment_scope",
        "n_screened_genomes": n_genomes,
        "n_original_proteins_in_stage1_scope": int(stage1["n_eq_thermoplasmatota_genes"]),
        "n_candidate_member_locus_rows_complete_membership": len(full_loci_df),
        "n_exact_family_ids": int(full_loci_df["family_id"].nunique()),
        "n_final_promotion_gate_false_rows_retained_as_membership_evidence": gate_false,
        "stage1_summary": str(stage1_summary_path), "stage1b_summary": str(stage1b_summary_path),
        "absence_rule": "Only exact family IDs in this complete assignment table have known absence over the declared 237-genome universe; failed gates never become absence",
    }


def build_family_assay_maps(
    members_df: pd.DataFrame, full_loci_df: pd.DataFrame,
    assay_provenance: dict[str, Any] | None = None,
) -> tuple[list[str], dict[str, set[str]], dict[str, set[str] | None], pd.DataFrame]:
    require_columns(members_df, ["family_id", "genome"], "47-family member evidence")
    require_columns(full_loci_df, ["family_id", "genome"], "full verified loci")
    families = sorted(members_df["family_id"].dropna().astype(str).unique())
    full_universe = set(full_loci_df["genome"].dropna().astype(str))
    member_maps: dict[str, set[str]] = {}
    universes: dict[str, set[str] | None] = {}
    coverage_rows = []
    for family in families:
        dossier_members = set(members_df.loc[members_df["family_id"].astype(str).eq(family), "genome"].dropna().astype(str))
        exact_full = full_loci_df.loc[full_loci_df["family_id"].astype(str).eq(family)]
        if len(exact_full):
            # Complete original membership, not just promoted/verified loci: failed
            # gates can still be genuine members and must never become absences.
            member_maps[family] = set(exact_full["genome"].dropna().astype(str)) | dossier_members
            universes[family] = full_universe
            status = "exact_family_in_receipt_validated_complete_assignment_universe" if assay_provenance else "exact_family_in_fixture_declared_complete_assignment_universe"
        else:
            member_maps[family] = dossier_members
            universes[family] = None
            status = "dossier_only_no_family_specific_absence_universe"
        coverage_rows.append({"family_id": family, "assay_coverage_status": status,
                              "n_full_locus_universe_genomes": len(full_universe) if universes[family] is not None else np.nan,
                              "n_known_member_genomes_all591": len(member_maps[family]),
                              "n_gate_false_member_loci_retained": int((~as_bool(exact_full["final_promotion_gate_pass"])).sum()) if len(exact_full) and "final_promotion_gate_pass" in exact_full else np.nan,
                              "assay_universe_provenance": (assay_provenance or {}).get("status", "unit_fixture_declared_scope_not_production_provenance"),
                              "is_design_lead": int(family in DESIGN_LEAD_FAMILIES)})
    return families, member_maps, universes, pd.DataFrame(coverage_rows)


def conditional_family_tests(
    partition_df: pd.DataFrame,
    members_df: pd.DataFrame,
    full_loci_df: pd.DataFrame,
    assay_provenance: dict[str, Any] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    families, member_maps, universes, coverage = build_family_assay_maps(members_df, full_loci_df, assay_provenance)
    m = len(families)
    eligible = partition_df.loc[partition_df["is_screen_positive_eligible"].eq(1)].copy()
    rows = []
    for outcome_set, outcome_col in [("primary", "primary_DETECTEDLINK"), ("provisional_sensitivity", "provisional_DETECTEDLINK")]:
        for family in families:
            members = member_maps[family]
            universe = universes[family]
            base = {
                "outcome_set": outcome_set, "family_id": family, "fixed_BH_family_size": m,
                "is_design_lead": int(family in DESIGN_LEAD_FAMILIES), "n_screen_positive_eligible_total": len(eligible),
                "n_eligible_detected_links_total": int(eligible[outcome_col].sum()),
                "n_known_member_genomes_eligible": int(eligible["genome"].isin(members).sum()),
            }
            if universe is None:
                rows.append(base | {"assay_coverage_status": "dossier_only_absence_not_assayed", "n_assay_covered_eligible": np.nan,
                                    "n_assay_covered_detected_links": np.nan, "present_positive": np.nan, "present_negative": np.nan,
                                    "absent_positive": np.nan, "absent_negative": np.nan, "odds_ratio": np.nan,
                                    "fisher_exact_pvalue": np.nan, "test_status": "non_estimable_no_verified_absence_denominator"})
                continue
            cohort = eligible.loc[eligible["genome"].isin(universe)].copy()
            y = cohort[outcome_col].to_numpy(dtype=float)
            present = cohort["genome"].isin(members).to_numpy(dtype=int)
            a = int(((present == 1) & (y == 1)).sum())
            b = int(((present == 1) & (y == 0)).sum())
            c = int(((present == 0) & (y == 1)).sum())
            d = int(((present == 0) & (y == 0)).sum())
            base |= {"assay_coverage_status": "receipt_validated_complete_assignment_exact_family" if assay_provenance else "fixture_declared_complete_assignment_exact_family", "n_assay_covered_eligible": len(cohort),
                     "n_assay_covered_detected_links": int(y.sum()), "present_positive": a, "present_negative": b,
                     "absent_positive": c, "absent_negative": d}
            if len(cohort) == 0 or len(np.unique(y)) < 2:
                rows.append(base | {"odds_ratio": np.nan, "fisher_exact_pvalue": np.nan,
                                    "test_status": "non_estimable_outcome_constant_or_zero_assay_denominator"})
                continue
            result = stats.fisher_exact([[a, b], [c, d]])
            odds = float(result.statistic)
            rows.append(base | {"odds_ratio": odds if np.isfinite(odds) else np.nan, "fisher_exact_pvalue": float(result.pvalue),
                                "test_status": "exploratory_conditional_fisher_not_phylogeny_controlled_limited_power"})
    out = pd.DataFrame(rows)
    out["bh_fdr_qvalue_fixed47"] = np.nan
    for outcome_set in out["outcome_set"].unique():
        mask = out["outcome_set"].eq(outcome_set)
        out.loc[mask, "bh_fdr_qvalue_fixed47"] = bh_adjust(out.loc[mask, "fisher_exact_pvalue"], m)
    out["interpretation_boundary"] = "DETECTEDLINK only; nonsignificance is not absence of association or infection"
    return out, coverage


def diagnose_binary_separation(y: np.ndarray, x: np.ndarray) -> dict[str, Any]:
    y, x = np.asarray(y, dtype=float), np.asarray(x, dtype=float)
    if len(y) != len(x) or not np.isin(x, [0, 1]).all() or not np.isin(y, [0, 1]).all():
        raise ValueError("Binary separation diagnostic requires aligned binary x/y")
    if np.ptp(x) == 0:
        return {"status": "constant_predictor_no_contrast_not_separation", "complete": False, "quasi_complete": False}
    if np.ptp(y) == 0:
        return {"status": "constant_outcome_not_predictor_separation", "complete": False, "quasi_complete": False}
    table = {f"x{xv}_y{yv}": int(((x == xv) & (y == yv)).sum()) for xv in [0, 1] for yv in [0, 1]}
    # A level is pure if it contains only one outcome class.
    level_pure = [len(np.unique(y[x == xv])) == 1 for xv in [0, 1]]
    if all(level_pure):
        status, complete, quasi = "complete_separation", True, True
    elif any(level_pure):
        status, complete, quasi = "quasi_complete_separation", False, True
    else:
        status, complete, quasi = "overlap_no_binary_separation", False, False
    return {"status": status, "complete": complete, "quasi_complete": quasi, "contingency": table}


def separation_lp_diagnostic(y: np.ndarray, X: np.ndarray) -> dict[str, Any]:
    """Actual LP separation witness; constants are reported separately by caller."""
    y, X = np.asarray(y, dtype=float), np.asarray(X, dtype=float)
    if not len(y) or np.ptp(y) == 0:
        return {"status": "constant_or_empty_outcome_separation_not_diagnosed", "complete": False, "quasi_complete": False}
    scores_sign = np.where(y == 1, 1.0, -1.0)
    # Complete separation: signed scores >= 1 for every observation.
    A_complete = -(X * scores_sign[:, None])
    b_complete = -np.ones(len(y))
    complete = optimize.linprog(np.zeros(X.shape[1]), A_ub=A_complete, b_ub=b_complete,
                                bounds=[(None, None)] * X.shape[1], method="highs")
    if complete.success:
        margins = X @ complete.x * scores_sign
        return {"status": "complete_separation_LP_witness", "complete": True, "quasi_complete": True,
                "minimum_signed_margin": float(margins.min()), "witness_max_abs_coefficient": float(np.abs(complete.x).max())}
    # Quasi separation: signed scores >= 0, with a nonzero separating vector.
    A_quasi = np.vstack([-(X * scores_sign[:, None]), -(scores_sign @ X)[None, :]])
    b_quasi = np.r_[np.zeros(len(y)), -1.0]
    quasi = optimize.linprog(np.zeros(X.shape[1]), A_ub=A_quasi, b_ub=b_quasi,
                             bounds=[(-1e4, 1e4)] * X.shape[1], method="highs")
    if quasi.success:
        margins = X @ quasi.x * scores_sign
        return {"status": "quasi_complete_separation_LP_witness", "complete": False, "quasi_complete": True,
                "minimum_signed_margin": float(margins.min()), "maximum_signed_margin": float(margins.max()),
                "witness_max_abs_coefficient": float(np.abs(quasi.x).max())}
    return {"status": "no_separation_witness_found_not_proof_of_power", "complete": False, "quasi_complete": False}


def build_design_matrix_audit(
    partition_df: pd.DataFrame,
    members_df: pd.DataFrame,
    full_loci_df: pd.DataFrame,
    outcome_col: str,
    cohort_label: str,
    assay_provenance: dict[str, Any] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    families, member_maps, universes, _ = build_family_assay_maps(members_df, full_loci_df, assay_provenance)
    for lead in DESIGN_LEAD_FAMILIES:
        if lead not in families or universes[lead] is None:
            raise ValueError(f"Design lead lacks exact full-locus assay universe: {lead}")
    df = partition_df.copy() if cohort_label == "full_591_diagnostic" else partition_df.loc[partition_df["is_screen_positive_eligible"].eq(1)].copy()
    phyla = sorted(df["canonical_phylum"].fillna("p__missing").astype(str).unique())
    reference = phyla[0]
    raw = pd.DataFrame({"genome": df["genome"].to_numpy(), "outcome_DETECTEDLINK": df[outcome_col].to_numpy(),
                        "canonical_phylum": df["canonical_phylum"].fillna("p__missing").to_numpy()})
    raw["Intercept"] = 1.0
    for phylum in phyla[1:]:
        raw[f"phylum_{phylum}"] = raw["canonical_phylum"].eq(phylum).astype(float)
    raw["canonical_completeness_frac"] = pd.to_numeric(df["canonical_completeness"], errors="coerce").to_numpy() / 100.0
    raw["log_valid_spacers_plus1_detection_covariate"] = np.log1p(df["n_valid_spacers"].to_numpy(dtype=float))
    full_universe = set(full_loci_df["genome"].astype(str))
    for i, lead in enumerate(DESIGN_LEAD_FAMILIES, 1):
        values = np.where(df["genome"].isin(full_universe), df["genome"].isin(member_maps[lead]).astype(float), np.nan)
        raw[f"lead{i}_{lead}"] = values
    predictor_cols = [c for c in raw.columns if c not in {"genome", "outcome_DETECTEDLINK", "canonical_phylum"}]
    complete = raw["outcome_DETECTEDLINK"].notna() & np.isfinite(raw[predictor_cols].to_numpy(dtype=float)).all(axis=1)
    complete_df = raw.loc[complete].copy()
    y = complete_df["outcome_DETECTEDLINK"].to_numpy(dtype=float)
    X = complete_df[predictor_cols].to_numpy(dtype=float)
    n, p = X.shape
    rank = int(np.linalg.matrix_rank(X)) if n else 0
    singular = np.linalg.svd(X, compute_uv=False) if n else np.array([])
    condition = float(singular[0] / singular[-1]) if len(singular) and singular[-1] > 0 else (float("inf") if len(singular) else np.nan)
    constants = [c for c in predictor_cols if c != "Intercept" and n and np.ptp(complete_df[c].to_numpy(dtype=float)) == 0]
    separation = separation_lp_diagnostic(y, X) if n else {"status": "empty_complete_case_matrix", "complete": False, "quasi_complete": False}
    lead_diagnostics = {}
    for i, lead in enumerate(DESIGN_LEAD_FAMILIES, 1):
        col = f"lead{i}_{lead}"
        lead_diagnostics[col] = diagnose_binary_separation(y, complete_df[col].to_numpy(dtype=float)) if n else {"status": "empty_complete_case_matrix"}
    n_events = int(y.sum()) if n else 0
    n_slopes = p - 1
    audit = {
        "cohort": cohort_label, "outcome_column": outcome_col, "phylum_reference_lexicographic": reference,
        "raw_matrix_rows": len(raw), "complete_case_rows_no_imputation": n,
        "raw_known_outcomes": int(raw["outcome_DETECTEDLINK"].notna().sum()),
        "raw_detected_link_events": int(raw["outcome_DETECTEDLINK"].sum()),
        "complete_case_detected_link_events": n_events, "complete_case_non_events": n - n_events,
        "n_parameters_including_intercept": p, "n_slope_predictors": n_slopes,
        "events_per_slope_predictor_heuristic": n_events / n_slopes if n_slopes else np.nan,
        "EPV_boundary": "EPV<10 is a limited-power heuristic, not mathematical non-estimability",
        "design_matrix_rank": rank, "rank_deficient": bool(rank < p),
        "singular_values": singular.tolist(), "condition_number": condition,
        "constant_predictors_no_contrast": constants,
        "missingness_by_raw_column": {c: int(raw[c].isna().sum()) for c in ["outcome_DETECTEDLINK"] + predictor_cols},
        "separation_LP": separation, "binary_lead_diagnostics": lead_diagnostics,
        "finite_MLE_diagnostic": "no_finite_unique_MLE_if_LP_witness_and_full_rank" if separation["quasi_complete"] else "not_rejected_by_separation_LP_rank_and_outcome_reported_separately",
        "model_fit_performed": False,
        "limited_power_warning": "Few events restrict power even when algebraic rank and separation do not prove non-estimability",
        "CRISPR_covariate_boundary": "Spacer counts are detection-opportunity covariates, not independent mechanistic evidence; Cas/spacer presence is circular for CRISPR-derived DETECTEDLINK",
    }
    return raw, audit


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    df = pd.read_csv(metadata_path, sep="\t", dtype={"sample_id": str, "site": str})
    require_columns(df, ["sample_id", "site", "compartment", "replicate", "trip", "lat", "lon"], "metadata")
    if df["sample_id"].duplicated().any():
        raise ValueError("Metadata contains duplicate sample_id entries")
    if df[["sample_id", "site", "compartment", "trip"]].isna().any().any():
        raise ValueError("Metadata has missing design identifiers")
    return df


def normalize_site(value: Any) -> str:
    text = str(value).strip()
    try:
        number = float(text)
        return str(int(number)) if number.is_integer() else text
    except ValueError:
        return text


def normalize_compartment(value: Any) -> str:
    text = str(value).strip()
    # Audited current Rund map synonym; no sample-ID-derived compartment inference.
    return "RootAdjacent" if text in {"PlantRoot", "RootAdjacent"} else text


def load_field_map(field_map_path: Path) -> pd.DataFrame:
    field = pd.read_csv(field_map_path, sep="\t", dtype={"sample_id": str})
    require_columns(field, ["sample_id", "dune_field", "site_mag", "trip_mag", "compartment_mag"], "audited dune-field map")
    field = field[["sample_id", "dune_field", "site_mag", "trip_mag", "compartment_mag"]].copy()
    if field["sample_id"].duplicated().any():
        raise ValueError("Dune-field map has duplicate sample_id rows")
    return field


def select_primary_environment_design(
    metadata_df: pd.DataFrame,
    field_map_df: pd.DataFrame,
    joint_matrix_samples: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Select design using metadata/availability only, before any outcome values."""
    joint = metadata_df.merge(field_map_df, on="sample_id", how="inner", validate="one_to_one")
    joint = joint.loc[joint["sample_id"].isin(joint_matrix_samples) & joint["dune_field"].notna()].copy()
    site_mismatch = joint["site"].map(normalize_site) != joint["site_mag"].map(normalize_site)
    trip_mismatch = joint["trip"].astype(str) != joint["trip_mag"].astype(str)
    compartment_mismatch = joint["compartment"].map(normalize_compartment) != joint["compartment_mag"].map(normalize_compartment)
    if (site_mismatch | trip_mismatch | compartment_mismatch).any():
        raise ValueError("Metadata and audited field map disagree on site/trip/compartment")
    joint["site"] = joint["site"].map(normalize_site)
    group_counts = joint.groupby(["trip", "compartment"], as_index=False).agg(
        n_samples=("sample_id", "size"), n_distinct_sites=("site", "nunique"), n_dune_fields=("dune_field", "nunique")
    )
    strata = joint.groupby(["trip", "compartment", "dune_field"])["site"].nunique().rename("n_sites").reset_index()
    for row in group_counts.itertuples(index=False):
        sub = strata.loc[strata["trip"].eq(row.trip) & strata["compartment"].eq(row.compartment)]
        group_counts.loc[group_counts["trip"].eq(row.trip) & group_counts["compartment"].eq(row.compartment), "n_permutable_fields"] = int(sub["n_sites"].ge(2).sum())
        group_counts.loc[group_counts["trip"].eq(row.trip) & group_counts["compartment"].eq(row.compartment), "n_singleton_fields"] = int(sub["n_sites"].eq(1).sum())
    if not len(group_counts):
        raise ValueError("No joint eligible metadata with valid dune-field strata")
    selected_row = group_counts.sort_values(["n_distinct_sites", "trip", "compartment"], ascending=[False, True, True], kind="mergesort").iloc[0]
    trip, compartment = selected_row["trip"], selected_row["compartment"]
    selected = joint.loc[joint["trip"].eq(trip) & joint["compartment"].eq(compartment)].sort_values(["site", "sample_id"]).copy()
    field_conflicts = selected.groupby("site")["dune_field"].nunique()
    if field_conflicts.gt(1).any():
        raise ValueError("Selected sites have conflicting dune-field assignments")
    selection = {
        "rule": "before outcomes: largest distinct-site count among joint metadata/matrix-available trip x compartment groups; ties lexicographic trip then compartment",
        "selected_trip": trip, "selected_compartment": compartment,
        "n_joint_metadata_samples": len(joint), "n_joint_metadata_sites": int(joint["site"].nunique()),
        "n_selected_samples": len(selected), "n_selected_sites": int(selected["site"].nunique()),
        "selected_sample_ids": selected["sample_id"].tolist(), "outcome_values_used_for_selection": False,
        "compartment_crosswalk": "current Rund compartment_mag RootAdjacent is the audited synonym of metadata PlantRoot; site/trip require exact normalized matches",
    }
    return group_counts, selected, selection


def aggregate_selected_sites(selected_metadata: pd.DataFrame, host: np.ndarray, virus: np.ndarray) -> pd.DataFrame:
    frame = selected_metadata[["sample_id", "site", "dune_field"]].copy().reset_index(drop=True)
    frame["host"] = np.asarray(host, dtype=float)
    frame["virus"] = np.asarray(virus, dtype=float)
    frame = frame.loc[np.isfinite(frame["host"]) & np.isfinite(frame["virus"])]
    if not len(frame):
        return pd.DataFrame(columns=["site", "dune_field", "host", "virus", "n_samples"])
    return frame.groupby(["site", "dune_field"], as_index=False, sort=True).agg(host=("host", "mean"), virus=("virus", "mean"), n_samples=("sample_id", "size"))


def stratified_permutation_indices(fields: np.ndarray, n_perms: int, seed: int) -> np.ndarray:
    fields = np.asarray(fields).astype(str)
    n = len(fields)
    rng = np.random.default_rng(seed)
    permutations = np.tile(np.arange(n, dtype=np.int32), (n_perms, 1))
    for field in sorted(set(fields)):
        idx = np.flatnonzero(fields == field)
        if len(idx) < 2:
            continue
        for b in range(n_perms):
            permutations[b, idx] = rng.permutation(idx)
    return permutations


def spearman_permutation_result(host: np.ndarray, virus: np.ndarray, permutations: np.ndarray) -> dict[str, Any]:
    host, virus = np.asarray(host, dtype=float), np.asarray(virus, dtype=float)
    if len(host) < 2 or np.ptp(host) == 0 or np.ptp(virus) == 0:
        return {"observed_rho": np.nan, "pvalue": np.nan, "n_extreme_including_ties": np.nan, "status": "non_estimable_zero_variance_or_insufficient_sites"}
    hr = stats.rankdata(host)
    vr = stats.rankdata(virus)
    hc = hr - hr.mean()
    vc = vr - vr.mean()
    observed = float(np.dot(hc, vc) / np.sqrt(np.dot(hc, hc) * np.dot(vc, vc)))
    perm_values = virus[permutations]
    ranks = stats.rankdata(perm_values, axis=1)
    centered = ranks - ranks.mean(axis=1, keepdims=True)
    denominator = np.sqrt(np.sum(centered ** 2, axis=1) * np.dot(hc, hc))
    null = np.divide(centered @ hc, denominator, out=np.full(len(permutations), np.nan), where=denominator > 0)
    extreme = np.isfinite(null) & (np.abs(null) >= abs(observed) - 1e-12)
    n_extreme = int(extreme.sum())
    return {"observed_rho": observed, "pvalue": (1 + n_extreme) / (len(permutations) + 1),
            "n_extreme_including_ties": n_extreme, "status": "exploratory_permutation_performed"}


def site_descriptives(site_df: pd.DataFrame) -> dict[str, Any]:
    if not len(site_df):
        return {"n_sites": 0, "site_mean_spearman_rho": np.nan, "site_jaccard": np.nan}
    h, v = site_df["host"].to_numpy(), site_df["virus"].to_numpy()
    rho = stats.spearmanr(h, v).statistic if np.ptp(h) > 0 and np.ptp(v) > 0 else np.nan
    hp, vp = h > 0, v > 0
    union = int((hp | vp).sum())
    return {"n_sites": len(site_df), "n_host_present_sites": int(hp.sum()), "n_virus_present_sites": int(vp.sum()),
            "n_codetected_sites": int((hp & vp).sum()), "site_mean_spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
            "site_jaccard": int((hp & vp).sum()) / union if union else np.nan}


def evaluate_environment_matched_associations(
    primary_pairs: pd.DataFrame,
    provisional_pairs: pd.DataFrame,
    controls_df: pd.DataFrame,
    mag_df: pd.DataFrame,
    vir_df: pd.DataFrame,
    selected_metadata: pd.DataFrame,
    n_perms: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    site_universe = selected_metadata.groupby(["site", "dune_field"], as_index=False).size().sort_values("site").reset_index(drop=True)
    fields = site_universe["dune_field"].to_numpy(dtype=str)
    field_sizes = pd.Series(fields).value_counts()
    samples = selected_metadata["sample_id"].tolist()
    rows, difference_rows = [], []
    key_counts = provisional_pairs["votu_id"].value_counts().to_dict()
    for analysis_set, pair_set in [("primary", primary_pairs), ("provisional_sensitivity", provisional_pairs)]:
        set_rows_start = len(rows)
        for _, pair in pair_set.sort_values(["canonical_591_host_mag", "target_seq_sha256"]).iterrows():
            host, target, votu = pair["canonical_591_host_mag"], pair["target_id"], pair["votu_id"]
            vkey = votu if votu in vir_df.index else (target if target in vir_df.index else None)
            base = {"analysis_set": analysis_set, "pair_evaluation_role": "candidate", "host_mag": host,
                    "control_for_host": pd.NA, "target_id": target, "votu_id": votu,
                    "target_seq_sha256": pair["target_seq_sha256"], "matched_abundance_key": vkey or pd.NA,
                    "fixed_BH_pair_family_size": len(pair_set), "shared_votu_key_pair_count": int(key_counts.get(votu, 0)),
                    "abundance_resolution_boundary": "shared_vOTU_key_may_represent_multiple_exact_targets_not_independent_viral_resolution"}
            if host not in mag_df.index or vkey is None:
                rows.append(base | {"estimability_status": "non_estimable_missing_host_or_viral_catalog_row_not_zero"})
                continue
            hvec = mag_df.loc[host, samples].to_numpy(dtype=float)
            vvec = vir_df.loc[vkey, samples].to_numpy(dtype=float)
            agg = aggregate_selected_sites(selected_metadata, hvec, vvec)
            agg = site_universe.merge(agg, on=["site", "dune_field"], how="left").dropna(subset=["host", "virus"])
            desc = site_descriptives(agg)
            agg_fields = agg["dune_field"].to_numpy(dtype=str)
            agg_perms = stratified_permutation_indices(agg_fields, n_perms, seed)
            counts = pd.Series(agg_fields).value_counts()
            if not len(agg) or not counts.ge(2).any():
                result = {"observed_rho": np.nan, "pvalue": np.nan, "n_extreme_including_ties": np.nan,
                          "status": "inference_not_performed_no_permutable_dune_field_stratum"}
            else:
                result = spearman_permutation_result(agg["host"].to_numpy(), agg["virus"].to_numpy(), agg_perms)
            # Draw the unblocked sensitivity from its own deterministic stream.
            ur = np.random.default_rng(seed + 1)
            un = np.vstack([ur.permutation(len(agg)) for _ in range(n_perms)])
            unblocked_result = spearman_permutation_result(agg["host"].to_numpy(), agg["virus"].to_numpy(), un)
            rows.append(base | desc | {
                "estimability_status": result["status"], "n_permutations": n_perms, "permutation_seed": seed,
                "field_stratified_permutation_pvalue": result["pvalue"], "n_extreme_including_observed_ties": result["n_extreme_including_ties"],
                "n_permutable_sites": int(sum(counts[counts >= 2])), "n_fixed_singleton_sites": int((counts == 1).sum()),
                "n_permutation_strata": int(len(counts)), "unblocked_site_sensitivity_pvalue": unblocked_result["pvalue"],
                "unblocked_status": "hypothesis_generating_only_not_field_or_spatial_controlled",
                "spatial_boundary": "within-field permutation reduces broad confounding; residual within-field spatial autocorrelation uncontrolled",
            })
            control_match = controls_df.loc[controls_df["candidate_host_mag"].eq(host)]
            if len(control_match) == 1 and pd.notna(control_match.iloc[0].get("matched_alternative_mag")):
                control = str(control_match.iloc[0]["matched_alternative_mag"])
                if control in mag_df.index:
                    cvec = mag_df.loc[control, samples].to_numpy(dtype=float)
                    cagg = aggregate_selected_sites(selected_metadata, cvec, vvec)
                    cagg = site_universe.merge(cagg, on=["site", "dune_field"], how="left").dropna(subset=["host", "virus"])
                    cdesc = site_descriptives(cagg)
                    cperms = stratified_permutation_indices(cagg["dune_field"].to_numpy(dtype=str), n_perms, seed)
                    ccounts = pd.Series(cagg["dune_field"]).value_counts()
                    if not len(cagg) or not ccounts.ge(2).any():
                        cresult = {"observed_rho": np.nan, "pvalue": np.nan, "n_extreme_including_ties": np.nan,
                                   "status": "inference_not_performed_no_permutable_dune_field_stratum"}
                    else:
                        cresult = spearman_permutation_result(cagg["host"].to_numpy(), cagg["virus"].to_numpy(), cperms)
                    cbase = base | {"pair_evaluation_role": "matched_alternative_MAG", "host_mag": control, "control_for_host": host}
                    rows.append(cbase | cdesc | {"estimability_status": cresult["status"], "n_permutations": n_perms,
                        "permutation_seed": seed, "field_stratified_permutation_pvalue": cresult["pvalue"],
                        "n_extreme_including_observed_ties": cresult["n_extreme_including_ties"],
                        "n_permutable_sites": int(sum(ccounts[ccounts >= 2])), "n_fixed_singleton_sites": int((ccounts == 1).sum()),
                        "n_permutation_strata": int(len(ccounts)), "spatial_boundary": "same selected design and seed as candidate; not an unrelated-taxon or infection-absence control"})
                    difference_base = {"analysis_set": analysis_set, "candidate_host_mag": host, "matched_alternative_mag": control,
                        "target_id": target, "votu_id": votu, "target_seq_sha256": pair["target_seq_sha256"],
                        "statistical_status": "descriptive_only_no_correlation_equality_or_host_specificity_test",
                        "null_boundary": "viral-value permutation tests joint independence from both host vectors, not equality of correlations"}
                    if agg[["site", "dune_field"]].equals(cagg[["site", "dune_field"]]) and np.isfinite(result["observed_rho"]) and np.isfinite(cresult["observed_rho"]):
                        difference_rows.append(difference_base | {
                            "same_site_support": True, "observed_candidate_rho": result["observed_rho"],
                            "observed_alternative_rho": cresult["observed_rho"],
                            "observed_rho_difference_candidate_minus_alternative": result["observed_rho"] - cresult["observed_rho"],
                            "interpretation": "observed direct difference only; no preferred-host or correlation-equality inference"})
                    else:
                        difference_rows.append(difference_base | {"same_site_support": False,
                            "observed_rho_difference_candidate_minus_alternative": np.nan,
                            "interpretation": "descriptive difference unavailable without identical site support and finite correlations"})
        # Fixed family includes every frozen pair; NA p-values still count in m.
        set_indices = list(range(set_rows_start, len(rows)))
        for role in ["candidate", "matched_alternative_MAG"]:
            idx = [i for i in set_indices if rows[i]["pair_evaluation_role"] == role]
            q = bh_adjust(pd.Series([rows[i].get("field_stratified_permutation_pvalue", np.nan) for i in idx]), len(pair_set))
            for i, value in zip(idx, q):
                rows[i]["field_stratified_bh_qvalue"] = value
        cand_idx = [i for i in set_indices if rows[i]["pair_evaluation_role"] == "candidate"]
        q = bh_adjust(pd.Series([rows[i].get("unblocked_site_sensitivity_pvalue", np.nan) for i in cand_idx]), len(pair_set))
        for i, value in zip(cand_idx, q):
            rows[i]["unblocked_sensitivity_bh_qvalue"] = value
    # Compute in log space to avoid overflow in the number of within-field permutations.
    permutation_space_log10 = sum(math.lgamma(int(size) + 1) / math.log(10) for size in field_sizes if size >= 2)
    design = {"n_sites_in_selected_design": len(site_universe), "n_permutable_sites": int(field_sizes[field_sizes >= 2].sum()),
              "n_fixed_singleton_sites": int((field_sizes == 1).sum()), "n_permutation_strata": int(len(field_sizes)),
              "n_permutations": n_perms, "permutation_seed": seed, "monte_carlo_p_resolution": 1 / (n_perms + 1),
              "within_field_permutation_space_log10": permutation_space_log10,
              "permutation_space_plus_one_resolution_log10": -math.log10(10 ** permutation_space_log10 + 1) if permutation_space_log10 < 6 else -permutation_space_log10,
              "field_site_counts": {str(k): int(v) for k, v in field_sizes.items()},
              "exchangeability_boundary": "Only whole-site viral values permuted within audited dune-field strata; singleton strata fixed; residual within-field spatial autocorrelation remains uncontrolled"}
    return pd.DataFrame(rows), pd.DataFrame(difference_rows), design


def descriptive_phylum_counts(partition_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for phylum, group in partition_df.groupby("canonical_phylum", dropna=False):
        rows.append({"canonical_phylum": phylum, "n_canonical_MAGs": len(group),
                     "n_zero_spacer_untested": int(group["n_valid_spacers"].eq(0).sum()),
                     "n_screen_positive_eligible": int(group["is_screen_positive_eligible"].sum()),
                     "n_primary_DETECTEDLINK_hosts_descriptive": int(group["is_primary_detected_link_host"].sum()),
                     "n_provisional_DETECTEDLINK_hosts_descriptive": int(group["is_provisional_detected_link_host"].sum()),
                     "model_type": "descriptive_only_not_phylogeny_controlled_test"})
    return pd.DataFrame(rows)


def write_json(path: Path, payload: Any) -> None:
    def convert(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): convert(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [convert(v) for v in value]
        if isinstance(value, (np.integer,)):
            return int(value)
        if isinstance(value, (np.floating, float)):
            return None if not np.isfinite(value) else float(value)
        if isinstance(value, (np.bool_,)):
            return bool(value)
        if isinstance(value, Path):
            return str(value)
        return value
    path.write_text(json.dumps(convert(payload), indent=2, allow_nan=False) + "\n")


def main() -> None:
    t0 = time.time()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exp3-root", type=Path, default=BRANCH_ROOT)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--mag-abund", type=Path)
    parser.add_argument("--vir-tpm", type=Path)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--screening", type=Path)
    parser.add_argument("--evidence", type=Path)
    parser.add_argument("--pairs", type=Path)
    parser.add_argument("--protein-members", type=Path)
    parser.add_argument("--full-loci", type=Path)
    parser.add_argument("--stage1-summary", type=Path, required=True)
    parser.add_argument("--stage1b-summary", type=Path, required=True)
    parser.add_argument("--field-map", type=Path, required=True)
    parser.add_argument("--root-disposition", type=Path, required=True)
    parser.add_argument("--provisional-disposition", type=Path, required=True)
    parser.add_argument("--n-perms", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20261109)
    args = parser.parse_args()
    if args.n_perms < 9999:
        raise ValueError("Root design requires >=9999 permutations")
    root = args.exp3_root
    in_dir, res_dir = root / "inputs", args.results_dir or root / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "metadata": args.metadata or (DEFAULT_METADATA if DEFAULT_METADATA.exists() else in_dir / "metadata.tsv"),
        "mag_abundance": args.mag_abund or (DEFAULT_MAG_ABUND if DEFAULT_MAG_ABUND.exists() else in_dir / "mag_abundance.tsv.gz"),
        "viral_tpm": args.vir_tpm or (DEFAULT_VIR_TPM if DEFAULT_VIR_TPM.exists() else in_dir / "viral_tpm.tsv.gz"),
        "taxonomy": args.taxonomy or in_dir / "taxonomy_reconciliation_all5438.tsv.gz",
        "screening": args.screening or res_dir / "stage1_mag_screening_status.tsv",
        "evidence": args.evidence or res_dir / "stage3_all591_archaeal_mag_host_evidence_v2.tsv",
        "pairs": args.pairs or res_dir / "stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv",
        "protein_members": args.protein_members or in_dir / "exp4_47_member_evidence.tsv",
        "full_loci": args.full_loci or in_dir / "task_a_stage1b_all_candidate_member_loci_verified.tsv.gz",
        "stage1_summary": args.stage1_summary, "stage1b_summary": args.stage1b_summary,
        "field_map": args.field_map, "root_disposition": args.root_disposition,
        "provisional_disposition": args.provisional_disposition,
    }
    for name, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Required input {name}: {path}")
    print("[v4] Loading compact ledgers and actual schemas", flush=True)
    pairs_raw = pd.read_csv(paths["pairs"], sep="\t")
    taxonomy = pd.read_csv(paths["taxonomy"], sep="\t")
    unresolved_targets = load_root_unresolved_targets(paths["root_disposition"])
    provisional_key, provisional_payload = load_provisional_candidate(paths["provisional_disposition"])
    ledger, primary, provisional = resolve_graded_candidate_pairs(pairs_raw, unresolved_targets, {provisional_key})
    primary = join_candidate_taxonomy(primary, taxonomy)
    provisional = join_candidate_taxonomy(provisional, taxonomy)
    ledger.to_csv(res_dir / "stage4_v4_candidate_pairs_graded_ledger.tsv", sep="\t", index=False)
    primary.to_csv(res_dir / "stage4_v4_primary_detected_link_pairs.tsv", sep="\t", index=False)
    provisional.to_csv(res_dir / "stage4_v4_provisional_sensitivity_pairs.tsv", sep="\t", index=False)

    screening = pd.read_csv(paths["screening"], sep="\t")
    evidence = pd.read_csv(paths["evidence"], sep="\t")
    partition = partition_detection_opportunity(screening, evidence, taxonomy, primary, provisional)
    partition.to_csv(res_dir / "stage4_v4_detection_opportunity_partition.tsv", sep="\t", index=False)
    descriptive_phylum_counts(partition).to_csv(res_dir / "stage4_v4_descriptive_phylum_counts.tsv", sep="\t", index=False)

    members = pd.read_csv(paths["protein_members"], sep="\t")
    full_loci = pd.read_csv(paths["full_loci"], sep="\t")
    assay_provenance = validate_assay_provenance(full_loci, paths["stage1_summary"], paths["stage1b_summary"])
    family_results, family_coverage = conditional_family_tests(partition, members, full_loci, assay_provenance)
    family_results.to_csv(res_dir / "stage4_v4_conditional_family_associations.tsv", sep="\t", index=False)
    family_coverage.to_csv(res_dir / "stage4_v4_family_assay_coverage.tsv", sep="\t", index=False)
    eligible_raw, eligible_audit = build_design_matrix_audit(partition, members, full_loci, "primary_DETECTEDLINK", "screen_positive_eligible_primary", assay_provenance)
    full_raw, full_audit = build_design_matrix_audit(partition, members, full_loci, "primary_DETECTEDLINK", "full_591_diagnostic", assay_provenance)
    eligible_raw.to_csv(res_dir / "stage4_v4_eligible_design_matrix_raw.tsv", sep="\t", index=False)
    full_raw.to_csv(res_dir / "stage4_v4_full591_design_matrix_raw_diagnostic.tsv", sep="\t", index=False)
    eligible_hash = compute_sha256(res_dir / "stage4_v4_eligible_design_matrix_raw.tsv")
    full_hash = compute_sha256(res_dir / "stage4_v4_full591_design_matrix_raw_diagnostic.tsv")
    eligible_audit["raw_matrix_sha256"], full_audit["raw_matrix_sha256"] = eligible_hash, full_hash
    eligible_audit["raw_matrix_columns"] = list(eligible_raw.columns)
    full_audit["raw_matrix_columns"] = list(full_raw.columns)
    write_json(res_dir / "stage4_v4_design_matrix_audit.json", {"eligible_primary": eligible_audit, "full591_diagnostic_only": full_audit})

    metadata = load_metadata(paths["metadata"])
    mag_df, vir_df = load_and_orient_abundance_matrices(paths["mag_abundance"], paths["viral_tpm"], set(metadata["sample_id"]))
    shared_samples = [s for s in metadata["sample_id"] if s in mag_df.columns and s in vir_df.columns]
    metrics = compute_abundance_metrics(mag_df, partition["genome"], shared_samples)
    metrics.to_csv(res_dir / "stage4_v4_host_abundance_assay_metrics.tsv", sep="\t", index=False)
    unresolved_hosts = set(ledger.loc[ledger["is_unresolved_root_disposition"].eq(1), "canonical_591_host_mag"])
    controls = find_matched_alternative_controls(
        set(provisional["canonical_591_host_mag"]), set(primary["canonical_591_host_mag"]),
        set(provisional["canonical_591_host_mag"]), unresolved_hosts, partition, metrics,
    )
    controls.to_csv(res_dir / "stage4_v4_matched_alternative_MAG_controls.tsv", sep="\t", index=False)
    field_map = load_field_map(paths["field_map"])
    group_counts, selected, selection = select_primary_environment_design(metadata, field_map, set(shared_samples))
    group_counts.to_csv(res_dir / "stage4_v4_environment_candidate_group_counts.tsv", sep="\t", index=False)
    selected[["sample_id", "site", "compartment", "replicate", "trip", "lat", "lon", "dune_field"]].to_csv(
        res_dir / "stage4_v4_selected_environment_samples.tsv", sep="\t", index=False)
    abundance, differences, env_design = evaluate_environment_matched_associations(
        primary, provisional, controls, mag_df, vir_df, selected, args.n_perms, args.seed)
    abundance.to_csv(res_dir / "stage4_v4_environment_matched_abundance.tsv", sep="\t", index=False)
    differences.to_csv(res_dir / "stage4_v4_candidate_control_correlation_differences.tsv", sep="\t", index=False)
    write_json(res_dir / "stage4_v4_environment_design.json", {"selection": selection, "permutation_design": env_design})

    tier1 = ledger["evidence_tier"].eq("Tier1_verified_direct")
    summary = {
        "stage": "stage4_statistical_audit_v4", "timestamp_UTC": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
        "elapsed_seconds": time.time() - t0, "runtime_command": " ".join(sys.argv), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "v3_disposition": "preserved_unchanged_invalid_inference_baseline_exploratory_descriptives_only_v4_supersedes_v3_p_and_q_values",
        "input_files": {name: {"path": str(path), "sha256": compute_sha256(path)} for name, path in paths.items()},
        "code": {"path": str(Path(__file__).resolve()), "sha256": compute_sha256(Path(__file__).resolve())},
        "candidate_counts_actual": {"raw_rows": len(ledger), "tier1_raw_rows": int(tier1.sum()),
            "alias_duplicate_rows": int(ledger["is_alias_duplicate"].sum()), "unresolved_rows": int(ledger["is_unresolved_root_disposition"].sum()),
            "boundary_sensitive_provisional_rows": int(ledger["is_boundary_sensitive_provisional"].sum()), "primary_pairs": len(primary),
            "primary_hosts": int(primary["canonical_591_host_mag"].nunique()), "provisional_pairs": len(provisional),
            "provisional_hosts": int(provisional["canonical_591_host_mag"].nunique())},
        "provisional_candidate_root_disposition": {key: provisional_payload[key] for key in [
            "host", "target", "target_sha256", "disposition", "canonical_spacer_sequence",
            "spacer_target_interval_1based_inclusive", "catalog_viral_interval", "fresh_checkv_viral_interval",
            "fresh_checkv_host_interval", "interpretation"]},
        "detection_counts_actual": partition["detection_opportunity_category"].value_counts().to_dict() | {
            "screen_positive_eligible": int(partition["is_screen_positive_eligible"].sum()),
            "primary_detected_hosts": int(partition["is_primary_detected_link_host"].sum()),
            "provisional_detected_hosts": int(partition["is_provisional_detected_link_host"].sum())},
        "family_tests_actual": {"frozen_families": int(family_coverage.shape[0]),
            "assay_provenance": assay_provenance,
            "full_locus_universe_genomes": int(full_loci["genome"].nunique()),
            "full_absence_universe_families": int(family_coverage["assay_coverage_status"].str.contains("complete_assignment_universe").sum()),
            "dossier_only_nonestimable_absence_families": int(family_coverage["assay_coverage_status"].str.startswith("dossier_only").sum()),
            "rows_by_outcome_set": family_results.groupby("outcome_set").size().to_dict(),
            "estimable_tests_by_outcome_set": family_results.groupby("outcome_set")["fisher_exact_pvalue"].count().to_dict()},
        "design_matrix_actual": {"eligible": {k: eligible_audit[k] for k in ["raw_matrix_rows", "complete_case_rows_no_imputation", "complete_case_detected_link_events", "design_matrix_rank", "n_parameters_including_intercept", "rank_deficient", "constant_predictors_no_contrast", "separation_LP"]},
                                  "full591_diagnostic": {k: full_audit[k] for k in ["raw_matrix_rows", "complete_case_rows_no_imputation", "complete_case_detected_link_events", "design_matrix_rank", "n_parameters_including_intercept", "rank_deficient"]}},
        "matched_controls_actual": {"candidate_hosts": len(controls), "matched": int(controls["matched_alternative_mag"].notna().sum()),
            "unmatched_NA": int(controls["matched_alternative_mag"].isna().sum())},
        "environment_design_actual": {"selection": selection, "permutations": env_design,
            "abundance_rows": len(abundance), "descriptive_candidate_alternative_difference_rows_no_pvalues": len(differences)},
        "scientific_boundaries": [
            "All v4 abundance and family tests are exploratory, not confirmatory or fully phylogeny/spatial controlled.",
            "DETECTEDLINK is conditional on CRISPR screening and assay coverage; zero spacers and missing membership assays are not biological negatives.",
            "Nonsignificant Fisher/permutation results do not establish absence of association or infection.",
            "Separate historical MAG/vOTU catalogs are not joint competitive coverage; shared vOTU keys do not resolve independent exact viruses.",
            "The second Thermoplasmatota candidate is boundary-sensitive provisional only: its exact spacer lies outside the catalog viral core but inside the longer CheckV interval; the whole 67,924-bp contig is not claimed viral.",
            "The older site_landforms.tsv was not used; the current Rund matched_profiles dune-field map is exact-joined and checked against metadata.",
            "Within-dune-field site permutation leaves residual within-field spatial autocorrelation uncontrolled; unblocked results are sensitivity only.",
            "Candidate-minus-alternative rho differences are descriptive only: viral permutation is a joint-independence null, not a general equality-of-correlations or host-preference test.",
            f"The {assay_provenance['n_screened_genomes']}-genome family assay universe is supported by independent Stage 1/1b complete-assignment receipts; all {assay_provenance['n_final_promotion_gate_false_rows_retained_as_membership_evidence']:,} gate-failed member rows remain membership evidence and never become absences.",
            "EPV is a power heuristic only; rank, constants and actual LP separation witnesses are reported separately; no overparameterized model was fit.",
        ],
    }
    write_json(res_dir / "stage4_v4_statistical_audit_summary.json", summary)
    output_paths = sorted(p for p in res_dir.glob("stage4_v4_*") if p.is_file() and p.name != "stage4_v4_output_manifest.tsv")
    manifest = pd.DataFrame([{"file": p.name, "bytes": p.stat().st_size, "sha256": compute_sha256(p)} for p in output_paths])
    manifest.to_csv(res_dir / "stage4_v4_output_manifest.tsv", sep="\t", index=False)
    print("[v4] Complete", json.dumps(summary["candidate_counts_actual"]), flush=True)
    print("[v4] Detection", json.dumps(summary["detection_counts_actual"]), flush=True)
    print("[v4] Environment", selection["selected_trip"], selection["selected_compartment"], selection["n_selected_sites"], flush=True)


if __name__ == "__main__":
    main()
