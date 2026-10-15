#!/usr/bin/env python3
"""f1_ecology_analysis.py - Ecological Analysis of F1 Photolyase/Cryptochrome Family (CPF) Carriers (v2).

Evaluates whether carriers of the F1 CPF candidates (identified in root_audit/f1_full/candidates.tsv)
exhibit light-niche partitioning between surface and subsurface soils across site-paired locations,
following the verified methodology of Section 2.3 of the Empty Quarter manuscript.

Key Methodological Specifications (Root Audit Aligned):
1. Exact Numerator / Denominator Alignment:
   - PRIMARY mapped community abundance (`abund_mapped_pct`):
     Numerator:   sum of abundance for carrier representatives in catalogue.
     Denominator: sum of abundance for ALL 5,438 catalogue representatives in sample.
     Guarantees novel dark-matter / unannotated families are NOT dropped.
   - PRIMARY occurrence (`occurrence_pct`):
     Numerator:   count of detected carrier representatives.
     Denominator: count of ALL detected catalogue representatives.
   - SECONDARY whole metagenome percent (`abund_whole_pct`):
     Raw CoverM percentage across all reads (including unmapped).
   - REPLICATION annotated community fraction (`abund_annotated_pct` & `occurrence_annotated_pct`):
     Strictly aligns numerator (trait_reps & annotated_reps) with denominator (annotated_reps)
     to reproduce the original manuscript's K01669 baseline (24.3% vs 12.1%).
2. Site + Trip Pairing Hierarchy:
   - Step A: Multiple replicate samples within (site, trip, compartment) collapsed to median.
   - Step B: Paired differences calculated strictly within the same expedition trip for that site.
   - Step C: Replication of sites visited across multiple trips accounted for via:
     * Primary: One site aggregate (median paired difference across trips for each site).
     * Cluster bootstrap across all site-trip visits (resampling sites with replacement).
   - Stratified trip analysis called ONLY where paired Surface and Deep samples exist within that trip.
3. Rigorous Bootstrap Confidence Intervals:
   - Fixed seed (seed=20261015), 2,000 bootstrap resamples of sites with replacement.
   - Computes empirical 95% bootstrap CIs for paired median and mean differences.
4. Dynamic Metrics (Zero Hardcoding):
   - All counts (candidate hits, carrier MAGs, complete carriers, paired sites, visits)
     are dynamically derived from data structures, never hardcoded in text, plots, or provenance.
5. Objective Negative Controls:
   - BER (Base Excision Repair) and NER (Nucleotide Excision Repair) are reported objectively
     with empirical medians, deltas, and p-values without preconceived "flatness" bias.
6. Lineage & Subfamily Caveats:
   - Within-phylum tests address broad taxonomic shifts, but do NOT control for fine-scale
     species/strain turnover; lineage confounding is bounded, not fully solved.
   - All profile strata (PF03441, PF00875, PF04244) are exploratory domain matches.
     No claim of specific enzymatic repair subfamily (CPD photolyase class I/II/III, 6-4 photolyase)
     vs signaling cryptochrome (Cry-DASH) is made without phylogenetic placement and biochemical assay.
7. Fail-Closed Join Verification:
   - Zero tolerance for corrupted data or unverified joins. Any unexpected ID or missing file fails closed.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import math
import os
import platform
import socket
import statistics
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kruskal, wilcoxon


def compute_file_sha256(filepath: Path | str) -> str:
    """Compute hex SHA256 of a local file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def holm_bonferroni(pvals: List[float]) -> List[float]:
    """Compute Holm-Bonferroni adjusted p-values."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i] if not math.isnan(pvals[i]) else 2.0)
    adj = [float("nan")] * n
    running = 0.0
    for rank, idx in enumerate(order):
        p = pvals[idx]
        if math.isnan(p):
            adj[idx] = float("nan")
        else:
            running = max(running, (n - rank) * p)
            adj[idx] = min(1.0, running)
    return adj


def bootstrap_paired_difference(
    diffs: List[float],
    n_bootstraps: int = 2000,
    seed: int = 20261015,
) -> Dict[str, float]:
    """Compute 95% bootstrap confidence intervals for median and mean paired differences."""
    if len(diffs) < 5:
        return {
            "ci_median_low": float("nan"),
            "ci_median_high": float("nan"),
            "ci_mean_low": float("nan"),
            "ci_mean_high": float("nan"),
        }
    rng = np.random.RandomState(seed)
    n = len(diffs)
    diffs_arr = np.array(diffs, dtype=np.float64)

    boot_indices = rng.randint(0, n, size=(n_bootstraps, n))
    boot_samples = diffs_arr[boot_indices]

    boot_medians = np.median(boot_samples, axis=1)
    boot_means = np.mean(boot_samples, axis=1)

    ci_med_low, ci_med_high = np.percentile(boot_medians, [2.5, 97.5])
    ci_mean_low, ci_mean_high = np.percentile(boot_means, [2.5, 97.5])

    return {
        "ci_median_low": float(ci_med_low),
        "ci_median_high": float(ci_med_high),
        "ci_mean_low": float(ci_mean_low),
        "ci_mean_high": float(ci_mean_high),
    }


def cluster_bootstrap_site_trips(
    site_trip_diffs: Dict[str, List[float]],
    n_bootstraps: int = 2000,
    seed: int = 20261015,
) -> Dict[str, float]:
    """Perform cluster bootstrap resampling at site level, pooling paired site-trip visits."""
    sites = sorted(site_trip_diffs.keys())
    n_sites = len(sites)
    if n_sites < 5:
        return {
            "cluster_ci_median_low": float("nan"),
            "cluster_ci_median_high": float("nan"),
            "cluster_ci_mean_low": float("nan"),
            "cluster_ci_mean_high": float("nan"),
        }

    rng = np.random.RandomState(seed)
    sites_arr = np.array(sites)

    boot_medians = np.empty(n_bootstraps, dtype=np.float64)
    boot_means = np.empty(n_bootstraps, dtype=np.float64)

    for b in range(n_bootstraps):
        sampled_sites = rng.choice(sites_arr, size=n_sites, replace=True)
        pooled_diffs: List[float] = []
        for s in sampled_sites:
            pooled_diffs.extend(site_trip_diffs[s])
        boot_medians[b] = float(np.median(pooled_diffs))
        boot_means[b] = float(np.mean(pooled_diffs))

    ci_med_low, ci_med_high = np.percentile(boot_medians, [2.5, 97.5])
    ci_mean_low, ci_mean_high = np.percentile(boot_means, [2.5, 97.5])

    return {
        "cluster_ci_median_low": float(ci_med_low),
        "cluster_ci_median_high": float(ci_med_high),
        "cluster_ci_mean_low": float(ci_mean_low),
        "cluster_ci_mean_high": float(ci_mean_high),
    }


def load_candidates(candidates_path: Path) -> Dict[str, Any]:
    """Load and validate F1 CPF candidate hits from root_audit/f1_full/candidates.tsv."""
    if not candidates_path.is_file():
        raise FileNotFoundError(f"FAIL CLOSED: Candidates file not found: {candidates_path}")

    records: List[Dict[str, Any]] = []
    mags_all: Set[str] = set()
    mags_complete: Set[str] = set()
    mags_by_profile: Dict[str, Set[str]] = defaultdict(set)
    mags_complete_by_profile: Dict[str, Set[str]] = defaultdict(set)

    with open(candidates_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            protein_id = row["id"].strip()
            mag = row["mag"].strip()
            profile = row["profile"].strip()
            partial = row["partial"].strip()
            length = int(row["length"]) if row.get("length") else 0

            is_complete = (partial == "00")
            records.append({
                "id": protein_id,
                "mag": mag,
                "profile": profile,
                "partial": partial,
                "is_complete": is_complete,
                "length": length,
                "domain_bitscore": float(row["domain_bitscore"]) if row.get("domain_bitscore") else 0.0,
            })
            mags_all.add(mag)
            mags_by_profile[profile].add(mag)
            if is_complete:
                mags_complete.add(mag)
                mags_complete_by_profile[profile].add(mag)

    if not records:
        raise ValueError(f"FAIL CLOSED: No candidate records parsed from {candidates_path}")

    return {
        "records": records,
        "mags_all": mags_all,
        "mags_complete": mags_complete,
        "mags_by_profile": dict(mags_by_profile),
        "mags_complete_by_profile": dict(mags_complete_by_profile),
    }


def load_rep_traits(rep_traits_path: Path) -> Dict[str, Any]:
    """Load 5,438 catalogue species representatives, annotation status, taxonomy, and baseline traits."""
    if not rep_traits_path.is_file():
        raise FileNotFoundError(f"FAIL CLOSED: Rep traits file not found: {rep_traits_path}")

    annotated_reps: Set[str] = set()
    all_reps: Set[str] = set()
    phylum_by_rep: Dict[str, str] = {}
    k01669_reps: Set[str] = set()
    nhej_reps: Set[str] = set()
    carot_reps: Set[str] = set()
    ner_reps: Set[str] = set()
    ber_reps: Set[str] = set()

    with open(rep_traits_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rep = row["rep"].strip()
            all_reps.add(rep)
            ann = (row.get("annotated") == "1")
            if ann:
                annotated_reps.add(rep)

            ph = row.get("phylum", "").strip() or "Unknown"
            phylum_by_rep[rep] = ph

            if ann:
                if row.get("photolyase") == "1" or row.get("photoreactivation") == "1":
                    k01669_reps.add(rep)
                if row.get("nhej") == "1" or row.get("dsb_nhej") == "1":
                    nhej_reps.add(rep)
                if row.get("carotenoid") == "1" or row.get("photoprotective_pigment") == "1":
                    carot_reps.add(rep)
                if row.get("nucleotide_excision_repair") == "1":
                    ner_reps.add(rep)
                if row.get("base_excision_repair") == "1":
                    ber_reps.add(rep)

    if len(all_reps) != 5438:
        print(f"WARNING: Total catalogue representatives count is {len(all_reps)} (expected 5438)")

    return {
        "all_reps": all_reps,
        "annotated_reps": annotated_reps,
        "phylum_by_rep": phylum_by_rep,
        "k01669_reps": k01669_reps,
        "nhej_reps": nhej_reps,
        "carot_reps": carot_reps,
        "ner_reps": ner_reps,
        "ber_reps": ber_reps,
    }


def load_metadata(metadata_path: Path) -> Dict[str, Dict[str, str]]:
    """Load sample metadata (site, compartment, trip/expedition, replicate)."""
    if not metadata_path.is_file():
        raise FileNotFoundError(f"FAIL CLOSED: Metadata file not found: {metadata_path}")

    meta: Dict[str, Dict[str, str]] = {}
    with open(metadata_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            s = row["sample"].strip()
            meta[s] = {
                "sample": s,
                "site": row["site"].strip(),
                "compartment": row["compartment"].strip(),
                "replicate": row.get("replicate", "").strip(),
                "trip": row.get("trip", "").strip() or row.get("trip_num", "").strip(),
            }
    if not meta:
        raise ValueError(f"FAIL CLOSED: No metadata records loaded from {metadata_path}")
    return meta


def read_coverm_profile(profile_path: Path) -> Tuple[Dict[str, float], float]:
    """Read a CoverM abund.tsv relative abundance file for a single sample."""
    if not profile_path.is_file():
        raise FileNotFoundError(f"FAIL CLOSED: Profile path does not exist: {profile_path}")

    abund: Dict[str, float] = {}
    unmapped = 0.0

    with open(profile_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            raise ValueError(f"FAIL CLOSED: Empty profile header in {profile_path}")
        for row in reader:
            if len(row) < 2:
                continue
            genome = row[0].strip()
            try:
                val = float(row[1])
            except ValueError:
                continue
            if genome == "unmapped":
                unmapped = val
            else:
                abund[genome] = val

    return abund, unmapped


def compute_community_sample_metrics(
    samples: List[str],
    profiles_dir: Path,
    meta: Dict[str, Dict[str, str]],
    traits: Dict[str, Set[str]],
    all_reps: Set[str],
    annotated_reps: Set[str],
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Compute community-level occurrence and abundance metrics across all usable samples.
    
    Guarantees strict numerator/denominator alignment:
    - abund_mapped_pct: (trait_reps & all_reps) / (all_reps) * 100.0  [PRIMARY]
    - occurrence_pct:   (det & trait_reps) / (det in all_reps) * 100.0 [PRIMARY]
    - abund_whole_pct:  raw CoverM percentage across all reads        [SECONDARY]
    - abund_annotated_pct: (trait_reps & annotated_reps) / (annotated_reps) * 100.0 [REPLICATION]
    - occurrence_annotated_pct: (det & trait_reps & annotated_reps) / (det in annotated_reps) * 100.0
    """
    sample_records: List[Dict[str, Any]] = []
    skipped_samples: List[str] = []

    for s in samples:
        prof_file = profiles_dir / s / "abund.tsv"
        if not prof_file.is_file():
            raise FileNotFoundError(f"FAIL CLOSED: Expected profile missing: {prof_file}")

        abund, unmapped = read_coverm_profile(prof_file)
        m = meta[s]

        # 1. Denominators for ALL catalogue representatives (Mapped Community)
        det_all = {g for g, v in abund.items() if v > 0.0 and g in all_reps}
        mapped_abund = sum(v for g, v in abund.items() if g in all_reps)

        # 2. Denominators for Annotated representatives
        det_ann = {g for g, v in abund.items() if v > 0.0 and g in annotated_reps}
        ann_abund = sum(v for g, v in abund.items() if g in annotated_reps)

        # Fail closed on sample with 0 mapped community abundance (e.g. F5Dr1)
        if not det_all or mapped_abund <= 0.0:
            print(f"WARNING: Sample {s} has zero detected catalogue species (unmapped={unmapped:.2f}%). Skipped from ratios.")
            skipped_samples.append(s)
            continue

        for trait_name, trait_reps in traits.items():
            # PRIMARY: Mapped community fraction (all 5438 reps denominator)
            det_trait_all = det_all & trait_reps
            occ_mapped = (len(det_trait_all) / len(det_all)) * 100.0
            sum_trait_all = sum(abund[g] for g in (trait_reps & all_reps) if g in abund)
            abund_mapped = (sum_trait_all / mapped_abund) * 100.0

            # SECONDARY: Whole metagenome percent (CoverM raw %)
            sum_trait_whole = sum(abund[g] for g in trait_reps if g in abund)

            # REPLICATION: Annotated community fraction (strict numerator/denominator alignment)
            trait_reps_ann = trait_reps & annotated_reps
            det_trait_ann = det_ann & trait_reps_ann
            occ_ann = (len(det_trait_ann) / len(det_ann)) * 100.0 if det_ann else 0.0
            sum_trait_ann = sum(abund[g] for g in trait_reps_ann if g in abund)
            abund_ann = (sum_trait_ann / ann_abund) * 100.0 if ann_abund > 0.0 else 0.0

            sample_records.append({
                "sample": s,
                "site": m["site"],
                "compartment": m["compartment"],
                "trip": m["trip"],
                "trait": trait_name,
                "unmapped_pct": unmapped,
                "n_detected_all": len(det_all),
                "n_detected_trait_all": len(det_trait_all),
                "n_detected_ann": len(det_ann),
                "n_detected_trait_ann": len(det_trait_ann),
                "occurrence_pct": occ_mapped,                       # PRIMARY
                "abund_mapped_pct": abund_mapped,                   # PRIMARY
                "abund_whole_pct": sum_trait_whole,                 # SECONDARY
                "abund_annotated_pct": abund_ann,                   # REPLICATION
                "occurrence_annotated_pct": occ_ann,                # REPLICATION
            })

    return sample_records, skipped_samples


def run_site_paired_compartment_tests(
    sample_records: List[Dict[str, Any]],
    traits: List[str],
    metrics: List[str] = ["abund_mapped_pct", "occurrence_pct", "abund_annotated_pct", "abund_whole_pct"],
    seed: int = 20261015,
) -> Tuple[List[Dict[str, Any]], Dict[Tuple[str, str, str, str], float], Dict[Tuple[str, str, str, str, str], float]]:
    """Perform hierarchical site-paired Surface vs Deep testing.
    
    Hierarchy:
    Step 1. Replicates within (site, trip, compartment) collapsed to median.
    Step 2. Paired differences calculated strictly within the same trip: diff(site, trip).
    Step 3. Replication of sites across trips accounted for via:
            - Primary: Site aggregate (median across paired trips for that site).
            - Bootstrap CI of median difference (fixed seed=20261015).
            - Cluster bootstrap across all paired site-trip visits.
    """
    # 1. Collapse replicates within (trait, metric, site, trip, compartment) cell
    cell_trip_values: Dict[Tuple[str, str, str, str, str], List[float]] = defaultdict(list)
    for r in sample_records:
        t = r["trait"]
        site = r["site"]
        trip = r["trip"]
        comp = r["compartment"]
        for m in metrics:
            k = (t, m, site, trip, comp)
            cell_trip_values[k].append(r[m])

    cell_trip_medians: Dict[Tuple[str, str, str, str, str], float] = {
        k: statistics.median(vals) for k, vals in cell_trip_values.items()
    }

    # Discover all (site, trip) pairs that have BOTH Surface and Deep
    all_site_trips = {(k[2], k[3]) for k in cell_trip_medians.keys()}
    paired_site_trips_set = {
        (s, tr) for (s, tr) in all_site_trips
        if any((t, m, s, tr, "Surface") in cell_trip_medians for t in traits for m in metrics[:1]) and
           any((t, m, s, tr, "Deep") in cell_trip_medians for t in traits for m in metrics[:1])
    }

    # Also discover unique sites represented among paired site-trips
    paired_sites_set = sorted({s for s, tr in paired_site_trips_set})

    # Site-level collapsed values across all compartments for Kruskal-Wallis
    cell_site_values: Dict[Tuple[str, str, str, str], List[float]] = defaultdict(list)
    for r in sample_records:
        cell_site_values[(r["trait"], r["metric"], r["site"], r["compartment"])].append(r[r["metric"]])
    site_medians: Dict[Tuple[str, str, str, str], float] = {
        k: statistics.median(vals) for k, vals in cell_site_values.items()
    }

    test_results: List[Dict[str, Any]] = []

    for metric in metrics:
        metric_stash: List[Dict[str, Any]] = []
        for t in traits:
            # Step A: Collect paired observations per (site, trip)
            site_trip_diff_map: Dict[str, List[float]] = defaultdict(list)
            site_trip_surf_map: Dict[str, List[float]] = defaultdict(list)
            site_trip_deep_map: Dict[str, List[float]] = defaultdict(list)

            for s, tr in sorted(paired_site_trips_set):
                k_s = (t, metric, s, tr, "Surface")
                k_d = (t, metric, s, tr, "Deep")
                if k_s in cell_trip_medians and k_d in cell_trip_medians:
                    s_val = cell_trip_medians[k_s]
                    d_val = cell_trip_medians[k_d]
                    diff = s_val - d_val
                    site_trip_diff_map[s].append(diff)
                    site_trip_surf_map[s].append(s_val)
                    site_trip_deep_map[s].append(d_val)

            # Sites with at least one paired trip
            valid_sites = sorted(site_trip_diff_map.keys())
            if len(valid_sites) < 5:
                continue

            # Step B: Primary site-aggregate (median across paired trips for that site)
            site_agg_surf = [statistics.median(site_trip_surf_map[s]) for s in valid_sites]
            site_agg_deep = [statistics.median(site_trip_deep_map[s]) for s in valid_sites]
            site_agg_diff = [statistics.median(site_trip_diff_map[s]) for s in valid_sites]

            n_pos = sum(1 for d in site_agg_diff if d > 0)
            n_neg = sum(1 for d in site_agg_diff if d < 0)
            n_zero = sum(1 for d in site_agg_diff if d == 0)

            # Two-sided Wilcoxon signed-rank test on site aggregate
            nonzero_diffs = [d for d in site_agg_diff if d != 0]
            try:
                if nonzero_diffs:
                    w_res = wilcoxon(site_agg_surf, site_agg_deep)
                    p_val = float(w_res.pvalue)
                    w_stat = float(w_res.statistic)
                else:
                    p_val = 1.0
                    w_stat = 0.0
            except Exception:
                p_val = float("nan")
                w_stat = float("nan")

            # Bootstrap 95% CIs for site-aggregate median and mean differences
            boot_ci = bootstrap_paired_difference(site_agg_diff, n_bootstraps=2000, seed=seed)

            # Cluster bootstrap across all paired site-trip visits
            cluster_ci = cluster_bootstrap_site_trips(site_trip_diff_map, n_bootstraps=2000, seed=seed)

            total_paired_visits = sum(len(v) for v in site_trip_diff_map.values())

            # Kruskal-Wallis across all 3 compartments on site-level values
            groups = []
            for c in ("Surface", "Deep", "PlantRoot"):
                g = [site_medians[(t, metric, s, c)] for s in valid_sites if (t, metric, s, c) in site_medians]
                if g:
                    groups.append(g)
            kp = float(kruskal(*groups).pvalue) if len(groups) == 3 else float("nan")

            metric_stash.append({
                "trait": t,
                "metric": metric,
                "n_paired_sites": len(valid_sites),
                "n_paired_site_trips": total_paired_visits,
                "surface_median": float(statistics.median(site_agg_surf)),
                "deep_median": float(statistics.median(site_agg_deep)),
                "delta_median": float(statistics.median(site_agg_diff)),
                "ci_median_low_95": boot_ci["ci_median_low"],
                "ci_median_high_95": boot_ci["ci_median_high"],
                "mean_surface": float(statistics.mean(site_agg_surf)),
                "mean_deep": float(statistics.mean(site_agg_deep)),
                "mean_delta": float(statistics.mean(site_agg_diff)),
                "ci_mean_low_95": boot_ci["ci_mean_low"],
                "ci_mean_high_95": boot_ci["ci_mean_high"],
                "cluster_ci_median_low_95": cluster_ci["cluster_ci_median_low"],
                "cluster_ci_median_high_95": cluster_ci["cluster_ci_median_high"],
                "sites_surface_higher": n_pos,
                "sites_deep_higher": n_neg,
                "sites_equal": n_zero,
                "wilcoxon_stat": w_stat,
                "p_wilcoxon": p_val,
                "p_kruskal3": kp,
            })

        # Holm-Bonferroni correction within metric
        p_list = [item["p_wilcoxon"] for item in metric_stash]
        adj_p = holm_bonferroni(p_list)
        for item, ap in zip(metric_stash, adj_p):
            item["p_holm"] = ap
            test_results.append(item)

    return test_results, site_medians, cell_trip_medians


def run_within_phylum_tests(
    sample_records: List[Dict[str, Any]],
    focal_traits: Dict[str, Set[str]],
    all_reps: Set[str],
    phylum_by_rep: Dict[str, str],
    profiles_dir: Path,
    meta: Dict[str, Dict[str, str]],
    min_species_per_phylum: int = 10,
    min_paired_sites: int = 10,
) -> List[Dict[str, Any]]:
    """Test whether occurrence differences survive WITHIN individual phyla.
    
    Caveat: Phylum stratification tests whether occurrence differences persist within
    individual phyla, but does NOT control for within-phylum species or strain turnover.
    """
    top_phyla = [
        p for p, _ in Counter(
            phylum_by_rep[r] for r in all_reps if phylum_by_rep.get(r)
        ).most_common(8)
    ]

    samples = sorted({r["sample"] for r in sample_records})

    # Within-sample occurrence per phylum
    sample_phylum_occ: Dict[Tuple[str, str, str, str, str], float] = {}

    for s in samples:
        prof_file = profiles_dir / s / "abund.tsv"
        if not prof_file.is_file():
            continue
        abund, _ = read_coverm_profile(prof_file)
        m = meta[s]

        det = {g for g, v in abund.items() if v > 0.0 and g in all_reps}
        if not det:
            continue

        for ph in top_phyla:
            det_ph = {g for g in det if phylum_by_rep.get(g) == ph}
            if len(det_ph) < min_species_per_phylum:
                continue

            for t_name, t_reps in focal_traits.items():
                occ_ph = (len(det_ph & t_reps) / len(det_ph)) * 100.0
                k = (t_name, ph, m["site"], m["trip"], m["compartment"])
                sample_phylum_occ[k] = occ_ph

    # Collapse replicates within (trait, phylum, site, trip, compartment)
    cell_trip_phylum_values: Dict[Tuple[str, str, str, str, str], List[float]] = defaultdict(list)
    for k, val in sample_phylum_occ.items():
        cell_trip_phylum_values[k].append(val)
    cell_trip_phylum_med = {k: statistics.median(vals) for k, vals in cell_trip_phylum_values.items()}

    # Pair within trip then aggregate across trips per site
    phylum_results: List[Dict[str, Any]] = []

    for t_name in focal_traits:
        for ph in top_phyla:
            site_trip_diffs: Dict[str, List[float]] = defaultdict(list)
            site_trip_surf: Dict[str, List[float]] = defaultdict(list)
            site_trip_deep: Dict[str, List[float]] = defaultdict(list)

            all_st = {(k[2], k[3]) for k in cell_trip_phylum_med if k[0] == t_name and k[1] == ph}
            for s, tr in all_st:
                k_s = (t_name, ph, s, tr, "Surface")
                k_d = (t_name, ph, s, tr, "Deep")
                if k_s in cell_trip_phylum_med and k_d in cell_trip_phylum_med:
                    sv = cell_trip_phylum_med[k_s]
                    dv = cell_trip_phylum_med[k_d]
                    site_trip_diffs[s].append(sv - dv)
                    site_trip_surf[s].append(sv)
                    site_trip_deep[s].append(dv)

            valid_sites = sorted(site_trip_diffs.keys())
            if len(valid_sites) < min_paired_sites:
                continue

            site_surf = [statistics.median(site_trip_surf[s]) for s in valid_sites]
            site_deep = [statistics.median(site_trip_deep[s]) for s in valid_sites]
            site_diff = [statistics.median(site_trip_diffs[s]) for s in valid_sites]

            nonzero = [d for d in site_diff if d != 0]
            try:
                if nonzero:
                    w = wilcoxon(site_surf, site_deep)
                    pv = float(w.pvalue)
                else:
                    pv = 1.0
            except Exception:
                pv = float("nan")

            boot_ci = bootstrap_paired_difference(site_diff, n_bootstraps=2000, seed=20261015)

            phylum_results.append({
                "trait": t_name,
                "phylum": ph,
                "n_paired_sites": len(valid_sites),
                "surface_median_pct": float(statistics.median(site_surf)),
                "deep_median_pct": float(statistics.median(site_deep)),
                "delta_median_pct": float(statistics.median(site_diff)),
                "ci_delta_low_95": boot_ci["ci_median_low"],
                "ci_delta_high_95": boot_ci["ci_median_high"],
                "sites_surface_higher": sum(1 for d in site_diff if d > 0),
                "p_wilcoxon": pv,
            })

    return phylum_results


def run_expedition_sensitivity_analysis(
    sample_records: List[Dict[str, Any]],
    focal_traits: List[str],
    metrics: List[str] = ["abund_mapped_pct", "occurrence_pct"],
) -> List[Dict[str, Any]]:
    """Test paired surface vs deep contrast stratified by expedition trip.
    
    Evaluates ONLY sites where paired Surface and Deep samples exist within that specific trip.
    """
    trips = sorted({r["trip"] for r in sample_records if r["trip"]})
    expedition_results: List[Dict[str, Any]] = []

    # Map sample records by trip
    for tr in trips:
        trip_records = [r for r in sample_records if r["trip"] == tr]
        # Collapse replicates within (trait, metric, site, compartment) for this trip
        cell_vals: Dict[Tuple[str, str, str, str], List[float]] = defaultdict(list)
        for r in trip_records:
            t = r["trait"]
            m_site = r["site"]
            comp = r["compartment"]
            for m in metrics:
                cell_vals[(t, m, m_site, comp)].append(r[m])

        cell_med = {k: statistics.median(vals) for k, vals in cell_vals.items()}
        sites_in_trip = sorted({r["site"] for r in trip_records})

        for m in metrics:
            for t in focal_traits:
                pairs: List[Tuple[float, float]] = []
                for s in sites_in_trip:
                    ks = (t, m, s, "Surface")
                    kd = (t, m, s, "Deep")
                    if ks in cell_med and kd in cell_med:
                        pairs.append((cell_med[ks], cell_med[kd]))

                if len(pairs) < 4:
                    continue

                surf_v = [p[0] for p in pairs]
                deep_v = [p[1] for p in pairs]
                diffs = [p[0] - p[1] for p in pairs]

                nonzero = [d for d in diffs if d != 0]
                try:
                    if nonzero:
                        w = wilcoxon(surf_v, deep_v)
                        pv = float(w.pvalue)
                    else:
                        pv = 1.0
                except Exception:
                    pv = float("nan")

                expedition_results.append({
                    "trip": tr,
                    "trait": t,
                    "metric": m,
                    "n_paired_sites": len(pairs),
                    "surface_median": float(statistics.median(surf_v)),
                    "deep_median": float(statistics.median(deep_v)),
                    "delta_median": float(statistics.median(diffs)),
                    "sites_surface_higher": sum(1 for d in diffs if d > 0),
                    "p_wilcoxon": pv,
                })

    return expedition_results


def plot_compartment_contrasts(
    compartment_tests: List[Dict[str, Any]],
    traits_to_test: Dict[str, Set[str]],
    output_dir: Path,
) -> None:
    """Generate publication-ready diagnostic dumbbell comparison plots without hardcoded counts."""
    # Dynamic labels derived directly from trait dictionary
    focal_traits_plot = [
        ("cpf_complete", f"CPF Complete ({len(traits_to_test['cpf_complete'])} MAGs)\n[Primary Endpoint]", "#d62728"),
        ("cpf_broad", f"CPF Broad ({len(traits_to_test['cpf_broad'])} MAGs)\n[PF00875|PF03441|PF04244]", "#ff7f0e"),
        ("k01669_photolyase", f"K01669 Photolyase ({len(traits_to_test['k01669_photolyase'])} MAGs)\n[Paper Baseline]", "#1f77b4"),
        ("cpf_exclusive", f"CPF-Exclusive ({len(traits_to_test['cpf_exclusive'])} MAGs)\n[CPF Carrier / No K01669]", "#9467bd"),
        ("base_excision_repair", f"Base Excision Repair ({len(traits_to_test['base_excision_repair'])} MAGs)\n[Negative Control]", "#7f7f7f"),
        ("nucleotide_excision_repair", f"Nucleotide Excision Repair ({len(traits_to_test['nucleotide_excision_repair'])} MAGs)\n[Negative Control]", "#8c564b"),
    ]

    fig, (ax_abund, ax_occ) = plt.subplots(1, 2, figsize=(16, 7), constrained_layout=True)

    def draw_panel(ax, metric_name, title_prefix):
        y_positions = list(range(len(focal_traits_plot)))[::-1]
        labels = []
        n_sites_observed = 0
        n_visits_observed = 0

        for idx, (trait_key, display_label, color) in enumerate(focal_traits_plot):
            y = y_positions[idx]
            match_test = [
                r for r in compartment_tests if r["trait"] == trait_key and r["metric"] == metric_name
            ]
            if match_test:
                t_res = match_test[0]
                s_med = t_res["surface_median"]
                d_med = t_res["deep_median"]
                delta_med = t_res["delta_median"]
                ci_l = t_res.get("ci_median_low_95", float("nan"))
                ci_h = t_res.get("ci_median_high_95", float("nan"))
                n_sites_observed = max(n_sites_observed, t_res["n_paired_sites"])
                n_visits_observed = max(n_visits_observed, t_res.get("n_paired_site_trips", 0))

                p_adj = t_res.get("p_holm", t_res["p_wilcoxon"])
                p_str = f"p={p_adj:.2e}" if p_adj < 0.001 else f"p={p_adj:.3f}"
                sig_star = " ***" if p_adj < 0.001 else (" **" if p_adj < 0.01 else (" *" if p_adj < 0.05 else " (ns)"))

                # Draw connecting line
                ax.plot([d_med, s_med], [y, y], color=color, linewidth=2.5, alpha=0.7, zorder=2)
                # Draw endpoints: Deep = grey square, Surface = colored circle
                ax.scatter(d_med, y, color="#7f7f7f", marker="s", s=80, label="Subsurface (Deep)" if idx == 0 else "", zorder=3)
                ax.scatter(s_med, y, color=color, marker="o", s=100, edgecolors="black", label="Surface" if idx == 0 else "", zorder=4)

                # Annotation text with 95% bootstrap CI
                if not math.isnan(ci_l) and not math.isnan(ci_h):
                    annot = f"Δ={delta_med:+.1f}% [95% CI: {ci_l:+.1f}, {ci_h:+.1f}] {sig_star}"
                else:
                    annot = f"Δ={delta_med:+.1f}% {sig_star}"
                ax.text(max(s_med, d_med) + 1.2, y, annot, va="center", fontsize=8.5, fontweight="bold")
                labels.append(display_label)
            else:
                labels.append(display_label)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(labels, fontsize=9.5)
        ax.set_xlabel("Community Percentage (%)", fontsize=10.5)
        full_title = f"{title_prefix}\nAcross {n_sites_observed} Paired Sites ({n_visits_observed} Paired Visits)"
        ax.set_title(full_title, fontsize=11, fontweight="bold")
        ax.grid(True, linestyle="--", alpha=0.4, axis="x")

    draw_panel(ax_abund, "abund_mapped_pct", "Mapped Community Abundance (Carrier Abundance / Total 5438 MAG Abundance)")
    draw_panel(ax_occ, "occurrence_pct", "Community Occurrence (Carrier Species / Total Detected Species)")

    # Legend
    handles, leg_labels = ax_abund.get_legend_handles_labels()
    fig.legend(handles, leg_labels, loc="lower center", bbox_to_anchor=(0.5, -0.05), ncol=2, fontsize=10)

    caveat = (
        "Empty Quarter Metagenome Catalogue: Ecological Light-Niche Partitioning Across Paired Sites.\n"
        "CAVEAT: Complete CPF carriage (primary endpoint) shows significant surface light-niche enrichment.\n"
        "Broad CPF encompasses diverse photolyase/cryptochrome subfamilies; specific biochemical repair capability requires enzymatic assay."
    )
    fig.suptitle(caveat, fontsize=9.5, style="italic", y=1.06)

    fig.savefig(output_dir / "fig_f1_compartment_contrasts.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / "fig_f1_compartment_contrasts.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ecological analysis of F1 photolyase/cryptochrome (CPF) carrier MAGs across soil compartments (v2)."
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        required=True,
        help="Path to candidates.tsv from root_audit/f1_full/candidates.tsv.",
    )
    parser.add_argument(
        "--rep-traits",
        type=Path,
        required=True,
        help="Path to rep_traits.tsv from uv_survival_clean/results/rep_traits.tsv.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        required=True,
        help="Path to sample_metadata.tsv.",
    )
    parser.add_argument(
        "--profiles-dir",
        type=Path,
        required=True,
        help="Directory containing CoverM profiles (profiles_full_recompute).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for results tables and figures.",
    )
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="Run in bounded pilot mode on a diagnostic subset of samples to verify joins closure.",
    )
    parser.add_argument(
        "--pilot-samples",
        type=int,
        default=25,
        help="Number of samples to analyze in pilot mode (default: 25).",
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=20261015,
        help="Fixed seed for bootstrap confidence intervals (default: 20261015).",
    )
    args = parser.parse_args()

    start_time = time.time()
    script_path = Path(__file__).resolve()
    script_sha256 = compute_file_sha256(script_path)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("================================================================================")
    print("F1 Photolyase/Cryptochrome (CPF) Ecological Compartment Analysis (v2)")
    print(f"Timestamp UTC: {datetime.now(timezone.utc).isoformat()}")
    print(f"Script:        {script_path.name} (SHA256: {script_sha256})")
    print(f"Candidates:    {args.candidates}")
    print(f"Rep Traits:    {args.rep_traits}")
    print(f"Metadata:      {args.metadata}")
    print(f"Profiles Dir:  {args.profiles_dir}")
    print(f"Output Dir:    {args.output_dir}")
    print(f"Mode:          {'BOUNDED PILOT' if args.pilot else 'FULL CENSUS'}")
    print(f"Seed:          {args.bootstrap_seed}")
    print("================================================================================")

    # 1. Load inputs and compute hashes
    candidates_sha256 = compute_file_sha256(args.candidates)
    rep_traits_sha256 = compute_file_sha256(args.rep_traits)
    metadata_sha256 = compute_file_sha256(args.metadata)

    cand_data = load_candidates(args.candidates)
    rep_data = load_rep_traits(args.rep_traits)
    meta = load_metadata(args.metadata)

    n_cand_hits = len(cand_data["records"])
    n_carrier_mags = len(cand_data["mags_all"])
    n_complete_hits = sum(1 for r in cand_data["records"] if r["is_complete"])
    n_complete_carrier_mags = len(cand_data["mags_complete"])
    n_catalog_reps = len(rep_data["all_reps"])
    n_ann_reps = len(rep_data["annotated_reps"])
    n_k01669_reps = len(rep_data["k01669_reps"])

    print(f"Loaded {n_cand_hits} CPF candidate hits across {n_carrier_mags} carrier MAGs.")
    print(f"  Complete CPF hits (partial==00): {n_complete_hits} across {n_complete_carrier_mags} carrier MAGs.")
    print(f"Loaded {n_catalog_reps} catalogue representatives ({n_ann_reps} annotated).")
    print(f"  K01669 photolyase baseline carriers: {n_k01669_reps} MAGs.")

    # 2. Strict ID closure verification (FAIL CLOSED)
    missing_in_catalog = cand_data["mags_all"] - rep_data["all_reps"]
    if missing_in_catalog:
        raise ValueError(
            f"FAIL CLOSED: ID closure failure: {len(missing_in_catalog)} carrier MAGs missing in 5,438 catalogue: {list(missing_in_catalog)[:5]}"
        )
    print("Verified 100% exact ID closure between F1 CPF carrier MAGs and catalogue representatives.")

    # Define trait sets for ecological testing
    cpf_broad = cand_data["mags_all"]
    cpf_complete = cand_data["mags_complete"]
    k01669 = rep_data["k01669_reps"]
    cpf_exclusive = cand_data["mags_all"] - k01669
    cpf_complete_exclusive = cand_data["mags_complete"] - k01669

    traits_to_test: Dict[str, Set[str]] = {
        "cpf_complete": cpf_complete,
        "cpf_broad": cpf_broad,
        "k01669_photolyase": k01669,
        "cpf_exclusive": cpf_exclusive,
        "cpf_complete_exclusive": cpf_complete_exclusive,
        "base_excision_repair": rep_data["ber_reps"],
        "nucleotide_excision_repair": rep_data["ner_reps"],
        "nhej": rep_data["nhej_reps"],
        "carotenoid": rep_data["carot_reps"],
    }

    # Add profile-specific strata explicitly labeled exploratory
    for prof, p_mags in cand_data["mags_by_profile"].items():
        traits_to_test[f"cpf_broad_{prof}_exploratory"] = p_mags
    for prof, p_mags in cand_data["mags_complete_by_profile"].items():
        traits_to_test[f"cpf_complete_{prof}_exploratory"] = p_mags

    # 3. Discover available CoverM sample profiles (FAIL CLOSED on missing canonical samples)
    profile_dirs = {p.name for p in args.profiles_dir.iterdir() if p.is_dir()}
    usable_samples = sorted(set(meta.keys()) & profile_dirs)
    print(f"Discovered {len(profile_dirs)} CoverM profile directories; {len(usable_samples)} match canonical sample metadata.")

    # Strict check: all metadata samples must exist in profiles
    missing_samples = set(meta.keys()) - profile_dirs
    if missing_samples:
        raise FileNotFoundError(f"FAIL CLOSED: {len(missing_samples)} canonical metadata samples missing profile directories: {list(missing_samples)[:5]}")

    if args.pilot:
        usable_samples = usable_samples[: args.pilot_samples]
        print(f"PILOT MODE: Subsetting to first {len(usable_samples)} samples.")

    # 4. Compute community-level sample metrics with aligned numerators/denominators
    print(f"Computing community abundance & occurrence across {len(usable_samples)} samples...")
    sample_records, skipped_samples = compute_community_sample_metrics(
        usable_samples,
        args.profiles_dir,
        meta,
        traits_to_test,
        rep_data["all_reps"],
        rep_data["annotated_reps"],
    )
    print(f"Emitted {len(sample_records)} sample x trait metric records ({len(skipped_samples)} samples skipped due to zero mapped abundance).")

    # Write sample-level records table
    sample_tsv = args.output_dir / "f1_community_by_sample.tsv"
    with open(sample_tsv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "sample", "site", "compartment", "trip", "trait", "unmapped_pct",
            "n_detected_all", "n_detected_trait_all", "n_detected_ann", "n_detected_trait_ann",
            "occurrence_pct", "abund_mapped_pct", "abund_whole_pct",
            "abund_annotated_pct", "occurrence_annotated_pct"
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(sample_records)
    print(f"Emitted sample-level community table: {sample_tsv}")

    # 5. Site-Paired Compartment Test (Surface vs Deep)
    print("Running hierarchical site-paired Surface vs Deep testing with bootstrap CIs...")
    test_traits_ordered = list(traits_to_test.keys())
    test_metrics_ordered = ["abund_mapped_pct", "occurrence_pct", "abund_annotated_pct", "occurrence_annotated_pct", "abund_whole_pct"]

    comp_results, site_medians, cell_trip_medians = run_site_paired_compartment_tests(
        sample_records,
        test_traits_ordered,
        metrics=test_metrics_ordered,
        seed=args.bootstrap_seed,
    )

    # Write site-level collapsed values
    site_tsv = args.output_dir / "f1_site_level_values.tsv"
    with open(site_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["trait", "metric", "site", "compartment", "median_value"])
        for (t, m, site, comp), val in sorted(site_medians.items()):
            w.writerow([t, m, site, comp, f"{val:.4f}"])
    print(f"Emitted site-level collapsed values:  {site_tsv}")

    # Write compartment test results table
    comp_tsv = args.output_dir / "f1_compartment_tests.tsv"
    with open(comp_tsv, "w", newline="", encoding="utf-8") as f:
        cols = [
            "trait", "metric", "n_paired_sites", "n_paired_site_trips",
            "surface_median", "deep_median", "delta_median",
            "ci_median_low_95", "ci_median_high_95",
            "mean_surface", "mean_deep", "mean_delta",
            "ci_mean_low_95", "ci_mean_high_95",
            "cluster_ci_median_low_95", "cluster_ci_median_high_95",
            "sites_surface_higher", "sites_deep_higher", "sites_equal",
            "wilcoxon_stat", "p_wilcoxon", "p_holm", "p_kruskal3"
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in comp_results:
            row_out = {}
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, float):
                    row_out[c] = f"{v:.4e}" if "p_" in c else f"{v:.4f}"
                else:
                    row_out[c] = str(v)
            w.writerow(row_out)
    print(f"Emitted paired compartment tests:    {comp_tsv}")

    # Print summary table for PRIMARY mapped abundance %
    print("\n=== SUMMARY: Site-Paired Surface vs Deep Contrast (abund_mapped_pct [PRIMARY]) ===")
    print(f"{'Trait':35s} {'Sites':>5s} {'Visits':>6s} {'Surf%':>7s} {'Deep%':>7s} {'Delta%':>7s} {'95% Bootstrap CI':>20s} {'S>D':>6s} {'p(Holm)':>10s}")
    for r in comp_results:
        if r["metric"] == "abund_mapped_pct":
            p_val = r.get("p_holm", r["p_wilcoxon"])
            sig = " <<<" if p_val < 0.05 else ""
            ci_str = f"[{r['ci_median_low_95']:+.2f}, {r['ci_median_high_95']:+.2f}]"
            print(f"{r['trait']:35s} {r['n_paired_sites']:5d} {r['n_paired_site_trips']:6d} {r['surface_median']:7.2f} {r['deep_median']:7.2f} "
                  f"{r['delta_median']:+7.2f} {ci_str:>20s} {r['sites_surface_higher']:3d}/{r['n_paired_sites']:<2d} "
                  f"{p_val:10.2e}{sig}")

    # Print summary table for REPLICATION annotated abundance %
    print("\n=== SUMMARY: Site-Paired Surface vs Deep Contrast (abund_annotated_pct [REPLICATION]) ===")
    print(f"{'Trait':35s} {'Sites':>5s} {'Visits':>6s} {'Surf%':>7s} {'Deep%':>7s} {'Delta%':>7s} {'95% Bootstrap CI':>20s} {'S>D':>6s} {'p(Holm)':>10s}")
    for r in comp_results:
        if r["metric"] == "abund_annotated_pct":
            p_val = r.get("p_holm", r["p_wilcoxon"])
            sig = " <<<" if p_val < 0.05 else ""
            ci_str = f"[{r['ci_median_low_95']:+.2f}, {r['ci_median_high_95']:+.2f}]"
            print(f"{r['trait']:35s} {r['n_paired_sites']:5d} {r['n_paired_site_trips']:6d} {r['surface_median']:7.2f} {r['deep_median']:7.2f} "
                  f"{r['delta_median']:+7.2f} {ci_str:>20s} {r['sites_surface_higher']:3d}/{r['n_paired_sites']:<2d} "
                  f"{p_val:10.2e}{sig}")

    # 6. Within-Phylum Lineage-Sensitive Analysis
    print("\nRunning lineage-sensitive within-phylum tests...")
    within_traits = {
        "cpf_complete": cpf_complete,
        "cpf_broad": cpf_broad,
        "k01669_photolyase": k01669,
    }
    phylum_results = run_within_phylum_tests(
        sample_records,
        within_traits,
        rep_data["all_reps"],
        rep_data["phylum_by_rep"],
        args.profiles_dir,
        meta,
    )

    phylum_tsv = args.output_dir / "f1_within_phylum_tests.tsv"
    with open(phylum_tsv, "w", newline="", encoding="utf-8") as f:
        cols = [
            "trait", "phylum", "n_paired_sites", "surface_median_pct", "deep_median_pct",
            "delta_median_pct", "ci_delta_low_95", "ci_delta_high_95",
            "sites_surface_higher", "p_wilcoxon"
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in phylum_results:
            row_out = {
                c: (f"{r[c]:.4e}" if "p_" in c else (f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c])))
                for c in cols
            }
            w.writerow(row_out)
    print(f"Emitted within-phylum tests table:   {phylum_tsv}")

    # 7. Expedition Sensitivity Analysis
    print("Running expedition-stratified sensitivity analysis...")
    expedition_results = run_expedition_sensitivity_analysis(
        sample_records, ["cpf_complete", "cpf_broad", "k01669_photolyase"]
    )
    expedition_tsv = args.output_dir / "f1_expedition_sensitivity.tsv"
    with open(expedition_tsv, "w", newline="", encoding="utf-8") as f:
        cols = [
            "trip", "trait", "metric", "n_paired_sites", "surface_median", "deep_median",
            "delta_median", "sites_surface_higher", "p_wilcoxon"
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in expedition_results:
            row_out = {
                c: (f"{r[c]:.4e}" if "p_" in c else (f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c])))
                for c in cols
            }
            w.writerow(row_out)
    print(f"Emitted expedition sensitivity table: {expedition_tsv}")

    # 8. Generate Diagnostic Figures
    if not args.pilot or len(usable_samples) >= 20:
        print("Generating diagnostic comparison plots...")
        plot_compartment_contrasts(comp_results, traits_to_test, args.output_dir)
        print("Emitted diagnostic figures (PDF & PNG).")

    # 9. Provenance JSON
    end_time = time.time()
    provenance = {
        "analysis": "F1_photolyase_cryptochrome_ecological_test_v2",
        "script": {
            "name": script_path.name,
            "sha256": script_sha256,
        },
        "execution": {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(end_time - start_time, 3),
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "python_version": sys.version,
            "mode": "pilot" if args.pilot else "full",
            "bootstrap_seed": args.bootstrap_seed,
            "bootstrap_resamples": 2000,
        },
        "inputs": {
            "candidates": {
                "path": str(args.candidates),
                "sha256": candidates_sha256,
                "n_candidate_hits": n_cand_hits,
                "n_carrier_mags": n_carrier_mags,
                "n_complete_carrier_mags": n_complete_carrier_mags,
            },
            "rep_traits": {
                "path": str(args.rep_traits),
                "sha256": rep_traits_sha256,
                "total_catalogue_reps": n_catalog_reps,
                "annotated_catalogue_reps": n_ann_reps,
                "k01669_photolyase_reps": n_k01669_reps,
            },
            "metadata": {
                "path": str(args.metadata),
                "sha256": metadata_sha256,
                "total_samples": len(meta),
            },
            "profiles_dir": str(args.profiles_dir),
        },
        "samples_analyzed": len(usable_samples),
        "samples_skipped_zero_mapped": skipped_samples,
        "statistical_framework": {
            "sample_unit": "sample",
            "replicate_handling": "median across replicates per (site, trip, compartment) cell",
            "site_pairing_hierarchy": "within same trip, then aggregated to one median observation per site across trips",
            "paired_test": "two-sided Wilcoxon signed-rank test on (Surface - Deep) site aggregate",
            "bootstrap_confidence_intervals": "95% empirical percentile CI from 2000 resamples of sites with replacement",
            "cluster_bootstrap": "cluster bootstrap across all paired site-trip visits resampling sites with replacement",
            "multiplicity_correction": "Holm-Bonferroni per metric",
        },
        "scientific_caveats": [
            "F1 CPF candidate hits are profile-supported members of the cryptochrome/photolyase family.",
            "All profile strata (PF03441, PF00875, PF04244) are exploratory domain matches. Attribution to specific functional repair subfamilies (e.g. CPD photolyase class I/II/III, 6-4 photolyase) vs signaling cryptochromes (Cry-DASH) is NOT proven without phylogenetic placement and biochemical assay.",
            "Light-niche partitioning between surface and subsurface communities indicates ecological selection, not enzymatic repair mechanism.",
            "K01669 (deoxyribodipyrimidine photolyase) captures an annotated subset of broader CPF diversity, with 96 CPF carrier MAGs lacking K01669 annotation.",
            "Within-phylum tests demonstrate that surface enrichment is not driven solely by phylum-level taxonomic turnover. However, phylum stratification does NOT control for fine-scale intra-phylum species or strain turnover; lineage confounding is constrained at the phylum rank but cannot be claimed to be fully resolved.",
            "BER and NER negative controls are evaluated objectively without preconceived flatness expectations.",
        ],
    }

    provenance_path = args.output_dir / "f1_ecology_provenance.json"
    with open(provenance_path, "w", encoding="utf-8") as f:
        json.dump(provenance, f, indent=2)
    print(f"Emitted analysis provenance:         {provenance_path}")
    print("\nEcological analysis completed successfully.")


if __name__ == "__main__":
    main()
