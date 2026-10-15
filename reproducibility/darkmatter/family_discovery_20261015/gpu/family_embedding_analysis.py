#!/usr/bin/env python3
"""family_embedding_analysis.py - Standalone Reproducible CPU Analysis & Plotting for Family Embeddings.

Evaluates high-dimensional ProstT5 protein embeddings against curated references
and performs multi-seed manifold projection (PCA and multi-seed t-SNE).

Key Design Principles & Verification Rules:
1. Strict ID closure and sequence hash verification between --embeddings NPZ and --metadata TSV.
2. Unrelated folds/classes kept strictly separate; never pool distinct families into a shared plot.
3. Multi-seed t-SNE (seeds 42, 123, 20261015) with PCA initialization and valid perplexities subject to N.
   No cherry-picking of "prettiest" seed; all seeds displayed side-by-side.
4. High-dimensional cosine nearest-reference distances computed in native 1024-d space (leave-one-out for references).
5. Explicit duplicate tracking: duplicate sequence counts and weights recorded in tables and reflected in point size.
6. Diagnostic plots: distance vs sequence length difference to detect truncation artifacts.
7. Explicit caveats: descriptive manifold projections and cosine metrics do NOT prove functional novelty or classification.
8. Machine-readable analysis provenance recording code SHA, input hashes, library versions, and graph trustworthiness.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import socket
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
import sklearn
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, trustworthiness


def compute_file_sha256(filepath: Path | str) -> str:
    """Compute hex SHA256 of a local file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def compute_string_sha256(s: str) -> str:
    """Compute hex SHA256 of a string."""
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def load_embeddings(embeddings_path: Path) -> Dict[str, Any]:
    """Load and validate the embeddings NPZ archive produced by prostt5_embed.py."""
    if not embeddings_path.is_file():
        raise FileNotFoundError(f"Embeddings file not found: {embeddings_path}")

    # allow_pickle=True required for string/object arrays stored in NPZ
    npz_data = np.load(embeddings_path, allow_pickle=True)
    available_keys = list(npz_data.keys())

    if "embeddings" not in available_keys:
        raise KeyError(f"NPZ archive missing 'embeddings' array. Found keys: {available_keys}")
    if "ids" not in available_keys:
        raise KeyError(f"NPZ archive missing 'ids' array. Found keys: {available_keys}")

    embeddings = npz_data["embeddings"]
    raw_ids = npz_data["ids"]

    if not isinstance(embeddings, np.ndarray):
        raise TypeError(f"Embeddings array is not a numpy ndarray: {type(embeddings)}")
    if embeddings.ndim != 2:
        raise ValueError(f"Embeddings array must be 2D (N, D), found shape {embeddings.shape}")
    if embeddings.shape[1] != 1024:
        raise ValueError(f"Expected embedding dimension 1024, found {embeddings.shape[1]}")

    n_samples = embeddings.shape[0]
    ids = [str(x).strip() for x in raw_ids]

    if len(ids) != n_samples:
        raise ValueError(f"Length mismatch: {n_samples} embeddings but {len(ids)} IDs")

    if len(set(ids)) != len(ids):
        duplicates = [x for x in ids if ids.count(x) > 1]
        raise ValueError(f"Duplicate IDs detected within embeddings NPZ: {set(duplicates)}")

    # Check for NaN / Inf
    if not np.isfinite(embeddings).all():
        nan_count = int(np.isnan(embeddings).sum())
        inf_count = int(np.isinf(embeddings).sum())
        raise ValueError(f"Embeddings contain non-finite values (NaN: {nan_count}, Inf: {inf_count})")

    lengths = None
    if "lengths" in available_keys:
        lengths = [int(x) for x in npz_data["lengths"]]

    raw_sha256 = None
    if "raw_sequence_sha256" in available_keys:
        raw_sha256 = [str(x).strip().lower() for x in npz_data["raw_sequence_sha256"]]

    clean_sha256 = None
    if "clean_sequence_sha256" in available_keys:
        clean_sha256 = [str(x).strip().lower() for x in npz_data["clean_sequence_sha256"]]
    elif "sequence_sha256" in available_keys:
        clean_sha256 = [str(x).strip().lower() for x in npz_data["sequence_sha256"]]

    return {
        "embeddings": embeddings.astype(np.float32),
        "ids": ids,
        "lengths": lengths,
        "raw_sha256": raw_sha256,
        "clean_sha256": clean_sha256,
    }


def load_metadata(metadata_path: Path) -> List[Dict[str, str]]:
    """Load and validate caller-provided metadata TSV."""
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    records: List[Dict[str, str]] = []
    seen_ids: Set[str] = set()

    with open(metadata_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Empty or headerless metadata file: {metadata_path}")

        fieldnames = [c.strip() for c in reader.fieldnames]
        field_map = {c.lower(): c for c in fieldnames}

        # Check required columns
        for req in ["id", "source", "family_class"]:
            if req not in field_map:
                raise ValueError(
                    f"Metadata TSV missing required column '{req}'. Present columns: {fieldnames}"
                )

        id_col = field_map["id"]
        source_col = field_map["source"]
        family_col = field_map["family_class"]
        length_col = field_map.get("length") or field_map.get("len")
        phylum_col = field_map.get("phylum")
        genus_col = field_map.get("genus")
        clean_sha_col = field_map.get("clean_sha256") or field_map.get("sequence_sha256")
        raw_sha_col = field_map.get("raw_sha256")

        for line_no, row in enumerate(reader, start=2):
            seq_id = row[id_col].strip()
            if not seq_id:
                raise ValueError(f"Empty ID at row {line_no} in metadata: {metadata_path}")
            if seq_id in seen_ids:
                raise ValueError(f"Duplicate ID '{seq_id}' in metadata TSV at row {line_no}")
            seen_ids.add(seq_id)

            source = row[source_col].strip()
            if not source:
                raise ValueError(f"Empty source for ID '{seq_id}' at row {line_no}")

            family_class = row[family_col].strip()
            if not family_class:
                raise ValueError(f"Empty family_class for ID '{seq_id}' at row {line_no}")

            length_val = row[length_col].strip() if length_col and row.get(length_col) else ""
            phylum_val = row[phylum_col].strip() if phylum_col and row.get(phylum_col) else "Unknown"
            genus_val = row[genus_col].strip() if genus_col and row.get(genus_col) else "Unknown"
            clean_sha = row[clean_sha_col].strip().lower() if clean_sha_col and row.get(clean_sha_col) else ""
            raw_sha = row[raw_sha_col].strip().lower() if raw_sha_col and row.get(raw_sha_col) else ""

            records.append({
                "id": seq_id,
                "source": source,
                "family_class": family_class,
                "length": length_val,
                "phylum": phylum_val or "Unknown",
                "genus": genus_val or "Unknown",
                "clean_sha256": clean_sha,
                "raw_sha256": raw_sha,
            })

    if not records:
        raise ValueError(f"Metadata file contains no data rows: {metadata_path}")

    return records


def verify_id_closure_and_joins(
    emb_data: Dict[str, Any], meta_records: List[Dict[str, str]]
) -> Dict[str, Dict[str, Any]]:
    """Enforce exact 1-to-1 ID closure and verify sequence hash/length consistency."""
    emb_ids = emb_data["ids"]
    emb_id_set = set(emb_ids)
    meta_id_map = {r["id"]: r for r in meta_records}
    meta_id_set = set(meta_id_map.keys())

    missing_in_meta = emb_id_set - meta_id_set
    extra_in_meta = meta_id_set - emb_id_set

    if missing_in_meta or extra_in_meta:
        err_msg = ["ID closure validation FAILED:"]
        if missing_in_meta:
            err_msg.append(f"  IDs in embeddings but missing in metadata ({len(missing_in_meta)}): {sorted(missing_in_meta)[:10]}")
        if extra_in_meta:
            err_msg.append(f"  IDs in metadata but missing in embeddings ({len(extra_in_meta)}): {sorted(extra_in_meta)[:10]}")
        raise ValueError("\n".join(err_msg))

    joined_data: Dict[str, Dict[str, Any]] = {}
    for i, seq_id in enumerate(emb_ids):
        meta = meta_id_map[seq_id]
        emb_vec = emb_data["embeddings"][i]

        # Length verification
        emb_len = emb_data["lengths"][i] if emb_data["lengths"] is not None else None
        meta_len = int(meta["length"]) if meta["length"] else None
        if emb_len is not None and meta_len is not None:
            if emb_len != meta_len:
                raise ValueError(
                    f"Length mismatch for ID '{seq_id}': embeddings length {emb_len} vs metadata length {meta_len}"
                )
        final_len = emb_len if emb_len is not None else (meta_len if meta_len is not None else 0)

        # Hash verification
        emb_clean_sha = emb_data["clean_sha256"][i] if emb_data["clean_sha256"] is not None else ""
        meta_clean_sha = meta["clean_sha256"]
        if emb_clean_sha and meta_clean_sha:
            if emb_clean_sha != meta_clean_sha:
                raise ValueError(
                    f"Clean SHA256 mismatch for ID '{seq_id}': embeddings {emb_clean_sha} vs metadata {meta_clean_sha}"
                )
        final_clean_sha = emb_clean_sha or meta_clean_sha

        emb_raw_sha = emb_data["raw_sha256"][i] if emb_data["raw_sha256"] is not None else ""
        meta_raw_sha = meta["raw_sha256"]
        if emb_raw_sha and meta_raw_sha:
            if emb_raw_sha != meta_raw_sha:
                raise ValueError(
                    f"Raw SHA256 mismatch for ID '{seq_id}': embeddings {emb_raw_sha} vs metadata {meta_raw_sha}"
                )
        final_raw_sha = emb_raw_sha or meta_raw_sha

        joined_data[seq_id] = {
            "id": seq_id,
            "embedding": emb_vec,
            "source": meta["source"],
            "family_class": meta["family_class"],
            "length": final_len,
            "phylum": meta["phylum"],
            "genus": meta["genus"],
            "clean_sha256": final_clean_sha,
            "raw_sha256": final_raw_sha,
        }

    return joined_data


def detect_exact_duplicates(
    items: List[Dict[str, Any]]
) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Identify exact sequence duplicate clusters based on sequence SHA or embedding distance."""
    # First group by clean sequence sha if available
    sha_groups: Dict[str, List[str]] = {}
    for item in items:
        sha = item.get("clean_sha256")
        if not sha:
            # Fallback to embedding SHA if no clean sequence sha
            sha = hashlib.sha256(item["embedding"].tobytes()).hexdigest()
        sha_groups.setdefault(sha, []).append(item["id"])

    dup_group_map: Dict[str, str] = {}
    dup_count_map: Dict[str, int] = {}

    grp_idx = 1
    for sha, group_ids in sha_groups.items():
        count = len(group_ids)
        grp_name = f"grp_{grp_idx:04d}" if count > 1 else "unique"
        for seq_id in group_ids:
            dup_group_map[seq_id] = grp_name
            dup_count_map[seq_id] = count
        if count > 1:
            grp_idx += 1

    return dup_group_map, dup_count_map


def compute_cosine_distances(
    family_items: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Compute high-dimensional cosine nearest-reference distances.

    Queries: min cosine distance to any curated reference.
    References: leave-one-out min cosine distance to any OTHER curated reference.
    """
    n = len(family_items)
    embeddings = np.array([item["embedding"] for item in family_items], dtype=np.float32)

    # Normalize vectors for cosine metric: d_cos(u, v) = 1.0 - (u_norm . v_norm)
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    normalized_emb = embeddings / norms

    # Pairwise cosine distance matrix: shape (N, N)
    sim_matrix = np.clip(np.matmul(normalized_emb, normalized_emb.T), -1.0, 1.0)
    dist_matrix = np.clip(1.0 - sim_matrix, 0.0, 2.0)

    # Distinguish references from queries
    is_ref = [item["source"].strip().lower() == "curated_reference" for item in family_items]
    ref_indices = [i for i, r in enumerate(is_ref) if r]
    n_refs = len(ref_indices)

    results: List[Dict[str, Any]] = []

    for i, item in enumerate(family_items):
        item_is_ref = is_ref[i]

        if item_is_ref:
            # Leave-one-out nearest reference
            other_ref_indices = [idx for idx in ref_indices if idx != i]
            if not other_ref_indices:
                nearest_ref_id = "N/A (single reference)"
                nearest_ref_source = "N/A"
                nearest_ref_len = None
                len_diff = None
                min_dist = float("nan")
            else:
                sub_dists = dist_matrix[i, other_ref_indices]
                best_sub_idx = int(np.argmin(sub_dists))
                nearest_idx = other_ref_indices[best_sub_idx]
                min_dist = float(sub_dists[best_sub_idx])
                nearest_ref_id = family_items[nearest_idx]["id"]
                nearest_ref_source = family_items[nearest_idx]["source"]
                nearest_ref_len = family_items[nearest_idx]["length"]
                len_diff = abs(item["length"] - nearest_ref_len) if (item["length"] and nearest_ref_len) else 0
        else:
            # Query vs all references
            if n_refs == 0:
                nearest_ref_id = "N/A (no references in family)"
                nearest_ref_source = "N/A"
                nearest_ref_len = None
                len_diff = None
                min_dist = float("nan")
            else:
                sub_dists = dist_matrix[i, ref_indices]
                best_sub_idx = int(np.argmin(sub_dists))
                nearest_idx = ref_indices[best_sub_idx]
                min_dist = float(sub_dists[best_sub_idx])
                nearest_ref_id = family_items[nearest_idx]["id"]
                nearest_ref_source = family_items[nearest_idx]["source"]
                nearest_ref_len = family_items[nearest_idx]["length"]
                len_diff = abs(item["length"] - nearest_ref_len) if (item["length"] and nearest_ref_len) else 0

        results.append({
            "family_class": item["family_class"],
            "id": item["id"],
            "source": item["source"],
            "length": item["length"],
            "phylum": item["phylum"],
            "genus": item["genus"],
            "nearest_reference_id": nearest_ref_id,
            "nearest_reference_source": nearest_ref_source,
            "nearest_reference_length": nearest_ref_len if nearest_ref_len is not None else "N/A",
            "length_difference": len_diff if len_diff is not None else "N/A",
            "cosine_distance": min_dist,
            "is_reference": item_is_ref,
        })

    return results


def select_valid_perplexities(n_samples: int, candidate_perplexities: List[float]) -> List[float]:
    """Select valid perplexities subject to N for sklearn t-SNE."""
    if n_samples <= 1:
        return []

    # In sklearn TSNE, perplexity must be strictly < n_samples.
    # Preferred bound: p <= (n_samples - 1) / 3.0
    valid = [p for p in candidate_perplexities if p < n_samples and p <= (n_samples - 1) / 3.0]
    if not valid:
        # Check if any candidate is strictly < n_samples
        valid = [p for p in candidate_perplexities if p < n_samples]

    if not valid:
        # Clamped fallback for tiny sample sizes (e.g. N <= 5 in smoke test)
        adjusted_p = max(1.0, float(min(5.0, (n_samples - 1) / 2.0)))
        return [round(adjusted_p, 1)]

    return valid


def run_dimensionality_reduction(
    family_items: List[Dict[str, Any]],
    candidate_perplexities: List[float],
    seeds: List[int],
    tsne_max_iter: int = 1000,
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
    """Compute PCA and multi-seed t-SNE projections with graph trustworthiness."""
    n_samples = len(family_items)
    embeddings = np.array([item["embedding"] for item in family_items], dtype=np.float32)

    coords_by_id: Dict[str, Dict[str, Any]] = {
        item["id"]: {
            "family_class": item["family_class"],
            "id": item["id"],
            "source": item["source"],
            "length": item["length"],
            "phylum": item["phylum"],
            "genus": item["genus"],
        }
        for item in family_items
    }

    dim_reduction_meta: Dict[str, Any] = {
        "n_samples": n_samples,
        "pca_variance_ratio": [],
        "perplexities_run": [],
        "trustworthiness": {},
    }

    # PCA Baseline
    if n_samples >= 2:
        n_comp = min(2, n_samples)
        pca = PCA(n_components=n_comp, random_state=42)
        pca_coords = pca.fit_transform(embeddings)
        var_ratios = [float(v) for v in pca.explained_variance_ratio_]
        dim_reduction_meta["pca_variance_ratio"] = var_ratios

        for i, item in enumerate(family_items):
            coords_by_id[item["id"]]["pca_1"] = float(pca_coords[i, 0])
            coords_by_id[item["id"]]["pca_2"] = float(pca_coords[i, 1]) if n_comp > 1 else 0.0
            coords_by_id[item["id"]]["pca_var_ratio_1"] = var_ratios[0]
            coords_by_id[item["id"]]["pca_var_ratio_2"] = var_ratios[1] if n_comp > 1 else 0.0
    else:
        for item in family_items:
            coords_by_id[item["id"]]["pca_1"] = 0.0
            coords_by_id[item["id"]]["pca_2"] = 0.0
            coords_by_id[item["id"]]["pca_var_ratio_1"] = 0.0
            coords_by_id[item["id"]]["pca_var_ratio_2"] = 0.0

    # Multi-seed t-SNE
    valid_perplexities = select_valid_perplexities(n_samples, candidate_perplexities)
    dim_reduction_meta["perplexities_run"] = valid_perplexities

    # Use method="exact" for N < 500 to prevent Barnes-Hut quad-tree degeneration on small N or exact duplicates
    tsne_method = "exact" if n_samples < 500 else "barnes_hut"

    for perp in valid_perplexities:
        for seed in seeds:
            field_x = f"tsne_perp{perp}_seed{seed}_1"
            field_y = f"tsne_perp{perp}_seed{seed}_2"
            tw_key = f"perp_{perp}_seed_{seed}"

            if n_samples >= 2:
                # Always use PCA initialization with specified seed
                tsne = TSNE(
                    n_components=2,
                    perplexity=perp,
                    init="pca",
                    method=tsne_method,
                    random_state=seed,
                    max_iter=tsne_max_iter,
                )
                tsne_coords = tsne.fit_transform(embeddings)

                # Compute trustworthiness if N >= 4 and k_neighbors < n_samples / 2
                tw_score = None
                k_neighbors = min(5, (n_samples - 1) // 2)
                if n_samples >= 4 and k_neighbors >= 1:
                    try:
                        tw_score = float(
                            trustworthiness(embeddings, tsne_coords, n_neighbors=k_neighbors, metric="cosine")
                        )
                    except Exception as e:
                        tw_score = None
                dim_reduction_meta["trustworthiness"][tw_key] = tw_score

                for i, item in enumerate(family_items):
                    coords_by_id[item["id"]][field_x] = float(tsne_coords[i, 0])
                    coords_by_id[item["id"]][field_y] = float(tsne_coords[i, 1])
            else:
                dim_reduction_meta["trustworthiness"][tw_key] = None
                for item in family_items:
                    coords_by_id[item["id"]][field_x] = 0.0
                    coords_by_id[item["id"]][field_y] = 0.0

    return coords_by_id, dim_reduction_meta


def get_source_style(source: str) -> Tuple[str, str, str]:
    """Assign consistent color, marker, and label for each source category."""
    s_lower = source.strip().lower()
    if s_lower == "curated_reference":
        return "#1f77b4", "D", "Curated Reference"
    elif s_lower.startswith("eq"):
        return "#d62728", "o", "EQ Environmental"
    elif s_lower.startswith("smag"):
        return "#2ca02c", "s", "SMAG Environmental"
    elif s_lower.startswith("gem"):
        return "#9467bd", "^", "GEM Environmental"
    else:
        return "#7f7f7f", "v", source


def plot_family_tsne_multiseed(
    family_class: str,
    family_items: List[Dict[str, Any]],
    coords_by_id: Dict[str, Dict[str, Any]],
    dup_count_map: Dict[str, int],
    valid_perplexities: List[float],
    seeds: List[int],
    dim_reduction_meta: Dict[str, Any],
    output_dir: Path,
) -> None:
    """Generate multi-seed comparison plots (PCA + 3 t-SNE seeds side-by-side).

    No cherry-picked 'prettiest' seed: all 3 seeds are explicitly displayed.
    Duplicate counts explicitly weight point size.
    """
    if not valid_perplexities:
        return

    # Create a figure for each valid perplexity
    for perp in valid_perplexities:
        fig, axes = plt.subplots(1, 4, figsize=(22, 5.2), constrained_layout=True)

        # Panel 1: PCA Baseline
        ax_pca = axes[0]
        var_ratios = dim_reduction_meta.get("pca_variance_ratio", [0.0, 0.0])
        var1_str = f"{var_ratios[0]*100:.1f}%" if len(var_ratios) > 0 else "N/A"
        var2_str = f"{var_ratios[1]*100:.1f}%" if len(var_ratios) > 1 else "N/A"

        legend_handles: Dict[str, Any] = {}
        for item in family_items:
            seq_id = item["id"]
            c = coords_by_id[seq_id]
            color, marker, label = get_source_style(item["source"])
            dup_cnt = dup_count_map.get(seq_id, 1)
            pt_size = 45 + 35 * min(dup_cnt - 1, 5)

            sc = ax_pca.scatter(
                c["pca_1"],
                c["pca_2"],
                c=color,
                marker=marker,
                s=pt_size,
                alpha=0.85,
                edgecolors="black",
                linewidths=0.7,
            )
            if label not in legend_handles:
                legend_handles[label] = sc

        ax_pca.set_title(f"PCA Baseline\n(Variance: PC1={var1_str}, PC2={var2_str})", fontsize=11, fontweight="bold")
        ax_pca.set_xlabel("PC1", fontsize=10)
        ax_pca.set_ylabel("PC2", fontsize=10)
        ax_pca.grid(True, linestyle="--", alpha=0.4)

        # Panels 2-4: t-SNE across the 3 seeds
        for idx, seed in enumerate(seeds):
            ax = axes[idx + 1]
            field_x = f"tsne_perp{perp}_seed{seed}_1"
            field_y = f"tsne_perp{perp}_seed{seed}_2"
            tw_key = f"perp_{perp}_seed_{seed}"
            tw = dim_reduction_meta["trustworthiness"].get(tw_key)
            tw_str = f"Trustworthiness: {tw:.3f}" if tw is not None else "Trustworthiness: N/A"

            for item in family_items:
                seq_id = item["id"]
                c = coords_by_id[seq_id]
                color, marker, label = get_source_style(item["source"])
                dup_cnt = dup_count_map.get(seq_id, 1)
                pt_size = 45 + 35 * min(dup_cnt - 1, 5)

                ax.scatter(
                    c.get(field_x, 0.0),
                    c.get(field_y, 0.0),
                    c=color,
                    marker=marker,
                    s=pt_size,
                    alpha=0.85,
                    edgecolors="black",
                    linewidths=0.7,
                )

            ax.set_title(f"t-SNE Seed {seed} (Perp {perp})\n{tw_str}", fontsize=11, fontweight="bold")
            ax.set_xlabel("t-SNE Dimension 1", fontsize=10)
            ax.set_ylabel("t-SNE Dimension 2", fontsize=10)
            ax.grid(True, linestyle="--", alpha=0.4)

        # Add single unified legend
        fig.legend(
            legend_handles.values(),
            legend_handles.keys(),
            loc="lower center",
            bbox_to_anchor=(0.5, -0.08),
            ncol=min(len(legend_handles), 5),
            frameon=True,
            fontsize=10,
        )

        # Scientific caveat subtitle
        caveat_text = (
            f"Family Class: {family_class} (N={len(family_items)}). "
            "Point size reflects exact sequence duplicate count. PCA initialization used across all seeds.\n"
            "CAVEAT: Manifold projections are descriptive visualizations only. "
            "Do NOT use 2D cluster area or distance to quantify biological novelty."
        )
        fig.suptitle(caveat_text, fontsize=10, style="italic", y=1.06)

        safe_fam = family_class.replace("/", "_").replace(" ", "_")
        pdf_path = output_dir / f"tsne_multiseed_{safe_fam}_perp{perp}.pdf"
        png_path = output_dir / f"tsne_multiseed_{safe_fam}_perp{perp}.png"

        fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def plot_family_distance_diagnostics(
    family_class: str,
    distance_records: List[Dict[str, Any]],
    dup_count_map: Dict[str, int],
    output_dir: Path,
) -> None:
    """Plot high-dimensional nearest-reference cosine distances and length-disparity diagnostics."""
    fig, (ax_box, ax_diag) = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)

    # Separate references vs queries
    ref_dists = [r["cosine_distance"] for r in distance_records if r["is_reference"] and not math.isnan(r["cosine_distance"])]
    query_records = [r for r in distance_records if not r["is_reference"] and not math.isnan(r["cosine_distance"])]
    query_dists = [r["cosine_distance"] for r in query_records]

    # Panel A: Nearest-Reference Cosine Distance Distribution
    box_data = []
    box_labels = []
    if ref_dists:
        box_data.append(ref_dists)
        box_labels.append("Curated Ref\n(Leave-One-Out)")
    if query_dists:
        # Group queries by source
        queries_by_src: Dict[str, List[float]] = {}
        for r in query_records:
            queries_by_src.setdefault(r["source"], []).append(r["cosine_distance"])
        for src, vals in queries_by_src.items():
            box_data.append(vals)
            box_labels.append(f"{src}\n(vs Curated Ref)")

    if box_data:
        bp = ax_box.boxplot(box_data, tick_labels=box_labels, patch_artist=True, widths=0.5)
        colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
        for patch, color in zip(bp["boxes"], colors * 3):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        # Add jittered points
        for grp_idx, vals in enumerate(box_data, start=1):
            jitter = np.random.normal(0, 0.04, size=len(vals)) if len(vals) > 1 else [0]
            ax_box.scatter(
                [grp_idx + j for j in jitter],
                vals,
                color="black",
                alpha=0.7,
                s=35,
                zorder=3,
            )

        ax_box.set_ylabel("High-Dimensional Cosine Distance", fontsize=10)
        ax_box.set_title("Nearest-Reference Cosine Distance\n(Baseline vs Query Groups)", fontsize=11, fontweight="bold")
        ax_box.grid(True, linestyle="--", alpha=0.4)
    else:
        ax_box.text(0.5, 0.5, "Insufficient distance pairs\nto construct boxplot", ha="center", va="center")
        ax_box.set_title("Nearest-Reference Cosine Distance", fontsize=11)

    # Panel B: Distance vs Length Difference Diagnostic
    diag_pts = [r for r in distance_records if not math.isnan(r["cosine_distance"]) and r["length_difference"] != "N/A"]
    if diag_pts:
        legend_handles: Dict[str, Any] = {}
        for r in diag_pts:
            color, marker, label = get_source_style(r["source"])
            seq_id = r["id"]
            dup_cnt = dup_count_map.get(seq_id, 1)
            pt_size = 45 + 35 * min(dup_cnt - 1, 5)

            sc = ax_diag.scatter(
                r["length_difference"],
                r["cosine_distance"],
                c=color,
                marker=marker,
                s=pt_size,
                alpha=0.85,
                edgecolors="black",
                linewidths=0.7,
            )
            if label not in legend_handles:
                legend_handles[label] = sc

        ax_diag.set_xlabel("|Length - Nearest Reference Length| (aa)", fontsize=10)
        ax_diag.set_ylabel("High-Dimensional Cosine Distance", fontsize=10)
        ax_diag.set_title("Distance vs Length Disparity Diagnostic\n(Detects Truncation / Fragment Artifacts)", fontsize=11, fontweight="bold")
        ax_diag.grid(True, linestyle="--", alpha=0.4)
        ax_diag.legend(legend_handles.values(), legend_handles.keys(), loc="best", fontsize=9)
    else:
        ax_diag.text(0.5, 0.5, "No length delta pairs available", ha="center", va="center")
        ax_diag.set_title("Distance vs Length Disparity Diagnostic", fontsize=11)

    # Scientific Caveat Subtitle
    caveat_text = (
        f"Family Class: {family_class} (Total N={len(distance_records)}). "
        "High-dimensional cosine distance computed directly in 1024-d ProstT5 space.\n"
        "CAVEAT: Quantitative distances do NOT prove functional novelty or classification. "
        "Distances are strictly bounded by reference database density, length disparities, and taxonomy sampling."
    )
    fig.suptitle(caveat_text, fontsize=9.5, style="italic", y=1.06)

    safe_fam = family_class.replace("/", "_").replace(" ", "_")
    pdf_path = output_dir / f"distance_diagnostics_{safe_fam}.pdf"
    png_path = output_dir / f"distance_diagnostics_{safe_fam}.png"

    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_distance_table(all_distances: List[Dict[str, Any]], output_path: Path) -> None:
    """Write comprehensive distance table with sequence duplicate details."""
    fieldnames = [
        "family_class",
        "id",
        "source",
        "length",
        "phylum",
        "genus",
        "nearest_reference_id",
        "nearest_reference_source",
        "nearest_reference_length",
        "length_difference",
        "cosine_distance",
        "is_reference",
        "is_duplicate",
        "duplicate_group_id",
        "duplicate_count",
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in all_distances:
            writer.writerow({
                "family_class": row["family_class"],
                "id": row["id"],
                "source": row["source"],
                "length": row["length"],
                "phylum": row["phylum"],
                "genus": row["genus"],
                "nearest_reference_id": row["nearest_reference_id"],
                "nearest_reference_source": row["nearest_reference_source"],
                "nearest_reference_length": row["nearest_reference_length"],
                "length_difference": row["length_difference"],
                "cosine_distance": f"{row['cosine_distance']:.6f}" if not math.isnan(row["cosine_distance"]) else "NaN",
                "is_reference": row["is_reference"],
                "is_duplicate": row["is_duplicate"],
                "duplicate_group_id": row["duplicate_group_id"],
                "duplicate_count": row["duplicate_count"],
            })


def write_coordinates_table(
    all_coords: List[Dict[str, Any]], fieldnames: List[str], output_path: Path
) -> None:
    """Write coordinates TSV containing PCA and multi-seed t-SNE projections."""
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in all_coords:
            out_row = {}
            for k in fieldnames:
                v = row.get(k, "")
                if isinstance(v, float):
                    out_row[k] = f"{v:.6f}"
                else:
                    out_row[k] = str(v)
            writer.writerow(out_row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Standalone reproducible CPU analysis and plotting for protein family embeddings."
    )
    parser.add_argument(
        "--embeddings",
        type=Path,
        required=True,
        help="Path to NPZ embeddings archive generated by prostt5_embed.py.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        required=True,
        help="Path to TSV metadata file with exact unique ID rows and taxonomic/source columns.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Path to directory where tables, plots, and provenance will be written.",
    )
    parser.add_argument(
        "--seeds",
        type=str,
        default="42,123,20261015",
        help="Comma-separated random seeds for t-SNE (default: 42,123,20261015).",
    )
    parser.add_argument(
        "--perplexities",
        type=str,
        default="5,30,50",
        help="Comma-separated candidate perplexities for t-SNE (default: 5,30,50).",
    )
    parser.add_argument(
        "--tsne-max-iter",
        type=int,
        default=1000,
        help="Maximum iterations for t-SNE optimization (default: 1000).",
    )
    args = parser.parse_args()

    start_time = time.time()
    script_path = Path(__file__).resolve()
    script_sha256 = compute_file_sha256(script_path)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    seeds = [int(s.strip()) for s in args.seeds.split(",") if s.strip()]
    candidate_perplexities = [float(p.strip()) for p in args.perplexities.split(",") if p.strip()]

    print("================================================================================")
    print("Family Embedding Analysis & Diagnostic Plotting (CPU)")
    print(f"Timestamp UTC: {datetime.now(timezone.utc).isoformat()}")
    print(f"Script:        {script_path} (SHA256: {script_sha256})")
    print(f"Embeddings:    {args.embeddings}")
    print(f"Metadata:      {args.metadata}")
    print(f"Output Dir:    {args.output_dir}")
    print(f"Seeds:         {seeds}")
    print(f"Perplexities:  {candidate_perplexities}")
    print("================================================================================")

    # 1. Load inputs and compute input hashes
    embeddings_sha256 = compute_file_sha256(args.embeddings)
    metadata_sha256 = compute_file_sha256(args.metadata)

    emb_data = load_embeddings(args.embeddings)
    meta_records = load_metadata(args.metadata)

    print(f"Loaded {len(emb_data['ids'])} embeddings (1024-d) from {args.embeddings.name}")
    print(f"Loaded {len(meta_records)} metadata records from {args.metadata.name}")

    # 2. Strict ID closure and join verification
    joined_data = verify_id_closure_and_joins(emb_data, meta_records)
    print(f"Verified 1-to-1 ID closure and sequence hash joins for all {len(joined_data)} items.")

    # 3. Group by family_class to keep unrelated folds/classes strictly separate
    families: Dict[str, List[Dict[str, Any]]] = {}
    for item in joined_data.values():
        fam = item["family_class"]
        families.setdefault(fam, []).append(item)

    print(f"Identified {len(families)} distinct family class(es): {list(families.keys())}")
    print("Unrelated folds/classes will be analyzed and plotted independently (no mixed manifold artifacts).")

    all_distance_records: List[Dict[str, Any]] = []
    all_coord_records: List[Dict[str, Any]] = []
    provenance_families_meta: Dict[str, Any] = {}

    all_coord_fieldnames_set: Set[str] = {
        "family_class",
        "id",
        "source",
        "length",
        "phylum",
        "genus",
        "pca_1",
        "pca_2",
        "pca_var_ratio_1",
        "pca_var_ratio_2",
    }

    # 4. Analyze each family class independently
    for fam_name, fam_items in families.items():
        print(f"\n--- Processing Family Class: {fam_name} (N={len(fam_items)}) ---")

        # Detect exact duplicates
        dup_group_map, dup_count_map = detect_exact_duplicates(fam_items)
        n_unique_seqs = len(set(dup_group_map.values()))
        n_dup_seqs = sum(1 for c in dup_count_map.values() if c > 1)
        print(f"  Exact duplicates check: {n_unique_seqs} unique sequence clusters, {n_dup_seqs} duplicate instances.")

        # Compute high-dimensional cosine nearest-reference distances
        dist_records = compute_cosine_distances(fam_items)
        for r in dist_records:
            seq_id = r["id"]
            r["is_duplicate"] = dup_count_map[seq_id] > 1
            r["duplicate_group_id"] = dup_group_map[seq_id]
            r["duplicate_count"] = dup_count_map[seq_id]
            all_distance_records.append(r)

        n_refs = sum(1 for r in dist_records if r["is_reference"])
        n_queries = len(dist_records) - n_refs
        print(f"  Sequences: {n_refs} curated references, {n_queries} queries.")

        # Dimensionality reduction (PCA + multi-seed t-SNE)
        coords_by_id, dim_meta = run_dimensionality_reduction(
            fam_items, candidate_perplexities, seeds, tsne_max_iter=args.tsne_max_iter
        )

        for c_dict in coords_by_id.values():
            all_coord_records.append(c_dict)
            all_coord_fieldnames_set.update(c_dict.keys())

        # Generate per-family plots
        print("  Generating multi-seed manifold projection plots...")
        plot_family_tsne_multiseed(
            fam_name,
            fam_items,
            coords_by_id,
            dup_count_map,
            dim_meta["perplexities_run"],
            seeds,
            dim_meta,
            args.output_dir,
        )

        print("  Generating distance and length-disparity diagnostic plots...")
        plot_family_distance_diagnostics(
            fam_name, dist_records, dup_count_map, args.output_dir
        )

        # Summary statistics for provenance
        ref_dists = [r["cosine_distance"] for r in dist_records if r["is_reference"] and not math.isnan(r["cosine_distance"])]
        query_dists = [r["cosine_distance"] for r in dist_records if not r["is_reference"] and not math.isnan(r["cosine_distance"])]

        provenance_families_meta[fam_name] = {
            "n_samples": len(fam_items),
            "n_curated_references": n_refs,
            "n_queries": n_queries,
            "n_duplicate_instances": n_dup_seqs,
            "perplexities_run": dim_meta["perplexities_run"],
            "trustworthiness": dim_meta["trustworthiness"],
            "mean_reference_leave_one_out_cosine_distance": float(np.mean(ref_dists)) if ref_dists else None,
            "mean_query_nearest_reference_cosine_distance": float(np.mean(query_dists)) if query_dists else None,
        }

    # 5. Write distance table
    distances_file = args.output_dir / "distances.tsv"
    write_distance_table(all_distance_records, distances_file)
    print(f"\nEmitted raw distance table: {distances_file}")

    # 6. Write coordinates table
    coord_fieldnames = [
        "family_class",
        "id",
        "source",
        "length",
        "phylum",
        "genus",
        "pca_1",
        "pca_2",
        "pca_var_ratio_1",
        "pca_var_ratio_2",
    ]
    # Add tsne columns in sorted order
    tsne_cols = sorted([c for c in all_coord_fieldnames_set if c.startswith("tsne_")])
    coord_fieldnames.extend(tsne_cols)

    coords_file = args.output_dir / "coordinates.tsv"
    write_coordinates_table(all_coord_records, coord_fieldnames, coords_file)
    print(f"Emitted coordinates table:    {coords_file}")

    # 7. Write analysis provenance
    end_time = time.time()
    provenance = {
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
            "library_versions": {
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "sklearn": sklearn.__version__,
                "matplotlib": matplotlib.__version__,
            },
        },
        "inputs": {
            "embeddings": {
                "path": str(args.embeddings),
                "sha256": embeddings_sha256,
                "n_samples": len(emb_data["ids"]),
                "dimensions": 1024,
            },
            "metadata": {
                "path": str(args.metadata),
                "sha256": metadata_sha256,
                "n_records": len(meta_records),
            },
        },
        "parameters": {
            "seeds": seeds,
            "candidate_perplexities": candidate_perplexities,
            "tsne_max_iter": args.tsne_max_iter,
            "tsne_init": "pca",
            "distance_metric": "cosine",
        },
        "scientific_caveats": [
            "Quantitative high-dimensional cosine distances are descriptive metrics and do NOT prove functional novelty or classification.",
            "Distance distributions are strictly bounded by reference database density, sequence length disparities, and taxonomic sampling.",
            "t-SNE coordinates are descriptive manifold visualizations; do NOT use 2D cluster area or spacing to quantify biological novelty.",
            "Unrelated protein folds/classes must not be pooled together to avoid manufactured separation artifacts.",
        ],
        "families": provenance_families_meta,
    }

    provenance_file = args.output_dir / "analysis_provenance.json"
    with open(provenance_file, "w", encoding="utf-8") as f:
        json.dump(provenance, f, indent=2)
    print(f"Emitted analysis provenance:  {provenance_file}")
    print("\nAnalysis and plotting completed successfully.")


if __name__ == "__main__":
    main()
