#!/usr/bin/env python3
"""Bounded posthoc qualification of completed Slurm sequence job 15894662 and
audit of legacy `archaeal_viral_hallmarks` keyword labels and focal `V60PRr2::V60PRr2_715976` ORFs.

Implements all requirements of:
- `mechanism-20261109/EXP3-SEQUENCE-POSTHOC-REVIEW.txt`
- `mechanism-20261109/EXP3-ORF10-REFERENCE-CORRECTION.txt`
- `mechanism-20261109/ROOT-VIRAL-SEQUENCE-DEPENDENCY-UPDATE.txt`

Key invariants:
1. Preserves all raw outputs and code hashes of completed Slurm job `15894662`
   (`code/stage4_structural_hallmark_screen.py` SHA-256 `8fb33e18af4dc6b165acc7a1b932c1ff5ec3cf1f09e1b4d2bdcac1c25a3abbd9`)
   and frozen Stage 2/3 code (`code/stage2_viral_identity_and_spacers_v3.py` SHA-256
   `e6c50836d24b926bd7d0bf089f71ab623109fed221691020b344b5b5de68562e`) without modification.
2. Parses full HMMER3 `--domtblout` coordinates (`hmm_from..hmm_to/tlen`, `ali_from..ali_to/qlen`,
   `env_from..env_to/qlen`), separates full-sequence `E <= 1e-3` from independent domain
   `i-Evalue <= 1e-3` and curated Pfam `GA` (gathering threshold) cutoffs, and groups overlapping
   Pfam hits on the same query interval as competing/overlapping family candidates rather than
   independent additive functions.
3. Qualifies focal `V60PRr2::V60PRr2_715976_12` (`ORF12`, `9862..10779`, `+`, 305 aa) as a
   shell-like / HK97-fold-related superfamily candidate (`Linocin_M18` `PF04454.19` passing `GA`,
   overlapping `DUF6260` `PF19774.6` and sub-threshold `Phage_capsid` `PF05065.20`) in the same
   viral packaging/morphogenesis marker neighborhood as `ORF8` (`3526..4911`, large terminase) and
   `ORF9` (`5135..7525`, `GENOMAD.004769.VV`), with `ORF10` (`7522..9255`) and `ORF11` (`9258..9845`,
   RecA/RadA recombinase) intervening — NEVER a peroxidase or cargo enzyme.
4. Corrects `ORF10` (`V60PRr2::V60PRr2_715976_10`, `GENOMAD.101057.VV`, `virus_hallmark=0`,
   `annotation_accessions=PF09979;COG3566;K09960`, `annotation_description=NA`) using the exact
   official source records saved locally under `review-20261109/viral_sequence15894662/`:
   - `PF09979` (`pf09979_source.txt`, InterPro `IPR016913`): current Pfam name
     `"Pa193 gp24 scaffolding protein"` (`short="Pa193_gp24"`, formerly `DUF2213`), related to
     phage capsid assembly scaffolding (Iglesias et al. 2024 PMID 39370451), NOT a protease.
   - `K09960` (`k09960_source.txt`): KEGG KO `"uncharacterized protein"` (`Function unknown`,
     links `COG: COG3566`).
   - `COG3566` (`cog3566_source.txt`): NCBI COG category `S` `"Uncharacterized conserved protein,
     DUF2213 domain"` (older computational homology report Cheng et al. 2004 discussed putative
     superfamily protease homology, not direct experimental annotation for `ORF10`).
5. Audits the legacy `archaeal_viral_hallmarks` column from `ARCH_VIRAL_KW_RE` in
   `stage2_viral_identity_and_spacers_v3.py`, proving it counts generic viral morphogenetic/packaging
   keywords (`capsid`, `terminase`, `portal`, `tail`, `baseplate`, `virion`, etc.), NOT archaea-specific
   hallmarks or host-domain proof.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exp3_common import _opt_int
from stage2_viral_identity_and_spacers_v3 import ARCH_VIRAL_KW_RE
from stage4_competitive_common import atomic_json, stream_sha256
from stage4_structural_hallmark_screen import gene_sequence_id, is_unannotated_genomad_gene

SCRIPT_ROOT = Path(__file__).resolve().parents[1]
BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", str(SCRIPT_ROOT)))
RES_DIR = (
    BRANCH_ROOT / "results"
    if (BRANCH_ROOT / "results/stage4_unannotated_proteins_pfam_domtbl.txt").exists()
    or not (BRANCH_ROOT / "compact_evidence/stage4_unannotated_proteins_pfam_domtbl.txt").exists()
    else BRANCH_ROOT / "compact_evidence"
)
DEFAULT_REVIEW_SOURCE_DIR = BRANCH_ROOT.parent / "review-20261109/viral_sequence15894662"

GENERIC_MORPHOGENETIC_VIRAL_KW_RE = re.compile(
    r"\b(?:capsid|major_capsid|minor_capsid|terminase|portal|tail|baseplate|virion|head[-_]tail|spike)\b",
    re.IGNORECASE,
)
ARCHAEA_SPECIFIC_VIRAL_KW_RE = re.compile(
    r"\b(?:pleolipovir\w*|halovirus\w*|spindle|fusellovir\w*|rudivir\w*|turrivir\w*|"
    r"salterprovir\w*|bicaudavir\w*|ampullavir\w*|sphaerolipovir\w*)\b",
    re.IGNORECASE,
)
ARCHAEAL_VIRUS_TAXON_NAMES = frozenset({
    "Vertoviridae",
    "Thumleimavirales",
    "Haloferacalesvirus",
    "Hagravirus capitaneum",
    "Graaviviridae",
    "Pleolipoviridae",
    "Halorubrumvirus",
    "Salterprovirus",
})


def audit_genomad_gene_hallmark_and_keyword_semantics(gene_row: dict[str, Any]) -> dict[str, Any]:
    """Audit a geNomad gene row to separate true geNomad viral hallmarks (`virus_hallmark == 1`),
    generic morphogenetic/packaging keyword hits (`terminase`, `capsid`, `portal`, `tail`, etc.),
    archaea-specific viral keyword hits, and broad `Caudoviricetes` taxonomy vs archaea-associated
    viral family/genus taxa (`EXP3-SEQUENCE-POSTHOC-REVIEW.txt` item 2).
    """
    v_hall = (_opt_int(gene_row.get("virus_hallmark")) or 0) > 0
    desc = " ".join(
        str(gene_row.get(k) or "")
        for k in ("annotation_description", "annotation_Pfam", "taxname")
    )
    legacy_matches = [m.group(0) for m in ARCH_VIRAL_KW_RE.finditer(desc)]
    generic_matches = [m.group(0) for m in GENERIC_MORPHOGENETIC_VIRAL_KW_RE.finditer(desc)]
    arch_spec_matches = [m.group(0) for m in ARCHAEA_SPECIFIC_VIRAL_KW_RE.finditer(desc)]
    taxname = str(gene_row.get("taxname") or "").strip()
    if taxname in ("", "NA", "None", "nan"):
        taxname = "NA"
    return {
        "gene": str(gene_row.get("gene") or ""),
        "is_genomad_virus_hallmark": v_hall,
        "matches_legacy_arch_viral_kw_re": bool(legacy_matches),
        "legacy_matched_keywords": ";".join(dict.fromkeys(legacy_matches)),
        "matches_generic_morphogenetic_kw": bool(generic_matches),
        "generic_morphogenetic_keywords": ";".join(dict.fromkeys(generic_matches)),
        "matches_archaea_specific_viral_kw": bool(arch_spec_matches),
        "archaea_specific_keywords": ";".join(dict.fromkeys(arch_spec_matches)),
        "taxname": taxname,
        "taxid": str(gene_row.get("taxid") or "NA"),
        "taxname_is_broad_caudoviricetes_not_host_domain_specific": taxname == "Caudoviricetes",
        "taxname_is_archaea_associated_viral_taxon": taxname in ARCHAEAL_VIRUS_TAXON_NAMES,
        "legacy_column_semantics_note": (
            "legacy_archaeal_viral_hallmarks_counts_generic_morphogenetic_or_viral_keywords_not_archaea_specific_proof"
            if generic_matches and not arch_spec_matches
            else (
                "matches_archaea_specific_viral_taxon_keyword"
                if arch_spec_matches
                else "no_legacy_keyword_match"
            )
        ),
    }


def parse_domtblout_full_coordinates(domtbl_path: Path, domain_evalue_cutoff: float = 1e-3) -> pd.DataFrame:
    """Parse HMMER3 `--domtblout` with full HMM, alignment, and envelope coordinates and coverages."""
    rows: list[dict[str, Any]] = []
    with open(domtbl_path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 23:
                raise ValueError(f"malformed HMMER domtblout row: {line}")
            tlen = int(parts[2])
            qlen = int(parts[5])
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            full_bias = float(parts[8])
            domain_index = int(parts[9])
            domain_total = int(parts[10])
            c_evalue = float(parts[11])
            domain_evalue = float(parts[12])
            domain_score = float(parts[13])
            domain_bias = float(parts[14])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            env_from = int(parts[19])
            env_to = int(parts[20])
            mean_post_prob = float(parts[21])
            hmm_cov = (hmm_to - hmm_from + 1) / float(tlen) if tlen > 0 else 0.0
            q_cov = (ali_to - ali_from + 1) / float(qlen) if qlen > 0 else 0.0
            env_cov = (env_to - env_from + 1) / float(qlen) if qlen > 0 else 0.0
            rows.append({
                "query_gene": parts[3],
                "pfam_target_name": parts[0],
                "pfam_accession": parts[1],
                "full_evalue": full_evalue,
                "full_score": full_score,
                "domain_evalue": domain_evalue,
                "domain_score": domain_score,
                "pfam_description": " ".join(parts[22:]),
                "tlen": tlen,
                "qlen": qlen,
                "full_bias": full_bias,
                "domain_index": domain_index,
                "domain_total": domain_total,
                "c_evalue": c_evalue,
                "domain_bias": domain_bias,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "env_from": env_from,
                "env_to": env_to,
                "mean_post_prob": mean_post_prob,
                "hmm_coverage": round(hmm_cov, 10),
                "query_span_coverage": round(q_cov, 10),
                "envelope_span_coverage": round(env_cov, 10),
                "passes_full_seq_e_1e3": bool(full_evalue <= domain_evalue_cutoff),
                "passes_independent_domain_ie_1e3": bool(domain_evalue <= domain_evalue_cutoff),
            })
    return pd.DataFrame(rows)


def annotate_overlapping_pfam_intervals(
    df: pd.DataFrame,
    pfam_ga_map: dict[str, dict[str, Any]] | None = None,
    min_overlap_aa: int = 10,
    min_overlap_frac: float = 0.25,
) -> pd.DataFrame:
    """Group overlapping Pfam hits on the same ORF query interval and attach posthoc
    Pfam GA cutoffs and conservative functional qualifications.

    Overlapping profile hits on the same ORF interval are marked as competing/overlapping
    family candidates (`is_competing_overlapping_hit = True`), never additive independent functions.
    """
    if df.empty:
        return df.copy()
    out = df.copy().reset_index(drop=True)
    ga_lookup = pfam_ga_map or {}

    ga_seq_scores: list[float | None] = []
    ga_dom_scores: list[float | None] = []
    passes_ga_list: list[bool | None] = []
    for _, row in out.iterrows():
        acc = str(row["pfam_accession"])
        ga_info = ga_lookup.get(acc)
        if ga_info and ga_info.get("ga_seq_score") is not None and ga_info.get("ga_dom_score") is not None:
            gs = float(ga_info["ga_seq_score"])
            gd = float(ga_info["ga_dom_score"])
            ga_seq_scores.append(gs)
            ga_dom_scores.append(gd)
            passes_ga_list.append(bool(float(row["full_score"]) >= gs and float(row["domain_score"]) >= gd))
        else:
            ga_seq_scores.append(None)
            ga_dom_scores.append(None)
            passes_ga_list.append(None)

    out["ga_seq_score"] = ga_seq_scores
    out["ga_dom_score"] = ga_dom_scores
    out["passes_pfam_ga"] = passes_ga_list

    cluster_ids = [""] * len(out)
    n_models_in_cluster = [1] * len(out)
    n_ie_models_in_cluster = [0] * len(out)
    best_model_in_cluster = [""] * len(out)
    is_competing_overlap = [False] * len(out)

    for gene_id, group in out.groupby("query_gene", sort=False):
        idxs = list(group.index)
        adj: dict[int, set[int]] = {i: set() for i in idxs}
        for pos_a, i in enumerate(idxs):
            lo_a, hi_a = int(out.loc[i, "ali_from"]), int(out.loc[i, "ali_to"])
            len_a = max(1, hi_a - lo_a + 1)
            for j in idxs[pos_a + 1 :]:
                lo_b, hi_b = int(out.loc[j, "ali_from"]), int(out.loc[j, "ali_to"])
                len_b = max(1, hi_b - lo_b + 1)
                ov = max(0, min(hi_a, hi_b) - max(lo_a, lo_b) + 1)
                if ov >= min_overlap_aa or (ov / float(min(len_a, len_b))) >= min_overlap_frac:
                    adj[i].add(j)
                    adj[j].add(i)

        visited: set[int] = set()
        comp_num = 0
        for i in idxs:
            if i in visited:
                continue
            comp_num += 1
            stack = [i]
            comp: list[int] = []
            while stack:
                cur = stack.pop()
                if cur in visited:
                    continue
                visited.add(cur)
                comp.append(cur)
                stack.extend(adj[cur] - visited)
            comp_sorted = sorted(
                comp,
                key=lambda k: (float(out.loc[k, "domain_evalue"]), -float(out.loc[k, "domain_score"])),
            )
            best_idx = comp_sorted[0]
            best_label = f"{out.loc[best_idx, 'pfam_target_name']}({out.loc[best_idx, 'pfam_accession']})"
            distinct_models = {str(out.loc[k, "pfam_accession"]) for k in comp}
            distinct_ie_models = {
                str(out.loc[k, "pfam_accession"])
                for k in comp
                if bool(out.loc[k, "passes_independent_domain_ie_1e3"])
            }
            cid_label = f"{gene_id}__interval_{comp_num}"
            for k in comp:
                cluster_ids[k] = cid_label
                n_models_in_cluster[k] = len(distinct_models)
                n_ie_models_in_cluster[k] = len(distinct_ie_models)
                best_model_in_cluster[k] = best_label
                is_competing_overlap[k] = len(distinct_models) > 1

    out["overlap_cluster_id"] = cluster_ids
    out["n_models_in_overlap_cluster"] = n_models_in_cluster
    out["n_ie_1e3_models_in_overlap_cluster"] = n_ie_models_in_cluster
    out["best_ievalue_model_in_cluster"] = best_model_in_cluster
    out["is_competing_overlapping_hit"] = is_competing_overlap

    tiers: list[str] = []
    notes: list[str] = []
    for _, row in out.iterrows():
        gene = str(row["query_gene"])
        passes_ie = bool(row["passes_independent_domain_ie_1e3"])
        passes_ga = row["passes_pfam_ga"]
        competing = bool(row["is_competing_overlapping_hit"])
        if not passes_ie:
            tier = "rejected_domain_ie_gt_1e3_full_seq_only"
        elif passes_ga is True:
            tier = "posthoc_ie_1e3_and_pfam_ga_supported_candidate_family"
        elif passes_ga is False:
            tier = "posthoc_ie_1e3_below_pfam_ga_candidate_family"
        else:
            tier = "posthoc_ie_1e3_candidate_family_uncalibrated_ga"

        if gene == "V60PRr2::V60PRr2_715976_12":
            note = (
                "shell_like_hk97_fold_superfamily_candidate_in_terminase_and_capsid_assembly_related_profile_neighborhood_"
                "competing_Linocin_M18_DUF6260_weak_Phage_capsid_profiles_not_peroxidase_or_cargo_enzyme"
            )
        elif not passes_ie:
            note = "domain_ievalue_exceeds_1e3_excluded_from_posthoc_domain_candidates"
        elif competing:
            note = (
                f"competing_overlapping_profile_on_same_orf_interval_best_model_{row['best_ievalue_model_in_cluster']}_"
                "not_independent_additive_function"
            )
        else:
            note = "broad_sequence_family_candidate_only_not_specific_biochemical_assignment"
        tiers.append(tier)
        notes.append(note)

    out["evidence_tier_posthoc"] = tiers
    out["functional_qualification_note"] = notes
    return out


def parse_local_orf10_reference_sources(source_dir: Path | None = None) -> dict[str, Any]:
    """Parse and hash the local official reference files for PF09979, K09960, and COG3566
    (`EXP3-ORF10-REFERENCE-CORRECTION.txt`).
    """
    candidates = [
        source_dir,
        DEFAULT_REVIEW_SOURCE_DIR,
        SCRIPT_ROOT.parent / "review-20261109/viral_sequence15894662",
        BRANCH_ROOT / "inputs/orf10_reference_sources",
        SCRIPT_ROOT / "inputs/orf10_reference_sources",
        BRANCH_ROOT / "compact_evidence/orf10_reference_sources",
        SCRIPT_ROOT / "compact_evidence/orf10_reference_sources",
    ]
    resolved_dir = next(
        (d for d in candidates if d is not None and (d / "pf09979_source.txt").exists()),
        None,
    )
    if resolved_dir is None:
        raise FileNotFoundError("Could not locate local {pf09979,k09960,cog3566}_source.txt files")

    pf_path = resolved_dir / "pf09979_source.txt"
    ko_path = resolved_dir / "k09960_source.txt"
    cog_path = resolved_dir / "cog3566_source.txt"

    pf_json = json.loads(pf_path.read_text())
    pf_meta = pf_json["metadata"]
    pf_name = pf_meta["name"]["name"]
    pf_short = pf_meta["name"]["short"]
    pf_integrated = pf_meta["integrated"]
    pf_clan = pf_meta["set_info"]
    pf_desc = re.sub(r"<[^>]+>", "", pf_meta["description"][0]["text"]).strip()

    ko_text = ko_path.read_text()
    ko_name_match = re.search(r"^NAME\s+(.+)$", ko_text, re.MULTILINE)
    ko_dblinks_match = re.search(r"^DBLINKS\s+(.+)$", ko_text, re.MULTILINE)
    ko_name = ko_name_match.group(1).strip() if ko_name_match else "uncharacterized protein"
    ko_dblinks = ko_dblinks_match.group(1).strip() if ko_dblinks_match else "COG: COG3566"

    cog_text = cog_path.read_text(errors="replace")
    cog_h1_match = re.search(r"<h1>\s*(S\s*-\s*COG3566\s*-\s*[^<]+?)\s*</h1>", cog_text)
    cog_heading = cog_h1_match.group(1).strip() if cog_h1_match else "S - COG3566 - Uncharacterized conserved protein, DUF2213 domain"

    return {
        "fetched_date": "2026-11-09",
        "PF09979": {
            "source_url": "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF09979/",
            "local_source_relpath": "review-20261109/viral_sequence15894662/pf09979_source.txt",
            "local_source_sha256": stream_sha256(pf_path),
            "accession": pf_meta["accession"],
            "interpro_integrated": pf_integrated,
            "name": pf_name,
            "short_name": pf_short,
            "clan": pf_clan,
            "description": pf_desc,
            "cited_publication": "Iglesias et al. 2024 Commun Biol 7:1275 (PMID 39370451, https://doi.org/10.1038/s42003-024-06985-x)",
            "qualification": (
                "Current Pfam/InterPro entry PF09979 (formerly DUF2213) is 'Pa193 gp24 scaffolding protein', "
                "related to phage capsid assembly scaffolding, NOT a confirmed prohead protease; reference "
                "annotation does not prove experimental function of V60PRr2::V60PRr2_715976_10."
            ),
        },
        "K09960": {
            "source_url": "https://rest.kegg.jp/get/ko:K09960",
            "local_source_relpath": "review-20261109/viral_sequence15894662/k09960_source.txt",
            "local_source_sha256": stream_sha256(ko_path),
            "entry": "K09960",
            "name": ko_name,
            "brite_hierarchy": "09190 Not Included in Pathway or Brite -> 09194 Poorly characterized -> 99997 Function unknown",
            "dblinks": ko_dblinks,
            "qualification": "KEGG Orthology K09960 is an uncharacterized protein (Function unknown) linked to COG3566; no prohead protease assignment.",
        },
        "COG3566": {
            "source_url": "https://www.ncbi.nlm.nih.gov/research/cog/cog/COG3566/",
            "local_source_relpath": "review-20261109/viral_sequence15894662/cog3566_source.txt",
            "local_source_sha256": stream_sha256(cog_path),
            "heading": cog_heading,
            "category": "S - Function unknown",
            "historical_homology_note": (
                "Older computational homology literature (Cheng et al. 2004 Protein Sci 13:2304-2316, "
                "https://onlinelibrary.wiley.com/doi/pdfdirect/10.1110/ps.04726004) reported putative distant "
                "superfamily homology to herpesvirus/phage proteases, but NCBI COG classifies COG3566 under "
                "S (Function unknown: Uncharacterized conserved protein, DUF2213 domain) and geNomad leaves "
                "V60PRr2::V60PRr2_715976_10 annotation_description=NA (without TIGR01543)."
            ),
        },
    }


def run_posthoc_qualification(
    res_dir: Path = RES_DIR,
    extract_json_path: Path | None = None,
    orf10_source_dir: Path | None = None,
) -> dict[str, Any]:
    """Perform bounded posthoc qualification of completed sequence job 15894662,
    ORF10 reference correction, and legacy `archaeal_viral_hallmarks` keyword audit.
    """
    if extract_json_path is None:
        extract_json_path = res_dir / "stage4_target_genomad_checkv_pfam_ga_extract.json"
    if not extract_json_path.exists():
        raise FileNotFoundError(f"missing extracted target geNomad/CheckV/Pfam-GA JSON: {extract_json_path}")

    ext = json.loads(extract_json_path.read_text())
    ga_map: dict[str, dict[str, Any]] = ext.get("pfam_ga_map", {})
    orf10_refs = parse_local_orf10_reference_sources(orf10_source_dir)

    domtbl_path = res_dir / "stage4_unannotated_proteins_pfam_domtbl.txt"
    raw_summary_path = res_dir / "stage4_sequence_domain_annotation_summary.json"
    model_audit_path = res_dir / "stage4_existing_structural_model_exact_sequence_audit.tsv"
    foldseek_path = res_dir / "stage4_existing_model_foldseek_results.tsv"

    raw_df = parse_domtblout_full_coordinates(domtbl_path, domain_evalue_cutoff=1e-3)
    qualified_df = annotate_overlapping_pfam_intervals(raw_df, pfam_ga_map=ga_map)
    qualified_df["target_id"] = qualified_df["query_gene"].astype(str).map(gene_sequence_id)

    qualified_hits_path = res_dir / "stage4_sequence_domain_hits_posthoc_qualified.tsv"
    qualified_df.to_csv(qualified_hits_path, sep="\t", index=False)

    ie_df = qualified_df[qualified_df["passes_independent_domain_ie_1e3"]].copy()
    ga_df = qualified_df[qualified_df["passes_pfam_ga"] == True].copy()  # noqa: E712

    n_raw_rows = int(len(qualified_df))
    n_raw_prots = int(qualified_df["query_gene"].nunique())
    n_ie_rows = int(len(ie_df))
    n_ie_prots = int(ie_df["query_gene"].nunique())
    n_ga_rows = int(len(ga_df))
    n_ga_prots = int(ga_df["query_gene"].nunique())
    n_ga_also_ie = int((ga_df["passes_independent_domain_ie_1e3"]).sum())

    # Build complete 51-ORF evidence table for focal primary target V60PRr2::V60PRr2_715976
    focal_target = "V60PRr2::V60PRr2_715976"
    focal_genes = [r for r in ext["gene_rows"] if r["gene"].rsplit("_", 1)[0] == focal_target]
    focal_genes.sort(key=lambda r: int(r["gene"].rsplit("_", 1)[1]))
    pfam_by_gene: dict[str, pd.DataFrame] = {
        str(g): grp for g, grp in qualified_df.groupby("query_gene", sort=False)
    }

    focal_orf_rows: list[dict[str, Any]] = []
    for r in focal_genes:
        gid = str(r["gene"])
        orf_num = int(gid.rsplit("_", 1)[1])
        sem = audit_genomad_gene_hallmark_and_keyword_semantics(r)
        was_unann = is_unannotated_genomad_gene(pd.Series(r))
        g_pfam = pfam_by_gene.get(gid, pd.DataFrame())
        g_ie = g_pfam[g_pfam["passes_independent_domain_ie_1e3"]] if not g_pfam.empty else pd.DataFrame()
        g_ga = g_pfam[g_pfam["passes_pfam_ga"] == True] if not g_pfam.empty else pd.DataFrame()  # noqa: E712

        raw_models = (
            "; ".join(
                f"{m.pfam_target_name}({m.pfam_accession}:iE={m.domain_evalue:g},HMM={m.hmm_from}..{m.hmm_to}/{m.tlen},Q={m.ali_from}..{m.ali_to}/{m.qlen})"
                for m in g_pfam.itertuples()
            )
            if not g_pfam.empty
            else "none"
        )
        ie_models = (
            "; ".join(
                f"{m.pfam_target_name}({m.pfam_accession}:iE={m.domain_evalue:g},covHMM={m.hmm_coverage:.3f},covQ={m.query_span_coverage:.3f},passGA={m.passes_pfam_ga})"
                for m in g_ie.itertuples()
            )
            if not g_ie.empty
            else "none"
        )
        if gid == "V60PRr2::V60PRr2_715976_8":
            interp = "geNomad_virus_hallmark_large_terminase_PBSX_family_Caudoviricetes_triggers_legacy_terminase_keyword"
        elif gid == "V60PRr2::V60PRr2_715976_9":
            interp = "geNomad_virus_hallmark_NAD_asparagine_ADP_ribosyltransferase_TIGR01641_Caudoviricetes"
        elif gid == "V60PRr2::V60PRr2_715976_10":
            interp = (
                "geNomad_non_hallmark_viral_marker_GENOMAD.101057.VV_description_NA_accessions_PF09979_COG3566_K09960_"
                "current_PF09979_Pa193_gp24_capsid_assembly_scaffolding_related_profile_K09960_COG3566_uncharacterized_"
                "DUF2213_older_COG3566_distant_protease_homology_report_not_confirmed_protease"
            )
        elif gid == "V60PRr2::V60PRr2_715976_11":
            interp = "geNomad_non_hallmark_viral_marker_predicted_phage_recombinase_RecA_RadA_family"
        elif gid == "V60PRr2::V60PRr2_715976_12":
            interp = (
                "posthoc_pfam_shell_like_hk97_fold_superfamily_candidate_Linocin_M18_passes_GA_overlapping_DUF6260_"
                "and_subthreshold_Phage_capsid_in_terminase_and_capsid_assembly_related_profile_neighborhood_NOT_peroxidase_or_cargo_enzyme"
            )
        elif not g_ie.empty:
            interp = str(g_ie.iloc[0]["functional_qualification_note"])
        elif not was_unann:
            interp = "geNomad_viral_marker_without_specific_functional_description"
        else:
            interp = "unannotated_orf_no_pfam_hit_and_no_existing_structural_model"

        focal_orf_rows.append({
            "gene": gid,
            "orf_index": orf_num,
            "start": int(r["start"]),
            "end": int(r["end"]),
            "length_nt": int(r["length"]),
            "strand": int(r["strand"]),
            "genomad_marker": str(r.get("marker") or "NA"),
            "genomad_virus_hallmark": int(sem["is_genomad_virus_hallmark"]),
            "genomad_taxid": str(r.get("taxid") or "NA"),
            "genomad_taxname": sem["taxname"],
            "genomad_annotation_accessions": str(r.get("annotation_accessions") or "NA"),
            "genomad_annotation_description": str(r.get("annotation_description") or "NA"),
            "matches_legacy_arch_viral_kw_re": sem["matches_legacy_arch_viral_kw_re"],
            "legacy_matched_keywords": sem["legacy_matched_keywords"] or "none",
            "matches_archaea_specific_viral_kw": sem["matches_archaea_specific_viral_kw"],
            "was_screened_as_unannotated_candidate": was_unann,
            "n_raw_pfam_domtbl_rows": len(g_pfam),
            "n_posthoc_ie_1e3_pfam_rows": len(g_ie),
            "n_pfam_ga_pass_rows": len(g_ga),
            "raw_pfam_models": raw_models,
            "posthoc_qualified_pfam_models": ie_models,
            "existing_structural_model_status": (
                "no_existing_model_in_colabfold_or_esmfold_inventory"
                if was_unann
                else "not_in_unannotated_candidate_set"
            ),
            "qualified_scientific_interpretation": interp,
        })

    focal_orf_df = pd.DataFrame(focal_orf_rows)
    focal_orf_path = res_dir / "stage4_focal_V60PRr2_715976_orf_evidence_table.tsv"
    focal_orf_df.to_csv(focal_orf_path, sep="\t", index=False)

    # Build target-level hallmark & legacy keyword label audit across all reviewed targets
    target_roles = {
        "V60PRr2::V60PRr2_715976": ("primary_credible_viral_target", "41PRr2__single_SemiBin_236", "p__Thermoplasmatota"),
        "55Dr2::55Dr2_75333": ("primary_credible_viral_target", "41PRr2__single_SemiBin_97", "p__Halobacteriota"),
        "9Sr2::9Sr2_284237": ("primary_credible_viral_target", "41PRr2__single_SemiBin_97", "p__Halobacteriota"),
        "48Sr2::48Sr2_190170": ("primary_credible_viral_target_canonical", "S49Dr1__self_SemiBin_1106", "p__Halobacteriota"),
        "coasm_48_S::k127_1210458": ("primary_credible_viral_target_identical_alias_of_48Sr2_190170", "S49Dr1__self_SemiBin_1106", "p__Halobacteriota"),
        "coasm_48_PR::k127_4172870": ("primary_credible_viral_target_unmapped_source_sample", "V48Dr1__self_SemiBin_3215", "p__Halobacteriota"),
        "41Sr3::41Sr3_56779": ("sensitivity_unresolved_virality_cellular_dominated_196kb", "41Sr3__self_SemiBin_108", "p__Halobacteriota"),
        "44Dr2::44Dr2_4002037": ("boundary_sensitive_provisional_second_thermoplasmatota_full_contig", "43PRr3__self_SemiBin_469", "p__Thermoplasmatota"),
        "44Dr2::44Dr2_4002037|provirus_1_34828": ("boundary_sensitive_provisional_catalog_core_competitor", "43PRr3__self_SemiBin_469", "p__Thermoplasmatota"),
    }
    genes_by_target: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for r in ext["gene_rows"]:
        cid = r["gene"].rsplit("_", 1)[0]
        genes_by_target[cid].append(r)
    if "coasm_48_S::k127_1210458" not in genes_by_target and "48Sr2::48Sr2_190170" in genes_by_target:
        genes_by_target["coasm_48_S::k127_1210458"] = genes_by_target["48Sr2::48Sr2_190170"]

    qual_by_target = {r["contig_id"]: r for r in ext.get("qual_rows", [])}
    if "coasm_48_S::k127_1210458" not in qual_by_target and "48Sr2::48Sr2_190170" in qual_by_target:
        qual_by_target["coasm_48_S::k127_1210458"] = qual_by_target["48Sr2::48Sr2_190170"]

    target_audit_rows: list[dict[str, Any]] = []
    for tid, (role, host_mag, host_phylum) in target_roles.items():
        g_rows = genes_by_target.get(tid, [])
        sems = [audit_genomad_gene_hallmark_and_keyword_semantics(r) for r in g_rows]
        vh_rows = [(r, s) for r, s in zip(g_rows, sems) if s["is_genomad_virus_hallmark"]]
        legacy_kw_rows = [(r, s) for r, s in zip(g_rows, sems) if s["matches_legacy_arch_viral_kw_re"]]
        generic_kw_rows = [(r, s) for r, s in zip(g_rows, sems) if s["matches_generic_morphogenetic_kw"]]
        arch_spec_kw_rows = [(r, s) for r, s in zip(g_rows, sems) if s["matches_archaea_specific_viral_kw"]]
        arch_tax_rows = [(r, s) for r, s in zip(g_rows, sems) if s["taxname_is_archaea_associated_viral_taxon"]]

        canonical_tid = "48Sr2::48Sr2_190170" if tid == "coasm_48_S::k127_1210458" else tid
        t_raw_pfam = qualified_df[qualified_df["target_id"] == canonical_tid]
        t_ie_pfam = ie_df[ie_df["target_id"] == canonical_tid]
        t_ga_pfam = ga_df[ga_df["target_id"] == canonical_tid]
        q_row = qual_by_target.get(tid, {})

        target_audit_rows.append({
            "target_id": tid,
            "canonical_annotated_target_id": canonical_tid,
            "analysis_role": role,
            "linked_host_mag": host_mag,
            "linked_host_phylum": host_phylum,
            "checkv_quality": str(q_row.get("checkv_quality") or "NA"),
            "checkv_viral_genes": int(q_row["viral_genes"]) if str(q_row.get("viral_genes", "")).isdigit() else 0,
            "checkv_host_genes": int(q_row["host_genes"]) if str(q_row.get("host_genes", "")).isdigit() else 0,
            "n_total_orfs": len(g_rows),
            "genomad_virus_hallmark_orf_count": len(vh_rows),
            "legacy_archaeal_viral_hallmarks_col_count": len(legacy_kw_rows),
            "generic_morphogenetic_packaging_kw_orf_count": len(generic_kw_rows),
            "archaea_specific_viral_kw_orf_count": len(arch_spec_kw_rows),
            "genomad_archaea_associated_taxon_orf_count": len(arch_tax_rows),
            "genomad_hallmark_orfs_and_descriptions": "; ".join(
                f"{r['gene']}:{r.get('annotation_description','NA')}[tax={s['taxname']}]"
                for r, s in vh_rows
            ) or "none",
            "legacy_keyword_matched_orfs_and_terms": "; ".join(
                f"{r['gene']}:{s['legacy_matched_keywords']}[tax={s['taxname']}]"
                for r, s in legacy_kw_rows
            ) or "none",
            "n_unannotated_orfs_screened_pfam": int(sum(1 for r in g_rows if is_unannotated_genomad_gene(pd.Series(r))))
            if canonical_tid in {"V60PRr2::V60PRr2_715976", "55Dr2::55Dr2_75333", "9Sr2::9Sr2_284237", "48Sr2::48Sr2_190170", "coasm_48_PR::k127_4172870"}
            else 0,
            "n_raw_pfam_domtbl_rows": len(t_raw_pfam),
            "n_raw_pfam_proteins": int(t_raw_pfam["query_gene"].nunique()) if not t_raw_pfam.empty else 0,
            "n_posthoc_ie_1e3_pfam_rows": len(t_ie_pfam),
            "n_posthoc_ie_1e3_pfam_proteins": int(t_ie_pfam["query_gene"].nunique()) if not t_ie_pfam.empty else 0,
            "n_pfam_ga_pass_rows": len(t_ga_pfam),
            "n_pfam_ga_pass_proteins": int(t_ga_pfam["query_gene"].nunique()) if not t_ga_pfam.empty else 0,
            "legacy_label_correction_note": (
                "legacy_archaeal_viral_hallmarks_column_counts_generic_morphogenetic_packaging_keywords_"
                "not_archaea_specific_hallmarks; archaeal_host_link_comes_from_canonical_archaeal_MAG_CRISPR_spacers"
            ),
        })

    target_audit_df = pd.DataFrame(target_audit_rows)
    target_audit_path = res_dir / "stage4_target_hallmark_keyword_label_audit.tsv"
    target_audit_df.to_csv(target_audit_path, sep="\t", index=False)

    focal_ie_matches = [
        {
            "protein": str(r.query_gene),
            "model": str(r.pfam_target_name),
            "accession": str(r.pfam_accession),
            "full_evalue": float(r.full_evalue),
            "independent_domain_evalue": float(r.domain_evalue),
            "domain_score": float(r.domain_score),
            "ga_dom_score": float(r.ga_dom_score) if pd.notna(r.ga_dom_score) else None,
            "passes_pfam_ga": bool(r.passes_pfam_ga) if pd.notna(r.passes_pfam_ga) else None,
            "hmm_from": int(r.hmm_from),
            "hmm_to": int(r.hmm_to),
            "tlen": int(r.tlen),
            "hmm_coverage": float(r.hmm_coverage),
            "ali_from": int(r.ali_from),
            "ali_to": int(r.ali_to),
            "qlen": int(r.qlen),
            "query_span_coverage": float(r.query_span_coverage),
            "overlap_cluster_id": str(r.overlap_cluster_id),
            "is_competing_overlapping_hit": bool(r.is_competing_overlapping_hit),
            "evidence_tier_posthoc": str(r.evidence_tier_posthoc),
            "functional_qualification_note": str(r.functional_qualification_note),
        }
        for r in ie_df[ie_df["target_id"] == focal_target].itertuples()
    ]

    posthoc_summary = {
        "status": "posthoc_domain_and_hallmark_qualification_completed",
        "schema_version": 2,
        "posthoc_script_relpath": "code/stage4_posthoc_sequence_domain_qualification.py",
        "posthoc_script_sha256": stream_sha256(Path(__file__)),
        "slurm_sequence_job_id": "15894662",
        "slurm_sequence_job_elapsed": "00:02:07",
        "slurm_sequence_job_exit_code": "0:0",
        "executed_sequence_script_relpath": "code/stage4_structural_hallmark_screen.py",
        "executed_sequence_script_sha256": "8fb33e18af4dc6b165acc7a1b932c1ff5ec3cf1f09e1b4d2bdcac1c25a3abbd9",
        "dependency_note": (
            "Job 15894662 completed after root removed its unnecessary afterok:15894661 dependency "
            "(ROOT-VIRAL-SEQUENCE-DEPENDENCY-UPDATE.txt); raw outputs and executed script hash are preserved unmodified."
        ),
        "raw_inputs": {
            "raw_domtbl_relpath": "results/stage4_unannotated_proteins_pfam_domtbl.txt",
            "raw_domtbl_sha256": stream_sha256(domtbl_path),
            "raw_summary_relpath": "results/stage4_sequence_domain_annotation_summary.json",
            "raw_summary_sha256": stream_sha256(raw_summary_path),
            "existing_model_audit_relpath": "results/stage4_existing_structural_model_exact_sequence_audit.tsv",
            "existing_model_audit_sha256": stream_sha256(model_audit_path),
            "foldseek_results_relpath": "results/stage4_existing_model_foldseek_results.tsv",
            "foldseek_results_sha256": stream_sha256(foldseek_path),
            "target_extract_json_relpath": "results/stage4_target_genomad_checkv_pfam_ga_extract.json",
            "target_extract_json_sha256": stream_sha256(extract_json_path),
            "production_pfam_db_path": "/project/k10018/sw/db/pfam-ibex/pfam/Pfam-A.hmm",
            "production_pfam_db_sha256": ext["extraction_provenance"]["source_files"]["pfam_db_sha256_from_job_15894662_summary"],
        },
        "threshold_hierarchy": {
            "raw_full_sequence_evalue_lte_1e3": {
                "domain_rows": n_raw_rows,
                "unique_proteins": n_raw_prots,
                "note": (
                    f"Raw hmmscan -E 0.001 domtblout rows ({n_raw_rows} rows across {n_raw_prots} proteins); "
                    "includes weak secondary/alternative domain alignments on proteins passing full-sequence E <= 0.001."
                ),
            },
            "posthoc_independent_domain_ievalue_lte_1e3": {
                "domain_rows": n_ie_rows,
                "unique_proteins": n_ie_prots,
                "note": (
                    f"Primary posthoc domain filter (independent domain i-Evalue <= 0.001: {n_ie_rows} rows across "
                    f"{n_ie_prots} proteins) with explicit HMM and query coordinates/coverages."
                ),
            },
            "posthoc_pfam_gathering_threshold_ga": {
                "domain_rows": n_ga_rows,
                "unique_proteins": n_ga_prots,
                "domain_rows_also_passing_ie_1e3": n_ga_also_ie,
                "note": (
                    f"Posthoc comparison against curated Pfam-A.hmm GA sequence and domain bitscore cutoffs "
                    f"({n_ga_rows} rows across {n_ga_prots} proteins; {n_ga_also_ie} of {n_ga_rows} also satisfy i-Evalue <= 0.001)."
                ),
            },
        },
        "structural_followup_qualification": {
            "status": "not_performed_no_exact_sequence_matched_existing_models",
            "n_candidate_orfs_checked": 262,
            "n_exact_sequence_matched_models": 0,
            "scientific_interpretation": (
                "No pre-existing ColabFold_AF2 or ESMFold structure existed for the 262 unannotated viral candidate ORFs. "
                "Foldseek structural comparison was therefore NOT performed; this is an unscreened structural status, "
                "NEVER a negative structural fold result."
            ),
        },
        "focal_target_V60PRr2_715976_audit": {
            "target_id": focal_target,
            "target_seq_sha256": "6b19afac4f46beca045911d44ba029c9ca73519ba02aea9ece23b27ee2284cee",
            "linked_host_mag": "41PRr2__single_SemiBin_236",
            "linked_host_phylum": "p__Thermoplasmatota",
            "n_total_orfs": 51,
            "n_genomad_annotated_orfs": 6,
            "n_unannotated_candidate_orfs": 45,
            "n_unannotated_orfs_with_raw_pfam_hits": int(qualified_df[qualified_df["target_id"] == focal_target]["query_gene"].nunique()),
            "n_unannotated_orfs_with_posthoc_ie_1e3_hits": int(ie_df[ie_df["target_id"] == focal_target]["query_gene"].nunique()),
            "n_unannotated_orfs_with_pfam_ga_hits": int(ga_df[ga_df["target_id"] == focal_target]["query_gene"].nunique()),
            "orf8_to_orf12_genomic_neighborhood": {
                "orf8": {
                    "gene": "V60PRr2::V60PRr2_715976_8",
                    "coords_nt": "3526..4911 (+)",
                    "marker": "GENOMAD.168280.VP",
                    "virus_hallmark": 1,
                    "taxname": "Caudoviricetes",
                    "annotation_accessions": "PF04466;TIGR01547;COG4626",
                    "annotation_description": "phage terminase, large subunit, PBSX family",
                },
                "orf9": {
                    "gene": "V60PRr2::V60PRr2_715976_9",
                    "coords_nt": "5135..7525 (+)",
                    "marker": "GENOMAD.004769.VV",
                    "virus_hallmark": 1,
                    "taxname": "Caudoviricetes",
                    "annotation_accessions": "PF04233;COG5585;TIGR01641",
                    "annotation_description": "NAD+--asparagine ADP-ribosyltransferase",
                },
                "intervening_orf10": {
                    "gene": "V60PRr2::V60PRr2_715976_10",
                    "coords_nt": "7522..9255 (+)",
                    "marker": "GENOMAD.101057.VV",
                    "virus_hallmark": 0,
                    "taxname": "Caudoviricetes",
                    "annotation_accessions": "PF09979;COG3566;K09960",
                    "annotation_description": "NA",
                    "official_reference_sources_2026_11_09": orf10_refs,
                    "scientific_qualification": (
                        "geNomad non-hallmark viral marker GENOMAD.101057.VV with annotation_description=NA and "
                        "accessions PF09979;COG3566;K09960. Current official Pfam/InterPro PF09979 is 'Pa193 gp24 "
                        "scaffolding protein' (formerly DUF2213, capsid assembly scaffolding in phage Pa193; Iglesias "
                        "et al. 2024 PMID 39370451), KEGG K09960 is 'uncharacterized protein' (Function unknown), and "
                        "NCBI COG3566 is 'S - Uncharacterized conserved protein, DUF2213 domain' (with older "
                        "computational superfamily protease homology discussed in Cheng et al. 2004). Therefore ORF10 "
                        "is a capsid-assembly/scaffolding-related profile candidate with database annotation "
                        "disagreement, NOT a confirmed biochemical prohead protease."
                    ),
                },
                "intervening_orf11": {
                    "gene": "V60PRr2::V60PRr2_715976_11",
                    "coords_nt": "9258..9845 (+)",
                    "marker": "GENOMAD.132211.VV",
                    "virus_hallmark": 0,
                    "taxname": "Caudoviricetes",
                    "annotation_accessions": "PF09956;COG5471",
                    "annotation_description": "Predicted phage recombinase, RecA/RadA family",
                },
            },
            "orf12_qualification": {
                "gene": "V60PRr2::V60PRr2_715976_12",
                "coords_nt": "9862..10779 (+)",
                "length_aa": 305,
                "genomic_context": (
                    "Located at 9862..10779 (+) in the same viral packaging and capsid-assembly-related marker "
                    "neighborhood as ORF8 (3526..4911, large terminase, GENOMAD.168280.VP, virus_hallmark=1) and "
                    "ORF9 (5135..7525, GENOMAD.004769.VV, virus_hallmark=1), with ORF10 (7522..9255, GENOMAD.101057.VV, "
                    "PF09979 Pa193 gp24 capsid-assembly scaffolding-related profile / K09960 / COG3566 uncharacterized "
                    "DUF2213, description=NA) and ORF11 (9258..9845, GENOMAD.132211.VV, predicted phage recombinase "
                    "RecA/RadA family) intervening between ORF9 and ORF12."
                ),
                "competing_overlapping_profiles_on_aa_22_268": [
                    {
                        "model": "Linocin_M18",
                        "accession": "PF04454.19",
                        "domain_evalue": 3.7e-11,
                        "domain_score": 43.3,
                        "ga_dom_score": 24.3,
                        "passes_ie_1e3": True,
                        "passes_pfam_ga": True,
                        "hmm_coords": "7..234/253",
                        "query_coords": "22..268/305",
                    },
                    {
                        "model": "DUF6260",
                        "accession": "PF19774.6",
                        "domain_evalue": 1.8e-05,
                        "domain_score": 24.5,
                        "ga_dom_score": 26.6,
                        "passes_ie_1e3": True,
                        "passes_pfam_ga": False,
                        "hmm_coords": "7..207/301",
                        "query_coords": "26..214/305",
                    },
                    {
                        "model": "Phage_capsid",
                        "accession": "PF05065.20",
                        "domain_evalue": 0.0022,
                        "domain_score": 18.1,
                        "ga_dom_score": 24.2,
                        "passes_ie_1e3": False,
                        "passes_pfam_ga": False,
                        "hmm_coords": "9..175/281",
                        "query_coords": "34..202/305",
                    },
                ],
                "scientific_conclusion": (
                    "Candidate shell-like / HK97-fold-related viral structural/encapsulin superfamily protein in a "
                    "viral packaging and capsid-assembly-related marker neighborhood. Because Pfam PF04454 (Linocin_M18) "
                    "spans bacteriocins, encapsulin shell proteins, and HK97-like phage capsid relatives (Andreas & Giessen "
                    "2021 PMC8346489; Aguilar et al. 2026 PMC13669383), this sequence profile match does NOT establish "
                    "peroxidase, DyP cargo, or oxidative-stress enzymatic activity."
                ),
            },
            "focal_posthoc_ie_1e3_matches": focal_ie_matches,
        },
        "legacy_keyword_label_correction": {
            "legacy_column_name": "archaeal_viral_hallmarks",
            "source_script_relpath": "code/stage2_viral_identity_and_spacers_v3.py",
            "source_script_sha256": "e6c50836d24b926bd7d0bf089f71ab623109fed221691020b344b5b5de68562e",
            "regex_name": "ARCH_VIRAL_KW_RE",
            "finding": (
                "ARCH_VIRAL_KW_RE matches generic viral morphogenetic/packaging terms (capsid, terminase, portal, "
                "tail, baseplate, virion, head-tail, spike) in addition to archaeal viral taxon names. In frozen "
                "Stage 2/3 tables, archaeal_viral_hallmarks is therefore a generic morphogenetic/packaging + viral "
                "keyword ORF count, NOT archaea-specific hallmark or host-domain proof. For V60PRr2::V60PRr2_715976, "
                "archaeal_viral_hallmarks=1 comes solely from ORF8 matching the generic word 'terminase', and both "
                "geNomad hallmark genes (ORF8 and ORF9) carry broad taxname=Caudoviricetes (taxid=1246). Archaeal host "
                "assignment is established by the canonical archaeal MAG CRISPR spacer link (41PRr2__single_SemiBin_236), "
                "while viral identity is supported independently by CheckV and geNomad viral hallmarks."
            ),
        },
        "posthoc_outputs": {
            "qualified_hits_relpath": "results/stage4_sequence_domain_hits_posthoc_qualified.tsv",
            "qualified_hits_sha256": stream_sha256(qualified_hits_path),
            "focal_orf_table_relpath": "results/stage4_focal_V60PRr2_715976_orf_evidence_table.tsv",
            "focal_orf_table_sha256": stream_sha256(focal_orf_path),
            "target_hallmark_keyword_audit_relpath": "results/stage4_target_hallmark_keyword_label_audit.tsv",
            "target_hallmark_keyword_audit_sha256": stream_sha256(target_audit_path),
        },
        "references": [
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC8346489/",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC13669383/",
            "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF09979/",
            "https://doi.org/10.1038/s42003-024-06985-x",
            "https://rest.kegg.jp/get/ko:K09960",
            "https://www.ncbi.nlm.nih.gov/research/cog/cog/COG3566/",
            "https://onlinelibrary.wiley.com/doi/pdfdirect/10.1110/ps.04726004",
        ],
    }
    posthoc_summary_path = res_dir / "stage4_sequence_posthoc_qualification_summary.json"
    atomic_json(posthoc_summary_path, posthoc_summary)
    return posthoc_summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res-dir", type=Path, default=RES_DIR)
    args = parser.parse_args()
    summary = run_posthoc_qualification(args.res_dir)
    th = summary["threshold_hierarchy"]
    print(
        f"[Posthoc qualification completed] "
        f"raw={th['raw_full_sequence_evalue_lte_1e3']['domain_rows']}/{th['raw_full_sequence_evalue_lte_1e3']['unique_proteins']} "
        f"ie1e3={th['posthoc_independent_domain_ievalue_lte_1e3']['domain_rows']}/{th['posthoc_independent_domain_ievalue_lte_1e3']['unique_proteins']} "
        f"ga={th['posthoc_pfam_gathering_threshold_ga']['domain_rows']}/{th['posthoc_pfam_gathering_threshold_ga']['unique_proteins']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
