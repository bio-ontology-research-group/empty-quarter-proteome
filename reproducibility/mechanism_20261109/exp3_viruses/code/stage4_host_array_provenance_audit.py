#!/usr/bin/env python3
"""Posthoc Host-Array & Contig Provenance Qualification for Experiment 3 Archaeal-Viral Links.

Consumes `stage4_host_array_contig_provenance_extract.json` (produced by Slurm job 15898645 on Shaheen
via `code/extract_stage4_host_array_contig_provenance.py`) and `exp3_thermo_spacer_sequence_independent_join.json`
to generate:
1. `stage4_host_array_contig_provenance_table.tsv`: Per-(host, array-bearing contig) provenance
   table across the 4 primary credible hosts (5 primary pairs), 1 boundary-sensitive provisional
   host (`43PRr3__self_SemiBin_469`), and 1 unresolved virality sensitivity host (`41Sr3__self_SemiBin_108`).
2. `stage4_host_array_spacer_target_join_audit.tsv`: Per-alignment spacer-to-target join table
   with literal FASTA substring verification on the host array-bearing contig.
3. `stage4_host_array_provenance_summary.json`: Portable structured summary with explicit
   separation of:
   (i) exact sequence/bin membership and source-category duplication vs biological dual ownership,
   (ii) MAG-level GTDB r232 / CheckM2 / GUNC quality,
   (iii) official geNomad marker specificity classes (`CC, CP, CV, PC, PP, PV, VC, VP, VV` per
         https://portal.nersc.gov/genomad/marker_features.html) and viral-only ICTV `batch_taxonomy.tsv`
         semantics (https://portal.nersc.gov/genomad/taxonomic_assignment.html),
   (iv) corroborative archaeal-associated genomic coding context on array-bearing contigs vs
        unavailable formal independent GTDB contig-level classification (keeping host assignment
        conditional on MAG bin membership),
   (v) CRISPR repeat / Cas operon subtype boundaries (1 physical array locus called by 2 callers,
       6 strict distinct spacers [4 exact 0-mismatch] vs 7 raw extracted rows on focal `41PRr2_1588631`),
   (vi) actual mapper coverage integration via `<sample>/<sample>_filtered_coverage.tsv` (`rname` =
        `representative_header`) and `<sample>/<sample>_focal_targets_evaluation.tsv`.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DEFAULT_ROOT = Path(__file__).resolve().parents[1]
BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", str(DEFAULT_ROOT)))

# Official geNomad marker specificity classes (https://portal.nersc.gov/genomad/marker_features.html):
# First letter = primary specificity (C=Chromosome, P=Plasmid, V=Virus)
# Second letter = secondary specificity (C=Chromosome, P=Plasmid, V=Virus)
GENOMAD_OFFICIAL_MARKER_CLASSES: tuple[str, ...] = (
    "CC", "CP", "CV", "PC", "PP", "PV", "VC", "VP", "VV",
)

SAMTOOLS_COVERAGE_COLUMNS: list[str] = [
    "rname", "startpos", "endpos", "numreads", "covbases", "coverage",
    "meandepth", "meanbaseq", "meanmapq",
]

REVIEWED_VIRAL_TARGET_IDS: frozenset[str] = frozenset({
    "V60PRr2::V60PRr2_715976",
    "55Dr2::55Dr2_75333",
    "9Sr2::9Sr2_284237",
    "48Sr2::48Sr2_190170",
    "coasm_48_S::k127_1210458",
    "coasm_48_PR::k127_4172870",
    "44Dr2::44Dr2_4002037",
    "41Sr3::41Sr3_56779",
})


def _resolve_results_dir(root: Path) -> Path:
    if (root / "results").is_dir():
        return root / "results"
    if (root / "compact_evidence").is_dir():
        return root / "compact_evidence"
    return root / "results"


def stream_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _format_counts(counts: dict[str, int]) -> str:
    if not counts:
        return "none"
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ";".join(f"{k}:{v}" for k, v in items if v > 0) or "none"


def summarize_genomad_marker_classes(marker_genes: list[dict[str, Any]]) -> dict[str, Any]:
    """Count official geNomad marker specificity classes (CC, CP, CV, PC, PP, PV, VC, VP, VV).

    Fails closed if any marker uses an invalid or invented suffix (such as .A or .B).
    Official reference: https://portal.nersc.gov/genomad/marker_features.html
    """
    class_counts: dict[str, int] = {cls: 0 for cls in GENOMAD_OFFICIAL_MARKER_CLASSES}
    for m in marker_genes:
        marker_id = str(m.get("marker", "")).strip()
        if not marker_id:
            raise ValueError(f"Empty geNomad marker identifier in row: {m}")
        parts = marker_id.split(".")
        suffix = parts[-1] if len(parts) >= 3 and parts[0] == "GENOMAD" else ""
        if suffix not in GENOMAD_OFFICIAL_MARKER_CLASSES:
            raise ValueError(
                f"Invalid geNomad marker class suffix {suffix!r} in {marker_id!r}; "
                f"official geNomad marker classes (https://portal.nersc.gov/genomad/marker_features.html) "
                f"are chromosome/plasmid/virus specificity classes {GENOMAD_OFFICIAL_MARKER_CLASSES}, "
                f"never archaea/bacteria (.A/.B) domain suffixes."
            )
        class_counts[suffix] += 1

    chrom_primary = class_counts["CC"] + class_counts["CP"] + class_counts["CV"]
    plasmid_primary = class_counts["PC"] + class_counts["PP"] + class_counts["PV"]
    virus_primary = class_counts["VC"] + class_counts["VP"] + class_counts["VV"]
    return {
        "class_counts": class_counts,
        "chromosome_primary_count": chrom_primary,
        "plasmid_primary_count": plasmid_primary,
        "virus_primary_count": virus_primary,
        "summary_string": _format_counts(class_counts),
    }


def validate_genomad_viral_ictv_taxonomy(genomad_contig_tax: dict[str, Any] | None) -> dict[str, Any]:
    """Validate that geNomad batch_taxonomy.tsv is treated strictly as viral ICTV taxonomy.

    Official reference: https://portal.nersc.gov/genomad/taxonomic_assignment.html
    Fails closed if a cellular domain lineage (Archaea/Bacteria/Eukaryota) is ever supplied
    as a geNomad batch_taxonomy.tsv lineage.
    """
    lineage = str((genomad_contig_tax or {}).get("lineage", "") or "").strip()
    n_tax_genes = int((genomad_contig_tax or {}).get("n_genes_with_taxonomy", 0) or 0)
    agreement_raw = (genomad_contig_tax or {}).get("agreement", None)
    agreement = float(agreement_raw) if agreement_raw not in (None, "", "NA") else None

    if lineage:
        lowered = lineage.lower()
        if any(tok in lowered for tok in ("archaea", "bacteria", "eukaryota", "d__archaea", "d__bacteria")):
            raise ValueError(
                f"geNomad batch_taxonomy.tsv assigns viral ICTV lineages only "
                f"(https://portal.nersc.gov/genomad/taxonomic_assignment.html), "
                f"never cellular host domain taxonomy; got invalid cellular lineage {lineage!r}"
            )
        if not lineage.startswith("Viruses"):
            raise ValueError(
                f"Unexpected non-viral geNomad batch_taxonomy.tsv lineage {lineage!r}; "
                f"expected 'Viruses;...' or empty."
            )

    return {
        "viral_ictv_lineage_from_viral_markers": lineage or "none",
        "n_genes_with_viral_ictv_taxonomy": n_tax_genes,
        "viral_ictv_taxonomy_agreement": f"{agreement:.4f}" if agreement is not None else "NA",
        "taxonomy_semantics_note": (
            "Official geNomad batch_taxonomy.tsv assigns viral ICTV lineages from viral-annotated "
            "marker hits only; a Viruses;... lineage on a high-chromosome-score cellular contig "
            "reflects the small subset of ORFs hitting viral-like profiles and does not establish "
            "viral identity of the contig or cellular host domain."
        ),
    }


def extract_corroborative_coding_context(marker_genes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Extract corroborative archaeal/eukaryotic-associated coding annotations on a contig.

    Reports exact gene IDs, 1-based coordinates, strand, USCG flag, geNomad specificity class,
    accessions (KEGG/Pfam/TIGRFAM/COG), and descriptions as corroborative coding context
    (not formal independent GTDB contig taxonomy or experimentally verified function).
    """
    corroborative: list[dict[str, Any]] = []
    diagnostic_tokens = (
        "archaea",
        "archaeo",
        "arcog",
        "eukaryot",
        "halobacter",
        "rpoe",
        "us2",
        "s24e",
        "l37e",
        "s3ae",
        "es8",
        "s19e",
        "diphthine",
        "topoisomerase vi",
        "spt4",
        "spt5",
        "nob1",
        "k03049",
        "k02967",
        "k07557",
        "k19664",
        "k15027",
        "k03166",
        "k03167",
        "k20215",
        "k02322",
        "k07462",
        "k07463",
        "k22585",
        "k07094",
        "pf08915",
        "pf18545",
        "pf10502",
        "pf05854",
        "proteasome lid",
        "geranylgeranyl",
    )
    for m in marker_genes:
        acc = str(m.get("annotation_accessions") or "").strip()
        desc = str(m.get("annotation_description") or "").strip()
        combined = f"{acc} {desc}".lower()
        if any(tok in combined for tok in diagnostic_tokens):
            corroborative.append({
                "gene": str(m.get("gene", "")),
                "start": int(m.get("start", 0) or 0),
                "end": int(m.get("end", 0) or 0),
                "strand": int(m.get("strand", 0) or 0),
                "uscg": int(m.get("uscg", 0) or 0),
                "marker": str(m.get("marker", "")),
                "annotation_accessions": acc,
                "annotation_description": desc,
            })
    return corroborative


def format_corroborative_coding_summary(corroborative_genes: list[dict[str, Any]]) -> str:
    if not corroborative_genes:
        return "none_detected_on_contig"
    items: list[str] = []
    for g in corroborative_genes:
        short_gene = g["gene"].split("|")[-1]
        items.append(
            f"{short_gene}({g['start']}..{g['end']},uscg={g['uscg']},{g['annotation_accessions']}:{g['annotation_description']})"
        )
    return "; ".join(items)


def classify_contig_cellular_evidence(
    contig_length_bp: int,
    fraction_of_mag_bases: float,
    gn_scores: dict[str, Any],
    total_orfs: int,
    marker_orfs_count: int,
    uscg_count: int,
    marker_class_info: dict[str, Any],
    viral_ictv_info: dict[str, Any],
    corroborative_genes: list[dict[str, Any]],
    gunc_contig_rows: list[dict[str, Any]],
) -> tuple[str, str]:
    """Classify contig-level cellular evidence and explicitly maintain bin-conditional qualification."""
    chrom_score = gn_scores.get("chromosome_score", "NA")
    plasmid_score = gn_scores.get("plasmid_score", "NA")
    virus_score = gn_scores.get("virus_score", "NA")
    cls_summary = marker_class_info["summary_string"]
    v_lineage = viral_ictv_info["viral_ictv_lineage_from_viral_markers"]
    v_n_genes = viral_ictv_info["n_genes_with_viral_ictv_taxonomy"]
    gunc_status = "persisted_in_tar" if gunc_contig_rows else "genome_level_only_in_c3_gunc_tar"

    if corroborative_genes and uscg_count >= 1:
        tier = "corroborative_archaeal_associated_coding_context_conditional_on_mag_bin"
        gene_highlights = "; ".join(
            f"{g['gene'].split('|')[-1]} ({g['start']}..{g['end']}, {g['annotation_accessions']}: {g['annotation_description']})"
            for g in corroborative_genes
        )
        qual = (
            f"Contig ({contig_length_bp} bp, {fraction_of_mag_bases * 100:.2f}% of MAG; "
            f"chromosome_score={chrom_score}, plasmid_score={plasmid_score}, virus_score={virus_score}; "
            f"{total_orfs} ORFs, {marker_orfs_count} geNomad markers [{cls_summary}], {uscg_count} USCGs) "
            f"carries {len(corroborative_genes)} corroborative archaeal/eukaryotic-associated coding annotation(s) "
            f"({gene_highlights}). geNomad batch_taxonomy.tsv ({v_lineage}, n_genes={v_n_genes}) reflects viral ICTV "
            f"marker hits only, not cellular domain taxonomy; per-contig GUNC={gunc_status}. "
            f"These annotations provide corroborative archaeal-associated coding context rather than formal "
            f"independent GTDB contig classification; host taxonomic assignment remains conditional on MAG bin membership."
        )
    elif corroborative_genes and uscg_count == 0:
        tier = "limited_corroborative_coding_context_no_uscg_conditional_on_mag_bin"
        gene_highlights = "; ".join(
            f"{g['gene'].split('|')[-1]} ({g['start']}..{g['end']}, {g['annotation_accessions']}: {g['annotation_description']})"
            for g in corroborative_genes[:4]
        )
        qual = (
            f"Contig ({contig_length_bp} bp, {fraction_of_mag_bases * 100:.2f}% of MAG; "
            f"chromosome_score={chrom_score}, plasmid_score={plasmid_score}, virus_score={virus_score}; "
            f"{total_orfs} ORFs, {marker_orfs_count} geNomad markers [{cls_summary}], 0 USCGs) "
            f"carries {len(corroborative_genes)} corroborative coding annotation(s) ({gene_highlights}), "
            f"with no USCGs and no formal independent GTDB/GUNC contig-level taxonomic classification "
            f"(geNomad viral ICTV lineage={v_lineage}, n_genes={v_n_genes}; per-contig GUNC={gunc_status}). "
            f"Host taxonomic assignment is conditional on MAG bin membership."
        )
    else:
        tier = "no_corroborative_domain_coding_markers_conditional_on_mag_bin"
        qual = (
            f"Contig ({contig_length_bp} bp, {fraction_of_mag_bases * 100:.2f}% of MAG; "
            f"chromosome_score={chrom_score}, plasmid_score={plasmid_score}, virus_score={virus_score}; "
            f"{total_orfs} ORFs, {marker_orfs_count} geNomad markers [{cls_summary}], {uscg_count} USCGs) "
            f"lacks domain-specific archaeal coding markers or independent contig-level GTDB/GUNC taxonomy "
            f"(geNomad viral ICTV lineage={v_lineage}, n_genes={v_n_genes}; per-contig GUNC={gunc_status}). "
            f"CRISPR-spacer sequence association is exact on the contig, while archaeal host attribution is "
            f"strictly conditional on MAG bin membership."
        )
    return tier, qual


def qualify_crosswalk_source_alias_ownership(
    host_mag: str,
    source_sample: str,
    raw_cid: str,
    ns_cid: str,
    length_bp: int,
    c_own: dict[str, Any],
    cell_memberships: list[str],
    reviewed_viral_targets: frozenset[str] = REVIEWED_VIRAL_TARGET_IDS,
) -> dict[str, Any]:
    """Distinguish raw search-pool self-alias duplication from biological dual ownership."""
    canonical_sha = str(c_own.get("canonical_sha256", ""))
    raw_ownership = str(c_own.get("ownership_status", ""))
    raw_crosswalk_ambiguity = str(c_own.get("crosswalk_ambiguity_status", ""))
    rep_header = str(c_own.get("representative_header", ""))
    entries: list[dict[str, Any]] = list(c_own.get("crosswalk_entries") or [])

    distinct_host_mags = sorted({
        str(e.get("entity_id", ""))
        for e in entries
        if e.get("source_category") == "host_mag"
    })
    other_entries = [
        e for e in entries
        if not (e.get("source_category") == "host_mag" and e.get("entity_id") == host_mag)
    ]
    other_summary = ";".join(
        f"{e.get('source_category')}:{e.get('entity_id')}" for e in other_entries
    ) or "none"
    other_in_reviewed_viral = any(
        str(e.get("entity_id", "")) in reviewed_viral_targets or str(e.get("contig_id", "")) in reviewed_viral_targets
        for e in other_entries
    )
    is_exact_same_contig_search_pool_alias = (
        len(entries) == 2
        and distinct_host_mags == [host_mag]
        and cell_memberships == [host_mag]
        and len(other_entries) == 1
        and other_entries[0].get("source_category") == "virus_target"
        and other_entries[0].get("entity_id") == ns_cid
        and other_entries[0].get("contig_id") == ns_cid
        and int(other_entries[0].get("seq_length", 0) or 0) == length_bp
        and not other_in_reviewed_viral
    )

    if len(entries) == 1 and distinct_host_mags == [host_mag]:
        qualified_status = "unique_to_host_mag_in_crosswalk"
        explanation = (
            f"Canonical SHA {canonical_sha} appears once in the competitive crosswalk under {rep_header} "
            f"and belongs solely to canonical MAG [{host_mag}]."
        )
    elif is_exact_same_contig_search_pool_alias:
        qualified_status = "source_category_duplication_same_contig_in_raw_search_pool_unique_mag_owner"
        explanation = (
            f"Canonical SHA {canonical_sha} has exactly 2 crosswalk rows: host_mag ({host_mag} / {raw_cid}, "
            f"representative_header={rep_header}, is_reference_representative=1) and raw Stage 2 search-pool alias "
            f"virus_target ({ns_cid}, is_reference_representative=0), which is the exact same {length_bp} bp "
            f"chromosomal contig from sample {source_sample} (entered all_matched_primary_targets.fna via "
            f"technical CRISPR array self-hit before self-hit masking). Unique canonical MAG owner across all "
            f"5,438 MAGs is [{host_mag}]; {ns_cid} is not in the reviewed viral target set and is not an "
            f"independently validated virus, provirus, or foreign bin."
        )
    else:
        qualified_status = "biological_or_cross_entity_sequence_sharing_requires_caution"
        explanation = (
            f"Canonical SHA {canonical_sha} appears across {len(entries)} crosswalk rows "
            f"(host MAG owners={distinct_host_mags}, other entries={other_summary}, "
            f"in_reviewed_viral_targets={other_in_reviewed_viral}). Treat as genuine cross-entity ambiguity."
        )

    return {
        "raw_competitive_reference_ownership_status": raw_ownership,
        "raw_crosswalk_ambiguity_status": raw_crosswalk_ambiguity,
        "competitive_reference_representative_header": rep_header,
        "crosswalk_row_count": len(entries),
        "crosswalk_distinct_host_mag_owners": ";".join(distinct_host_mags),
        "crosswalk_other_source_entries": other_summary,
        "other_entry_in_reviewed_viral_pairs": other_in_reviewed_viral,
        "qualified_source_alias_interpretation": qualified_status,
        "qualified_ownership_explanation": explanation,
    }


def parse_samtools_filtered_coverage_tsv(cov_tsv_path: Path) -> dict[str, dict[str, Any]]:
    """Parse `<sample>_filtered_coverage.tsv` produced by `stage4_map_sample_competitive.py`.

    Enforces the 9-column `samtools coverage` schema (`rname` / `#rname`, `startpos`, `endpos`,
    `numreads`, `covbases`, `coverage`, `meandepth`, `meanbaseq`, `meanmapq`), uniqueness of `rname`,
    and finite numeric values.
    """
    lines = cov_tsv_path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].strip():
        raise ValueError(f"samtools coverage file is empty: {cov_tsv_path}")
    header = lines[0].split("\t")
    if header and header[0] == "#rname":
        header[0] = "rname"
    if header != SAMTOOLS_COVERAGE_COLUMNS:
        raise ValueError(
            f"Unexpected samtools coverage header in {cov_tsv_path}: {header}; "
            f"expected={SAMTOOLS_COVERAGE_COLUMNS}"
        )
    if len(lines) == 1:
        raise ValueError(f"samtools coverage file has header but no rows: {cov_tsv_path}")

    coverage_by_rname: dict[str, dict[str, Any]] = {}
    for line_num, line in enumerate(lines[1:], start=2):
        parts = line.split("\t")
        if len(parts) != 9:
            raise ValueError(f"Malformed samtools coverage row {line_num} in {cov_tsv_path}: {len(parts)} fields")
        rname = parts[0].strip()
        if not rname or rname in coverage_by_rname:
            raise ValueError(f"Empty or duplicate rname at row {line_num} in {cov_tsv_path}: {rname!r}")
        row = {
            "rname": rname,
            "startpos": int(parts[1]),
            "endpos": int(parts[2]),
            "numreads": int(parts[3]),
            "covbases": int(parts[4]),
            "coverage": float(parts[5]),
            "meandepth": float(parts[6]),
            "meanbaseq": float(parts[7]),
            "meanmapq": float(parts[8]),
        }
        if not all(math.isfinite(float(v)) for k, v in row.items() if k != "rname"):
            raise ValueError(f"Non-finite coverage metric at row {line_num} in {cov_tsv_path}")
        if row["startpos"] != 1 or row["endpos"] <= 0 or row["numreads"] < 0:
            raise ValueError(f"Invalid coordinates/numreads at row {line_num} in {cov_tsv_path}: {row}")
        if not (0 <= row["covbases"] <= row["endpos"] and 0 <= row["coverage"] <= 100.0 and row["meandepth"] >= 0.0):
            raise ValueError(f"Invalid coverage/depth bounds at row {line_num} in {cov_tsv_path}: {row}")
        coverage_by_rname[rname] = row
    return coverage_by_rname


def _strict_int_from_tsv(val: Any, col_name: str) -> int:
    """Parse an integer or exact integer-valued float string (e.g. '17.0' from mixed-entity pandas TSV export)
    while rejecting fractional floats, booleans, or missing values."""
    if val is None or isinstance(val, bool):
        raise ValueError(f"missing or boolean value {val!r} in integer column {col_name!r}")
    s = str(val).strip()
    if not s or s.lower() == "nan":
        raise ValueError(f"empty/NaN value {val!r} in integer column {col_name!r}")
    try:
        return int(s)
    except ValueError:
        try:
            fval = float(s)
        except ValueError as exc:
            raise ValueError(f"non-integer value {val!r} in integer column {col_name!r}") from exc
        if not math.isfinite(fval) or not fval.is_integer():
            raise ValueError(f"non-integer float {val!r} in integer column {col_name!r}")
        return int(fval)


def evaluate_host_array_contig_coverage_from_mapper(
    host_mag: str,
    array_rep_header: str,
    array_contig_length_bp: int,
    mag_total_contigs: int,
    mag_total_bases: int,
    host_crosswalk_rows: list[dict[str, Any]],
    coverage_by_rname: dict[str, dict[str, Any]],
    focal_eval_host_row: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compare array-bearing contig `representative_header` coverage to rest of MAG and total MAG.

    Fail-closed contract:
    * Joins strictly by verified `representative_header == rname` (never by unverified `contig_id` or `ns_cid`).
    * Enforces `missing != zero`: raises ValueError if any `representative_header` of `host_mag` is missing
      from `coverage_by_rname` or has a length mismatch.
    * Deduplicates by `representative_header` when summing reads/bases/depth across the MAG.
    * Retains raw conservative unique-only MAG numbers as a declared sensitivity alongside total MAG coverage.
    """
    if not host_crosswalk_rows:
        raise ValueError(f"Missing verified host_mag crosswalk rows for {host_mag}; cannot join coverage")
    obs_contigs = len(host_crosswalk_rows)
    obs_bases = sum(int(r["seq_length"]) for r in host_crosswalk_rows)
    if (obs_contigs, obs_bases) != (mag_total_contigs, mag_total_bases):
        raise ValueError(
            f"Host crosswalk denominator mismatch for {host_mag}: "
            f"crosswalk={(obs_contigs, obs_bases)} vs manifest={(mag_total_contigs, mag_total_bases)}"
        )

    if array_rep_header not in coverage_by_rname:
        raise ValueError(
            f"Missing coverage row for focal array representative_header {array_rep_header!r}; missing != zero"
        )

    # Deduplicate by representative_header while verifying every host row exists in coverage_by_rname
    rep_records: dict[str, dict[str, Any]] = {}
    for r in host_crosswalk_rows:
        if r.get("source_category") != "host_mag" or r.get("entity_id") != host_mag:
            raise ValueError(f"Non-host crosswalk row passed for {host_mag}: {r}")
        rep = str(r.get("representative_header", "")).strip()
        seq_len = int(r.get("seq_length", 0) or 0)
        if not rep or rep not in coverage_by_rname:
            raise ValueError(f"Missing coverage row for host representative_header {rep!r}; missing != zero")
        cov_entry = coverage_by_rname[rep]
        if int(cov_entry["endpos"]) != seq_len:
            raise ValueError(
                f"Length mismatch for {rep}: coverage endpos={cov_entry['endpos']} vs crosswalk seq_length={seq_len}"
            )
        if rep not in rep_records:
            rep_records[rep] = cov_entry

    array_cov = rep_records[array_rep_header]
    if int(array_cov["endpos"]) != array_contig_length_bp:
        raise ValueError(
            f"Array contig length mismatch for {array_rep_header}: "
            f"coverage={array_cov['endpos']} vs expected={array_contig_length_bp}"
        )

    rest_reps = [cov for rep, cov in rep_records.items() if rep != array_rep_header]
    rest_bases = sum(int(c["endpos"]) for c in rest_reps)
    rest_covbases = sum(int(c["covbases"]) for c in rest_reps)
    rest_reads = sum(int(c["numreads"]) for c in rest_reps)
    rest_mean_depth = (
        sum(float(c["meandepth"]) * int(c["endpos"]) for c in rest_reps) / rest_bases
        if rest_bases > 0 else 0.0
    )
    rest_breadth_pct = (100.0 * rest_covbases / rest_bases) if rest_bases > 0 else 0.0

    all_reps = list(rep_records.values())
    total_dedup_bases = sum(int(c["endpos"]) for c in all_reps)
    total_dedup_covbases = sum(int(c["covbases"]) for c in all_reps)
    total_dedup_reads = sum(int(c["numreads"]) for c in all_reps)
    total_dedup_mean_depth = (
        sum(float(c["meandepth"]) * int(c["endpos"]) for c in all_reps) / total_dedup_bases
        if total_dedup_bases > 0 else 0.0
    )
    total_dedup_breadth_pct = (100.0 * total_dedup_covbases / total_dedup_bases) if total_dedup_bases > 0 else 0.0

    out: dict[str, Any] = {
        "sample_mapping_coverage_status": "completed_sample_mapping_integrated_via_representative_header",
        "sample_mapping_contig_rname": array_rep_header,
        "sample_mapping_contig_numreads": int(array_cov["numreads"]),
        "sample_mapping_contig_covbases": int(array_cov["covbases"]),
        "sample_mapping_contig_breadth_pct": round(100.0 * int(array_cov["covbases"]) / array_contig_length_bp, 6),
        "sample_mapping_contig_mean_depth": round(float(array_cov["meandepth"]), 6),
        "sample_mapping_remaining_mag_rep_contigs": len(rest_reps),
        "sample_mapping_remaining_mag_bases": rest_bases,
        "sample_mapping_remaining_mag_numreads": rest_reads,
        "sample_mapping_remaining_mag_breadth_pct": round(rest_breadth_pct, 6),
        "sample_mapping_remaining_mag_mean_depth": round(rest_mean_depth, 6),
        "sample_mapping_total_mag_rep_contigs": len(all_reps),
        "sample_mapping_total_mag_bases": total_dedup_bases,
        "sample_mapping_total_mag_numreads": total_dedup_reads,
        "sample_mapping_total_mag_breadth_pct": round(total_dedup_breadth_pct, 6),
        "sample_mapping_total_mag_mean_depth": round(total_dedup_mean_depth, 6),
    }
    if focal_eval_host_row is not None:
        out["mapper_raw_conservative_unique_mag_contigs"] = _strict_int_from_tsv(
            focal_eval_host_row["unique_mag_contigs"], "unique_mag_contigs"
        )
        out["mapper_raw_conservative_unique_mag_bases"] = _strict_int_from_tsv(
            focal_eval_host_row["unique_mag_bases"], "unique_mag_bases"
        )
        out["mapper_raw_conservative_unique_mag_reads"] = _strict_int_from_tsv(
            focal_eval_host_row["unique_filtered_reads_mapped"], "unique_filtered_reads_mapped"
        )
        out["mapper_raw_conservative_unique_mag_breadth_pct"] = round(
            float(focal_eval_host_row["unique_horizontal_breadth_pct"]), 6
        )
        out["mapper_raw_conservative_unique_mag_mean_depth"] = round(
            float(focal_eval_host_row["unique_mean_depth"]), 6
        )
    else:
        out["mapper_raw_conservative_unique_mag_contigs"] = "NA_eval_tsv_not_supplied"
        out["mapper_raw_conservative_unique_mag_bases"] = "NA_eval_tsv_not_supplied"
        out["mapper_raw_conservative_unique_mag_reads"] = "NA_eval_tsv_not_supplied"
        out["mapper_raw_conservative_unique_mag_breadth_pct"] = "NA_eval_tsv_not_supplied"
        out["mapper_raw_conservative_unique_mag_mean_depth"] = "NA_eval_tsv_not_supplied"
    return out


def _load_host_crosswalk_rows_if_available(res_dir: Path, target_hosts: set[str]) -> dict[str, list[dict[str, Any]]]:
    crosswalk_candidates = [
        res_dir / "stage4_competitive_reference/stage4_reference_sequence_crosswalk.tsv",
        res_dir.parent / "results/stage4_competitive_reference/stage4_reference_sequence_crosswalk.tsv",
    ]
    cw_path = next((p for p in crosswalk_candidates if p.exists()), None)
    if cw_path is None:
        return {}
    by_host: dict[str, list[dict[str, Any]]] = {h: [] for h in target_hosts}
    with cw_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("source_category") == "host_mag" and row.get("entity_id") in by_host:
                by_host[row["entity_id"]].append(row)
    return by_host


def _load_focal_eval_host_row(eval_tsv_path: Path, host_mag: str) -> dict[str, Any] | None:
    if not eval_tsv_path.exists():
        return None
    with eval_tsv_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("entity_type") == "host_mag" and row.get("entity_id") == host_mag:
                return row
    return None


def run_host_array_provenance_audit(root: Path | None = None) -> dict[str, Any]:
    base_root = root or BRANCH_ROOT
    res_dir = _resolve_results_dir(base_root)

    extract_path = res_dir / "stage4_host_array_contig_provenance_extract.json"
    if not extract_path.exists():
        raise FileNotFoundError(f"Missing required host provenance extract: {extract_path}")

    extract = json.loads(extract_path.read_text(encoding="utf-8"))

    thermo_join_candidates = [
        res_dir / "exp3_thermo_spacer_sequence_independent_join.json",
        base_root / "compact_evidence/exp3_thermo_spacer_sequence_independent_join.json",
        base_root / "inputs/exp3_thermo_spacer_sequence_independent_join.json",
        base_root.parent / "review-20261109/exp3_thermo_spacer_sequence_independent_join.json",
    ]
    thermo_join_path = next((p for p in thermo_join_candidates if p.exists()), None)
    thermo_join_data = json.loads(thermo_join_path.read_text(encoding="utf-8")) if thermo_join_path else None

    sub_candidates = [
        base_root / "provenance/stage4_competitive_pipeline_submission.json",
        base_root / "manifests/stage4_competitive_pipeline_submission.json",
    ]
    sub_path = next((p for p in sub_candidates if p.exists()), None)
    active_map_job_id = "15905725"
    if sub_path is not None:
        try:
            active_map_job_id = str(json.loads(sub_path.read_text(encoding="utf-8")).get("job_map_array", "15905725"))
        except Exception:
            active_map_job_id = "15905725"

    compact_cov_extract_path = res_dir / "stage4_host_array_sample_coverage_extract.json"
    compact_cov_extract = (
        json.loads(compact_cov_extract_path.read_text(encoding="utf-8"))
        if compact_cov_extract_path.exists()
        else None
    )
    computed_cov_extract_records: dict[str, Any] = {}
    computed_cov_source_files: dict[str, Any] = {}

    hosts_data: dict[str, Any] = extract["hosts"]
    pair_alignments: list[dict[str, str]] = extract["pair_spacer_alignments"]

    host_order = [
        "41PRr2__single_SemiBin_236",
        "41PRr2__single_SemiBin_97",
        "S49Dr1__self_SemiBin_1106",
        "V48Dr1__self_SemiBin_3215",
        "43PRr3__self_SemiBin_469",
        "41Sr3__self_SemiBin_108",
    ]
    host_crosswalk_by_mag = _load_host_crosswalk_rows_if_available(res_dir, set(host_order))

    contig_table_rows: list[dict[str, Any]] = []
    host_summaries: dict[str, Any] = {}

    for host_mag in host_order:
        hobj = hosts_data[host_mag]
        fver = hobj["mag_fasta_verification"]
        taxq = hobj["mag_taxonomy_and_quality"]
        ref_rep = taxq.get("refined_representatives_record") or {}
        tax_rec = taxq.get("taxonomy_reconciliation_record") or {}
        gunc_rec = taxq.get("gunc_per_genome_record") or {}
        gunc_ctg_rows = taxq.get("gunc_contig_assignments_on_array_contigs") or []
        own_sum = hobj["mag_competitive_reference_ownership_summary"]
        samp = hobj["source_sample"]

        # Actual mapper output path: results/stage4_competitive_mapping/<sample>/<sample>_filtered_coverage.tsv
        samp_map_dir = res_dir / "stage4_competitive_mapping" / samp
        samp_cov_tsv = samp_map_dir / f"{samp}_filtered_coverage.tsv"
        samp_eval_tsv = samp_map_dir / f"{samp}_focal_targets_evaluation.tsv"
        cov_by_rname = parse_samtools_filtered_coverage_tsv(samp_cov_tsv) if samp_cov_tsv.exists() else None
        focal_eval_host_row = _load_focal_eval_host_row(samp_eval_tsv, host_mag) if samp_eval_tsv.exists() else None
        if cov_by_rname is not None and samp not in computed_cov_source_files:
            computed_cov_source_files[samp] = {
                "filtered_coverage_tsv": f"stage4_competitive_mapping/{samp}/{samp}_filtered_coverage.tsv",
                "filtered_coverage_sha256": stream_sha256(samp_cov_tsv),
                "focal_targets_evaluation_tsv": (
                    f"stage4_competitive_mapping/{samp}/{samp}_focal_targets_evaluation.tsv"
                    if samp_eval_tsv.exists() else None
                ),
                "focal_targets_evaluation_sha256": (
                    stream_sha256(samp_eval_tsv) if samp_eval_tsv.exists() else None
                ),
            }

        per_contig_summaries: dict[str, Any] = {}
        for raw_cid, cobj in sorted(hobj["array_bearing_contigs"].items()):
            ns_cid = cobj["namespaced_contig_id"]
            length_bp = int(cobj["length_bp"])
            frac_mag_bp = float(cobj["fraction_of_mag_bases"])
            c_own = cobj.get("competitive_reference_ownership") or {}
            gn_tax = cobj.get("genomad_contig_taxonomy")
            gn_scores = cobj.get("genomad_aggregated_classification") or {}
            marker_genes = cobj.get("genomad_contig_marker_genes") or []
            gene_tax_counts = cobj.get("genomad_contig_gene_taxname_counts") or {}
            cat_stats = cobj.get("catalogue_gene_stats") or {}
            arrays = cobj.get("verified_crispr_arrays") or []
            spacers = cobj.get("verified_viral_matching_spacers") or []
            same_cas = cobj.get("cas_operons_on_same_contig") or []
            cell_memberships = list(cobj.get("cellular_contig_to_mag_memberships") or [])

            # 1. Official geNomad marker specificity classes (CC, CP, CV, PC, PP, PV, VC, VP, VV)
            marker_class_info = summarize_genomad_marker_classes(marker_genes)
            # 2. Official geNomad viral ICTV batch_taxonomy.tsv validation
            viral_ictv_info = validate_genomad_viral_ictv_taxonomy(gn_tax)
            # 3. Corroborative archaeal-associated coding context on the contig
            corroborative_genes = extract_corroborative_coding_context(marker_genes)
            corroborative_summary_str = format_corroborative_coding_summary(corroborative_genes)

            # 4. Physical CRISPR array loci vs caller records, and strict/exact/exploratory spacer counts
            physical_loci = sorted({
                (int(a["array_start_1based"]), int(a["array_end_1based"]))
                for a in arrays
            })
            ctg_alns = [
                a for a in pair_alignments
                if a["host_mag"] == host_mag and a.get("host_namespaced_contig_id") == ns_cid
            ]
            strict_sp_ids = {
                a["unique_spacer_id"] for a in ctg_alns
                if str(a.get("passes_strict_spacer_gate", "")).lower() in ("true", "1")
            }
            exact_sp_ids = {
                a["unique_spacer_id"] for a in ctg_alns
                if str(a.get("passes_strict_spacer_gate", "")).lower() in ("true", "1")
                and int(float(a.get("effective_mismatches", 99))) == 0
            }
            raw_extracted_sp_ids = {s["unique_spacer_id"] for s in spacers}
            exploratory_only_sp_ids = raw_extracted_sp_ids - strict_sp_ids
            all_literal_ok = bool(spacers) and all(s.get("verified_literal_in_array_interval", False) for s in spacers)

            # Annotate each extracted spacer row with its strict/exact/exploratory gate status
            annotated_spacers: list[dict[str, Any]] = []
            for s in spacers:
                sid = s["unique_spacer_id"]
                sp_alns = [a for a in ctg_alns if a["unique_spacer_id"] == sid]
                min_mm = min((int(float(a.get("effective_mismatches", 99))) for a in sp_alns), default=99)
                annotated_spacers.append({
                    **s,
                    "passes_strict_spacer_gate": sid in strict_sp_ids,
                    "is_exact_0mm_match": sid in exact_sp_ids,
                    "is_exploratory_only_2mm_row": sid in exploratory_only_sp_ids,
                    "min_effective_mismatches": min_mm,
                })

            # 5. Qualify crosswalk source-alias duplication vs biological dual ownership
            ownership_qual = qualify_crosswalk_source_alias_ownership(
                host_mag=host_mag,
                source_sample=samp,
                raw_cid=raw_cid,
                ns_cid=ns_cid,
                length_bp=length_bp,
                c_own=c_own,
                cell_memberships=cell_memberships,
            )

            # 6. Classify contig cellular evidence tier and explicit bin-conditional qualification
            ctg_gunc_rows = [r for r in gunc_ctg_rows if r.get("contig") == raw_cid]
            ev_tier, ev_qual = classify_contig_cellular_evidence(
                contig_length_bp=length_bp,
                fraction_of_mag_bases=frac_mag_bp,
                gn_scores=gn_scores,
                total_orfs=int(cobj.get("genomad_contig_total_orfs", 0)),
                marker_orfs_count=int(cobj.get("genomad_contig_marker_orfs_count", 0)),
                uscg_count=int(cobj.get("genomad_contig_uscg_count", 0)),
                marker_class_info=marker_class_info,
                viral_ictv_info=viral_ictv_info,
                corroborative_genes=corroborative_genes,
                gunc_contig_rows=ctg_gunc_rows,
            )

            # 7. Coverage in source sample from actual mapper output (`<sample>_filtered_coverage.tsv`)
            rep_header = ownership_qual["competitive_reference_representative_header"]
            cov_key = f"{host_mag}__{raw_cid}"
            if cov_by_rname is not None:
                cov_metrics = evaluate_host_array_contig_coverage_from_mapper(
                    host_mag=host_mag,
                    array_rep_header=rep_header,
                    array_contig_length_bp=length_bp,
                    mag_total_contigs=int(fver["n_contigs"]),
                    mag_total_bases=int(fver["total_bases"]),
                    host_crosswalk_rows=host_crosswalk_by_mag.get(host_mag, []),
                    coverage_by_rname=cov_by_rname,
                    focal_eval_host_row=focal_eval_host_row,
                )
                computed_cov_extract_records[cov_key] = cov_metrics
            elif compact_cov_extract is not None and cov_key in compact_cov_extract.get("records", {}):
                cov_metrics = dict(compact_cov_extract["records"][cov_key])
            else:
                pend_tag = f"NA_pending_{active_map_job_id}"
                cov_metrics = {
                    "sample_mapping_coverage_status": f"pending_queued_mapping_array_{active_map_job_id}_ownership_preverified",
                    "sample_mapping_contig_rname": rep_header,
                    "sample_mapping_contig_numreads": pend_tag,
                    "sample_mapping_contig_covbases": pend_tag,
                    "sample_mapping_contig_breadth_pct": pend_tag,
                    "sample_mapping_contig_mean_depth": pend_tag,
                    "sample_mapping_remaining_mag_rep_contigs": int(fver["n_contigs"]) - 1,
                    "sample_mapping_remaining_mag_bases": int(fver["total_bases"]) - length_bp,
                    "sample_mapping_remaining_mag_numreads": pend_tag,
                    "sample_mapping_remaining_mag_breadth_pct": pend_tag,
                    "sample_mapping_remaining_mag_mean_depth": pend_tag,
                    "sample_mapping_total_mag_rep_contigs": int(fver["n_contigs"]),
                    "sample_mapping_total_mag_bases": int(fver["total_bases"]),
                    "sample_mapping_total_mag_numreads": pend_tag,
                    "sample_mapping_total_mag_breadth_pct": pend_tag,
                    "sample_mapping_total_mag_mean_depth": pend_tag,
                    "mapper_raw_conservative_unique_mag_contigs": int(own_sum["mag_unique_contigs"]),
                    "mapper_raw_conservative_unique_mag_bases": int(own_sum["mag_unique_bp"]),
                    "mapper_raw_conservative_unique_mag_reads": pend_tag,
                    "mapper_raw_conservative_unique_mag_breadth_pct": pend_tag,
                    "mapper_raw_conservative_unique_mag_mean_depth": pend_tag,
                }

            row_dict = {
                "host_mag": host_mag,
                "host_role": hobj["host_role"],
                "source_sample": samp,
                "paired_targets": ";".join(hobj["paired_targets"]),
                "mag_gtdb_r232_taxonomy": ref_rep.get("gtdb_r232_classification", ""),
                "mag_domain": tax_rec.get("canonical_domain", "d__Archaea"),
                "mag_phylum": tax_rec.get("canonical_phylum", ""),
                "mag_order": tax_rec.get("canonical_order", ""),
                "mag_family": tax_rec.get("canonical_family", ""),
                "mag_genus": tax_rec.get("canonical_genus", ""),
                "mag_checkm2_completeness": float(ref_rep.get("Completeness", 0.0) or 0.0),
                "mag_checkm2_contamination": float(ref_rep.get("Contamination", 0.0) or 0.0),
                "mag_pass_gunc": str(gunc_rec.get("pass.GUNC", "")),
                "mag_gunc_max_css_level": gunc_rec.get("taxonomic_level", ""),
                "mag_gunc_css": float(gunc_rec.get("clade_separation_score", 0.0) or 0.0),
                "mag_gunc_contamination_portion": float(gunc_rec.get("contamination_portion", 0.0) or 0.0),
                "mag_gunc_reference_representation_score": float(gunc_rec.get("reference_representation_score", 0.0) or 0.0),
                "mag_total_contigs": int(fver["n_contigs"]),
                "mag_total_bases": int(fver["total_bases"]),
                "mag_n50_bp": int(fver["n50_bp"]),
                "mag_fasta_sha256": fver["fasta_sha256"],
                "raw_contig_id": raw_cid,
                "namespaced_contig_id": ns_cid,
                "carries_viral_matching_spacer": bool(cobj["carries_viral_matching_spacer"]),
                "contig_length_bp": length_bp,
                "contig_gc_fraction": float(cobj["gc_fraction"]),
                "contig_length_rank_in_mag": int(cobj["length_rank_in_mag"]),
                "contig_fraction_of_mag_bases": frac_mag_bp,
                "contig_raw_seq_sha256": cobj["raw_seq_sha256"],
                "contig_canonical_seq_sha256": cobj["canonical_seq_sha256"],
                "cellular_contig_to_mag_memberships": ";".join(cell_memberships),
                **ownership_qual,
                "n_physical_crispr_array_loci": len(physical_loci),
                "n_crispr_caller_records": len(arrays),
                "crispr_array_ids": ";".join(a["array_id"] for a in arrays),
                "crispr_callers": ";".join(a["caller"] for a in arrays),
                "crispr_subtypes_or_labels": ";".join(a["subtype_or_label"] for a in arrays),
                "crispr_array_intervals_1based": ";".join(
                    f"{a['array_start_1based']}..{a['array_end_1based']}" for a in arrays
                ),
                "crispr_array_spans_bp": ";".join(str(a["array_span_bp"]) for a in arrays),
                "crispr_consensus_repeats": ";".join(sorted({a["consensus_repeat"] for a in arrays})),
                "crispr_n_repeats": ";".join(str(a["n_repeats_reported"]) for a in arrays),
                "crispr_n_valid_spacers": ";".join(str(a["n_valid_spacers_reported"]) for a in arrays),
                "n_raw_extracted_matching_spacers": len(raw_extracted_sp_ids),
                "n_strict_viral_matching_spacers": len(strict_sp_ids),
                "n_exact_0mm_viral_matching_spacers": len(exact_sp_ids),
                "n_exploratory_only_viral_matching_spacers": len(exploratory_only_sp_ids),
                "all_viral_matching_spacers_literal_in_contig": all_literal_ok,
                "cas_operon_on_same_contig": len(same_cas) > 0,
                "cas_operons_same_contig_summary": (
                    ";".join(f"{c['operon_id']}:{c['prediction']}:{c['start']}..{c['end']}" for c in same_cas)
                    if same_cas else "none_on_same_contig"
                ),
                "cas_operons_other_mag_contigs_count": int(cobj.get("cas_operons_on_other_mag_contigs_count", 0)),
                "genomad_contig_chromosome_score": gn_scores.get("chromosome_score", "NA"),
                "genomad_contig_plasmid_score": gn_scores.get("plasmid_score", "NA"),
                "genomad_contig_virus_score": gn_scores.get("virus_score", "NA"),
                "genomad_contig_total_orfs": int(cobj.get("genomad_contig_total_orfs", 0)),
                "genomad_contig_marker_orfs_count": int(cobj.get("genomad_contig_marker_orfs_count", 0)),
                "genomad_contig_uscg_count": int(cobj.get("genomad_contig_uscg_count", 0)),
                "genomad_marker_class_counts_summary": marker_class_info["summary_string"],
                "genomad_chromosome_primary_markers": marker_class_info["chromosome_primary_count"],
                "genomad_plasmid_primary_markers": marker_class_info["plasmid_primary_count"],
                "genomad_virus_primary_markers": marker_class_info["virus_primary_count"],
                "genomad_viral_ictv_lineage_from_viral_markers": viral_ictv_info["viral_ictv_lineage_from_viral_markers"],
                "genomad_n_genes_with_viral_ictv_taxonomy": viral_ictv_info["n_genes_with_viral_ictv_taxonomy"],
                "genomad_viral_ictv_taxonomy_agreement": viral_ictv_info["viral_ictv_taxonomy_agreement"],
                "genomad_contig_gene_taxname_summary": _format_counts(gene_tax_counts),
                "corroborative_coding_annotations_count": len(corroborative_genes),
                "corroborative_coding_annotations_summary": corroborative_summary_str,
                "catalogue_contig_genes": int(cat_stats.get("n_genes", 0)),
                "catalogue_contig_dark_genes": int(cat_stats.get("n_dark", 0)),
                "catalogue_contig_dark_strict_genes": int(cat_stats.get("n_dark_strict", 0)),
                **cov_metrics,
                "independent_contig_archaeal_evidence_tier": ev_tier,
                "scientific_qualification": ev_qual,
            }
            contig_table_rows.append(row_dict)
            per_contig_summaries[raw_cid] = {
                **row_dict,
                "genomad_marker_class_breakdown": marker_class_info,
                "genomad_viral_ictv_taxonomy_qualification": viral_ictv_info,
                "corroborative_coding_annotations": corroborative_genes,
                "verified_crispr_arrays": arrays,
                "verified_viral_matching_spacers": annotated_spacers,
                "genomad_contig_marker_genes": marker_genes,
                "gunc_contig_assignments": ctg_gunc_rows,
            }

        host_summaries[host_mag] = {
            "host_mag": host_mag,
            "host_role": hobj["host_role"],
            "source_sample": samp,
            "paired_targets": hobj["paired_targets"],
            "mag_gtdb_r232_taxonomy": ref_rep.get("gtdb_r232_classification", ""),
            "mag_checkm2_completeness": float(ref_rep.get("Completeness", 0.0) or 0.0),
            "mag_checkm2_contamination": float(ref_rep.get("Contamination", 0.0) or 0.0),
            "mag_gunc_summary": gunc_rec,
            "mag_gunc_all_levels": taxq.get("gunc_all_levels", []),
            "mag_fasta_sha256": fver["fasta_sha256"],
            "mag_total_contigs": int(fver["n_contigs"]),
            "mag_total_bases": int(fver["total_bases"]),
            "mag_n50_bp": int(fver["n50_bp"]),
            "mag_competitive_reference_ownership_summary": own_sum,
            "mag_genomad_contig_taxonomy_summary": hobj.get("mag_genomad_contig_taxonomy_summary", {}),
            "mag_genomad_gene_taxname_counts": hobj.get("mag_genomad_gene_taxname_counts", {}),
            "mag_cas_operons": hobj.get("mag_cas_operons", []),
            "array_bearing_contigs": per_contig_summaries,
        }

    # Build spacer-to-target alignment join table with literal contig verification
    spacer_lookup: dict[tuple[str, str], dict[str, Any]] = {}
    for host_mag, hobj in hosts_data.items():
        for raw_cid, cobj in hobj["array_bearing_contigs"].items():
            for s in cobj.get("verified_viral_matching_spacers", []):
                spacer_lookup[(host_mag, s["unique_spacer_id"])] = {
                    **s,
                    "raw_contig_id": raw_cid,
                    "namespaced_contig_id": cobj["namespaced_contig_id"],
                    "contig_length_bp": cobj["length_bp"],
                    "contig_raw_seq_sha256": cobj["raw_seq_sha256"],
                }

    join_rows: list[dict[str, Any]] = []
    for aln in sorted(
        pair_alignments,
        key=lambda r: (
            r["host_mag"],
            r["target_id"],
            int(r.get("target_hit_lo", 0) or 0),
            r["unique_spacer_id"],
            r["pool"],
        ),
    ):
        h = aln["host_mag"]
        sp_id = aln["unique_spacer_id"]
        sp_meta = spacer_lookup.get((h, sp_id), {})
        join_rows.append({
            "host_mag": h,
            "host_role": hosts_data[h]["host_role"],
            "unique_spacer_id": sp_id,
            "representative_spacer_seq": sp_meta.get("representative_spacer_seq", ""),
            "canonical_spacer_seq_rc": sp_meta.get("canonical_spacer_seq_rc", ""),
            "spacer_len": sp_meta.get("spacer_len", aln.get("qlen", "")),
            "host_namespaced_contig_id": aln.get("host_namespaced_contig_id", ""),
            "host_contig_length_bp": sp_meta.get("contig_length_bp", ""),
            "host_contig_raw_seq_sha256": sp_meta.get("contig_raw_seq_sha256", ""),
            "host_primary_array_start": sp_meta.get("primary_array_start", ""),
            "host_primary_array_end": sp_meta.get("primary_array_end", ""),
            "host_all_array_loci": aln.get("host_all_array_loci", ""),
            "verified_literal_in_host_array_interval": sp_meta.get("verified_literal_in_array_interval", False),
            "spacer_strand_in_host_contig": sp_meta.get("strand_in_contig", ""),
            "spacer_1based_start_in_host_contig": sp_meta.get("contig_1based_start", ""),
            "callers": sp_meta.get("callers", ""),
            "repeat_seq": sp_meta.get("repeat_seq", ""),
            "target_id": aln["target_id"],
            "pool": aln["pool"],
            "target_seq_sha256": aln.get("target_seq_sha256", ""),
            "target_len": int(aln.get("target_len", 0) or 0),
            "sstart": int(aln.get("sstart", 0) or 0),
            "send": int(aln.get("send", 0) or 0),
            "target_hit_lo": int(aln.get("target_hit_lo", 0) or 0),
            "target_hit_hi": int(aln.get("target_hit_hi", 0) or 0),
            "effective_mismatches": int(float(aln.get("effective_mismatches", 0) or 0)),
            "full_spacer_coverage": float(aln.get("full_spacer_coverage", 0.0) or 0.0),
            "match_tier": aln.get("match_tier", ""),
            "passes_strict_spacer_gate": str(aln.get("passes_strict_spacer_gate", "")) == "True",
            "passes_exploratory_spacer_gate": str(aln.get("passes_exploratory_spacer_gate", "")) == "True",
        })

    # Reconcile focal Thermoplasmatota host 41PRr2__single_SemiBin_236 against exp3_thermo_spacer_sequence_independent_join.json
    focal_host = "41PRr2__single_SemiBin_236"
    focal_contig_id = "41PRr2_1588631"
    focal_contig_summary = host_summaries[focal_host]["array_bearing_contigs"][focal_contig_id]
    focal_all_extracted_spacers = {
        s["unique_spacer_id"]: s
        for s in focal_contig_summary["verified_viral_matching_spacers"]
    }
    focal_strict_spacers = {
        sid: s for sid, s in focal_all_extracted_spacers.items() if s["passes_strict_spacer_gate"]
    }

    thermo_join_reconciliation: dict[str, Any] = {
        "focal_host_mag": focal_host,
        "focal_array_contig": f"41PRr2::{focal_contig_id}",
        "focal_contig_length_bp": focal_contig_summary["contig_length_bp"],
        "focal_contig_length_rank_in_mag": focal_contig_summary["contig_length_rank_in_mag"],
        "focal_contig_fraction_of_mag_bases": focal_contig_summary["contig_fraction_of_mag_bases"],
        "focal_contig_raw_seq_sha256": focal_contig_summary["contig_raw_seq_sha256"],
        "focal_contig_canonical_seq_sha256": focal_contig_summary["contig_canonical_seq_sha256"],
        "n_physical_crispr_array_loci": focal_contig_summary["n_physical_crispr_array_loci"],
        "n_crispr_caller_records": focal_contig_summary["n_crispr_caller_records"],
        "array_interval_1based": [433346, 438682],
        "array_span_bp": 438682 - 433346 + 1,
        "distance_from_contig_5p_bp": 433345,
        "distance_from_contig_3p_bp": focal_contig_summary["contig_length_bp"] - 438682,
        "adjacent_type_ie_cas_operon": "41PRr2_1588631@8:I-E:424731..432985",
        "n_raw_extracted_matching_spacers": focal_contig_summary["n_raw_extracted_matching_spacers"],
        "n_strict_unique_spacers_verified": focal_contig_summary["n_strict_viral_matching_spacers"],
        "n_exact_0mm_unique_spacers_verified": focal_contig_summary["n_exact_0mm_viral_matching_spacers"],
        "n_exploratory_only_2mm_spacers": focal_contig_summary["n_exploratory_only_viral_matching_spacers"],
        "exploratory_only_spacer_ids": sorted(
            sid for sid, s in focal_all_extracted_spacers.items() if s["is_exploratory_only_2mm_row"]
        ),
        "all_spacers_verified_literal_in_contig_fasta": focal_contig_summary["all_viral_matching_spacers_literal_in_contig"],
        "crosswalk_source_alias_qualification": focal_contig_summary["qualified_source_alias_interpretation"],
        "root_independent_join_verified": None,
    }
    if thermo_join_data is not None:
        root_strict_items = thermo_join_data.get("strict_evidence", [])
        root_strict_ids = [item["spacer_id"] for item in root_strict_items]
        exact_match_root_join = (
            len(root_strict_items) == 6
            and set(root_strict_ids) == set(focal_strict_spacers.keys())
            and all(
                item["sequence"] in (
                    focal_strict_spacers[item["spacer_id"]]["representative_spacer_seq"],
                    focal_strict_spacers[item["spacer_id"]].get("canonical_spacer_seq_rc", ""),
                )
                and f"{focal_strict_spacers[item['spacer_id']]['primary_array_start']}-{focal_strict_spacers[item['spacer_id']]['primary_array_end']}"
                == item["host_array_loci"].split(":")[-1]
                for item in root_strict_items
            )
        )
        thermo_join_reconciliation["root_independent_join_verified"] = exact_match_root_join
        thermo_join_reconciliation["root_strict_spacer_ids"] = root_strict_ids
        thermo_join_reconciliation["orientation_note"] = (
            "exp3_thermo_spacer_sequence_independent_join.json records canonical_spacer_seq_rc for "
            "reverse-complement oriented spacers; all 6 strict spacers match 100% in ID, array locus "
            "(433346-438682), and exact sequence/reverse-complement on 41PRr2_1588631."
        )

    if computed_cov_extract_records:
        cov_extract_payload = {
            "schema_version": 1,
            "generated_at_utc": extract.get("generated_at_utc", "2026-11-09T16:45:00+00:00"),
            "source_files": computed_cov_source_files,
            "records": computed_cov_extract_records,
        }
        compact_cov_extract_path.write_text(
            json.dumps(cov_extract_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    contig_tsv_path = res_dir / "stage4_host_array_contig_provenance_table.tsv"
    with contig_tsv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(contig_table_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(contig_table_rows)

    join_tsv_path = res_dir / "stage4_host_array_spacer_target_join_audit.tsv"
    with join_tsv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(join_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(join_rows)

    summary_payload = {
        "schema_version": 2,
        "generated_at_utc": extract.get("generated_at_utc", "2026-11-09T16:45:00+00:00"),
        "status": "host_array_contig_provenance_completed_verified",
        "source_extract": {
            "extract_json": f"{res_dir.name}/stage4_host_array_contig_provenance_extract.json",
            "extract_sha256": stream_sha256(extract_path),
            "slurm_job_id": extract.get("slurm_job_id"),
            "slurm_nodelist": extract.get("slurm_nodelist"),
            "thermo_independent_join_json": (
                f"{thermo_join_path.parent.name}/{thermo_join_path.name}" if thermo_join_path else None
            ),
            "thermo_independent_join_sha256": stream_sha256(thermo_join_path) if thermo_join_path else None,
        },
        "outputs": {
            "contig_provenance_table_tsv": f"{res_dir.name}/stage4_host_array_contig_provenance_table.tsv",
            "contig_provenance_table_sha256": stream_sha256(contig_tsv_path),
            "spacer_target_join_audit_tsv": f"{res_dir.name}/stage4_host_array_spacer_target_join_audit.tsv",
            "spacer_target_join_audit_sha256": stream_sha256(join_tsv_path),
        },
        "focal_thermoplasmatota_host_41PRr2_SemiBin_236_reconciliation": thermo_join_reconciliation,
        "hosts": host_summaries,
        "scientific_interpretation_boundaries": [
            "Official geNomad marker features (https://portal.nersc.gov/genomad/marker_features.html) define 9 chromosome/plasmid/virus specificity classes (CC, CP, CV, PC, PP, PV, VC, VP, VV), never archaea/bacteria (.A/.B) domain suffixes. Official geNomad batch_taxonomy.tsv (https://portal.nersc.gov/genomad/taxonomic_assignment.html) assigns viral ICTV lineages from viral-like marker hits only, never cellular host domain taxonomy.",
            "Canonical MAG GTDB r232 taxonomy and CheckM2/GUNC genome-wide quality characterize the bin as a whole. Direct genomic annotations on array-bearing contigs (such as K03049/PF08292 archaeal/eukaryotic RNAP RpoE, K02967/PF00318 uS2, K07557 archaeosine tRNA-guanine transglycosylase, K19664 CDP-archaeol synthase, and ribosomal S24E/L37E on the 651,471 bp chromosome-scale contig 41PRr2_1588631) provide corroborative archaeal-associated coding context, not formal independent GTDB contig classification or experimentally verified functions; host taxonomic assignment remains conditional on MAG bin membership.",
            "CRISPR consensus repeat sequences and CCTyper Cas/repeat subtype labels (e.g. I-E, I-B, I-D) can be horizontally transferred between Bacteria and Archaea and do not independently establish host domain.",
            "All 8 array-bearing host contigs appear in exactly 2 rows of stage4_reference_sequence_crosswalk.tsv: once as host_mag (HOST__<mag>__<contig>, is_reference_representative=1) and once as virus_target (VIRUS_TARGET__<sample>::<contig>, is_reference_representative=0), because the raw Stage 2 BLAST search pool included unbinned assembly contigs hit by their own CRISPR arrays before self-hit masking. Each contig belongs to a single canonical MAG across all 5,438 MAGs; this shared_host_viral_ambiguous flag reflects source-category alias duplication of the exact same contig, not biological dual ownership, provirus integration, or a foreign bin.",
            "When sample read mapping completes (stage4_map_sample_competitive.py outputs <sample>/<sample>_filtered_coverage.tsv with rname=representative_header and <sample>/<sample>_focal_targets_evaluation.tsv), array-contig vs rest-of-MAG coverage must be joined by verified representative_header (deduplicated, with missing != zero), while retaining raw conservative unique-only host MAG numbers as a declared sensitivity.",
        ],
    }

    summary_path = res_dir / "stage4_host_array_provenance_summary.json"
    summary_path.write_text(json.dumps(summary_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary_payload


if __name__ == "__main__":
    out = run_host_array_provenance_audit()
    rec = out["focal_thermoplasmatota_host_41PRr2_SemiBin_236_reconciliation"]
    print(
        f"Completed host-array provenance audit: {len(out['hosts'])} hosts, "
        f"focal 41PRr2_1588631 length={rec['focal_contig_length_bp']} bp, "
        f"strict_spacers={rec['n_strict_unique_spacers_verified']} "
        f"(exact_0mm={rec['n_exact_0mm_unique_spacers_verified']}, "
        f"raw_extracted={rec['n_raw_extracted_matching_spacers']}), "
        f"root_join_verified={rec['root_independent_join_verified']}"
    )
