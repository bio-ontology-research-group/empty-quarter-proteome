#!/usr/bin/env python3
"""Reproducibility verifier for the Empty Quarter MAGs / proteome paper.

For a mixed historical and current registry, this script
extracts or recomputes values from the curated result tables under ../results/ and checksum-pinned compact
inputs under this directory, and checks it against the value printed in the paper.

The active October manuscript driver is darkmatter/run_current.sh; the older registry
retains historical analyses. H1-H4 are active FROZEN entries and are also checked by
darkmatter/check_round7_registered.py.

Four tiers:
  COMPUTE     -- recomputed here from a local analysis table; PASS/FAIL on numeric agreement.
  FROZEN      -- read or derived from a compact SHA-256-pinned result artifact; PASS/FAIL on numeric
                 agreement after the artifact-manifest gate. This reproduces manuscript extraction,
                 not the upstream heavy computation.
  TABLE       -- the value is the recorded statistic of a curated source-of-truth table that cannot be
                 cheaply recomputed locally (e.g. a p-value from a 16-sample GEM run); only source-file
                 presence is checked. Documented, not independently recomputed or value-validated.
  PROVENANCE  -- produced by heavy IBEX/ORIX pipelines (51.9M-protein clustering, GTDB, ESMFold). Not
                 locally re-runnable; a SHA-256-pinned ledger checks each registered value and upstream
                 source without claiming to rerun or revalidate the remote computation.

Exit code is non-zero if any COMPUTE/FROZEN claim fails, a TABLE source is absent, a frozen artifact
changes, or a manuscript boundary guard fails, so this doubles as a CI gate.

Usage:  python3 reproducibility/verify.py [--quiet]
"""
import csv
import hashlib
import json
import math
import os
import re
import statistics as st
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
REPRO = os.path.join(ROOT, "reproducibility")
sys.path.insert(0, os.path.join(REPRO, "darkmatter", "round7_final"))
import verify_round7 as v7
sys.path.insert(0, os.path.join(REPRO, "darkmatter", "family_followup_20261015"))
import verify_family


def rows(path, delim="\t"):
    with open(path) as f:
        return list(csv.reader(f, delimiter=delim))


def drows(path, delim="\t"):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter=delim))


def jread(path):
    with open(path) as f:
        return json.load(f)


# catalogue size N (denominator for prevalence %): unique MAGs in the dark-fraction table
def catalog_N():
    return len({r["mag"] for r in drows(os.path.join(RES, "picrust2_bias", "dark_fraction.tsv"))})


# ---------------------------------------------------------------------------
# COMPUTE recipes -- each returns a float recomputed from a local table.
# ---------------------------------------------------------------------------
def c_mags():
    return float(catalog_N())


def c_catalogue(field):
    data = jread(os.path.join(REPRO, "canonical_catalogue_summary.json"))
    return float(data["counts"][field])


def c_drep_reconciliation(metric):
    values = {
        row["metric"]: row["value"]
        for row in drows(os.path.join(REPRO, "drep_membership_reconciliation.tsv"))
    }
    return float(values[metric])


def c_r232(field):
    data = jread(os.path.join(REPRO, "r232_taxonomic_novelty_summary.json"))
    return float(data[field])


def c_r232_gap(rank, field="all_quality_tiers", cohort="legacy_hq_mq"):
    for row in drows(os.path.join(REPRO, "r232_first_unassigned_census.tsv")):
        if row["cohort"] == cohort and row["first_unassigned_rank"] == rank:
            return float(row[field])
    raise KeyError((cohort, rank, field))


def c_r232_named(rank, field="all_legacy_representatives"):
    for row in drows(os.path.join(REPRO, "r232_named_taxon_census.tsv")):
        if row["rank"] == rank:
            return float(row[field])
    raise KeyError((rank, field))


def c_phylogeny(field):
    data = jread(os.path.join(REPRO, "r232_phylogeny_adjudication_summary.json"))
    return float(data[field])


def c_structure_metric(*keys):
    data = jread(
        os.path.join(REPRO, "current_structure_metric_sensitivity_summary.json")
    )
    for key in keys:
        data = data[key]
    return float(data)


def c_structure_metric_wilson(index):
    """Wilson 95% interval from the frozen qTM control numerator/denominator."""
    k = c_structure_metric(
        "known_fold_sensitivity", "query_normalized",
        "confidence_length_qualified", "false_dual_nohits",
    )
    n = c_structure_metric(
        "known_fold_sensitivity", "query_normalized",
        "confidence_length_qualified", "controls",
    )
    z = 1.959963984540054
    p = k / n
    denominator = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / denominator
    margin = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denominator
    return round(100.0 * (centre - margin, centre + margin)[index], 2)


def c_tier_snapshot_nonclosure():
    data = {
        row["quantity"]: int(row["count"])
        for row in drows(
            os.path.join(REPRO, "historical_tier_snapshot_accounting.tsv")
        )
    }
    return float(
        data["tier1_historical_input"]
        + data["tier2_historical_input"]
        - data["assembly_catalogue_uniref50_dark"]
    )


def c_orthogroups():
    data = drows(os.path.join(RES, "darkmatter_ecology", "darkfam_meta.tsv"))
    ids = [row["family_id"] for row in data]
    if len(ids) != len(set(ids)):
        raise ValueError("dark-family metadata contains duplicate family IDs")
    return float(len(ids))


def _e2_summary(metric):
    path = os.path.join(ROOT, "RQ", "E2", "table_positive_families_summary.tsv")
    for row in drows(path):
        if row["metric"] == metric:
            return float(row["value"])
    raise KeyError(metric)


def _e2_tier(tier):
    path = os.path.join(ROOT, "RQ", "E2", "table_structural_tiers.tsv")
    for row in drows(path):
        if row["tier"] == tier:
            return float(row["n_families"])
    raise KeyError(tier)


def _specificity_intersection(field):
    path = os.path.join(RES, "global_specificity", "cross_catalogue_intersection.tsv")
    data = drows(path)
    if len(data) != 1:
        raise ValueError("cross-catalogue intersection table must have exactly one row")
    row = data[0]
    total = int(row["n_total"])
    hits = int(row["n_with_hit_in_at_least_one"])
    nohits = int(row["n_no_detected_homolog"])
    if hits + nohits != total:
        raise ValueError("cross-catalogue intersection does not close")
    if field == "pct_no_detected_homolog":
        return round(100.0 * nohits / total, 1)
    return float(row[field])


def _hydB(group):
    rs = rows(os.path.join(RES, "atmospheric", "hydB_class.tsv"))
    prot = [r for r in rs if len(r) > 1 and r[1] == f"[NiFe]_Group_{group}"]
    mags = {r[0].split("|")[0] for r in prot}
    return len(prot), len(mags)


def c_hydB_1l_prot():
    return float(_hydB("1l")[0])


def c_hydB_1l_mag():
    return float(_hydB("1l")[1])


def c_hydB_1h_prot():
    return float(_hydB("1h")[0])


def c_hydB_1h_mag():
    return float(_hydB("1h")[1])


def _coxL_formI():
    rs = drows(os.path.join(RES, "atmospheric", "coxL_class.tsv"))
    prot = [r for r in rs if r["form_I_call"] == "form_I"]
    mags = {r["query"].split("|")[0] for r in prot}
    return len(prot), len(mags)


def c_coxL_prot():
    return float(_coxL_formI()[0])


def c_coxL_mag():
    return float(_coxL_formI()[1])


def c_guild_core():
    with open(os.path.join(RES, "atmospheric", "guild_core_mags.txt")) as f:
        return float(sum(1 for line in f if line.strip()))


def _aridity_q1(field):
    path = os.path.join(RES, "atmospheric", "aridity", "aridity_summary.tsv")
    for r in drows(path):
        if r["question"] == "Q1_1l_share_vs_precip":
            return float(r[field])
    raise KeyError("Q1_1l_share_vs_precip")


def c_aridity_rho():
    return round(_aridity_q1("value"), 3)


def c_aridity_p():
    return round(_aridity_q1("p_permutation"), 3)


def c_aridity_n():
    return _aridity_q1("n_sites")


def _uv_group_pct(group):
    for r in drows(os.path.join(RES, "uv_survival", "group_prevalence.tsv")):
        if r["functional_group"] == group:
            return float(r["pct_of_catalog"])
    raise KeyError(group)


def c_ner():
    return _uv_group_pct("nucleotide_excision_repair")


def c_ber():
    return _uv_group_pct("base_excision_repair")


def c_nhej():
    return _uv_group_pct("dsb_nhej")


def _compartment(field, comp):
    for r in drows(os.path.join(RES, "uv_survival", "compartment.tsv")):
        if r["compartment"] == comp:
            return float(r[field])
    raise KeyError(comp)


def c_photo_surface():
    return _compartment("pct_photolyase", "Surface")


def c_photo_deep():
    return _compartment("pct_photolyase", "Deep")


def _bgc_pct(product):
    for r in drows(os.path.join(RES, "bgc_product_summary.tsv")):
        if r["product_type"] == product:
            return round(100.0 * int(r["n_MAGs"]) / catalog_N(), 1)
    raise KeyError(product)


def c_bgc_terpene():
    return _bgc_pct("terpene")


def c_bgc_nrps():
    return _bgc_pct("NRPS")


def c_bgc_ectoine():
    return _bgc_pct("ectoine")


def _metal_mags(metal):
    for r in drows(os.path.join(RES, "amr_metal", "metal_census.tsv")):
        if r["metal_class"] == metal:
            return float(r["n_MAGs"])
    raise KeyError(metal)


def c_arsenic():
    return _metal_mags("ARSENIC")


def c_mercury():
    return _metal_mags("MERCURY")


def c_betalactam():
    for r in drows(os.path.join(RES, "amr_class_summary.tsv")):
        if r["amr_class"] == "BETA-LACTAM":
            return round(100.0 * int(r["n_MAGs"]) / catalog_N(), 1)
    raise KeyError("BETA-LACTAM")


def c_pathogen_pct():
    tot = sum(int(r["n_MAGs"]) for r in drows(os.path.join(RES, "pathogens", "named_pathogen_species.tsv")))
    return round(100.0 * tot / catalog_N(), 2)


def _dark_dnds_values():
    values = []
    for r in drows(os.path.join(RES, "darkmatter_ecology", "dnds.tsv")):
        if r["set"] == "dark" and r["median_dNdS"]:
            values.append(float(r["median_dNdS"]))
    if not values:
        raise ValueError("no evaluable dark-family dN/dS values")
    return values


def c_dnds_dark_median():
    return round(st.median(_dark_dnds_values()), 3)


def _dnds_dark_floor(floor):
    sizes = {
        row["family_id"]: int(row["n_members"])
        for row in drows(
            os.path.join(RES, "darkmatter_ecology", "darkfam_meta.tsv")
        )
    }
    values = [
        float(row["median_dNdS"])
        for row in drows(os.path.join(RES, "darkmatter_ecology", "dnds.tsv"))
        if row["set"] == "dark"
        and row["median_dNdS"]
        and sizes.get(row["family_id"], -1) >= floor
    ]
    if not values:
        raise ValueError(f"no evaluable dark families at size floor {floor}")
    return len(values), st.median(values)


def _abundance_stats():
    ab = os.path.join(RES, "darkmatter_ecology", "darkfam_abundance.tsv")
    npres, peaks = [], []
    with open(ab) as f:
        nsamp = len(f.readline().rstrip("\n").split("\t")) - 1
        if nsamp != 551:
            raise ValueError(f"expected 551 abundance columns, observed {nsamp}")
        for line in f:
            v = [float(x) for x in line.rstrip("\n").split("\t")[1:]]
            if len(v) != nsamp:
                raise ValueError(
                    f"abundance row has {len(v)} values; expected {nsamp}"
                )
            if any(x < 0 for x in v):
                raise ValueError("abundance matrix contains a negative value")
            npres.append(sum(1 for x in v if x > 0))
            peaks.append(max(v) if v else 0.0)
    n = len(npres)
    return dict(
        single=round(100.0 * sum(1 for x in npres if x == 1) / n, 1),
        peak10=round(100.0 * sum(1 for x in peaks if x >= 10) / n, 1),
        medpeak=round(st.median(peaks), 0),
        medprev=round(st.median(npres), 0),
    )


_AB = None


def _ab():
    global _AB
    if _AB is None:
        _AB = _abundance_stats()
    return _AB


def c_b3d_single():
    return _ab()["single"]


def c_b3d_peak10():
    return _ab()["peak10"]


def c_b3d_medpeak():
    return _ab()["medpeak"]


def c_b3d_medprev():
    return _ab()["medprev"]


def _conc():
    return drows(os.path.join(RES, "picrust2_bias", "per_sample_concordance.tsv"))


def c_nsti():
    return round(st.median(float(r["nsti"]) for r in _conc()), 2)


def c_pic_ko():
    values = {float(r["n_pic_ko"]) for r in _conc()}
    if len(values) != 1:
        raise ValueError(f"n_pic_ko is not constant: {len(values)} values")
    return values.pop()


def c_mag_ko():
    return round(st.median(float(r["n_mag_ko"]) for r in _conc()), 0)


def c_recall():
    return round(100.0 * st.median(float(r["recall"]) for r in _conc()), 0)


def c_no_ko():
    d = drows(os.path.join(RES, "picrust2_bias", "dark_fraction.tsv"))
    tot = sum(int(r["total_cds"]) for r in d)
    ko = sum(int(r["n_ko"]) for r in d)
    return round(100.0 * (1 - ko / tot), 1)


def c_rebuild_n():
    path = os.path.join(RES, "picrust2_bias", "rebuild", "per_sample_rebuilt.tsv")
    return len(drows(path))


def _c_richness_ratio(threshold, column):
    path = os.path.join(RES, "picrust2_bias", "rebuild", "threshold_sweep.tsv")
    for r in drows(path):
        if float(r["threshold"]) == threshold:
            return round(float(r[column]), 2)
    raise KeyError(f"threshold {threshold:g}")


def c_presence_ratio_zero():
    return _c_richness_ratio(0.0, "ratio_picrust_over_mag_presence")


def c_presence_ratio_1e5():
    return _c_richness_ratio(1e-5, "ratio_picrust_over_mag_presence")


def c_shared_relative_ratio_zero():
    return _c_richness_ratio(0.0, "ratio_picrust_over_mag_relab")


def c_shared_relative_ratio_1e5():
    return _c_richness_ratio(1e-5, "ratio_picrust_over_mag_relab")


def _signature_rhos():
    path = os.path.join(RES, "picrust2_bias", "marquee_functions.tsv")
    return [
        float(r["spearman_abund"])
        for r in drows(path)
        if float(r["mag_prevalence"]) >= 0.5
    ]


def c_signature_rho_min():
    return round(min(_signature_rhos()), 2)


def c_signature_rho_med():
    return round(st.median(_signature_rhos()), 2)


def c_signature_rho_max():
    return round(max(_signature_rhos()), 2)


def _gem_endemic():
    for r in drows(os.path.join(RES, "global_specificity", "global_specificity_summary.tsv")):
        if r["catalog"] == "gem":
            return float(r["pct_endemic"])
    raise KeyError("gem")


def c_gem():
    return _gem_endemic()


def c_deserts_union():
    for r in drows(os.path.join(RES, "global_specificity", "desert_specificity_summary.tsv")):
        if r["catalog"] == "deserts_combined":
            return float(r["pct_endemic"])
    raise KeyError("deserts_combined")


def c_mgnify():
    for row in drows(os.path.join(RES, "global_specificity", "global_specificity_summary.tsv")):
        if row["catalog"] == "mgnify":
            if int(row["n_with_hit"]) + int(row["n_endemic"]) != int(row["n_total"]):
                raise ValueError("MGnify summary does not close")
            return float(row["pct_endemic"])
    raise KeyError("mgnify")


def verify_artifact_manifest():
    """Fail closed if a compact manuscript input differs from its frozen hash."""
    manifest_path = os.path.join(REPRO, "current_manuscript_artifacts.tsv")
    errors = []
    if not os.path.isfile(manifest_path):
        return 0, ["artifact manifest missing: current_manuscript_artifacts.tsv"]
    manifest_rows = drows(manifest_path)
    for row in manifest_rows:
        path = os.path.join(REPRO, row["artifact"])
        if not os.path.isfile(path):
            errors.append(f"artifact missing: {row['artifact']}")
            continue
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        observed = digest.hexdigest()
        if observed != row["sha256"]:
            errors.append(
                f"artifact hash mismatch: {row['artifact']} "
                f"expected={row['sha256']} observed={observed}"
            )
    return len(manifest_rows), errors


def verify_dnds_legacy_provenance():
    """Enforce every locally deposited or Git-addressed legacy dN/dS hash."""
    path = os.path.join(REPRO, "dnds_legacy_provenance.tsv")
    if not os.path.isfile(path):
        return ["dN/dS legacy provenance ledger missing"]
    errors = []
    for row in drows(path):
        source = row["path_or_git_object"]
        try:
            if re.match(r"^[0-9a-f]{40}:.+", source):
                payload = subprocess.check_output(["git", "show", source], cwd=ROOT)
            else:
                with open(os.path.join(ROOT, source), "rb") as handle:
                    payload = handle.read()
        except (OSError, subprocess.CalledProcessError) as error:
            errors.append(f"dN/dS provenance source unavailable: {source}: {error}")
            continue
        observed = hashlib.sha256(payload).hexdigest()
        if observed != row["sha256"]:
            errors.append(
                f"dN/dS provenance hash mismatch: {source} "
                f"expected={row['sha256']} observed={observed}"
            )
    return errors


def verify_manuscript_boundaries():
    """Guard claim-bearing prose against superseded or unbounded language."""
    manuscript_path = os.path.join(ROOT, "manuscript", "main.tex")
    paths = [
        manuscript_path,
        os.path.join(ROOT, "paper", "main.tex"),
        os.path.join(ROOT, "README.md"),
        os.path.join(ROOT, "handover.md"),
        os.path.join(ROOT, "RQ.md"),
        os.path.join(ROOT, "CLAUDE.md"),
        os.path.join(ROOT, "reproducibility", "README.md"),
        os.path.join(ROOT, "manuscript", "README.md"),
        os.path.join(ROOT, "code", "phase_d_global.py"),
        os.path.join(ROOT, "code", "phase_c_adaptive.py"),
        os.path.join(ROOT, "code", "phase_d_endemism.py"),
        os.path.join(ROOT, "code", "ibex", "96_darkfam_dnds_aggregate.sh"),
        os.path.join(ROOT, "RQ", "C2", "140_pangenome_openness.py"),
        os.path.join(ROOT, "RQ", "C2", "140_pangenome_openness.sbatch"),
        os.path.join(ROOT, "RQ", "C2", "finalize.py"),
        os.path.join(ROOT, "RQ", "F2", "make_fig_bgc_novelty.py"),
    ]
    for tree in (RES, os.path.join(ROOT, "RQ"), os.path.join(ROOT, "docs")):
        for directory, _, filenames in os.walk(tree):
            for filename in filenames:
                if filename.endswith(".md") or filename == "_gen_answered.py":
                    paths.append(os.path.join(directory, filename))
    texts = {}
    for path in paths:
        if os.path.isfile(path):
            with open(path) as handle:
                texts[path] = handle.read()
    text = texts.get(manuscript_path, "")
    forbidden = {
        "SeqCode-eligible": "R232 empty-rank counts are not SeqCode eligibility calls",
        "no homolog anywhere sequenced": "database searches cannot establish universal absence",
        "no homolog anywhere)": "database searches cannot establish universal absence",
        "genuine endemism, not undersampling": "specificity is relative to named comparators",
        "genuine desert-endemic reservoir": "specificity is relative to named comparators",
        "no detectable homolog anywhere": "named catalogue searches cannot establish universal absence",
        "no homolog at all": "named catalogue searches cannot establish universal absence",
        "has not yet been run": "the current known-fold calibration has been run",
        "Domain-unplaced rhizosphere bins": "the active R232 HQ/MQ catalogue has no domain gap",
        "validated floor of 10,007": "the superseded RQ/A3 count is not the current calibrated result",
        "b1a.novelfloor": "the superseded verifier claim no longer exists",
        "68,808 strict novel": "predicted-3Di families are candidates, not strict novel folds",
        "globally endemic": "specificity is relative to named comparators",
        "Genuinely endemic:": "specificity is relative to named comparators",
        "active arms-race selection": "the B4 analysis is an exploratory screen",
        "selection localises to": "the site analysis cannot localise substitutions to interfaces",
        "molecular signature of an arms race": "the B4 analysis is an exploratory screen",
        "No taxonomic novelty (44 phyla, 0 novel)": "the current R232 catalogue has extensive lower-rank novelty",
        "Pervasive purifying selection with a diversifying minority": "pairwise NG86 and the historical site screen do not support this selection claim",
        "per-site diversifying selection": "the historical FUBAR screen is exploratory and not an interface-level conclusion",
        "under purifying selection and skew": "the osmolyte screen is an exploratory association",
        "Endemic genes are MORE constrained": "the unclustered occupancy comparison is descriptive",
        "pangenomes are wide open": "the historical pangenome workflow is not in the active provenance chain",
        "dark proteome lives in the open accessory": "the historical pangenome association is not a mechanism claim",
        "Almost none of the Empty Quarter's secondary-metabolite": "the fragmented BGC analysis is an upper bound",
        "EQ dark matter is not generic desert dark matter": "desert comparisons are relative to named catalogues",
        "EQ/desert-ENDEMIC": "an MGnify90 no-hit is not evidence of geographic endemism",
        "globally SHARED": "an MGnify90 hit is not a universal biogeographic category",
        "which catalog the paper uses": "the active manuscript uses multiple explicitly labelled universes",
        "The paper does not currently use the refined+dereplicated": "the 5,438 representatives are current for taxonomy and recurrent-family analyses",
        "the paper = single_sample 15,469-MAG": "the 15,469-MAG universe is a historical functional snapshot",
        "structured MORE constrained": "the environment-association comparison is an unclustered diagnostic",
        "ADAPTIVE if structured families": "the historical screen does not establish adaptation",
        "under purifying selection => real": "pairwise within-population NG86 does not authenticate every ORF",
        "absent even from MGnify": "a named-database no-hit is not universal absence",
        "EQ/desert-endemic": "an MGnify90 no-hit does not establish geographic endemism",
        "globally-shared environmental": "a named-database hit is not a universal biogeographic class",
        "endemic to a SINGLE site": "single-site is an occupancy description, not endemism",
        "selection: endemic median": "the historical unclustered diagnostic does not establish selection",
        "globally unprecedented": "a BiG-FAM no-hit would still be database-relative",
        "novel family with no characterized": "the fragmented GCF result is only an operational no-MIBiG class",
    }
    required = {
        "14{,}426": "current refined-bin count",
        "5{,}438": "current representative count",
        "GTDB-Tk v2.7.2": "current taxonomy software",
        "PRJEB125153": "retained ENA study accession",
        "237 Thermoplasmatota": "lineage analysis cohort",
        "24.16\\%": "genomic-context ascertainment bound",
        "4{,}210{,}592": "equal-depth reference count",
        "18.54": "primary matched-reference effect",
        "714 of these": "paired crop recovery numerator",
        "964 models": "same-run whole-chain-unmatched denominator",
        "Cropping jointly changes": "joint query-extent perturbation bound",
        "exact correction requires rescoring": "alignment-score correction limitation",
        "public release or reviewer access remains": "public data availability boundary",
    }

    errors = []
    for phrase, reason in forbidden.items():
        for path, document in texts.items():
            if phrase in document:
                errors.append(
                    f"forbidden phrase {phrase!r} in {os.path.relpath(path, ROOT)}: {reason}"
                )
    bounded_patterns = {
        r"novel(?:[-\s]+)fold(?:[-\s]+)famil": (
            re.compile(r"candidate|historical|legacy|reference[- ]relative", re.I),
            "novel-fold-family wording must be explicitly candidate-, historical-, "
            "legacy- or reference-relative",
        ),
    }
    for pattern, (allowance, reason) in bounded_patterns.items():
        matcher = re.compile(pattern, re.I)
        for path, document in texts.items():
            for match in matcher.finditer(document):
                context = document[max(0, match.start() - 120):match.end() + 120]
                if not allowance.search(context):
                    lineno = document.count("\n", 0, match.start()) + 1
                    errors.append(
                        f"unbounded pattern {pattern!r} in "
                        f"{os.path.relpath(path, ROOT)}:{lineno}: {reason}"
                    )
    for phrase, reason in required.items():
        if phrase not in text:
            errors.append(f"required manuscript phrase missing {phrase!r}: {reason}")
    return errors


def verify_provenance_ledger():
    """Pin every offline claim to its exact recorded value and upstream source."""
    path = os.path.join(REPRO, "historical_provenance_ledger.tsv")
    if not os.path.isfile(path):
        return ["historical provenance ledger missing"]
    ledger_rows = drows(path)
    ledger = {row["claim_id"]: row for row in ledger_rows}
    errors = []
    if len(ledger) != len(ledger_rows):
        errors.append("historical provenance ledger contains duplicate claim IDs")
    registered = {
        cid: (str(paper), source)
        for cid, _, _, paper, tier, source, _, _ in CLAIMS
        if tier == P
    }
    if set(ledger) != set(registered):
        errors.append(
            "historical provenance ledger IDs differ from registered PROVENANCE claims: "
            f"missing={sorted(set(registered) - set(ledger))} "
            f"extra={sorted(set(ledger) - set(registered))}"
        )
    for cid in sorted(set(ledger) & set(registered)):
        paper, source = registered[cid]
        row = ledger[cid]
        if row["value"] != paper:
            errors.append(
                f"provenance ledger value mismatch for {cid}: ledger={row['value']!r} registry={paper!r}"
            )
        if row["upstream_source"] != source:
            errors.append(
                f"provenance ledger source mismatch for {cid}: "
                f"ledger={row['upstream_source']!r} registry={source!r}"
            )
        if row["local_scope"] != "transcription_and_provenance_only":
            errors.append(f"unexpected provenance scope for {cid}: {row['local_scope']!r}")
    return errors


# ---------------------------------------------------------------------------
# Registry. (id, section, description, paper_value, tier, source, fn, tol)
#   COMPUTE: fn recomputes from a local analysis table.
#   FROZEN: fn extracts/derives from a SHA-256-pinned compact artifact.
#   TABLE: fn is None; only local source presence is checked (not the recorded value).
#   PROVENANCE: fn is None; upstream location is documented but not checked locally.
# ---------------------------------------------------------------------------
C = "COMPUTE"; F = "FROZEN"; T = "TABLE"; P = "PROVENANCE"
CLAIMS = [
    # --- current canonical scale / R232 taxonomy (A0, Table 1) ---
    ("scale.read_pairs", "A0", "paired-read sample universe", 573, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("paired_read_samples"), 0),
    ("scale.metadata", "A0", "samples with canonical metadata", 551, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("metadata_samples"), 0),
    ("scale.binnable", "A0", "binnable samples", 550, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("binnable_samples"), 0),
    ("scale.coassembly_samples", "A0", "samples assigned to coassembly groups", 527, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("coassembly_samples"), 0),
    ("scale.coassembly_groups", "A0", "coassembly groups", 126, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("coassembly_groups"), 0),
    ("scale.refined_bins", "A0", "canonical refined bins", 14426, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("refined_bins"), 0),
    ("scale.refined_hq", "A0", "canonical refined HQ bins", 4225, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("refined_hq"), 0),
    ("scale.refined_mq", "A0", "canonical refined MQ bins", 8209, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("refined_mq"), 0),
    ("scale.refined_lq", "A0", "canonical refined LQ bins", 1992, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("refined_lq"), 0),
    ("scale.refined_hq_mq", "A0", "canonical refined HQ+MQ bins", 12434, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("refined_hq") + c_catalogue("refined_mq"), 0),
    ("scale.drep_inputs", "A0", "members in the frozen historical dRep Cdb table", 12433, F, "reproducibility/drep_membership_reconciliation.tsv", lambda: c_drep_reconciliation("historical_drep_members"), 0),
    ("scale.drep_current_only", "A0", "current HQ/MQ identifiers absent from historical Cdb membership", 1, F, "reproducibility/drep_membership_reconciliation.tsv", lambda: c_drep_reconciliation("current_only_count"), 0),
    ("scale.drep_historical_only", "A0", "historical Cdb identifiers absent from the current HQ/MQ set", 0, F, "reproducibility/drep_membership_reconciliation.tsv", lambda: c_drep_reconciliation("historical_only_count"), 0),
    ("scale.salvage", "A0a", "coassembly LQ salvage bins retained as a separate prospective cohort", 304, F, "reproducibility/canonical_catalogue_summary.json", lambda: c_catalogue("coassembly_salvage_bins"), 0),
    ("scale.representatives", "A0", "ANI95 HQ/MQ representatives classified with GTDB R232", 5438, F, "reproducibility/r232_taxonomic_novelty_summary.json", lambda: c_r232("legacy_representatives"), 0),
    ("scale.rep_hq", "A0", "HQ representatives", 2259, F, "reproducibility/r232_taxonomic_novelty_summary.json", lambda: c_r232("legacy_hq"), 0),
    ("scale.rep_mq", "A0", "MQ representatives", 3179, F, "reproducibility/r232_taxonomic_novelty_summary.json", lambda: c_r232("legacy_mq"), 0),
    ("scale.named_phyla", "A0", "named phyla represented in the ANI95 catalogue", 51, F, "reproducibility/r232_named_taxon_census.tsv", lambda: c_r232_named("phylum"), 0),
    ("scale.gap_domain", "A0", "first assignment gaps at domain", 0, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("domain"), 0),
    ("scale.gap_phylum", "A0", "first assignment gaps at phylum", 0, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("phylum"), 0),
    ("scale.gap_class", "A0", "first assignment gaps at class", 0, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("class"), 0),
    ("scale.gap_order_all", "A0", "first assignment gaps at order", 11, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("order"), 0),
    ("scale.gap_order_hq", "A0", "HQ first assignment gaps at order", 6, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("order", "hq"), 0),
    ("scale.gap_family", "A0", "first assignment gaps at family", 142, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("family"), 0),
    ("scale.gap_family_hq", "A0", "HQ first assignment gaps at family", 64, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("family", "hq"), 0),
    ("scale.gap_genus", "A0", "first assignment gaps at genus", 1343, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("genus"), 0),
    ("scale.gap_genus_hq", "A0", "HQ first assignment gaps at genus", 530, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("genus", "hq"), 0),
    ("scale.gap_species", "A0", "first assignment gaps at species", 3755, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("species"), 0),
    ("scale.gap_species_hq", "A0", "HQ first assignment gaps at species", 1555, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("species", "hq"), 0),
    ("scale.fully_assigned", "A0", "representatives assigned through species", 187, F, "reproducibility/r232_first_unassigned_census.tsv", lambda: c_r232_gap("fully_assigned"), 0),
    ("scale.r226_changed", "A0", "R232 classifications changed from R226", 1898, F, "reproducibility/r232_taxonomic_novelty_summary.json", lambda: c_r232("legacy_classifications_changed_from_r226"), 0),
    ("scale.phylogeny_candidates", "A0a", "R232 empty-order candidates staged into de novo phylogenies", 10, F, "reproducibility/r232_phylogeny_adjudication_summary.json", lambda: c_phylogeny("candidate_mags"), 0),
    ("scale.phylogeny_present", "A0a", "staged candidates present in de novo trees", 10, F, "reproducibility/r232_phylogeny_adjudication_summary.json", lambda: c_phylogeny("tree_present"), 0),
    ("scale.order_assignments", "A0a", "automatic novel-rank assignments after phylogeny/RED", 0, F, "reproducibility/r232_phylogeny_adjudication_summary.json", lambda: c_phylogeny("automatic_novel_rank_assignments"), 0),
    ("scale.manual_proposals", "A0a", "candidates advanced to manual taxonomic proposal", 0, F, "reproducibility/r232_phylogeny_adjudication_summary.json", lambda: c_phylogeny("priority_manual_taxonomic_review"), 0),
    # Historical functional snapshot, deliberately not treated as the current taxonomy catalogue.
    ("scale.functional_mags", "B1", "MAGs in the historical R226-era functional/Bakta snapshot", 15469, C, "results/picrust2_bias/dark_fraction.tsv", c_mags, 0),
    ("scale.proteins", "B1", "predicted proteins in the historical Bakta snapshot (job 46528255)", "51920778", P, "workbench/CHANGELOG.md:2039", None, None),
    # --- A1 atmospheric (Table 2) ---
    ("a1.hydB1l.prot", "A1", "hydB group 1l proteins", 5292, C, "results/atmospheric/hydB_class.tsv", c_hydB_1l_prot, 0),
    ("a1.hydB1l.mag", "A1", "hydB group 1l MAGs", 2609, C, "results/atmospheric/hydB_class.tsv", c_hydB_1l_mag, 0),
    ("a1.hydB1h.prot", "A1", "hydB group 1h proteins", 986, C, "results/atmospheric/hydB_class.tsv", c_hydB_1h_prot, 0),
    ("a1.hydB1h.mag", "A1", "hydB group 1h MAGs", 893, C, "results/atmospheric/hydB_class.tsv", c_hydB_1h_mag, 0),
    ("a1.coxL.prot", "A1", "coxL form-I proteins", 5645, C, "results/atmospheric/coxL_class.tsv", c_coxL_prot, 0),
    ("a1.coxL.mag", "A1", "coxL form-I MAGs", 3828, C, "results/atmospheric/coxL_class.tsv", c_coxL_mag, 0),
    ("a1.guildcore", "A1", "guild core (H2+CO) MAGs", 1204, C, "results/atmospheric/guild_core_mags.txt", c_guild_core, 0),
    ("a1.aridity.rho", "A1", "group-1l share vs precipitation Spearman rho", -0.172, C, "results/atmospheric/aridity/aridity_summary.tsv", c_aridity_rho, 0.001),
    ("a1.aridity.p", "A1", "group-1l share vs precipitation permutation p", 0.196, C, "results/atmospheric/aridity/aridity_summary.tsv", c_aridity_p, 0.001),
    ("a1.aridity.n", "A1", "group-1l share vs precipitation sites", 58, C, "results/atmospheric/aridity/aridity_summary.tsv", c_aridity_n, 0),
    ("a1.aridity.cross", "A1", "cross-biome group-1l share has no monotonic aridity gradient", "NEGATIVE", T, "results/atmospheric/aridity/cross_biome_hydrogenase.tsv", None, None),
    # --- A2 UV/desiccation ---
    ("a2.ner", "A2", "NER prevalence %", 99.4, C, "results/uv_survival/group_prevalence.tsv", c_ner, 0.1),
    ("a2.ber", "A2", "BER prevalence %", 99.5, C, "results/uv_survival/group_prevalence.tsv", c_ber, 0.1),
    ("a2.nhej", "A2", "NHEJ prevalence %", 78.4, C, "results/uv_survival/group_prevalence.tsv", c_nhej, 0.1),
    ("a2.photo.surf", "A2", "photolyase surface %", 22.4, C, "results/uv_survival/compartment.tsv", c_photo_surface, 0.1),
    ("a2.photo.deep", "A2", "photolyase deep %", 13.1, C, "results/uv_survival/compartment.tsv", c_photo_deep, 0.1),
    # --- A3 osmolyte (table-backed stat) ---
    ("a3.osmolyte.or", "A3", "osmolyte class enrichment OR/p (null)", "OR 1.10, p=0.45", T, "results/darkmatter_ecology/OSMOLYTE.md", None, None),
    # --- A4 BGC ---
    ("a4.terpene", "A4", "terpene BGC % of MAGs", 64.5, C, "results/bgc_product_summary.tsv", c_bgc_terpene, 0.05),
    ("a4.nrps", "A4", "NRPS BGC % of MAGs", 25.3, C, "results/bgc_product_summary.tsv", c_bgc_nrps, 0.3),
    ("a4.ectoine", "A4", "ectoine BGC % of MAGs", 17.0, C, "results/bgc_product_summary.tsv", c_bgc_ectoine, 0.5),
    # --- A5 resistome / pathogens ---
    ("a5.arsenic", "A5", "arsenic-resistance MAGs", 2423, C, "results/amr_metal/metal_census.tsv", c_arsenic, 0),
    ("a5.mercury", "A5", "mercury-resistance MAGs", 1077, C, "results/amr_metal/metal_census.tsv", c_mercury, 0),
    ("a5.betalactam", "A5", "beta-lactam resistome % of MAGs", 24.7, C, "results/amr_class_summary.tsv", c_betalactam, 0.3),
    ("a5.pathogens", "A5", "named pathogen MAGs % (rounds to 0.05)", 0.05, C, "results/pathogens/named_pathogen_species.tsv", c_pathogen_pct, 0.02),
    # --- B1 MAG-Bakta catalogue chain (Table: tab:ladder) ---
    # Every source below is a FILE, not a repository. $M = /ibex/scratch/projects/c2014/rob/metagenomics.
    ("b1.reps", "B1", "linclust 90%id/80%cov representatives of the 51.9M Bakta proteins", "23416070", P, "workbench/CHANGELOG.md:2039", None, None),
    ("b1.uniref.hit", "B1", "UniRef50-hit representatives (of 23416070)", "20577040", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/seq_vs_uniref50/hit_ids.txt", None, None),
    ("b1.dark", "B1", "UniRef50 no-hit (sequence-dark) representatives = 12.1% of 23416070", "2839030 / 12.1%", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/seq_vs_uniref50/nohit_ids.txt", None, None),
    ("b1.structnovel", "B1", "no-AFDB-match representatives in the predicted-3Di screen (of 2839030 sequence-dark)", "2533446 / 89.2%", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/foldseek_results/struct_nohit_ids.txt", None, None),
    ("b1.afdbhit", "B1", "AFDB-hit representatives in the predicted-3Di screen (of 2839030 sequence-dark)", "305584 / 10.8%", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/prostt5_foldseek_orix/merged_afdb_hit_ids.txt", None, None),
    ("b1.clusters30", "B1", "30%-id clusters of ALL 2839030 sequence-dark reps (not of the 2533446 subset)", "1875540", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/novelty_analysis/novel_family_summary.tsv", None, None),
    ("b1.families", "B1", "multi-member sequence-dark families (>=2 members; of 1875540 clusters)", "334810", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/novelty_analysis/novel_family_summary.tsv", None, None),
    ("b1.family.members", "B1", "representatives inside the 334810 multi-member families (of 2839030)", "1298300", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/novelty_analysis/novel_family_summary.tsv", None, None),
    ("b1.singletons", "B1", "sequence-dark singleton clusters (of 1875540 clusters)", "1540730", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/novelty_analysis/novel_family_summary.tsv", None, None),
    # --- B1 assembly catalogue / tier analysis (Table: tab:tiers). DIFFERENT catalogue from the rows above. ---
    ("b1.assembly.catalog", "B1", "assembly protein catalogue (both binning strategies, incl. unbinned)", "753305533", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog/assembly_proteins.mmseqsdb.index", None, None),
    ("b1.assembly.clusters90", "B1", "assembly-catalogue 90%-id clusters (of 753305533 proteins)", "449265306", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog/catalog_ladder/catalog_ladder_summary.tsv", None, None),
    ("b1.assembly.dark", "B1", "UniRef50-dark assembly clusters (of 449265306 clusters)", "68646640 / 15.3%", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog/catalog_ladder/catalog_ladder_summary.tsv", None, None),
    ("b1.tier1.input", "B1", "historical tier1 MAG-associated input from its own pre-correction snapshot", "11418079", P, "workbench/results/tier1_structural_novelty/SUMMARY.md", None, None),
    ("b1.tier1.novel", "B1", "tier1 no-AFDB-match (raw rate; fragment-inflated)", "10534391 / 92.26%", P, "workbench/results/tier1_structural_novelty/SUMMARY.md", None, None),
    ("b1.tier1.fragment", "B1", "tier1 fragment caveat: 2339154 no-hits >=100 aa = 22.2% of raw no-hits; 74.0% no-hit among >=100-aa inputs", "2339154 / 22.2% / 74.0%", P, "workbench/results/tier1_structural_novelty/SUMMARY.md", None, None),
    ("b1.tier2.input", "B1", "historical tier2 non-MAG input from a separate snapshot", "57379961", P, "workbench/CHANGELOG.md:1307", None, None),
    ("b1.tier.nonclosure", "B1", "tier1 + tier2 exceed the assembly-catalogue dark count; snapshots do not close", 151400, F, "reproducibility/historical_tier_snapshot_accounting.tsv", c_tier_snapshot_nonclosure, 0),
    ("b1.tier1.superseded", "B1", "later re-tiering vs DAS_Tool-refined bins: 9515846 tier1 reps demoted to tier2", "11418079 -> 1902233", P, "workbench/CHANGELOG.md:1306", None, None),
    ("b1.tier2.structinput", "B1", "tier2 structural-analysis input (>=100 aa subset of 57379961)", "13088646", P, "workbench/results/tier2_structural_novelty/README.md", None, None),
    ("b1.tier2.novel", "B1", "tier2 (non-MAG, >=100 aa) no-AFDB-match of 13088646", "10113496 / 77.27%", P, "workbench/results/tier2_structural_novelty/README.md", None, None),
    ("b1.tier2.families", "B1", "historical candidate novel-fold families (separate tier2 clustering; not used elsewhere)", "156743", P, "workbench/results/tier2_structural_novelty/headline.txt", None, None),
    ("b1.tier1.clusters", "B1", "tier1 structural clusters / multi-member (of 11418079)", "10898521 / 220177", P, "workbench/results/tier2_structural_novelty/headline_tier1.txt", None, None),
    ("b1.folds", "B1", "candidate novel-fold families (>=2 members, 0 AFDB hit; of 220177 tier1 multi-member)", "68808", P, "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/tier1_cluster/novel_fold_analysis/family_summary.tsv", None, None),
    ("b1.folds.proteins", "B1", "proteins inside the 68808 historical candidate novel-fold families", "169651", P, "workbench/results/tier2_structural_novelty/headline_tier1.txt", None, None),
    ("b1a.structures", "B1a", "explicit ESMFold candidate structures in the current-reference integration", 68527, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_structures"), 0),
    ("b1a.qtm.all", "B1a", "query-normalised dual current-PDB/AFDB no-hits among all candidate models", 41127, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "query_normalized", "unfiltered_dual_nohits"), 0),
    ("b1a.qtm.all.no_rows", "B1a", "unfiltered qTM candidate no-hits with neither reference row", 6097, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "query_normalized", "unfiltered_dual_nohits_neither_reference_reported_a_row"), 0),
    ("b1a.qtm.all.low", "B1a", "unfiltered qTM candidate no-hits with a retrieved row below threshold", 35030, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "query_normalized", "unfiltered_dual_nohits_some_row_reported_but_both_scores_below_threshold"), 0),
    ("b1a.alntm.all", "B1a", "alignment-normalised dual current-PDB/AFDB no-hits among all candidate models", 16386, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "alignment_normalized", "unfiltered_dual_nohits"), 0),
    ("b1a.alntm.all.no_rows", "B1a", "unfiltered alignment-TM candidate no-hits with neither reference row", 6097, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "alignment_normalized", "unfiltered_dual_nohits_neither_reference_reported_a_row"), 0),
    ("b1a.alntm.all.low", "B1a", "unfiltered alignment-TM candidate no-hits with a retrieved row below threshold", 10289, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "alignment_normalized", "unfiltered_dual_nohits_some_row_reported_but_both_scores_below_threshold"), 0),
    ("b1a.candidate.qualified_n", "B1a", "candidate models passing length >=100 aa and mean pLDDT >=0.70", 1875, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "query_normalized", "confidence_length_eligible"), 0),
    ("b1a.qtm.qualified", "B1a", "query-normalised dual no-hits after the confidence/length gate", 608, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "query_normalized", "confidence_length_gated_dual_nohits"), 0),
    ("b1a.qtm.qualified_pct", "B1a", "query-normalised candidate dual-nohit rate %", 32.43, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: round(100.0 * c_structure_metric("candidate_sensitivity", "query_normalized", "confidence_length_gated_dual_nohit_rate"), 2), 0.01),
    ("b1a.qtm.qualified.no_rows", "B1a", "gated qTM candidate no-hits with neither reference row", 74, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "query_normalized", "confidence_length_gated_dual_nohits_neither_reference_reported_a_row"), 0),
    ("b1a.qtm.qualified.low", "B1a", "gated qTM candidate no-hits with a retrieved row below threshold", 534, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("candidate_sensitivity", "query_normalized", "confidence_length_gated_dual_nohits_some_row_reported_but_both_scores_below_threshold"), 0),

    ("b1a.control.qualified_n", "B1a", "confidence- and length-qualified known-fold controls", 3013, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "query_normalized", "confidence_length_qualified", "controls"), 0),
    ("b1a.control.qtm.qualified", "B1a", "query-normalised known-fold false dual no-hits", 101, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "query_normalized", "confidence_length_qualified", "false_dual_nohits"), 0),
    ("b1a.control.qtm.pct", "B1a", "query-normalised known-fold false dual-nohit rate %", 3.35, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: round(100.0 * c_structure_metric("known_fold_sensitivity", "query_normalized", "confidence_length_qualified", "false_dual_nohit_rate"), 2), 0.01),
    ("b1a.control.qtm.no_rows", "B1a", "query-normalised false no-hits with neither reference row", 66, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "query_normalized", "confidence_length_qualified", "false_dual_nohits_neither_reference_reported_a_row"), 0),
    ("b1a.control.qtm.low", "B1a", "query-normalised false no-hits with a row but both qTM scores below 0.5", 35, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "query_normalized", "confidence_length_qualified", "false_dual_nohits_some_row_reported_but_both_scores_below_threshold"), 0),

    ("b1a.control.wilson_lo", "B1a", "query-normalised known-fold false dual-nohit Wilson 95% lower bound %", 2.77, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric_wilson(0), 0.01),
    ("b1a.control.wilson_hi", "B1a", "query-normalised known-fold false dual-nohit Wilson 95% upper bound %", 4.06, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric_wilson(1), 0.01),
    ("b1a.control.all_n", "B1a", "all structure-defined known-fold controls", 3684, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "query_normalized", "all", "controls"), 0),
    ("b1a.control.qtm.all", "B1a", "unfiltered query-normalised known-fold false dual no-hits", 440, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "query_normalized", "all", "false_dual_nohits"), 0),
    ("b1a.control.qtm.all_pct", "B1a", "unfiltered query-normalised known-fold false dual-nohit rate %", 11.94, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: round(100.0 * c_structure_metric("known_fold_sensitivity", "query_normalized", "all", "false_dual_nohit_rate"), 2), 0.01),
    ("b1a.control.qtm.all.no_rows", "B1a", "unfiltered qTM control no-hits with neither reference row", 110, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "query_normalized", "all", "false_dual_nohits_neither_reference_reported_a_row"), 0),
    ("b1a.control.qtm.all.low", "B1a", "unfiltered qTM control no-hits with a retrieved row below threshold", 330, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "query_normalized", "all", "false_dual_nohits_some_row_reported_but_both_scores_below_threshold"), 0),
    ("b1a.control.alntm.all", "B1a", "unfiltered alignment-normalised known-fold false dual no-hits", 236, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "alignment_normalized", "all", "false_dual_nohits"), 0),
    ("b1a.control.alntm.all_pct", "B1a", "unfiltered alignment-normalised known-fold false dual-nohit rate %", 6.41, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: round(100.0 * c_structure_metric("known_fold_sensitivity", "alignment_normalized", "all", "false_dual_nohit_rate"), 2), 0.01),
    ("b1a.control.alntm.all.no_rows", "B1a", "unfiltered alignment-TM control no-hits with neither reference row", 110, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "alignment_normalized", "all", "false_dual_nohits_neither_reference_reported_a_row"), 0),
    ("b1a.control.alntm.all.low", "B1a", "unfiltered alignment-TM control no-hits with a retrieved row below threshold", 126, F, "reproducibility/current_structure_metric_sensitivity_summary.json", lambda: c_structure_metric("known_fold_sensitivity", "alignment_normalized", "all", "false_dual_nohits_some_row_reported_but_both_scores_below_threshold"), 0),
    # --- Object D: the 23,382 dark ORTHOGROUPS. B2/B3/B5/C3 are all measured on THIS object, ---
    # --- which is a 50%-id clustering over the 5,438 dRep MAGs -- NOT the 334,810 families above. ---
    ("b.orthogroups", "B2", "dark orthogroups (50% id over 5438 dRep MAGs, dark rep, >=3 members) = denominator of B2/B3/B5/C3", 23382, C, "results/darkmatter_ecology/darkfam_meta.tsv", c_orthogroups, 0),
    ("b.orthogroups.faa", "B2", "materialised query set for the endemism searches (23382 records)", "23382", P, "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology/endemism/query.faa", None, None),
    # --- B2 selection (on the 23,382 orthogroups) ---
    ("b2.dark.median", "B2", "filtered pairwise NG86 median dN/dS for recurrent dark orthogroups", 0.391, C, "results/darkmatter_ecology/dnds.tsv", c_dnds_dark_median, 0.001),
    ("b2.dark.n", "B2", "dark orthogroups with an evaluable codon alignment (of 23382)", 22767, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: float(_dnds_dark_floor(3)[0]), 0),
    ("b2.dark.floor5.med", "B2", "dark median dN/dS at >=5 total family members", 0.418, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: _dnds_dark_floor(5)[1], 0.001),
    ("b2.dark.floor5.n", "B2", "evaluable dark families at >=5 total members", 8769, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: float(_dnds_dark_floor(5)[0]), 0),
    ("b2.dark.floor10.med", "B2", "dark median dN/dS at >=10 total family members", 0.449, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: _dnds_dark_floor(10)[1], 0.001),
    ("b2.dark.floor10.n", "B2", "evaluable dark families at >=10 total members", 2479, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: float(_dnds_dark_floor(10)[0]), 0),
    ("b2.dark.floor20.med", "B2", "dark median dN/dS at >=20 total family members", 0.489, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: _dnds_dark_floor(20)[1], 0.001),
    ("b2.dark.floor20.n", "B2", "evaluable dark families at >=20 total members", 671, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: float(_dnds_dark_floor(20)[0]), 0),
    ("b2.dark.floor30.med", "B2", "dark median dN/dS at >=30 total family members", 0.499, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: _dnds_dark_floor(30)[1], 0.001),
    ("b2.dark.floor30.n", "B2", "evaluable dark families at >=30 total members", 302, C, "results/darkmatter_ecology/{dnds.tsv,darkfam_meta.tsv}", lambda: float(_dnds_dark_floor(30)[0]), 0),
    ("b2.site_restricted.med", "B2", "descriptive <=3-site vs EQ-wide(>=30) orthogroup pairwise NG86 dN/dS", "0.332 / 0.420", T, "results/darkmatter_ecology/ENDEMISM.md", None, None),
    # --- B3 ecology / abundance (on the 23,382 orthogroups) ---
    ("b3.occupancy", "B3", "orthogroup site occupancy: EQ-wide 10669/23382 vs single-site 543/23382", "45.6 / 2.3", T, "results/darkmatter_ecology/ENDEMISM.md", None, None),
    ("b3d.single", "B3", "dark orthogroups single-sample % (of 23382)", 1.3, C, "results/darkmatter_ecology/darkfam_abundance.tsv", c_b3d_single, 0.2),
    ("b3d.peak10", "B3", "dark orthogroups peak>=10 % (of 23382)", 99.5, C, "results/darkmatter_ecology/darkfam_abundance.tsv", c_b3d_peak10, 0.2),
    ("b3d.medpeak", "B3", "median peak abundance ~585 (23382 x 551 matrix)", 585, C, "results/darkmatter_ecology/darkfam_abundance.tsv", c_b3d_medpeak, 3),
    ("b3d.medprev", "B3", "median orthogroup in N samples (of 551)", 52, C, "results/darkmatter_ecology/darkfam_abundance.tsv", c_b3d_medprev, 1),
    # --- B4 exploratory candidate-tail screen ---
    ("b4.screen_tail", "B4", "historical pairwise-NG86 screen-tail families with mean dN/dS > 1", 148, C, "RQ/E2/table_positive_families_summary.tsv", lambda: _e2_summary("n_dark_families_dNdS_gt1"), 0),
    ("b4.screen_heuristic", "B4", "screen-tail families passing minimum-dS and pair-count heuristic", 123, C, "RQ/E2/table_positive_families_summary.tsv", lambda: _e2_summary("n_robust (mean_dS>0.05 & n_pairs>=3)"), 0),
    ("b4.confident_struct", "B4", "screen-tail families with a confident structural analogue", 8, C, "RQ/E2/table_structural_tiers.tsv", lambda: _e2_tier("confident_homolog (prob>=.9,TM>=.5,pLDDT>=.5)"), 0),
    ("b4.recognition_hyp", "B4", "confident analogues annotated as recognition/binding/surface/transport hypotheses", 5, C, "RQ/E2/table_structural_tiers.tsv", lambda: _e2_tier("  of confident: recognition/binding/surface/transport hypotheses"), 0),
    ("b4.sites", "B4", "nominal FUBAR-positive sites in predicted-pocket families", 332, C, "RQ/E2/table_structural_tiers.tsv", lambda: _e2_tier("nominal FUBAR-positive sites in predicted-pocket families"), 0),
    ("b4.sites_in_pocket", "B4", "nominal FUBAR-positive sites inside predicted pockets", 23, C, "RQ/E2/table_structural_tiers.tsv", lambda: _e2_tier("nominal FUBAR-positive sites inside predicted pockets"), 0),
    # --- B5 decode (on the 23,382 orthogroups) ---
    ("b5.enzymes", "B5", "well-supported specific-EC dark enzymes (from the 23382-orthogroup decode)", "10", T, "results/darkmatter_ecology/family_ec.tsv", None, None),
    # --- C1 reference bias ---
    ("c1.nsti", "C1", "weighted NSTI median", 0.16, C, "results/picrust2_bias/per_sample_concordance.tsv", c_nsti, 0.01),
    ("c1.pic_ko", "C1", "PICRUSt2 KOs/sample under the non-zero rule (constant, not a median)", 10543, C, "results/picrust2_bias/per_sample_concordance.tsv", c_pic_ko, 0),
    ("c1.mag_ko", "C1", "KOs encoded/sample (median)", 4482, C, "results/picrust2_bias/per_sample_concordance.tsv", c_mag_ko, 5),
    ("c1.recall", "C1", "KO presence-recall % (median)", 93, C, "results/picrust2_bias/per_sample_concordance.tsv", c_recall, 1),
    ("c1.no_ko", "C1", "MAG proteins with no KO %", 56.7, C, "results/picrust2_bias/dark_fraction.tsv", c_no_ko, 0.5),
    ("c1.rebuild_n", "C1", "samples in threshold rebuild", 440, C, "results/picrust2_bias/rebuild/per_sample_rebuilt.tsv", c_rebuild_n, 0),
    ("c1.presence_ratio_zero", "C1", "asymmetric PICRUSt2-threshold/MAG-any-presence ratio at zero", 2.40, C, "results/picrust2_bias/rebuild/threshold_sweep.tsv", c_presence_ratio_zero, 0.01),
    ("c1.presence_ratio_1e5", "C1", "asymmetric PICRUSt2-threshold/MAG-any-presence ratio at 1e-5", 1.11, C, "results/picrust2_bias/rebuild/threshold_sweep.tsv", c_presence_ratio_1e5, 0.01),
    ("c1.shared_rel_zero", "C1", "shared relative-abundance threshold ratio at zero", 31.41, C, "results/picrust2_bias/rebuild/threshold_sweep.tsv", c_shared_relative_ratio_zero, 0.01),
    ("c1.shared_rel_1e5", "C1", "shared relative-abundance threshold ratio at 1e-5", 14.63, C, "results/picrust2_bias/rebuild/threshold_sweep.tsv", c_shared_relative_ratio_1e5, 0.01),
    ("c1.signature_rho_min", "C1", "signature-function abundance Spearman minimum", -0.12, C, "results/picrust2_bias/marquee_functions.tsv", c_signature_rho_min, 0.01),
    ("c1.signature_rho_med", "C1", "signature-function abundance Spearman median", 0.09, C, "results/picrust2_bias/marquee_functions.tsv", c_signature_rho_med, 0.01),
    ("c1.signature_rho_max", "C1", "signature-function abundance Spearman maximum", 0.42, C, "results/picrust2_bias/marquee_functions.tsv", c_signature_rho_max, 0.01),
    # --- C2 cross-feeding (table-backed) ---
    ("c2.mro", "C2", "SMETANA global MRO ~0.7; guild p=0.69", "MRO 0.7 / p=0.69", T, "results/smetana/CROSSFEED.md", None, None),
    # --- C3 specificity ---
    ("c3.gem", "C3", "dark orthogroups absent from GEM % (of 23382)", 47.4, C, "results/global_specificity/global_specificity_summary.tsv", c_gem, 0.1),
    ("c3.deserts", "C3", "dark orthogroups absent from desert union % (of 23382)", 31.9, C, "results/global_specificity/desert_specificity_summary.tsv", c_deserts_union, 0.1),
    ("c3.mgnify", "C3", "dark orthogroups absent from MGnify90 % (of 23382)", 36.6, C, "results/global_specificity/global_specificity_summary.tsv", c_mgnify, 0.1),
    ("c3.intersection_n", "C3", "no detected homolog in the named MGnify90/GEM/five-desert comparison", 5011, C, "results/global_specificity/cross_catalogue_intersection.tsv", lambda: _specificity_intersection("n_no_detected_homolog"), 0),
    ("c3.intersection", "C3", "no-detected-homolog fraction in the named MGnify90/GEM/five-desert comparison", 21.4, C, "results/global_specificity/cross_catalogue_intersection.tsv", lambda: _specificity_intersection("pct_no_detected_homolog"), 0.1),
    # =========================================================================
    # H1: Within-Thermoplasmatota PGLS & Reference Proximity (Round 7 Final)
    # =========================================================================
    ("h1.reps", "H1", "within-Thermoplasmatota species representatives evaluated in PGLS", 237, F, "reproducibility/darkmatter/round7_final/h1/v3_closure_verification.json", lambda: v7.h1_closure_val("n_primary"), 0),
    ("h1.pgls.rbh_coef", "H1", "PGLS shared-protein coverage slope (frac_rbh_desert)", -0.1375, F, "reproducibility/darkmatter/round7_final/h1/v3_pgls_model_summary.tsv", lambda: round(v7.h1_pgls_model("frac_rbh_desert", "pgls_coef"), 4), 0.0001),
    ("h1.pgls.rbh_se", "H1", "PGLS shared-protein coverage standard error", 0.0195, F, "reproducibility/darkmatter/round7_final/h1/v3_pgls_model_summary.tsv", lambda: round(v7.h1_pgls_model("frac_rbh_desert", "pgls_se_unb"), 4), 0.0001),
    ("h1.pgls.rbh_p", "H1", "PGLS shared-protein coverage two-sided p-value", 1.88e-11, F, "reproducibility/darkmatter/round7_final/h1/v3_pgls_model_summary.tsv", lambda: float(f"{v7.h1_pgls_model('frac_rbh_desert', 'pgls_pval_unb'):.2e}"), 1e-12),
    ("h1.pgls.aai_coef", "H1", "PGLS amino-acid identity slope (per percentage point)", 0.00170, F, "reproducibility/darkmatter/round7_final/h1/v3_pgls_model_summary.tsv", lambda: round(v7.h1_pgls_model("aai", "pgls_coef"), 5), 0.00001),
    ("h1.pgls.aai_se", "H1", "PGLS amino-acid identity standard error", 0.00041, F, "reproducibility/darkmatter/round7_final/h1/v3_pgls_model_summary.tsv", lambda: round(v7.h1_pgls_model("aai", "pgls_se_unb"), 5), 0.00001),
    ("h1.pgls.aai_p", "H1", "PGLS amino-acid identity two-sided p-value", 4.07e-05, F, "reproducibility/darkmatter/round7_final/h1/v3_pgls_model_summary.tsv", lambda: float(f"{v7.h1_pgls_model('aai', 'pgls_pval_unb'):.2e}"), 1e-6),
    ("h1.pgls.lambda", "H1", "PGLS profile maximum likelihood Pagel's lambda", 0.9158, F, "reproducibility/darkmatter/round7_final/h1/v3_closure_verification.json", lambda: round(v7.h1_closure_val("primary_lambda_ml"), 4), 0.0001),
    ("h1.pgls.lambda_ci_lo", "H1", "Pagel's lambda profile likelihood 95% CI lower bound", 0.8463, F, "reproducibility/darkmatter/round7_final/h1/v3_closure_verification.json", lambda: round(v7.h1_closure_lambda_ci(0), 4), 0.0001),
    ("h1.pgls.lambda_ci_hi", "H1", "Pagel's lambda profile likelihood 95% CI upper bound", 0.9578, F, "reproducibility/darkmatter/round7_final/h1/v3_closure_verification.json", lambda: round(v7.h1_closure_lambda_ci(1), 4), 0.0001),
    ("h1.pgls.r2_gls", "H1", "whitened generalised-least-squares R-squared at lambda_ML", 0.4691, F, "reproducibility/darkmatter/round7_final/h1/v3_closure_verification.json", lambda: round(v7.h1_closure_val("primary_r2_gls"), 4), 0.0001),
    ("h1.meta.distinct_refs", "H1", "distinct reference genomes selected across 237 pairs", 32, F, "reproducibility/darkmatter/round7_final/h1/v3_comparator_metadata_summary.tsv", lambda: v7.h1_comparator_meta("n_unique_nearest_comparators"), 0),
    ("h1.meta.gtdb_selected", "H1", "pairs selecting a GTDB R232 reference", 229, F, "reproducibility/darkmatter/round7_final/h1/v3_comparator_metadata_summary.tsv", lambda: v7.h1_comparator_meta("n_nearest_gtdb_r232_rep"), 0),
    ("h1.meta.mgnify_selected", "H1", "pairs selecting an MGnify reference", 8, F, "reproducibility/darkmatter/round7_final/h1/v3_comparator_metadata_summary.tsv", lambda: v7.h1_comparator_meta("n_nearest_mgnify"), 0),
    ("h1.meta.env_soil", "H1", "pairs with soil environmental metadata", 142, F, "reproducibility/darkmatter/round7_final/h1/v3_comparator_metadata_summary.tsv", lambda: v7.h1_comparator_meta("n_biome_soil"), 0),
    ("h1.meta.env_marine", "H1", "pairs with marine sediment environmental metadata", 3, F, "reproducibility/darkmatter/round7_final/h1/v3_comparator_metadata_summary.tsv", lambda: v7.h1_comparator_meta("n_biome_marine_sediment"), 0),
    ("h1.meta.env_unknown", "H1", "pairs with unknown environmental metadata", 92, F, "reproducibility/darkmatter/round7_final/h1/v3_comparator_metadata_summary.tsv", lambda: v7.h1_comparator_meta("n_biome_unannotated_gtdb_only"), 0),
    ("h1.tree.total_taxa", "H1", "taxa in GTDB-Tk R232 archaeal marker tree", 1914, F, "reproducibility/darkmatter/round7_final/h1/thermo_tree_taxa.tsv", lambda: v7.h1_tree_taxa_count(), 0),
    ("h1.sens.caller_complete_n", "H1", "genomes retained in complete-gene-call proxy specification", 232, F, "reproducibility/darkmatter/round7_final/h1/v3_closure_verification.json", lambda: v7.h1_closure_val("n_caller_complete_sensitivity"), 0),
    ("h1.sens.missing_checkm2_faa", "H1", "genomes lacking archived CheckM2 protein calls", 5, F, "reproducibility/darkmatter/round7_final/h1/v3_closure_verification.json", lambda: float(len(v7._read_json(os.path.join(v7.H1_DIR, "v3_closure_verification.json"))["missing_checkm2_faa_genomes"])), 0),
    ("h1.closure.status", "H1", "closure verification pass for within-lineage PGLS pipeline", 1.0, F, "reproducibility/darkmatter/round7_final/h1/v3_closure_verification.json", lambda: v7.h1_closure_val("status"), 0),
    # =========================================================================
    # H2: Flanking-Orthologue Recurrence & Genomic Context (Round 7 Final)
    # =========================================================================
    ("h2.funnel.qualified", "H2", "qualified dark structural families in attrition universe", 2838, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_attrition("K", "gate0_confidence_qualified") + v7.h2_attrition("QnotK", "gate0_confidence_qualified"), 0),
    ("h2.funnel.binned", "H2", "dark families with a member on a refined-bin contig", 2426, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_attrition("K", "gate1_refined_mag_binned") + v7.h2_attrition("QnotK", "gate1_refined_mag_binned"), 0),
    ("h2.funnel.drep_species", "H2", "dark families with member in recognized 95% ANI species", 2378, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_attrition("K", "gate2_drep_species_ge1") + v7.h2_attrition("QnotK", "gate2_drep_species_ge1"), 0),
    ("h2.funnel.multispecies", "H2", "dark families occurring in at least two species", 1461, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_attrition("K", "gate3_multi_species_ge2") + v7.h2_attrition("QnotK", "gate3_multi_species_ge2"), 0),
    ("h2.funnel.multispecies_k", "H2", "low-query-score dark families in at least two species", 570, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_attrition("K", "gate3_multi_species_ge2"), 0),
    ("h2.funnel.multispecies_qnotk", "H2", "reference-matched dark families in at least two species", 891, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_attrition("QnotK", "gate3_multi_species_ge2"), 0),
    ("h2.admission.tier1_evaluable", "H2", "complete-window evaluable families admitted in Tier 1", 353, F, "reproducibility/darkmatter/round7_final/h2/pilot_family_conservation.tsv", lambda: v7.h2_family_table_stat("tier1_evaluable_total"), 0),
    ("h2.admission.tier1_k", "H2", "complete-window evaluable low-query-score families", 173, F, "reproducibility/darkmatter/round7_final/h2/pilot_family_conservation.tsv", lambda: v7.h2_family_table_stat("tier1_evaluable_k"), 0),
    ("h2.admission.tier1_qnotk", "H2", "complete-window evaluable reference-matched families", 180, F, "reproducibility/darkmatter/round7_final/h2/pilot_family_conservation.tsv", lambda: v7.h2_family_table_stat("tier1_evaluable_qnotk"), 0),
    ("h2.admission.tier1_gate_k_pct", "H2", "Tier 1 complete-window retention % in low-query-score group", 30.35, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_window_retention_pct("K"), 0.01),
    ("h2.admission.tier1_gate_qnotk_pct", "H2", "Tier 1 complete-window retention % in reference-matched group", 20.20, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_window_retention_pct("QnotK"), 0.01),
    ("h2.admission.tier1_multispecies_pct", "H2", "Tier 1 complete-window retention % of multispecies families", 24.16, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_multispecies_gate_retention_pct(), 0.01),
    ("h2.unadj.k_obs_pct", "H2", "low-query-score observed recurrence %", 40.46, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_recurrence_stat("K", "observed_pct"), 0.01),
    ("h2.unadj.qnotk_obs_pct", "H2", "reference-matched observed recurrence %", 32.22, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_recurrence_stat("QnotK", "observed_pct"), 0.01),
    ("h2.unadj.k_exp_pct", "H2", "low-query-score same-contig positional expectation %", 8.71, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_recurrence_stat("K", "expected_pct"), 0.01),
    ("h2.unadj.qnotk_exp_pct", "H2", "reference-matched same-contig positional expectation %", 6.56, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_recurrence_stat("QnotK", "expected_pct"), 0.01),
    ("h2.unadj.k_excess_pct", "H2", "low-query-score background-adjusted excess recurrence pp", 31.75, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_recurrence_stat("K", "excess_pp"), 0.01),
    ("h2.unadj.qnotk_excess_pct", "H2", "reference-matched background-adjusted excess recurrence pp", 25.67, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_recurrence_stat("QnotK", "excess_pp"), 0.01),
    ("h2.contrast.tier1_mh_diff", "H2", "primary covariate-adjusted difference in excess recurrence (MH weight)", 0.0354, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: round(v7.h2_contrast_stat("tier1_primary", "mh_excess_diff"), 4), 0.0001),
    ("h2.contrast.tier1_p_perm", "H2", "primary between-group stratified permutation p-value", 0.407, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: round(v7.h2_contrast_stat("tier1_primary", "p_perm"), 3), 0.001),
    ("h2.strict.pairs", "H2", "strict exact-phylum length-matched pairs", 96, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_strict_matching_stat("pairs"), 0),
    ("h2.strict.recurrence_pct", "H2", "strict matched pairs recurrence % in each group", 31.25, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_strict_matching_stat("recurrence_pct"), 0.01),
    ("h2.strict.excess_diff", "H2", "strict matched pairs difference in excess recurrence", -0.0054, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: round(v7.h2_strict_matching_stat("excess_diff"), 4), 0.0001),
    ("h2.strict.p_perm", "H2", "strict matched pairs sign-flip permutation p-value", 0.923, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: round(v7.h2_strict_matching_stat("p_perm"), 3), 0.001),
    ("h2.stop3.evaluable", "H2", "evaluable families under 3-prime stop match sensitivity", 637, F, "reproducibility/darkmatter/round7_final/h2/pilot_family_conservation.tsv", lambda: v7.h2_family_table_stat("tier2_evaluable_total"), 0),
    ("h2.stop3.k", "H2", "low-query-score families under 3-prime stop match sensitivity", 260, F, "reproducibility/darkmatter/round7_final/h2/pilot_family_conservation.tsv", lambda: v7.h2_family_table_stat("tier2_evaluable_k"), 0),
    ("h2.stop3.qnotk", "H2", "reference-matched families under 3-prime stop match sensitivity", 377, F, "reproducibility/darkmatter/round7_final/h2/pilot_family_conservation.tsv", lambda: v7.h2_family_table_stat("tier2_evaluable_qnotk"), 0),
    ("h2.stop3.mh_diff", "H2", "adjusted excess recurrence difference under 3-prime stop match", 0.0184, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: round(v7.h2_contrast_stat("tier2_stop3_sensitivity", "mh_excess_diff"), 4), 0.0001),
    ("h2.stop3.p_perm", "H2", "stratified permutation p-value under 3-prime stop match", 0.562, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: round(v7.h2_contrast_stat("tier2_stop3_sensitivity", "p_perm"), 3), 0.001),
    ("h2.disjoint.families", "H2", "families retained in species-disjoint sensitivity", 112, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_species_disjoint_stat("families"), 0),
    ("h2.disjoint.k", "H2", "low-query-score families in species-disjoint sensitivity", 51, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_species_disjoint_stat("k_families"), 0),
    ("h2.disjoint.qnotk", "H2", "reference-matched families in species-disjoint sensitivity", 61, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_species_disjoint_stat("qnotk_families"), 0),
    ("h2.disjoint.k_excess", "H2", "species-disjoint mean excess recurrence in low-query-score", 0.2619, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: round(v7.h2_species_disjoint_stat("k_excess"), 4), 0.0001),
    ("h2.disjoint.qnotk_excess", "H2", "species-disjoint mean excess recurrence in reference-matched", 0.1984, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: round(v7.h2_species_disjoint_stat("qnotk_excess"), 4), 0.0001),
    ("h2.partial.evaluable", "H2", "evaluable families retained in partial-window analysis", 1320, F, "reproducibility/darkmatter/round7_final/h2/pilot_family_conservation.tsv", lambda: v7.h2_family_table_stat("tier3_evaluable_total"), 0),
    ("h2.join.bins_verified", "H2", "refined bins verified by SHA-256 and Bakta logs", 2041, F, "reproducibility/darkmatter/round7_final/h2/pilot_join_metrics.json", lambda: v7.h2_bin_verification("verified"), 0),
    ("h2.join.bins_total", "H2", "total refined bins containing candidate loci", 2126, F, "reproducibility/darkmatter/round7_final/h2/pilot_join_metrics.json", lambda: v7.h2_bin_verification("total"), 0),
    ("h2.join.bins_missing", "H2", "bins lacking required Bakta annotation logs", 85, F, "reproducibility/darkmatter/round7_final/h2/pilot_join_metrics.json", lambda: v7.h2_bin_verification("missing"), 0),
    ("h2.loci.total", "H2", "selected loci across Tier 1 evaluable families", 877, F, "reproducibility/darkmatter/round7_final/h2/pilot_family_conservation.tsv", lambda: v7.h2_family_table_stat("tier1_selected_loci"), 0),
    ("h2.loci.bin_contigs", "H2", "unique bin-contig pairs occupied by selected loci", 869, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_locus_independence_stat("unique_bin_contigs"), 0),
    ("h2.loci.bins", "H2", "unique refined bins occupied by selected loci", 617, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_locus_independence_stat("unique_bins"), 0),
    ("h2.loci.species", "H2", "unique 95% ANI species groups occupied by selected loci", 568, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_locus_independence_stat("unique_drep_species"), 0),
    ("h2.loci.multi_contigs", "H2", "contigs contributing to multiple families", 8, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_locus_independence_stat("contigs_with_ge2_families"), 0),
    ("h2.loci.multi_bins", "H2", "bins contributing to multiple families", 145, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_locus_independence_stat("bins_with_ge2_families"), 0),
    ("h2.loci.multi_species", "H2", "species groups contributing to multiple families", 153, F, "reproducibility/darkmatter/round7_final/h2/pilot_conservation_summary.json", lambda: v7.h2_locus_independence_stat("species_with_ge2_families"), 0),
    # =========================================================================
    # H3: Matched-Reference Comparison (SMAG vs GEM) (Round 7 Final)
    # =========================================================================
    ("h3.matched.proteins_per_arm", "H3", "exact matched protein quota per reference arm", 4210592, F, "reproducibility/darkmatter/round7_final/h3/matching_closure_audit.tsv", lambda: v7.h3_closure_val("harmonized_seed20261014", "matched_proteins"), 0),
    ("h3.matched.residues_per_arm", "H3", "exact matched residue quota per reference arm", 1295970055, F, "reproducibility/darkmatter/round7_final/h3/matching_closure_audit.tsv", lambda: v7.h3_closure_val("harmonized_seed20261014", "matched_residues"), 0),
    ("h3.primary.dark_smag_hits", "H3", "dark queries with hit in SMAG matched reference", 30064, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "smag_hits"), 0),
    ("h3.primary.dark_gem_hits", "H3", "dark queries with hit in GEM matched reference", 10399, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "gem_hits"), 0),
    ("h3.primary.dark_smag_pct", "H3", "dark query hit rate % in SMAG matched reference", 28.35, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "smag_hit_pct"), 0.01),
    ("h3.primary.dark_gem_pct", "H3", "dark query hit rate % in GEM matched reference", 9.81, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "gem_hit_pct"), 0.01),
    ("h3.primary.dark_diff_pp", "H3", "paired dark query recovery difference percentage points", 18.54, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "paired_delta_pct"), 0.01),
    ("h3.primary.control_smag_pct", "H3", "control query hit rate % in SMAG matched reference", 93.99, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Control", "smag_hit_pct"), 0.01),
    ("h3.primary.control_gem_pct", "H3", "control query hit rate % in GEM matched reference", 90.89, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Control", "gem_hit_pct"), 0.01),
    ("h3.primary.control_diff_pp", "H3", "paired control query recovery difference percentage points", 3.10, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Control", "paired_delta_pct"), 0.01),
    ("h3.contingency.dark_both", "H3", "dark queries with hits in both SMAG and GEM", 8015, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "both_hits_a"), 0),
    ("h3.contingency.dark_smag_only", "H3", "dark queries with hits in SMAG only", 22049, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "smag_only_b"), 0),
    ("h3.contingency.dark_gem_only", "H3", "dark queries with hits in GEM only", 2384, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "gem_only_c"), 0),
    ("h3.contingency.dark_neither", "H3", "dark queries with hits in neither reference", 73594, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Dark", "neither_d"), 0),
    ("h3.contingency.ctrl_both", "H3", "control queries with hits in both SMAG and GEM", 94277, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Control", "both_hits_a"), 0),
    ("h3.contingency.ctrl_smag_only", "H3", "control queries with hits in SMAG only", 5396, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Control", "smag_only_b"), 0),
    ("h3.contingency.ctrl_gem_only", "H3", "control queries with hits in GEM only", 2108, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Control", "gem_only_c"), 0),
    ("h3.contingency.ctrl_neither", "H3", "control queries with hits in neither reference", 4261, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014", "Control", "neither_d"), 0),
    ("h3.draws.harmonized_draw2_diff", "H3", "harmonized GEM draw 2 dark recovery difference pp", 16.83, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261015", "Dark", "paired_delta_pct"), 0.01),
    ("h3.draws.harmonized_draw3_diff", "H3", "harmonized GEM draw 3 dark recovery difference pp", 16.91, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261016", "Dark", "paired_delta_pct"), 0.01),
    ("h3.draws.published_draw1_diff", "H3", "published GEM calls draw 1 dark recovery difference pp", 18.46, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("published_seed20261014", "Dark", "paired_delta_pct"), 0.01),
    ("h3.draws.published_draw2_diff", "H3", "published GEM calls draw 2 dark recovery difference pp", 16.87, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("published_seed20261015", "Dark", "paired_delta_pct"), 0.01),
    ("h3.draws.published_draw3_diff", "H3", "published GEM calls draw 3 dark recovery difference pp", 16.91, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("published_seed20261016", "Dark", "paired_delta_pct"), 0.01),
    ("h3.overlap.shared_exact_seqs", "H3", "exact shared sequences retained between primary arms", 22310, F, "reproducibility/darkmatter/round7_final/h3/matching_closure_audit.tsv", lambda: v7.h3_closure_val("harmonized_seed20261014", "post_matching_shared_seqs"), 0),
    ("h3.no_overlap.proteins_per_arm", "H3", "proteins per arm after removing shared sequences", 4185525, F, "reproducibility/darkmatter/round7_final/h3/matching_closure_audit.tsv", lambda: v7.h3_closure_val("harmonized_seed20261014_no_overlap", "matched_proteins"), 0),
    ("h3.no_overlap.residues_per_arm", "H3", "residues per arm after removing shared sequences", 1289129312, F, "reproducibility/darkmatter/round7_final/h3/matching_closure_audit.tsv", lambda: v7.h3_closure_val("harmonized_seed20261014_no_overlap", "matched_residues"), 0),
    ("h3.no_overlap.dark_diff_pp", "H3", "dark difference pp with zero shared sequence overlap", 18.62, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014_no_overlap", "Dark", "paired_delta_pct"), 0.01),
    ("h3.id30.dark_diff_pp", "H3", "dark difference pp under >=30% sequence identity filter", 18.14, F, "reproducibility/darkmatter/round7_final/h3/matched_full_census_summary.tsv", lambda: v7.h3_census_val("harmonized_seed20261014_id30", "Dark", "paired_delta_pct"), 0.01),
    ("h3.manifest.zenodo_parts", "H3", "verified archive parts from SMAG Zenodo deposit 8223844", 373, F, "reproducibility/darkmatter/round7_final/h3/smag_zenodo_parts_manifest.tsv", lambda: v7.h3_zenodo_parts_count(), 0),
    # =========================================================================
    # H4: Paired PDB Structural Searches & Retrieval Sensitivity (Round 7 Final)
    # =========================================================================
    ("h4.census.total_models", "H4", "finite population models evaluated in paired PDB census", 5851, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_total_models(), 0),
    ("h4.census.total_crops", "H4", "confidence-defined crops extracted across all models", 7773, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_total_crops(), 0),
    ("h4.census.models_no_region", "H4", "models lacking an eligible confidence-defined region", 161, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_models_no_region(), 0),
    ("h4.pooled_k.total_models", "H4", "pooled initial low-query-score dark models (K_local_pdb + K_local_afdb_only + K_dual_none)", 1031, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_pooled_k("model_count"), 0),
    ("h4.pooled_k.eligible_regions", "H4", "initial low-query-score dark models with an eligible region", 984, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_pooled_k("models_eligible"), 0),
    ("h4.cpu.pooled_k_wc_matched", "H4", "paired CPU PDB whole-chain matches in initial low-query dark group", 67, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_pooled_k("wc_matched"), 0),
    ("h4.cpu.pooled_k_wc_unmatched", "H4", "paired CPU PDB whole-chain absences in initial low-query dark group", 964, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_pooled_k("wc_unmatched"), 0),
    ("h4.cpu.pooled_k_crop_recovered", "H4", "whole-chain-unmatched dark models recovering a crop match", 714, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_pooled_k("crop_recovered"), 0),
    ("h4.cpu.pooled_k_crop_recovered_pct", "H4", "recovery rate % of whole-chain-unmatched dark models", 74.07, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_pooled_k("crop_recovery_pct"), 0.01),
    ("h4.cpu.pooled_k_crop_matched_total", "H4", "initial low-query dark models with at least one matched crop", 775, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_pooled_k("crop_matched_total"), 0),
    ("h4.cpu.ctrl_tier1_crop_matched", "H4", "control crop matches after excluding source PDB entry", 2889, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_control_crops("tier_1"), 0),
    ("h4.cpu.ctrl_tier1_crop_pct", "H4", "control crop match % under source-entry exclusion", 95.88, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_control_crops_pct("tier_1"), 0.01),
    ("h4.cpu.ctrl_tier2_crop_matched", "H4", "control crop matches under 30% sequence identity exclusion", 2842, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_control_crops("tier_2"), 0),
    ("h4.cpu.ctrl_tier3_crop_matched", "H4", "control crop matches under 25% sequence identity exclusion", 2823, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_control_crops("tier_3"), 0),
    ("h4.cpu.ctrl_tier4_crop_matched", "H4", "control crop matches under 20% sequence identity exclusion", 2796, F, "reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv", lambda: v7.h4_census_control_crops("tier_4"), 0),
    ("h4.adaptive.wc_expanded", "H4", "CPU whole chains with hit cap expanded to 5,000", 39, F, "reproducibility/darkmatter/round7_final/h4/v7_recovery_summary.json", lambda: v7.h4_cap_audit("wc_expanded"), 0),
    ("h4.adaptive.crop_expanded", "H4", "CPU crops with hit cap expanded to 5,000", 55, F, "reproducibility/darkmatter/round7_final/h4/v7_recovery_summary.json", lambda: v7.h4_cap_audit("crop_expanded"), 0),
    ("h4.adaptive.wc_saturated_unmatched", "H4", "whole chains remaining cap-saturated and unmatched", 3, F, "reproducibility/darkmatter/round7_final/h4/v7_recovery_summary.json", lambda: v7.h4_cap_audit("wc_saturated_unmatched"), 0),
    ("h4.adaptive.crop_saturated_unmatched", "H4", "crops remaining cap-saturated and unmatched", 1, F, "reproducibility/darkmatter/round7_final/h4/v7_recovery_summary.json", lambda: v7.h4_cap_audit("crop_saturated_unmatched"), 0),
    ("h4.gpu.wc_rows", "H4", "ORIX GPU whole-chain alignment rows returned", 2853102, F, "reproducibility/darkmatter/round7_final/h4/gpu_raw_audit.json", lambda: v7.h4_gpu_raw_stat("wc_rows"), 0),
    ("h4.gpu.crop_rows", "H4", "ORIX GPU crop alignment rows returned", 3228450, F, "reproducibility/darkmatter/round7_final/h4/gpu_raw_audit.json", lambda: v7.h4_gpu_raw_stat("crop_rows"), 0),
    ("h4.gpu.pooled_k_wc_unmatched", "H4", "GPU whole-chain absences in initial low-query dark group", 961, F, "reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv", lambda: v7.h4_gpu_census_pooled_k("gpu_wc_unmatched"), 0),
    ("h4.gpu.pooled_k_crop_recovered", "H4", "GPU whole-chain absences recovering a crop match", 731, F, "reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv", lambda: v7.h4_gpu_census_pooled_k("gpu_recovered"), 0),
    ("h4.gpu.pooled_k_crop_matched_total", "H4", "initial low-query dark models with matched crop on GPU", 794, F, "reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv", lambda: v7.h4_gpu_census_pooled_k("gpu_crop_matched"), 0),
    ("h4.gpu.wc_disagreements", "H4", "CPU vs GPU whole-chain binary classification differences", 41, F, "reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv", lambda: v7.h4_gpu_census_pooled_k("wc_disagreements"), 0),
    ("h4.gpu.crop_disagreements", "H4", "CPU vs GPU crop binary classification differences", 43, F, "reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv", lambda: v7.h4_gpu_census_pooled_k("crop_disagreements"), 0),
    ("h4.gpu.pooled_k_cap_unmatched", "H4", "dark queries cap-saturated while remaining unmatched on GPU", 0, F, "reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv", lambda: v7.h4_gpu_census_pooled_k("gpu_cap_unmatched"), 0),
    ("h4.gpu.ctrl_crop_matched", "H4", "control crop matches on GPU under source exclusion", 2896, F, "reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv", lambda: v7.h4_gpu_control_stat("crop_matched"), 0),
    ("h4.gpu.ctrl_wc_cap_unmatched", "H4", "control whole chains cap-saturated and unmatched on GPU", 2, F, "reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv", lambda: v7.h4_gpu_control_stat("wc_cap_unmatched"), 0),
    ("h4.anomaly.cpu_rows_audited", "H4", "CPU raw alignment rows audited across four hit files", 6845932, F, "reproducibility/darkmatter/round7_final/h4/tm_anomaly_census.json", lambda: v7.h4_cpu_anomaly_stat("rows_audited"), 0),
    ("h4.anomaly.cpu_alntm_gt1", "H4", "CPU raw alignment rows with alnTM > 1 (all controls)", 2194, F, "reproducibility/darkmatter/round7_final/h4/tm_anomaly_census.json", lambda: v7.h4_cpu_anomaly_stat("alntm_gt1"), 0),
    ("h4.anomaly.gpu_alntm_gt1", "H4", "GPU raw alignment rows with alnTM > 1 (all controls)", 1708, F, "reproducibility/darkmatter/round7_final/h4/gpu_raw_audit.json", lambda: v7.h4_gpu_raw_stat("alntm_gt1"), 0),
    ("h4.anomaly.cpu_exclusion_changed_models", "H4", "distinct control models with crop status changed by alnTM>1 exclusion", 9, F, "reproducibility/darkmatter/round7_final/h4/h4_exclusion_fact_audit.json", lambda: v7.h4_model_crop_status_changes(), 0),
    # Root additions after post-integration review round 1.
    ('h1.sens.same_class_n', 'H1', 'genomes with same-class comparator', 236, F, 'reproducibility/darkmatter/round7_final/h1/v3_pgls_sensitivity.tsv', lambda: v7.h1_same_class_n(), 0),
    ('h1.tree.gtdb_thermoplasmatota', 'H1', 'tree tips in gtdb_thermoplasmatota', 1647, F, 'reproducibility/darkmatter/round7_final/h1/thermo_tree_taxa.tsv', lambda: v7.h1_tree_taxa_count('gtdb_thermoplasmatota'), 0),
    ('h1.tree.outgroup', 'H1', 'tree tips in outgroup', 30, F, 'reproducibility/darkmatter/round7_final/h1/thermo_tree_taxa.tsv', lambda: v7.h1_tree_taxa_count('outgroup'), 0),
    ('h3.overlap.pre_matching', 'H3', 'shared exact sequences before quota matching', 25067, F, 'reproducibility/darkmatter/round7_final/h3/matching_closure_audit.tsv', lambda: float(next(r['pre_matching_shared_seqs'] for r in v7._read_tsv(os.path.join(v7.H3_DIR,'matching_closure_audit.tsv')) if r['label']=='harmonized_seed20261014')), 0),
    ('h4.cpu.other_wc', 'H4', 'other wc', 1408, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Q_not_K_all','wc_matched'), 0),
    ('h4.cpu.other_crop', 'H4', 'other crop', 1578, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Q_not_K_all','crop_matched'), 0),
    ('h4.cpu.initial_pdb_n', 'H4', 'initial pdb n', 1107, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Q_not_K_pdb','model_count'), 0),
    ('h4.cpu.initial_pdb_wc', 'H4', 'initial pdb wc', 1107, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Q_not_K_pdb','wc_matched'), 0),
    ('h4.cpu.initial_pdb_crop', 'H4', 'initial pdb crop', 1009, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Q_not_K_pdb','crop_matched'), 0),
    ('h4.cpu.initial_afdb_n', 'H4', 'initial afdb n', 700, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Q_not_K_afdb_only','model_count'), 0),
    ('h4.cpu.initial_afdb_wc', 'H4', 'initial afdb wc', 301, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Q_not_K_afdb_only','wc_matched'), 0),
    ('h4.cpu.initial_afdb_crop', 'H4', 'initial afdb crop', 569, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Q_not_K_afdb_only','crop_matched'), 0),
    ('h4.cpu.control_wc', 'H4', 'control wc', 2878, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_stratum_value('Controls','wc_matched'), 0),
    ('h4.k.trans_persistent_low', 'H4', 'pooled K trans_persistent_low', 167, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_k_transition('trans_persistent_low'), 0),
    ('h4.k.trans_raw_zero', 'H4', 'pooled K trans_raw_zero', 38, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_k_transition('trans_raw_zero'), 0),
    ('h4.k.trans_post_empty', 'H4', 'pooled K trans_post_empty', 0, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_k_transition('trans_post_empty'), 0),
    ('h4.k.trans_unm_no_reg', 'H4', 'pooled K trans_unm_no_reg', 45, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_k_transition('trans_unm_no_reg'), 0),
    ('h4.k.trans_retained', 'H4', 'pooled K trans_retained', 61, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_k_transition('trans_retained'), 0),
    ('h4.k.trans_lost_low', 'H4', 'pooled K trans_lost_low', 4, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_k_transition('trans_lost_low'), 0),
    ('h4.k.trans_m_no_reg', 'H4', 'pooled K trans_m_no_reg', 2, F, 'reproducibility/darkmatter/round7_final/h4/v7_stratum_census.tsv', lambda: v7.h4_k_transition('trans_m_no_reg'), 0),
    ('h4.gpu.control_tier2', 'H4', 'GPU finite-depth control crop matches', 2816, F, 'reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv', lambda: v7.h4_gpu_control_by_tier(2), 0),
    ('h4.gpu.control_tier3', 'H4', 'GPU finite-depth control crop matches', 2788, F, 'reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv', lambda: v7.h4_gpu_control_by_tier(3), 0),
    ('h4.gpu.control_tier4', 'H4', 'GPU finite-depth control crop matches', 2753, F, 'reproducibility/darkmatter/round7_final/h4/gpu_cpu_census.tsv', lambda: v7.h4_gpu_control_by_tier(4), 0),
    ('h4.anomaly.wc_changed', 'H4', 'unique control whole-chain classifications changed under strict tiers', 4, F, 'reproducibility/darkmatter/round7_final/h4/h4_exclusion_fact_audit.json', lambda: float(v7._read_json(os.path.join(v7.H4_DIR,'h4_exclusion_fact_audit.json'))['wc_status_changes']['unique_queries']), 0),
]


CLAIMS.extend(verify_family.claims())


def fmt(v):
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def manifest():
    """Dump the human-readable claims manifest as TSV (reproducibility/claims.tsv)."""
    out = ["claim_id\tsection\tdescription\tpaper_value\ttier\tsource"]
    for cid, sec, desc, paper, tier, source, fn, tol in CLAIMS:
        out.append(f"{cid}\t{sec}\t{desc}\t{fmt(paper)}\t{tier}\t{source}")
    print("\n".join(out))
    return 0


def main():
    if "--manifest" in sys.argv:
        return manifest()
    quiet = "--quiet" in sys.argv
    artifact_count, artifact_errors = verify_artifact_manifest()
    manuscript_errors = verify_manuscript_boundaries()
    provenance_errors = verify_provenance_ledger()
    dnds_errors = verify_dnds_legacy_provenance()
    gate_errors = artifact_errors + manuscript_errors + provenance_errors + dnds_errors
    npass_compute = npass_frozen = ntable = nprovenance = 0
    nfail = len(gate_errors)
    fails = [("gate", error) for error in gate_errors]
    if gate_errors:
        print(f"GATES: FAIL ({len(gate_errors)} issue(s))")
        for error in gate_errors:
            print(f"  - {error}")
    else:
        print(
            f"GATES: PASS ({artifact_count} frozen artifacts; "
            "manuscript boundary language current)"
        )
    if not quiet:
        print(f"{'STATUS':8} {'ID':22} {'PAPER':>14}  {'COMPUTED':>12}  DESCRIPTION")
        print("-" * 100)
    for cid, sec, desc, paper, tier, source, fn, tol in CLAIMS:
        if tier in (C, F):
            try:
                got = fn()
            except Exception as e:
                status, comp = "ERROR", str(e)[:30]
                fails.append((cid, comp)); nfail += 1
                if not quiet:
                    print(f"{status:8} {cid:22} {fmt(paper):>14}  {comp:>12}  {desc}")
                continue
            ok = abs(got - float(paper)) <= tol
            status = "PASS" if ok else "FAIL"
            if ok:
                if tier == C:
                    npass_compute += 1
                else:
                    npass_frozen += 1
            else:
                nfail += 1; fails.append((cid, f"paper={paper} got={fmt(got)}"))
            if not quiet:
                print(f"{status:8} {cid:22} {fmt(paper):>14}  {fmt(got):>12}  {desc}")
        elif tier == T:
            sp = os.path.join(ROOT, source.split(" ")[0]) if "/" in source else None
            exists = sp is not None and os.path.exists(sp)
            status = "SOURCE" if exists else "MISS"
            if status == "MISS":
                nfail += 1; fails.append((cid, f"source missing: {source}"))
            else:
                ntable += 1
            tag = f"[{tier}]"
            if not quiet:
                print(f"{status:8} {cid:22} {fmt(paper):>14}  {tag:>12}  {desc}  <- {source}")
        elif tier == P:
            nprovenance += 1
            if not quiet:
                print(f"{'LEDGER':8} {cid:22} {fmt(paper):>14}  {'[PROVENANCE]':>12}  {desc}  <- {source}")
        else:
            nfail += 1
            fails.append((cid, f"unknown tier: {tier}"))
    if not quiet:
        print("-" * 100)
    print(
        f"COMPUTE: {npass_compute} pass | FROZEN: {npass_frozen} pass | "
        f"TABLE sources: {ntable} present | PROVENANCE: {nprovenance} ledger-pinned | "
        f"FAIL: {nfail}"
    )
    if fails:
        print("\nFAILURES / MISSING:")
        for cid, why in fails:
            print(f"  - {cid}: {why}")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
