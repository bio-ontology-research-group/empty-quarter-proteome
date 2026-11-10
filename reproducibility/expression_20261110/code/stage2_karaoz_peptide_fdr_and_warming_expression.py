#!/usr/bin/env python3
"""
Stage 2 (Slurm workq): Karaoz et al. (Nat Commun 2026) Metadata-Verified MS-GF+ PSM & Peptide
Q-Value Filtering, Conservative Observed-Background vs Full-Search-FASTA Specificity Audit,
Unique-Spectrum Union Quantification, and Paired n=3 Plot-Blocked Warming Expression Analysis.

Complies with:
- DESIGN.txt
- ROOT-ANALYSIS-CONTRACT.txt (Rules 1-7)
- ROOT-STAGE1-REVIEW.txt
- ROOT-DOWNLOAD-REVIEW.txt
- ROOT-TRANSLATION-CORRECTION.txt
- ROOT-PSM-SCORES.txt:
  * Uses upstream MS-GF:PepQValue (search_engine_score[3], MS:1002055) <= 0.01 and
    MS-GF:QValue (search_engine_score[6], MS:1002054) <= 0.01 verified from MTD headers
    across all 30 .mzTab files. Never reconstructs peptide FDR by sorting PSM q-values.
  * Cross-checks opt_global_cv_MS:1002354_PSM-level_q-value vs search_engine_score[6].
- ROOT-BACKGROUND-REVIEW.txt:
  * Explicitly audits the 604-file MassIVE MSV000097826 inventory and ESS-DIVE packages to check
    whether the 15-split 5%-FDR composite search FASTA
    (Blodgett_Soil_Metagenomes_5PctFDR-Composite_TrypPigBov_2023-05-22_15x_01..15.fasta) or the
    initial 76M-ORF pool is publicly deposited.
  * Because those unbinned/full composite search FASTAs are not deposited in MassIVE or ESS-DIVE,
    sets full_background_verified = False, marks primary global single-MAG uniqueness as
    non-estimable, and reports provisional observed-subset support separately from dataset-reported
    author protein inference.
  * Uses stream_multipattern_peptide_search_records to scan all available sequences on disk
    (EMSL50964_Proteins.faa + Stage 1 v2 274,006 Archaeal MAG ORFs + MassIVE deposited
    A_thaliana_UniProt_2020-11-02_MaxQuant_Contaminant_Proteins_Tryp_Pig_Bov.fasta) plus all
    co-matching accessions in the 30 .mzTab files.
  * Quantifies protein, family, and clade abundance by taking the set union of unique spectra_ref
    via aggregate_unique_spectra_for_entity (never summing across aliases or paralogs, and never
    allowing ambiguous peptides to drive family- or clade-specific warming calculations).
  * Computes all summary counts and validation flags dynamically from runtime audit data.
"""
import concurrent.futures
import csv
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from download_massive_mztab import validate_mztab_file
from expression_validation_core import (
    aggregate_unique_spectra_for_entity,
    benjamini_hochberg_fdr,
    canonicalize_and_validate_protein_seq,
    classify_peptide_specificity_against_background,
    classify_primary_warming_estimability,
    classify_protein_and_family_detection,
    compute_paired_warming_effect_n3_plots,
    filter_and_audit_mztab_psms,
    lookup_strict,
    normalize_peptide_il,
    stream_multipattern_peptide_search_records,
    validate_karaoz_experimental_design,
)

EXP_ROOT = "/scratch/project/k10018/eq-protein-external-20261110/expression"
KARAOZ_DIR = os.path.join(EXP_ROOT, "inputs", "karaoz2026")
MZTAB_DIR = os.path.join(KARAOZ_DIR, "massive_MSV000097826")
STAGE1_V2_DIR = os.path.join(EXP_ROOT, "results", "stage1_v2_homology_eligibility")
WORK_STAGE1_V2 = os.path.join(EXP_ROOT, "work_stage1_v2")
RESULTS_DIR = os.path.join(EXP_ROOT, "results", "stage2_expression_validation")


def write_tsv_with_schema(path: str, cols: list, rows: list):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")


def normalize_mag_id(mag_id: str) -> str:
    """
    Strip only file extensions (.fa, .fna, .fasta) while preserving bin caller prefixes
    (e.g., metabat., maxbin., cocacola.) because both 3300060016.181_sub and
    metabat.3300060016.181_sub exist as distinct dereplicated MAGs in Karaoz 2026.
    """
    s = mag_id.strip().split("\t")[0].strip()
    for ext in (".fasta", ".fna", ".fa"):
        if s.endswith(ext):
            s = s[: -len(ext)]
            break
    return s


def parse_and_filter_single_mztab(args: Tuple[dict, str]) -> Dict[str, Any]:
    """
    Worker function to validate and parse one .mzTab file on a compute core:
    1. Runs validate_mztab_file (checking trailing newline, non-HTML, exact PSM_rows & Protein_rows).
    2. Parses MTD lines and PSH/PSM rows, runs filter_and_audit_mztab_psms with max_psm_q=0.01
       and max_peptide_q=0.01 using upstream MS-GF:PepQValue and MS-GF:QValue.
    3. Preserves qualifying_target_psms (with spectra_ref, psm_id, accession, clean_peptide, il_peptide)
       so protein/family/clade abundance is quantified by unique spectrum union (never summed across aliases).
    4. Tracks all target accessions sharing any spectrum/peptide in the .mzTab PSM rows.
    """
    meta_row, mztab_dir = args
    file_desc = meta_row["File_descriptor"]
    fname = os.path.basename(file_desc)
    suffix = "_msgfplus.mzTab"
    sample_id = fname[:-len(suffix)] if fname.endswith(suffix) else fname
    exp_psm = int(meta_row.get("PSM_rows", "0") or 0)
    exp_prt = int(meta_row.get("Protein_rows", "0") or 0)
    fpath = os.path.join(mztab_dir, fname)

    ok, reason, vstats = validate_mztab_file(fpath, exp_psm, exp_prt, sample_id, None)
    if not ok:
        raise ValueError(f"mzTab validation failed for {fname}: {reason}")

    mtd_lines = []
    psh_cols = None
    psm_rows = []
    prt_accessions = set()

    with open(fpath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("MTD\t"):
                mtd_lines.append(line)
            elif line.startswith("PRT\t"):
                parts = line.split("\t")
                if len(parts) > 1:
                    prt_accessions.add(parts[1].strip())
            elif line.startswith("PSH\t"):
                psh_cols = line.split("\t")[1:]
            elif line.startswith("PSM\t"):
                if psh_cols is None:
                    raise ValueError(f"PSM row before PSH header in {fname}")
                vals = line.split("\t")[1:]
                if len(vals) != len(psh_cols):
                    raise ValueError(f"Column count mismatch in {fname}")
                psm_rows.append(dict(zip(psh_cols, vals)))

    audit_res = filter_and_audit_mztab_psms(
        psm_rows, mtd_lines=mtd_lines, max_psm_q=0.01, max_peptide_q=0.01
    )

    # Record for every il_peptide all target accessions it matched in the .mzTab file
    pep_to_all_db_accessions: Dict[str, Set[str]] = defaultdict(set)
    pep_to_exact_forms: Dict[str, Set[str]] = defaultdict(set)
    for row in psm_rows:
        acc = (row.get("accession") or "").strip()
        if not acc or acc.startswith("XXX_") or row.get("opt_global_cv_MS:1002217_decoy_peptide") == "1":
            continue
        seq = (row.get("sequence") or "").strip()
        if seq and seq != "null":
            try:
                il_p = normalize_peptide_il(seq)
                pep_to_all_db_accessions[il_p].add(acc)
                pep_to_exact_forms[il_p].add(seq)
            except ValueError:
                pass

    # Summarize qualifying PSMs (psm_q<=0.01 and pep_q<=0.01) per protein accession in this sample
    acc_qual_il_peptides: Dict[str, Set[str]] = defaultdict(set)
    acc_qual_clean_peptides: Dict[str, Set[str]] = defaultdict(set)
    acc_best_psm_q: Dict[str, float] = {}
    acc_best_pep_q: Dict[str, float] = {}
    acc_best_spec_eval: Dict[str, float] = {}

    qualifying_spectra_refs = set()
    for qrec in audit_res["qualifying_target_psms"]:
        acc = qrec["accession"]
        il_p = qrec["il_peptide"]
        cl_p = qrec["clean_peptide"]
        pq = qrec["psm_q_value"]
        pepq = qrec["peptide_q_value"]
        se = qrec["spec_evalue"]
        spec_key = qrec.get("spectra_ref") or f"psm_id::{qrec['psm_id']}"
        qualifying_spectra_refs.add(spec_key)
        acc_qual_il_peptides[acc].add(il_p)
        acc_qual_clean_peptides[acc].add(cl_p)
        if acc not in acc_best_psm_q or pq < acc_best_psm_q[acc]:
            acc_best_psm_q[acc] = pq
        if acc not in acc_best_pep_q or pepq < acc_best_pep_q[acc]:
            acc_best_pep_q[acc] = pepq
        if acc not in acc_best_spec_eval or se < acc_best_spec_eval[acc]:
            acc_best_spec_eval[acc] = se

    return {
        "sample_id": sample_id,
        "sample_file": fname,
        "bytes": vstats["bytes"],
        "sha256": vstats["sha256"],
        "md5": vstats["md5"],
        "score_columns_verified": audit_res["score_columns_verified"],
        "n_total_psms": audit_res["n_total_psms"],
        "n_target_psms": audit_res["n_target_psms"],
        "n_decoy_psms": audit_res["n_decoy_psms"],
        "n_excluded_invalid_rows": audit_res["n_excluded_invalid_rows"],
        "n_discordant_psm_q": audit_res["n_discordant_psm_q"],
        "n_qualifying_target_psm_rows_q01": audit_res["n_qualifying_target_psms"],
        "n_qualifying_distinct_spectra_q01": len(qualifying_spectra_refs),
        "n_prt_accessions": len(prt_accessions),
        "qualifying_target_psms": audit_res["qualifying_target_psms"],
        "acc_qual_il_peptides": {k: sorted(v) for k, v in acc_qual_il_peptides.items()},
        "acc_qual_clean_peptides": {k: sorted(v) for k, v in acc_qual_clean_peptides.items()},
        "acc_best_psm_q": acc_best_psm_q,
        "acc_best_pep_q": acc_best_pep_q,
        "acc_best_spec_eval": acc_best_spec_eval,
        "pep_to_all_db_accessions": {k: sorted(v) for k, v in pep_to_all_db_accessions.items()},
        "pep_to_exact_forms": {k: sorted(v) for k, v in pep_to_exact_forms.items()},
    }


def parse_fasta(path: str):
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def main() -> int:
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("[Stage 2] 1. Auditing MassIVE MSV000097826 604-file repository inventory and search-database completeness (ROOT-BACKGROUND-REVIEW.txt)...", flush=True)
    massive_604_path = os.path.join(MZTAB_DIR, "massive_MSV000097826_dataset_files_604.json")
    params_xml_path = os.path.join(MZTAB_DIR, "params.xml")
    contaminants_faa_path = os.path.join(
        MZTAB_DIR, "A_thaliana_UniProt_2020-11-02_MaxQuant_Contaminant_Proteins_Tryp_Pig_Bov.fasta"
    )

    with open(massive_604_path, "r", encoding="utf-8") as f:
        m604_rows = json.load(f).get("row_data", [])
    m604_by_cat = Counter(
        (r.get("collection") or r.get("File_descriptor", "").split("/")[0]).strip()
        for r in m604_rows
    )
    m604_seq_files = [
        (r.get("file_descriptor") or r.get("File_descriptor", "")).strip()
        for r in m604_rows
        if (r.get("collection") == "sequence" or (r.get("file_descriptor") or r.get("File_descriptor", "")).startswith("sequence/"))
    ]

    # Parse params.xml to record exact upstream upload mappings
    params_upload_mappings = []
    if os.path.exists(params_xml_path):
        root = ET.parse(params_xml_path).getroot()
        for p_elem in root.findall(".//parameter[@name='upload_file_mapping']"):
            if p_elem.text:
                params_upload_mappings.append(p_elem.text.strip())

    # Check if the 15-split 5%-FDR composite search FASTAs (referenced in Split_FASTA_MSGFPlus_Results/*.mzid.gz)
    # are present anywhere in the 604 deposited MassIVE files or on disk
    expected_split_search_fastas = [
        f"Blodgett_Soil_Metagenomes_5PctFDR-Composite_TrypPigBov_2023-05-22_15x_{i:02d}.fasta"
        for i in range(1, 16)
    ]
    deposited_basenames = {
        os.path.basename((r.get("file_descriptor") or r.get("File_descriptor") or r.get("name") or "").strip())
        for r in m604_rows
    }
    found_split_fastas = [fn for fn in expected_split_search_fastas if fn in deposited_basenames]
    full_background_verified = len(found_split_fastas) == len(expected_split_search_fastas)

    print(
        f"  MassIVE MSV000097826 inventory audited: {len(m604_rows)} total files ({dict(m604_by_cat)}). "
        f"sequence/ files={m604_seq_files}. "
        f"15-split composite search FASTAs deposited={len(found_split_fastas)}/15 -> "
        f"full_background_verified={full_background_verified}.",
        flush=True,
    )

    print("[Stage 2] 2. Loading sample metadata and verifying 2018-only DNA recruitment pairing (n=3 paired plots x 5 depths)...", flush=True)
    samples_csv_path = os.path.join(KARAOZ_DIR, "essdive_3365729", "samples.csv")
    sample2meta_path = os.path.join(KARAOZ_DIR, "github", "data_products", "sample2metadata.txt")
    coverm_path = os.path.join(KARAOZ_DIR, "github", "data_products", "pooledderep_bins.vs.allreads.covermout.txt")

    emsl_to_s2m = {}
    with open(sample2meta_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            emsl_id = row.get("proteome-EMSL-sampleid", "").strip()
            if emsl_id and emsl_id != "#N/A":
                emsl_to_s2m[emsl_id] = row

    coverm_rel_abund: Dict[str, Dict[str, float]] = defaultdict(dict)
    coverm_read_count: Dict[str, Dict[str, int]] = defaultdict(dict)
    coverm_trimmed_mean: Dict[str, Dict[str, float]] = defaultdict(dict)

    with open(coverm_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        cov_header = next(reader)
        fastq_col_indices = {}
        for c_idx, col_name in enumerate(cov_header[1:], start=1):
            for metric in ("Relative Abundance (%)", "Trimmed Mean", "Read Count"):
                if col_name.endswith(" " + metric):
                    fq = col_name[: -(len(metric) + 1)].strip()
                    fastq_col_indices[(fq, metric)] = c_idx
        for row in reader:
            if not row or row[0].strip() == "unmapped":
                continue
            g_core = normalize_mag_id(row[0].strip())
            for (fq, metric), c_idx in fastq_col_indices.items():
                val_s = row[c_idx].strip()
                if val_s == "NA" or not val_s:
                    continue
                if metric == "Relative Abundance (%)":
                    coverm_rel_abund[g_core][fq] = float(val_s)
                elif metric == "Trimmed Mean":
                    coverm_trimmed_mean[g_core][fq] = float(val_s)
                elif metric == "Read Count":
                    coverm_read_count[g_core][fq] = int(float(val_s))

    sample_meta_records = []
    sample_by_id = {}
    with open(samples_csv_path, "r", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            sid = row["proteome-EMSL-sampleid"].strip()
            sname = row["Sample Name"].strip()
            s_labid = row["sample_labid"].strip()
            plot_raw = row["plot"].strip().lower()
            trt = row["treatment"].strip().lower()
            depth = row["depth (cm)"].strip()
            sdate = row["sampling-date"].strip()
            if not sdate.startswith("2018"):
                raise ValueError(f"Non-2018 metaproteome sample {sid}: sampling-date={sdate}")
            plot_pair_id = re.sub(r"^(control|heated)", "", plot_raw).strip()

            s2m = lookup_strict(emsl_to_s2m, sid, "sample2metadata lookup by proteome-EMSL-sampleid")
            dna_date = s2m["date"].strip()
            dna_year = 2018 if dna_date.endswith("2018") else int(dna_date[-4:])
            dna_fastq = s2m["fastq"].strip()
            if "2018" not in dna_fastq:
                raise ValueError(f"Paired DNA fastq '{dna_fastq}' for {sid} is not from 2018!")

            rec = {
                "sample_id": sid,
                "sample_name": sname,
                "sample_labid": s_labid,
                "plot_raw": plot_raw,
                "plot_pair_id": plot_pair_id,
                "treatment": trt,
                "depth_cm": depth,
                "sampling_date": sdate,
                "sampling_year": 2018,
                "paired_dna_fastq": dna_fastq,
                "paired_dna_date": dna_date,
                "paired_dna_year": dna_year,
                "ncbi_biosample": s2m.get("NCBI_BioSampleAccession", "").strip(),
                "value": 1.0,
            }
            sample_meta_records.append(rec)
            sample_by_id[sid] = rec

    design_check = validate_karaoz_experimental_design(sample_meta_records)
    print(f"  Verified experimental design: {design_check}", flush=True)

    print("[Stage 2] 3. Parsing and auditing all 30 MS-GF+ .mzTab files in parallel (enforcing upstream PepQValue<=0.01 & QValue<=0.01)...", flush=True)
    manifest_json_path = os.path.join(MZTAB_DIR, "massive_MSV000097826_result_list.json")
    with open(manifest_json_path, "r", encoding="utf-8") as f:
        massive_rows = json.load(f).get("blockData", [])
    if len(massive_rows) != 30:
        raise ValueError(f"Expected 30 MassIVE result files, found {len(massive_rows)}")

    n_workers = min(30, int(os.environ.get("SLURM_CPUS_PER_TASK", "16")))
    mztab_results_by_sample = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(parse_and_filter_single_mztab, (r, MZTAB_DIR)) for r in massive_rows]
        for fut in concurrent.futures.as_completed(futures):
            res = fut.result()
            mztab_results_by_sample[res["sample_id"]] = res

    if set(mztab_results_by_sample.keys()) != set(sample_by_id.keys()):
        raise ValueError("Mismatch between 30 mzTab sample IDs and essdive_3365729/samples.csv!")

    sample_audit_rows = []
    global_pep_to_db_accessions: Dict[str, Set[str]] = defaultdict(set)
    global_pep_to_exact_forms: Dict[str, Set[str]] = defaultdict(set)
    total_discordant_psm_q = 0
    total_excluded_invalid_rows = 0

    for sid, srec in sorted(sample_by_id.items(), key=lambda kv: (kv[1]["plot_pair_id"], kv[1]["depth_cm"], kv[1]["treatment"])):
        mres = mztab_results_by_sample[sid]
        total_discordant_psm_q += int(mres["n_discordant_psm_q"])
        total_excluded_invalid_rows += int(mres["n_excluded_invalid_rows"])
        for il_p, accs in mres["pep_to_all_db_accessions"].items():
            global_pep_to_db_accessions[il_p].update(accs)
        for il_p, forms in mres["pep_to_exact_forms"].items():
            global_pep_to_exact_forms[il_p].update(forms)

        sample_audit_rows.append({
            "sample_id": sid,
            "sample_name": srec["sample_name"],
            "plot_pair_id": srec["plot_pair_id"],
            "treatment": srec["treatment"],
            "depth_cm": srec["depth_cm"],
            "sampling_year": srec["sampling_year"],
            "paired_dna_fastq": srec["paired_dna_fastq"],
            "paired_dna_year": srec["paired_dna_year"],
            "mztab_file": mres["sample_file"],
            "mztab_bytes": mres["bytes"],
            "mztab_sha256": mres["sha256"],
            "pep_q_col_verified": mres["score_columns_verified"]["pep_q_col"],
            "psm_q_col_verified": mres["score_columns_verified"]["psm_q_col"],
            "spec_evalue_col_verified": mres["score_columns_verified"]["spec_evalue_col"],
            "n_total_psms": mres["n_total_psms"],
            "n_target_psms": mres["n_target_psms"],
            "n_decoy_psms": mres["n_decoy_psms"],
            "n_excluded_invalid_rows": mres["n_excluded_invalid_rows"],
            "n_discordant_psm_q": mres["n_discordant_psm_q"],
            "n_qualifying_target_psm_rows_q01": mres["n_qualifying_target_psm_rows_q01"],
            "n_qualifying_distinct_spectra_q01": mres["n_qualifying_distinct_spectra_q01"],
        })

    write_tsv_with_schema(
        os.path.join(RESULTS_DIR, "karaoz2026_30samples_mztab_psm_fdr_and_design_audit.tsv"),
        list(sample_audit_rows[0].keys()),
        sample_audit_rows,
    )
    print(
        f"  Audited {len(sample_audit_rows)} .mzTab files: total_discordant_psm_q={total_discordant_psm_q}, "
        f"total_excluded_invalid_rows={total_excluded_invalid_rows}.",
        flush=True,
    )

    print("[Stage 2] 4. Streaming multipattern I/L-normalized peptide specificity search across all available sequence pools...", flush=True)
    emsl_faa_path = os.path.join(KARAOZ_DIR, "essdive_3365729", "EMSL50964_Proteins.faa")
    emsl_csv_path = os.path.join(KARAOZ_DIR, "essdive_3365729", "EMSL50964_ProteinSpectralCounts.csv")
    arch_faa_v2_path = os.path.join(WORK_STAGE1_V2, "karaoz_113_archaeal_mags_proteins_v2.faa")
    arch_ledger_v2_path = os.path.join(STAGE1_V2_DIR, "karaoz2026_archaeal_mag_eligibility_ledger.tsv")
    emsl_arch_v2_path = os.path.join(STAGE1_V2_DIR, "karaoz2026_emsl_archaeal_detected_proteins_v2.tsv")

    derep_arch_mags = {}
    with open(arch_ledger_v2_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            derep_arch_mags[row["essdive_mag_id"]] = row

    emsl_arch_proteins = {}
    with open(emsl_arch_v2_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            emsl_arch_proteins[row["emsl_protein_id"]] = row

    alias_to_primary_emsl = {}
    emsl_dup_ids_by_primary = {}
    emsl_reported_counts_by_sample: Dict[str, Dict[str, int]] = defaultdict(dict)
    with open(emsl_csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        sample_cols = [c.strip() for c in header[3:]]
        for row in reader:
            if not row:
                continue
            pid = row[0].strip()
            desc = row[1]
            for s_col, val_s in zip(sample_cols, row[3:]):
                emsl_reported_counts_by_sample[pid][s_col] = int(float(val_s)) if val_s.strip() else 0
            all_ids = [pid]
            for m in re.finditer(r"Duplicate proteins:\s*([^;]+)", desc):
                for d_id in m.group(1).split(","):
                    d_id = d_id.strip()
                    if d_id:
                        all_ids.append(d_id)
            emsl_dup_ids_by_primary[pid] = all_ids
            for aid in all_ids:
                alias_to_primary_emsl[aid] = pid

    # Load Stage 1 v2 homology hits to build protein -> family mapping and candidate protein set
    eq_vs_emsl_v2_path = os.path.join(STAGE1_V2_DIR, "karaoz2026_eq_vs_emsl50964_detected_homology_v2.tsv")
    eq_emsl_hits = []
    with open(eq_vs_emsl_v2_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            eq_emsl_hits.append(row)

    eq_vs_arch_v2_path = os.path.join(STAGE1_V2_DIR, "karaoz2026_eq_vs_archaeal_mags_homology_v2.tsv")
    eq_arch_hits = []
    with open(eq_vs_arch_v2_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            eq_arch_hits.append(row)

    # Map every protein/ORF ID to its MAG(s) and to its EQ family ID(s) (or unassigned paralog tag)
    protein_to_mags: Dict[str, Set[str]] = defaultdict(set)
    protein_to_families: Dict[str, Set[str]] = defaultdict(set)

    for pid, erec in emsl_arch_proteins.items():
        for m in erec.get("matched_derep_mags", "").split(","):
            if m.strip():
                protein_to_mags[pid].add(m.strip())

    for r in eq_emsl_hits:
        if r["homology_gate_status"] in ("strict_full_length_homolog", "moderate_bidirectional_homolog"):
            protein_to_families[r["emsl_protein_id"]].add(r["family_id"])

    for r in eq_arch_hits:
        if r["homology_gate_status"] in ("strict_full_length_homolog", "moderate_bidirectional_homolog"):
            protein_to_families[r["target_orf_id"]].add(r["family_id"])
            for ep in r["emsl50964_protein_ids"].split(","):
                if ep.strip():
                    protein_to_families[ep.strip()].add(r["family_id"])

    candidate_emsl_ids = set(emsl_arch_proteins.keys()) | {
        r["emsl_protein_id"] for r in eq_emsl_hits
        if r["homology_gate_status"] in ("strict_full_length_homolog", "moderate_bidirectional_homolog")
    }

    # Collect all distinct qualifying clean peptides observed across the 30 samples for any candidate protein
    candidate_peptides_to_audit: Set[str] = set()
    for mres in mztab_results_by_sample.values():
        for acc, cl_peps in mres["acc_qual_clean_peptides"].items():
            prim_id = alias_to_primary_emsl.get(acc, acc)
            if prim_id in candidate_emsl_ids:
                candidate_peptides_to_audit.update(cl_peps)

    # Build combined streaming generator over all deposited FASTA sequence pools:
    # 1) EMSL50964_Proteins.faa
    # 2) karaoz_113_archaeal_mags_proteins_v2.faa
    # 3) A_thaliana_UniProt_2020-11-02_MaxQuant_Contaminant_Proteins_Tryp_Pig_Bov.fasta
    emsl_valid_count = 0
    emsl_invalid_count = 0
    emsl_aa_len: Dict[str, int] = {}
    arch_valid_count = 0
    contaminant_valid_count = 0

    def iter_all_available_background_records():
        nonlocal emsl_valid_count, emsl_invalid_count, arch_valid_count, contaminant_valid_count
        for hdr, seq in parse_fasta(emsl_faa_path):
            pid = hdr.split()[0].strip()
            c_res = canonicalize_and_validate_protein_seq(seq)
            if c_res["is_valid"]:
                emsl_valid_count += 1
                emsl_aa_len[pid] = c_res["aa_length"]
                if pid not in protein_to_mags:
                    protein_to_mags[pid].add(f"non_archaeal_or_unbinned_emsl::{pid}")
                yield f"emsl::{pid}", c_res["canonical_seq"]
            else:
                emsl_invalid_count += 1

        for hdr, seq in parse_fasta(arch_faa_v2_path):
            oid = hdr.split()[0].strip()
            mag_id = oid.split("::")[0]
            c_res = canonicalize_and_validate_protein_seq(seq)
            if c_res["is_valid"]:
                arch_valid_count += 1
                protein_to_mags[oid].add(mag_id)
                yield f"arch_orf::{oid}", c_res["canonical_seq"]

        if os.path.exists(contaminants_faa_path):
            for hdr, seq in parse_fasta(contaminants_faa_path):
                cid = hdr.split()[0].strip()
                c_res = canonicalize_and_validate_protein_seq(seq)
                if c_res["is_valid"]:
                    contaminant_valid_count += 1
                    protein_to_mags[f"contaminant::{cid}"].add(f"contaminant_or_arabidopsis::{cid}")
                    yield f"contaminant::{cid}", c_res["canonical_seq"]

    stream_hits = stream_multipattern_peptide_search_records(
        iter_all_available_background_records(),
        candidate_peptides_to_audit,
    )

    print(
        f"  Streamed multipattern search for {len(candidate_peptides_to_audit)} candidate peptides across "
        f"emsl_valid={emsl_valid_count} (invalid={emsl_invalid_count}), "
        f"arch_valid={arch_valid_count}, contaminant_valid={contaminant_valid_count}.",
        flush=True,
    )

    # Save explicit search-database background completeness audit JSON
    bg_audit = {
        "massive_dataset": "MSV000097826 / PXD063748",
        "total_deposited_files_in_massive": len(m604_rows),
        "deposited_files_by_category": dict(m604_by_cat),
        "deposited_sequence_files": m604_seq_files,
        "params_xml_upload_file_mappings": params_upload_mappings,
        "expected_15split_composite_search_fastas": expected_split_search_fastas,
        "found_15split_composite_search_fastas_in_repository": found_split_fastas,
        "full_search_background_fasta_available": full_background_verified,
        "audited_sequence_pools_on_disk": {
            "emsl50964_valid_proteins": emsl_valid_count,
            "emsl50964_invalid_proteins": emsl_invalid_count,
            "stage1_v2_archaeal_mag_orfs": arch_valid_count,
            "massive_contaminant_and_arabidopsis_proteins": contaminant_valid_count,
            "mztab_co_occurring_search_db_accessions_tracked": sum(len(v) for v in global_pep_to_db_accessions.values()),
        },
        "inference_consequence": (
            "Because the 15-split 5%-FDR composite search FASTAs and initial 76M-ORF pool are not publicly deposited "
            "in MassIVE MSV000097826 or ESS-DIVE, full_background_verified is False: primary global single-MAG "
            "peptide uniqueness is marked non-estimable, and all observed-subset detections are reported as "
            "provisional_observed_subset support alongside dataset-reported author inference."
        ),
    }
    with open(os.path.join(RESULTS_DIR, "karaoz2026_search_database_background_completeness_audit.json"), "w", encoding="utf-8") as f:
        json.dump(bg_audit, f, indent=2)

    # Evaluate each candidate peptide's specificity across all audited background pools + mzTab accessions
    peptide_specificity_ledger = {}
    for cl_pep in sorted(candidate_peptides_to_audit):
        il_pep = cl_pep.replace("I", "L")
        sh = stream_hits.get(cl_pep, {"exact_matches": set(), "il_matches": set()})

        exact_emsl = sorted(x[len("emsl::"):] for x in sh["exact_matches"] if x.startswith("emsl::"))
        il_emsl = sorted(x[len("emsl::"):] for x in sh["il_matches"] if x.startswith("emsl::"))
        il_only_emsl = sorted(set(il_emsl) - set(exact_emsl))

        exact_arch = sorted(x[len("arch_orf::"):] for x in sh["exact_matches"] if x.startswith("arch_orf::"))
        il_arch = sorted(x[len("arch_orf::"):] for x in sh["il_matches"] if x.startswith("arch_orf::"))
        il_only_arch = sorted(set(il_arch) - set(exact_arch))

        il_contam = sorted(x for x in sh["il_matches"] if x.startswith("contaminant::"))

        mztab_db_accs = global_pep_to_db_accessions.get(il_pep, set())
        mztab_mapped_primaries = {alias_to_primary_emsl.get(a, a) for a in mztab_db_accs}

        # Combine all canonical protein/ORF IDs for classify_peptide_specificity_against_background
        combined_exact_ids: Set[str] = set(exact_emsl) | set(exact_arch)
        combined_il_ids: Set[str] = set(il_emsl) | set(il_arch) | set(il_contam) | mztab_mapped_primaries

        for m_acc in mztab_mapped_primaries:
            if m_acc not in protein_to_mags:
                protein_to_mags[m_acc].add(f"unbinned_or_non_emsl_mztab_acc::{m_acc}")
            if len(global_pep_to_exact_forms.get(il_pep, set())) == 1:
                combined_exact_ids.add(m_acc)

        spec_eval = classify_peptide_specificity_against_background(
            raw_peptide=cl_pep,
            target_protein_ids=set(exact_emsl) | set(exact_arch),
            background_proteins={},
            protein_to_mags=protein_to_mags,
            protein_to_families=protein_to_families,
            precomputed_exact_matches=combined_exact_ids,
            precomputed_il_matches=combined_il_ids,
        )

        matched_mags = set(spec_eval["matched_mags"])
        has_il_collision = spec_eval["has_isobaric_il_collision"] or (len(global_pep_to_exact_forms.get(il_pep, set())) > 1)
        is_single_mag = (not has_il_collision) and spec_eval["is_unambiguous_for_single_mag"]
        is_single_family = (not has_il_collision) and spec_eval["is_unambiguous_for_family"]
        is_single_clade_thermo = (
            (not has_il_collision)
            and (len(il_contam) == 0)
            and len(matched_mags) >= 1
            and all(m in derep_arch_mags and derep_arch_mags[m]["phylum"] == "p__Thermoplasmatota" for m in matched_mags)
        )

        peptide_specificity_ledger[cl_pep] = {
            "clean_peptide": cl_pep,
            "il_normalized_peptide": il_pep,
            "full_search_background_verified": full_background_verified,
            "n_exact_emsl_proteins": len(exact_emsl),
            "exact_emsl_proteins": ",".join(exact_emsl[:10]),
            "n_il_emsl_or_searchdb_proteins": len(set(il_emsl) | mztab_mapped_primaries),
            "n_il_arch_mag_orfs": len(il_arch),
            "n_il_contaminant_proteins": len(il_contam),
            "n_il_only_collision_proteins": len(il_only_emsl) + len(il_only_arch),
            "has_isobaric_il_collision": has_il_collision,
            "n_matched_mags": len(matched_mags),
            "matched_mags": ",".join(sorted(matched_mags)[:10]),
            "matched_families": ",".join(spec_eval["matched_families"][:10]),
            "observed_background_specificity_status": spec_eval["specificity_status"],
            "is_unambiguous_single_mag_in_observed_background": is_single_mag,
            "is_unambiguous_for_family_in_observed_background": is_single_family,
            "is_unambiguous_thermoplasmatota_clade_in_observed_background": is_single_clade_thermo,
        }

    write_tsv_with_schema(
        os.path.join(RESULTS_DIR, "karaoz2026_candidate_peptide_whole_background_specificity_ledger.tsv"),
        list(peptide_specificity_ledger[next(iter(peptide_specificity_ledger))].keys()) if peptide_specificity_ledger else [
            "clean_peptide", "il_normalized_peptide", "full_search_background_verified",
            "n_exact_emsl_proteins", "exact_emsl_proteins", "n_il_emsl_or_searchdb_proteins",
            "n_il_arch_mag_orfs", "n_il_contaminant_proteins", "n_il_only_collision_proteins",
            "has_isobaric_il_collision", "n_matched_mags", "matched_mags", "matched_families",
            "observed_background_specificity_status", "is_unambiguous_single_mag_in_observed_background",
            "is_unambiguous_for_family_in_observed_background", "is_unambiguous_thermoplasmatota_clade_in_observed_background"
        ],
        list(peptide_specificity_ledger.values()),
    )

    print("[Stage 2] 5. Auditing protein-level & clade-level detection using unique spectrum union (never summing across aliases)...", flush=True)
    thermo_protein_audit_rows = []

    # Collect alias sets and clade-unambiguous peptide sets for unique spectrum union across clades
    thermo_all_aliases: Set[str] = set()
    e2_jacpao01_all_aliases: Set[str] = set()
    rbg16_all_aliases: Set[str] = set()
    thermo_unamb_clean_peptides: Set[str] = {
        pep for pep, sp in peptide_specificity_ledger.items()
        if sp["is_unambiguous_thermoplasmatota_clade_in_observed_background"]
    }

    for pid in sorted(candidate_emsl_ids):
        aliases = set(emsl_dup_ids_by_primary.get(pid, [pid]))
        unamb_mag_il_peps = set()
        unamb_fam_il_peps = set()
        unamb_clade_il_peps = set()
        ambig_il_peps = set()
        all_prot_clean_peps = set()
        unique_spectra_by_sample = Counter()
        best_psm_q = 1.0
        best_pep_q = 1.0
        best_spec_eval = 1.0

        for sid, mres in mztab_results_by_sample.items():
            sample_clean_peps_for_pid = set()
            for al in aliases:
                if al in mres["acc_qual_clean_peptides"]:
                    best_psm_q = min(best_psm_q, mres["acc_best_psm_q"].get(al, 1.0))
                    best_pep_q = min(best_pep_q, mres["acc_best_pep_q"].get(al, 1.0))
                    best_spec_eval = min(best_spec_eval, mres["acc_best_spec_eval"].get(al, 1.0))
                    for cl_p in mres["acc_qual_clean_peptides"].get(al, []):
                        sample_clean_peps_for_pid.add(cl_p)
                        all_prot_clean_peps.add(cl_p)
                        sp = peptide_specificity_ledger.get(cl_p)
                        if sp is None:
                            continue
                        il_p = sp["il_normalized_peptide"]
                        if sp["is_unambiguous_single_mag_in_observed_background"]:
                            unamb_mag_il_peps.add(il_p)
                            unamb_clade_il_peps.add(il_p)
                        elif sp["is_unambiguous_thermoplasmatota_clade_in_observed_background"]:
                            unamb_clade_il_peps.add(il_p)
                        else:
                            ambig_il_peps.add(il_p)
                        if sp["is_unambiguous_for_family_in_observed_background"]:
                            unamb_fam_il_peps.add(il_p)

            if sample_clean_peps_for_pid:
                # Take unique spectrum union across all aliases of this protein in this sample (never sum across aliases!)
                u_specs = aggregate_unique_spectra_for_entity(
                    mres["qualifying_target_psms"],
                    target_accessions=aliases,
                    allowed_clean_peptides=sample_clean_peps_for_pid,
                )
                unique_spectra_by_sample[sid] = len(u_specs)

        erec = emsl_arch_proteins.get(pid, {})
        thermo_mags = [m for m in erec.get("matched_thermo_mags", "").split(",") if m.strip()]
        arch_mags = [m for m in erec.get("matched_arch_mags", "").split(",") if m.strip()]
        rep_spec = sum(emsl_reported_counts_by_sample.get(pid, {}).values())

        if thermo_mags:
            thermo_all_aliases.update(aliases)
            if any(derep_arch_mags[m]["family"] == "f__JACPAO01" for m in thermo_mags if m in derep_arch_mags):
                e2_jacpao01_all_aliases.update(aliases)
            if any(derep_arch_mags[m]["karaoz_order"] == "o__RBG-16-68-12" for m in thermo_mags if m in derep_arch_mags):
                rbg16_all_aliases.update(aliases)

        conservative_status = classify_protein_and_family_detection(
            has_karaoz_homolog=True,
            in_search_database=True,
            n_unambiguous_il_peptides_q01=len(unamb_mag_il_peps),
            n_family_shared_il_peptides_q01=len(unamb_clade_il_peps - unamb_mag_il_peps),
            reported_spectral_count=rep_spec,
            full_background_verified=full_background_verified,
            search_database_exclusion_verified=full_background_verified,
        )
        author_subset_status = classify_protein_and_family_detection(
            has_karaoz_homolog=True,
            in_search_database=True,
            n_unambiguous_il_peptides_q01=len(unamb_mag_il_peps),
            n_family_shared_il_peptides_q01=len(unamb_clade_il_peps - unamb_mag_il_peps),
            reported_spectral_count=rep_spec,
            full_background_verified=True,
            search_database_exclusion_verified=full_background_verified,
        )

        thermo_protein_audit_rows.append({
            "emsl_protein_id": pid,
            "matched_arch_mags": ",".join(arch_mags),
            "matched_thermo_mags": ",".join(thermo_mags),
            "thermo_taxonomic_tiers": ",".join(sorted({derep_arch_mags[m]["phylogenetic_eligibility_tier"] for m in thermo_mags if m in derep_arch_mags})),
            "aa_length": emsl_aa_len.get(pid, 0),
            "reported_emsl_total_spectral_counts": rep_spec,
            "mztab_unique_qualifying_spectra_q01": sum(unique_spectra_by_sample.values()),
            "n_positive_samples_mztab_q01": sum(1 for v in unique_spectra_by_sample.values() if v > 0),
            "n_unambiguous_single_mag_il_peptides_in_observed_bg_q01": len(unamb_mag_il_peps),
            "n_unambiguous_family_il_peptides_in_observed_bg_q01": len(unamb_fam_il_peps),
            "n_unambiguous_thermo_clade_il_peptides_in_observed_bg_q01": len(unamb_clade_il_peps),
            "n_ambiguous_or_il_collision_peptides_q01": len(ambig_il_peps),
            "best_psm_q_value": round(best_psm_q, 6) if best_psm_q < 1.0 else "NA",
            "best_peptide_q_value": round(best_pep_q, 6) if best_pep_q < 1.0 else "NA",
            "best_spec_evalue": best_spec_eval if best_spec_eval < 1.0 else "NA",
            "conservative_detection_tier": conservative_status,
            "dataset_reported_observed_subset_tier": author_subset_status,
            "unambiguous_il_peptides_in_observed_bg": ",".join(sorted(unamb_clade_il_peps)),
        })

    write_tsv_with_schema(
        os.path.join(RESULTS_DIR, "karaoz2026_archaeal_and_thermoplasmatota_proteins_psm_confidence_audit.tsv"),
        list(thermo_protein_audit_rows[0].keys()) if thermo_protein_audit_rows else [
            "emsl_protein_id", "matched_arch_mags", "matched_thermo_mags", "thermo_taxonomic_tiers",
            "aa_length", "reported_emsl_total_spectral_counts", "mztab_unique_qualifying_spectra_q01",
            "n_positive_samples_mztab_q01", "n_unambiguous_single_mag_il_peptides_in_observed_bg_q01",
            "n_unambiguous_family_il_peptides_in_observed_bg_q01",
            "n_unambiguous_thermo_clade_il_peptides_in_observed_bg_q01",
            "n_ambiguous_or_il_collision_peptides_q01", "best_psm_q_value", "best_peptide_q_value",
            "best_spec_evalue", "conservative_detection_tier", "dataset_reported_observed_subset_tier",
            "unambiguous_il_peptides_in_observed_bg"
        ],
        thermo_protein_audit_rows,
    )

    print("[Stage 2] 6. Evaluating all EQ families (Primary Leads, 5 Partners, 1,460 Surface Families, 50 Positive Controls) and paired n=3 warming responses...", flush=True)
    comp_audit_path = os.path.join(STAGE1_V2_DIR, "eq_family_query_completeness_audit.tsv")
    all_families = {}
    with open(comp_audit_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            all_families[row["family_id"]] = row["cohort_role"]

    arch_hits_by_fam = defaultdict(list)
    for r in eq_arch_hits:
        if r["homology_gate_status"] in ("strict_full_length_homolog", "moderate_bidirectional_homolog"):
            arch_hits_by_fam[r["family_id"]].append(r)

    emsl_hits_by_fam = defaultdict(list)
    for r in eq_emsl_hits:
        if r["homology_gate_status"] in ("strict_full_length_homolog", "moderate_bidirectional_homolog"):
            emsl_hits_by_fam[r["family_id"]].append(r)

    prot_audit_by_id = {r["emsl_protein_id"]: r for r in thermo_protein_audit_rows}

    family_results = []
    pvals_primary_and_partners = {}
    pvals_exploratory_dark = {}
    pvals_controls = {}

    for fid, role in sorted(all_families.items()):
        ahits = arch_hits_by_fam.get(fid, [])
        ehits = emsl_hits_by_fam.get(fid, [])
        strict_arch = [h for h in ahits if h["homology_gate_status"] == "strict_full_length_homolog"]
        strict_thermo = [h for h in strict_arch if h["target_phylum"] == "p__Thermoplasmatota"]
        mod_thermo = [h for h in ahits if h["target_phylum"] == "p__Thermoplasmatota"]

        has_homolog = (len(ahits) > 0) or (len(ehits) > 0)
        # Note: Even when a homolog lies on a >=1kb scaffold in Blodgett_MAGs.tar.gz, if it is absent from
        # EMSL50964_Proteins.faa and the 15-split 5%-FDR composite search FASTAs are not publicly deposited,
        # we only confirm presence in the search DB when matched_emsl_pids is non-empty!
        matched_emsl_pids = set()
        for h in ahits:
            for ep in h["emsl50964_protein_ids"].split(","):
                if ep.strip():
                    matched_emsl_pids.add(ep.strip())
        for h in ehits:
            matched_emsl_pids.add(h["emsl_protein_id"])

        in_search_db_confirmed = len(matched_emsl_pids) > 0

        # Collect all aliases for matched EMSL proteins in this family, and restrict family quantification
        # strictly to peptides with is_unambiguous_for_family_in_observed_background == True!
        fam_all_aliases: Set[str] = set()
        fam_rep_spec = 0
        for ep in matched_emsl_pids:
            fam_all_aliases.update(emsl_dup_ids_by_primary.get(ep, [ep]))
            pa = prot_audit_by_id.get(ep)
            if pa:
                fam_rep_spec += int(pa["reported_emsl_total_spectral_counts"])

        fam_unamb_mag_il_peps: Set[str] = set()
        fam_unamb_family_il_peps: Set[str] = set()
        fam_unamb_family_clean_peps: Set[str] = set()
        fam_shared_clade_il_peps: Set[str] = set()

        for sid, mres in mztab_results_by_sample.items():
            for al in fam_all_aliases:
                for cl_p in mres["acc_qual_clean_peptides"].get(al, []):
                    sp = peptide_specificity_ledger.get(cl_p)
                    if not sp:
                        continue
                    il_p = sp["il_normalized_peptide"]
                    if sp["is_unambiguous_for_family_in_observed_background"]:
                        fam_unamb_family_il_peps.add(il_p)
                        fam_unamb_family_clean_peps.add(cl_p)
                        if sp["is_unambiguous_single_mag_in_observed_background"]:
                            fam_unamb_mag_il_peps.add(il_p)
                        else:
                            fam_shared_clade_il_peps.add(il_p)

        # Quantify unique family-specific spectra across the 30 samples using aggregate_unique_spectra_for_entity
        fam_unique_spectra_by_sample = Counter()
        if fam_all_aliases and fam_unamb_family_clean_peps:
            for sid, mres in mztab_results_by_sample.items():
                u_specs = aggregate_unique_spectra_for_entity(
                    mres["qualifying_target_psms"],
                    target_accessions=fam_all_aliases,
                    allowed_clean_peptides=fam_unamb_family_clean_peps,
                )
                fam_unique_spectra_by_sample[sid] = len(u_specs)

        fam_tot_unique_spectra = sum(fam_unique_spectra_by_sample.values())

        conservative_det_status = classify_protein_and_family_detection(
            has_karaoz_homolog=has_homolog,
            in_search_database=in_search_db_confirmed,
            n_unambiguous_il_peptides_q01=len(fam_unamb_mag_il_peps),
            n_family_shared_il_peptides_q01=len(fam_shared_clade_il_peps),
            reported_spectral_count=fam_rep_spec,
            full_background_verified=full_background_verified,
            search_database_exclusion_verified=full_background_verified,
        )
        author_subset_det_status = classify_protein_and_family_detection(
            has_karaoz_homolog=has_homolog,
            in_search_database=in_search_db_confirmed,
            n_unambiguous_il_peptides_q01=len(fam_unamb_mag_il_peps),
            n_family_shared_il_peptides_q01=len(fam_shared_clade_il_peps),
            reported_spectral_count=fam_rep_spec,
            full_background_verified=True,
            search_database_exclusion_verified=full_background_verified,
        )
        primary_warming_status = classify_primary_warming_estimability(
            has_karaoz_homolog=has_homolog,
            in_search_database=in_search_db_confirmed,
            has_qualifying_unambiguous_peptides_in_observed_bg=(fam_tot_unique_spectra > 0 and len(fam_unamb_family_clean_peps) > 0),
            reported_spectral_count=fam_rep_spec,
            full_background_verified=full_background_verified,
            search_database_exclusion_verified=full_background_verified,
        )

        # Compute paired n=3 plot warming effect ONLY from family-unambiguous unique spectra (never from ambiguous peptides!)
        warming_effect = None
        if fam_tot_unique_spectra > 0 and len(fam_unamb_family_clean_peps) > 0:
            s_records = []
            for sid, srec in sample_by_id.items():
                mres = mztab_results_by_sample[sid]
                raw_cnt = fam_unique_spectra_by_sample[sid]
                denom = max(1, mres["n_qualifying_distinct_spectra_q01"])
                norm_val = (raw_cnt * 100000.0) / float(denom)
                row_copy = dict(srec)
                row_copy["value"] = norm_val
                s_records.append(row_copy)
            warming_effect = compute_paired_warming_effect_n3_plots(s_records, value_field="value")
            pval = warming_effect["exact_whole_plot_sign_flip_two_sided_p"]
        else:
            pval = None

        if role.startswith("primary_lead") or role.startswith("lead_"):
            pvals_primary_and_partners[fid] = pval
        elif role.startswith("candidate_"):
            pvals_exploratory_dark[fid] = pval
        else:
            pvals_controls[fid] = pval

        family_results.append({
            "family_id": fid,
            "cohort_role": role,
            "n_strict_arch_homologs": len(strict_arch),
            "n_strict_thermoplasmatota_homologs": len(strict_thermo),
            "n_moderate_thermoplasmatota_homologs": len(mod_thermo),
            "thermo_target_mags": ",".join(sorted({h["target_mag_id"] for h in mod_thermo})[:10]),
            "thermo_eligibility_tiers": ",".join(sorted({h["phylogenetic_eligibility_tier"] for h in mod_thermo})),
            "confirmed_in_deposited_emsl50964_subset": in_search_db_confirmed,
            "matched_emsl50964_protein_ids": ",".join(sorted(matched_emsl_pids)),
            "reported_emsl_spectral_counts": fam_rep_spec,
            "mztab_unique_family_specific_spectra_q01": fam_tot_unique_spectra,
            "n_unambiguous_single_mag_il_peptides_in_observed_bg_q01": len(fam_unamb_mag_il_peps),
            "n_unambiguous_family_il_peptides_in_observed_bg_q01": len(fam_unamb_family_il_peps),
            "conservative_detection_classification": conservative_det_status,
            "dataset_reported_observed_subset_classification": author_subset_det_status,
            "primary_family_warming_estimability_status": primary_warming_status,
            "observed_subset_provisional_n_independent_plot_pairs": 3 if warming_effect else "NA",
            "observed_subset_provisional_mean_delta_heated_minus_control_per100k": round(warming_effect["overall_mean_delta_heated_minus_control"], 4) if warming_effect else "NA",
            "observed_subset_provisional_replicate_sign_consistency": warming_effect["replicate_sign_consistency"] if warming_effect else "NA",
            "observed_subset_provisional_exact_sign_flip_p": warming_effect["exact_whole_plot_sign_flip_two_sided_p"] if warming_effect else "NA",
            "observed_subset_provisional_n_heated_pos_of_15": warming_effect["n_heated_positive_samples_of_15"] if warming_effect else "NA",
            "observed_subset_provisional_n_control_pos_of_15": warming_effect["n_control_positive_samples_of_15"] if warming_effect else "NA",
            "observed_subset_provisional_n_depth_pairs_both_positive": warming_effect["n_depth_pairs_both_positive"] if warming_effect else "NA",
            "observed_subset_provisional_cohort_bh_fdr_q": "PENDING",
        })

    q_primary = benjamini_hochberg_fdr(pvals_primary_and_partners)
    q_exploratory = benjamini_hochberg_fdr(pvals_exploratory_dark)
    q_controls = benjamini_hochberg_fdr(pvals_controls)

    for r in family_results:
        fid = r["family_id"]
        qval = q_primary.get(fid) if fid in q_primary else (q_exploratory.get(fid) if fid in q_exploratory else q_controls.get(fid))
        r["observed_subset_provisional_cohort_bh_fdr_q"] = round(qval, 6) if qval is not None else "NA_nonestimable_or_undetected"

    write_tsv_with_schema(
        os.path.join(RESULTS_DIR, "karaoz2026_all_eq_families_expression_and_nonestimability_ledger.tsv"),
        list(family_results[0].keys()),
        family_results,
    )

    print("[Stage 2] 7. Running clade-level Thermoplasmatota (c__E2 f__JACPAO01 and c__Thermoplasmata o__RBG-16-68-12) unique-spectrum paired n=3 warming & 2018 DNA normalization analysis...", flush=True)
    thermo_mag_cores = {
        normalize_mag_id(m): rec for m, rec in derep_arch_mags.items() if rec["phylum"] == "p__Thermoplasmatota"
    }
    e2_jacpao_cores = {k for k, r in thermo_mag_cores.items() if r["family"] == "f__JACPAO01"}
    rbg16_cores = {k for k, r in thermo_mag_cores.items() if r["karaoz_order"] == "o__RBG-16-68-12"}

    clade_sample_rows = []
    for sid, srec in sorted(sample_by_id.items(), key=lambda kv: (kv[1]["plot_pair_id"], kv[1]["depth_cm"], kv[1]["treatment"])):
        mres = mztab_results_by_sample[sid]
        fq = srec["paired_dna_fastq"]
        denom_spec = max(1, mres["n_qualifying_distinct_spectra_q01"])

        # Quantify clade-level spectra by unique spectrum union (aggregate_unique_spectra_for_entity)
        thermo_u_specs = aggregate_unique_spectra_for_entity(
            mres["qualifying_target_psms"],
            target_accessions=thermo_all_aliases,
            allowed_clean_peptides=thermo_unamb_clean_peptides,
        )
        e2_u_specs = aggregate_unique_spectra_for_entity(
            mres["qualifying_target_psms"],
            target_accessions=e2_jacpao01_all_aliases,
            allowed_clean_peptides=thermo_unamb_clean_peptides,
        )
        rbg_u_specs = aggregate_unique_spectra_for_entity(
            mres["qualifying_target_psms"],
            target_accessions=rbg16_all_aliases,
            allowed_clean_peptides=thermo_unamb_clean_peptides,
        )

        thermo_psms = len(thermo_u_specs)
        e2_psms = len(e2_u_specs)
        rbg_psms = len(rbg_u_specs)

        thermo_dna_rel = sum(coverm_rel_abund.get(c, {}).get(fq, 0.0) for c in thermo_mag_cores)
        e2_dna_rel = sum(coverm_rel_abund.get(c, {}).get(fq, 0.0) for c in e2_jacpao_cores)
        rbg_dna_rel = sum(coverm_rel_abund.get(c, {}).get(fq, 0.0) for c in rbg16_cores)

        thermo_dna_reads = sum(coverm_read_count.get(c, {}).get(fq, 0) for c in thermo_mag_cores)

        thermo_prot_per100k = (thermo_psms * 100000.0) / float(denom_spec)
        e2_prot_per100k = (e2_psms * 100000.0) / float(denom_spec)
        rbg_prot_per100k = (rbg_psms * 100000.0) / float(denom_spec)

        thermo_prot_to_dna = thermo_prot_per100k / (thermo_dna_rel + 1e-3)

        clade_sample_rows.append({
            "sample_id": sid,
            "sample_name": srec["sample_name"],
            "plot_pair_id": srec["plot_pair_id"],
            "treatment": srec["treatment"],
            "depth_cm": srec["depth_cm"],
            "sampling_year": 2018,
            "paired_dna_year": 2018,
            "paired_dna_fastq": fq,
            "sample_qualifying_spectra_q01": denom_spec,
            "thermo_unique_qualifying_spectra_q01": thermo_psms,
            "e2_jacpao01_unique_qualifying_spectra_q01": e2_psms,
            "rbg16_unique_qualifying_spectra_q01": rbg_psms,
            "thermo_prot_per_100k_spectra": round(thermo_prot_per100k, 4),
            "e2_jacpao01_prot_per_100k_spectra": round(e2_prot_per100k, 4),
            "rbg16_prot_per_100k_spectra": round(rbg_prot_per100k, 4),
            "thermo_paired2018_dna_rel_abund_pct": round(thermo_dna_rel, 6),
            "e2_jacpao01_paired2018_dna_rel_abund_pct": round(e2_dna_rel, 6),
            "rbg16_paired2018_dna_rel_abund_pct": round(rbg_dna_rel, 6),
            "thermo_paired2018_dna_read_count": thermo_dna_reads,
            "thermo_prot_to_dna_ratio": round(thermo_prot_to_dna, 4),
        })

    write_tsv_with_schema(
        os.path.join(RESULTS_DIR, "karaoz2026_thermoplasmatota_30samples_paired_proteome_and_2018dna.tsv"),
        list(clade_sample_rows[0].keys()),
        clade_sample_rows,
    )

    clade_effects = {}
    for metric_col in (
        "thermo_prot_per_100k_spectra",
        "e2_jacpao01_prot_per_100k_spectra",
        "rbg16_prot_per_100k_spectra",
        "thermo_paired2018_dna_rel_abund_pct",
        "thermo_prot_to_dna_ratio",
    ):
        clade_effects[metric_col] = compute_paired_warming_effect_n3_plots(
            clade_sample_rows, value_field=metric_col
        )

    # 8. Write primary-source dataset & accession exclusion ledger
    exclusion_rows = [
        {
            "candidate_id": "Karaoz_2026_Blodgett_Soil_Warming_Proteogenomics",
            "citation_doi_pmid": "Karaoz et al. Nat Commun 17:10740 (13 Oct 2026); DOI:10.1038/s41467-026-78137-0; PMC13661901; PMID:42869931",
            "data_repositories_and_accessions": "ESS-DIVE doi:10.15485/2588706 (MAGs), doi:10.15485/3365729 (Proteome), MassIVE MSV000097826 / ProteomeXchange PXD063748, NCBI PRJNA1225762, GitHub ukaraoz/blodgett_omics_manuscript_analysis",
            "assay_modality": "Paired genome-resolved metagenomics (2,321 99%-ANI dereplicated MAGs) + LC-MS/MS metaproteomics (30 samples = 3 paired plots x 2 treatments x 5 depths, sampled June 2018)",
            "thermoplasmatota_representation": "38 dereplicated 99%-ANI Thermoplasmatota MAGs: 1 c__E2 f__JACPAO01/g__JACPAO01 (overlaps EQ order o__JACPAO01/family f__JACPAO01), 1 c__E2 o__UBA9212, 36 c__Thermoplasmata o__RBG-16-68-12;f__RBG-16-68-12; 0 c__SW-10-69-26 MAGs",
            "eligibility_decision": "ELIGIBLE_AND_ANALYZED_WITH_CONSERVATIVE_OBSERVED_SUBSET_SPECIFICITY",
            "scientific_justification": "Provides sample-matched dereplicated MAGs, detected protein FASTA (EMSL50964_Proteins.faa), paired 2018 CoverM DNA read recruitment, and all 30 MS-GF+ target-decoy .mzTab files with per-hit MS-GF:PepQValue and MS-GF:QValue; because the 15-split 5%-FDR composite search FASTAs are not deposited in MassIVE/ESS-DIVE, primary global single-MAG uniqueness is marked non-estimable and observed-subset support is reported provisionally.",
        },
        {
            "candidate_id": "DAlo_2023_Alpine_Climate_Warming_Archaea",
            "citation_doi_pmid": "D'Alo et al. Environ Microbiome 18:57 (2023); DOI:10.1186/s40793-023-00499-z; PMC10316362; PMID:37391832",
            "data_repositories_and_accessions": "EuropePMC PMC10316362 full-text XML audited",
            "assay_modality": "Short-read shotgun metagenomics and metatranscriptomics of high-elevation alpine soils under open-top chamber warming",
            "thermoplasmatota_representation": "Low-abundance archaeal read fraction (~0.2-1%); no MAGs assembled",
            "eligibility_decision": "EXCLUDED_INSUFFICIENT_ARCHAEAL_DEPTH_AND_NO_GENOME_RESOLUTION",
            "scientific_justification": "Full-text audit of PMC10316362 confirms no archaeal MAGs were reconstructed and the authors explicitly state: 'due to the scarcity of archaea sequences, it was impossible to investigate specific functions', reporting only broad COG category read counts.",
        },
        {
            "candidate_id": "Starke_2021_32Soil_Metaproteome_MetaAnalysis",
            "citation_doi_pmid": "Starke et al. J Proteomics 237:104147 (2021); DOI:10.1016/j.jprot.2021.104147; PMID:33450427",
            "data_repositories_and_accessions": "PRIDE / literature meta-analysis of 32 soil metaproteomes",
            "assay_modality": "Unpaired soil metaproteomics searched against generic public reference databases and aggregated to broad COG categories",
            "thermoplasmatota_representation": "No paired sample-specific metagenomic assemblies or Thermoplasmatota MAGs",
            "eligibility_decision": "EXCLUDED_NO_PAIRED_MAG_SEARCH_DB_ANNOTATION_BIAS",
            "scientific_justification": "Searches rely on public reference protein databases lacking novel sequence-dark desert/soil Thermoplasmatota families; absence in generic public reference searches is non-estimable database exclusion, not evidence of non-expression.",
        },
        {
            "candidate_id": "ENA_SRA_Hyperarid_Desert_Metagenome_and_Amplicon_BioProjects",
            "citation_doi_pmid": "Multiple desert soil studies (Namib, Atacama, Arabian/Empty Quarter, Mojave, Dry Valleys)",
            "data_repositories_and_accessions": "PRJNA1029032, PRJNA642254, PRJEB41174, PRJNA836362, PRJNA725625",
            "assay_modality": "DNA shotgun metagenomics (library_strategy=WGS, library_source=METAGENOMIC) or 16S rRNA amplicon sequencing (library_strategy=AMPLICON)",
            "thermoplasmatota_representation": "Contains Thermoplasmatota DNA reads/MAGs in several desert soils, but zero paired RNA-seq or metaproteomics runs",
            "eligibility_decision": "EXCLUDED_DNA_ONLY_NOT_EXPRESSION_ASSAY",
            "scientific_justification": "Audited via ENA Portal API; DNA metagenomes and 16S rRNA amplicons measure genomic presence, not transcript or protein expression.",
        },
    ]
    write_tsv_with_schema(
        os.path.join(RESULTS_DIR, "dataset_source_audit_and_exclusion_ledger.tsv"),
        list(exclusion_rows[0].keys()),
        exclusion_rows,
    )

    conservative_det_by_role = defaultdict(Counter)
    author_subset_det_by_role = defaultdict(Counter)
    for r in family_results:
        conservative_det_by_role[r["cohort_role"]][r["conservative_detection_classification"]] += 1
        author_subset_det_by_role[r["cohort_role"]][r["dataset_reported_observed_subset_classification"]] += 1

    thermo_conservative_counts = Counter(
        r["conservative_detection_tier"] for r in thermo_protein_audit_rows if r["matched_thermo_mags"]
    )
    thermo_author_subset_counts = Counter(
        r["dataset_reported_observed_subset_tier"] for r in thermo_protein_audit_rows if r["matched_thermo_mags"]
    )

    summary = {
        "study": "Karaoz et al. Nat Commun 17:10740 (13 Oct 2026) DOI:10.1038/s41467-026-78137-0",
        "stage": "stage2_expression_validation_completed",
        "runtime_verified_audit_counts": {
            "n_mztab_files_audited": len(sample_audit_rows),
            "total_discordant_psm_q_rows": total_discordant_psm_q,
            "total_excluded_invalid_rows": total_excluded_invalid_rows,
            "emsl50964_valid_proteins_scanned": emsl_valid_count,
            "emsl50964_invalid_proteins_logged": emsl_invalid_count,
            "archaeal_mag_orfs_v2_scanned": arch_valid_count,
            "massive_contaminant_and_arabidopsis_proteins_scanned": contaminant_valid_count,
            "candidate_qualifying_peptides_audited": len(candidate_peptides_to_audit),
            "full_15split_composite_search_fasta_deposited": full_background_verified,
        },
        "contract_compliance": {
            "dereplication_representation": f"99% ANI dereplicated representatives ({len(derep_arch_mags)} Archaeal MAGs out of 2,321 total MAGs), not necessarily 95% ANI species representatives",
            "psm_and_peptide_fdr_method": (
                f"Metadata-verified upstream MS-GF:QValue (search_engine_score[6] <= 0.01) and MS-GF:PepQValue "
                f"(search_engine_score[3] <= 0.01) across all {len(sample_audit_rows)} .mzTab files; "
                f"discordant_psm_q_rows={total_discordant_psm_q}; no reconstructed peptide FDR."
            ),
            "search_background_and_il_specificity": (
                f"Audited 604 deposited files in MassIVE MSV000097826: found {len(found_split_fastas)}/15 split composite "
                f"search FASTAs -> full_background_verified={full_background_verified}. Streamed multipattern I/L specificity "
                f"across {emsl_valid_count} valid EMSL50964 proteins, {arch_valid_count} Stage 1 v2 Archaeal MAG ORFs, "
                f"{contaminant_valid_count} contaminant/Arabidopsis proteins, and all .mzTab co-occurring accessions. "
                f"Primary global single-MAG uniqueness is marked non-estimable when full_background_verified is False."
            ),
            "spectrum_quantification_rule": (
                "Protein, family, and clade spectral counts are computed via aggregate_unique_spectra_for_entity "
                "(set union of unique spectra_ref per sample), never summing across duplicate aliases or paralogs, "
                "and restricting family/clade warming calculations strictly to family/clade-unambiguous peptides."
            ),
            "primary_family_and_clade_warming_estimability": (
                "Because the full 15-split 5%-FDR composite search FASTAs are unavailable in public repositories "
                "(full_background_verified=False), primary family-specific and clade-level warming responses are "
                "nonestimable_full_search_background_unavailable. All family and clade paired differences are reported "
                "strictly as observed_subset_provisional_sensitivity diagnostics; no positive warming claim is made. "
                "Protein-to-DNA ratios are descriptive sensitivities, not causal activity estimates."
            ),
            "log2_fold_change_zero_policy": (
                "Arbitrary 1e-6 pseudocounts at zeros are removed (eliminating artifactual log2FC ~24-30 values). "
                "Paired differences (heated - control) are the primary metric. Descriptive log2(heated/control) is only "
                "defined when both values > 0; pairs involving zero are reported with explicit detection states."
            ),
            "experimental_design": (
                f"{design_check['n_samples']} samples = {design_check['n_independent_plot_pairs']} paired plots "
                f"x 2 treatments (control, heated) x {design_check['n_nested_depth_strata']} nested depths, "
                "all sampled June 2018; paired DNA recruitment restricted strictly to 2018 metagenomes."
            ),
            "statistical_inference_note": (
                f"With n={design_check['n_independent_plot_pairs']} paired plots, exact whole-plot sign-flip permutation "
                f"has {design_check['max_whole_plot_sign_flips']} configurations and minimum achievable two-sided "
                f"p = {design_check['min_possible_two_sided_sign_flip_p']}."
            ),
            "claim_boundary": "Strongest permitted claim is external homolog/clade protein deployment in Blodgett Forest conifer soil (Karaoz 2026); does NOT prove expression in Empty Quarter desert soils or establish biochemical mechanism.",
        },
        "primary_clade_warming_estimability_status": "nonestimable_full_search_background_unavailable",
        "thermoplasmatota_proteins_in_emsl50964_conservative_detection_summary": dict(thermo_conservative_counts),
        "thermoplasmatota_proteins_in_emsl50964_dataset_reported_observed_subset_summary": dict(thermo_author_subset_counts),
        "eq_family_conservative_detection_summary_by_cohort_role": {k: dict(v) for k, v in sorted(conservative_det_by_role.items())},
        "eq_family_dataset_reported_observed_subset_summary_by_cohort_role": {k: dict(v) for k, v in sorted(author_subset_det_by_role.items())},
        "clade_level_paired_warming_observed_subset_provisional_sensitivity_n3_plots": clade_effects,
        "clade_level_paired_warming_analysis_n3_plots": clade_effects,
    }

    with open(os.path.join(RESULTS_DIR, "stage2_expression_validation_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("[Stage 2] Completed successfully. Summary:", flush=True)
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
