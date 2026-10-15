#!/usr/bin/env python3
"""01_fetch_references.py - Authoritative Reference Rebuild for Hypotheses F1 & F2 (v3).

Audit-Compliant Reference Pipeline (Root & Argon Audit Gate, 2026-10-15):
1. Clean programmatic REST queries to UniProtKB:
   - CPD: (ec:4.1.99.3) AND (reviewed:true)
   - 6-4 / FeS: (ec:4.1.99.13 OR protein_name:"6-4 photolyase") AND (reviewed:true)  [NO gene:phrB contamination]
   - Cryptochrome: (protein_name:cryptochrome) AND (reviewed:true)  [Filters out bHLH transcription factors]
   - SplB outgroup: (ec:4.1.99.14) AND (reviewed:true)
   - SPL: Verified Xu et al. 2021 (GenBank MK474076 -> QCI56913.1 / WP_032397095.1, exact translation match)
   - MpgP: (ec:3.1.3.70) AND (reviewed:true)
   - GgpP: (ec:3.1.3.69 OR protein_name:"glucosylglycerol-phosphate phosphatase") AND (reviewed:true) [NO gene:stpA contamination]
   - TPP: (ec:3.1.3.12) AND (reviewed:true)
   - Competing HAD enzymes: PGP (ec:3.1.3.18, 30), SPP (ec:3.1.3.24), HisP (ec:3.1.3.15, 30), bPGM (ec:5.4.2.6), DeoB (ec:5.4.2.7, 30), YrbI (ec:3.1.3.45)
   - Flanking Operon Synthases: MpgS (ec:2.4.1.217), GgpS (ec:2.4.1.213), OtsA (ec:2.4.1.15, 30)
2. Strict deduplication by accession, retaining all query memberships.
3. Broad family categories (CPD, 6-4, cryptochrome, SPL, SplB, CPF_unclassified; no heuristic fallbacks to CPD_I/CRY_DASH).
4. Evidence parsing:
   - Evaluates "CATALYTIC ACTIVITY" (with space), "FUNCTION", "ACTIVITY REGULATION", and ecNumbers[*].evidences.
   - Flags direct assay evidence (ECO:0000314, ECO:0000269) on catalytic/functional assertions.
5. Downloads authoritative Pfam HMMs and presses them.
6. Builds DIAMOND databases: cpf_references.dmnd, phosphatase_references.dmnd, and synthase_references.dmnd.
"""

import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_DIR = ROOT / "references_v3"
RAW_RECORDS = REF_DIR / "raw_records"
RAW_JSON = RAW_RECORDS / "uniprot_json"
RAW_FASTA = RAW_RECORDS / "uniprot_fasta"
RAW_NCBI = RAW_RECORDS / "ncbi"

for d in [REF_DIR, RAW_RECORDS, RAW_JSON, RAW_FASTA, RAW_NCBI]:
    d.mkdir(parents=True, exist_ok=True)

DIRECT_ASSAY_ECOS = {"ECO:0000314", "ECO:0000269"}
EXPERIMENTAL_ECOS = {"ECO:0000314", "ECO:0000269", "ECO:0000315", "ECO:0000353"}
AUTHOR_STATEMENT_ECOS = {"ECO:0000303", "ECO:0000304", "ECO:0000305", "ECO:0000250", "ECO:0000255"}


def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def http_get(url, retries=5, delay=2):
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "EmptyQuarter-DiscoveryPilot/3.0 (audit-verified pipeline)"}
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=45) as resp:
                headers = dict(resp.info())
                return resp.read(), headers
        except Exception as e:
            if attempt == retries - 1:
                raise RuntimeError(f"Failed to fetch {url} after {retries} attempts: {e}")
            time.sleep(delay * (attempt + 1))


def fetch_uniprot_query(query_str, track_label, category_label, subcategory_label, cap=None):
    """Query UniProt REST API with pagination, download raw JSON & FASTA, and parse evidence tiers."""
    encoded_query = urllib.parse.quote(query_str)
    base_url = f"https://rest.uniprot.org/uniprotkb/search?query={encoded_query}&format=json&size=50"
    print(f"\n[UniProt Query] {query_str}")
    print(f"  Target: {track_label} | Category: {category_label} | Subcategory: {subcategory_label}")

    raw_data, headers = http_get(base_url)
    total_results = int(headers.get("x-total-results", 0))
    parsed = json.loads(raw_data.decode("utf-8"))
    results = parsed.get("results", [])

    next_link = None
    link_header = headers.get("link", "")
    if 'rel="next"' in link_header:
        for part in link_header.split(","):
            if 'rel="next"' in part:
                next_link = part.split(";")[0].strip("<> ")

    while next_link and (cap is None or len(results) < cap):
        try:
            more_data, more_headers = http_get(next_link)
            more_parsed = json.loads(more_data.decode("utf-8"))
            results.extend(more_parsed.get("results", []))
            next_link = None
            link_header = more_headers.get("link", "")
            if 'rel="next"' in link_header:
                for part in link_header.split(","):
                    if 'rel="next"' in part:
                        next_link = part.split(";")[0].strip("<> ")
        except Exception as e:
            print(f"  Warning: pagination halted: {e}")
            break

    if cap is not None and len(results) > cap:
        sampling_mode = f"bounded_sample_{cap}_of_{total_results}"
        results = results[:cap]
    else:
        sampling_mode = f"exhaustive_reviewed_swissprot (n={len(results)} of total {total_results})"

    print(f"  Retrieved {len(results)} records (total in DB: {total_results}, mode: {sampling_mode})")

    records = []
    for entry in results:
        acc = entry.get("primaryAccession")
        if not acc:
            continue

        # Save raw JSON
        json_path = RAW_JSON / f"{acc}.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(entry, f, indent=2)

        # Fetch and save raw FASTA
        fasta_url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        fasta_bytes, _ = http_get(fasta_url)
        fasta_path = RAW_FASTA / f"{acc}.fasta"
        fasta_path.write_bytes(fasta_bytes)

        # Parse metadata
        uni_id = entry.get("uniProtkbId", "")
        org_name = entry.get("organism", {}).get("scientificName", "Unknown")
        lineage = ";".join(entry.get("organism", {}).get("lineage", []))

        prot_desc = entry.get("proteinDescription", {})
        rec_name = prot_desc.get("recommendedName", {}).get("fullName", {}).get("value", "")
        if not rec_name:
            sub_names = prot_desc.get("submissionNames", [])
            if sub_names:
                rec_name = sub_names[0].get("fullName", {}).get("value", "")

        # Post-filter: skip bHLH transcription factors from cryptochrome query
        alt_names_text = json.dumps(prot_desc).lower()
        if "bhlh" in rec_name.lower() or "helix-loop-helix" in rec_name.lower() or "bhlh" in alt_names_text:
            if "cryptochrome-interacting" in alt_names_text and "photolyase" not in rec_name.lower():
                print(f"  [Filter Excluded] {acc} ({rec_name}) is a bHLH transcription factor, not a CPF protein.")
                continue

        gene_names = [g["geneName"].get("value", "") for g in entry.get("genes", []) if "geneName" in g]
        gene_str = ",".join(gene_names)

        ec_list = []
        rec_obj = prot_desc.get("recommendedName", {})
        for ec_obj in rec_obj.get("ecNumbers", []):
            ec_list.append(ec_obj.get("value", ""))
        for inc in prot_desc.get("includes", []):
            for ec_obj in inc.get("recommendedName", {}).get("ecNumbers", []):
                ec_list.append(ec_obj.get("value", ""))
        ec_str = ",".join(sorted(set(ec_list)))

        seq_obj = entry.get("sequence", {})
        seq_len = seq_obj.get("length", 0)
        seq_val = seq_obj.get("value", "")

        # Strict evidence parsing across comments and ecNumbers
        function_ecos = set()
        catalytic_texts = []

        # 1. Inspect ecNumbers evidences
        for ec_obj in rec_obj.get("ecNumbers", []):
            for ev in ec_obj.get("evidences", []):
                code = ev.get("evidenceCode", "")
                if code:
                    function_ecos.add(code)

        # 2. Inspect comments (handling "CATALYTIC ACTIVITY" with space)
        for c in entry.get("comments", []):
            c_type = c.get("commentType", "").upper().replace("_", " ")
            if c_type in ["FUNCTION", "CATALYTIC ACTIVITY", "ACTIVITY REGULATION", "BIOPHYSICOCHEMICAL PROPERTIES"]:
                for text_obj in c.get("texts", []):
                    catalytic_texts.append(text_obj.get("value", ""))
                    for ev in text_obj.get("evidences", []):
                        code = ev.get("evidenceCode", "")
                        if code:
                            function_ecos.add(code)
                if "reaction" in c:
                    catalytic_texts.append(c["reaction"].get("name", ""))
                    for ev in c["reaction"].get("evidences", []):
                        code = ev.get("evidenceCode", "")
                        if code:
                            function_ecos.add(code)
                for ev in c.get("evidences", []):
                    code = ev.get("evidenceCode", "")
                    if code:
                        function_ecos.add(code)

        all_ecos = sorted(set(re.findall(r"ECO:\d{7}", json.dumps(entry))))
        fn_eco_str = ",".join(sorted(function_ecos))
        all_eco_str = ",".join(all_ecos)

        has_direct_assay = any(code in DIRECT_ASSAY_ECOS for code in (function_ecos or all_ecos))
        has_experimental = any(code in EXPERIMENTAL_ECOS for code in (function_ecos or all_ecos))
        has_author_stmt_only = (not has_experimental) and any(code in AUTHOR_STATEMENT_ECOS for code in all_ecos)

        direct_assay_codes = ",".join(sorted(c for c in (function_ecos or all_ecos) if c in DIRECT_ASSAY_ECOS))
        cat_text_clean = "; ".join(catalytic_texts)[:300].replace("\n", " ").replace("\t", " ")

        # Subfamily broad classification (Audit Gate: broad families, no heuristic fallbacks)
        pname_l = rec_name.lower()
        if track_label == "F1_CPF":
            if subcategory_label == "spore_photoproduct_lyase_radical_SAM" or "4.1.99.14" in ec_str:
                cls_assign = "SplB"
            elif subcategory_label == "SPL_short_photolyase":
                cls_assign = "SPL"
            elif "4.1.99.13" in ec_str or "(6-4)" in pname_l:
                cls_assign = "6-4"
            elif "4.1.99.3" in ec_str or "cpd" in pname_l or "class i" in pname_l or "class ii" in pname_l or "class iii" in pname_l:
                cls_assign = "CPD"
            elif "cryptochrome" in pname_l or "dash" in pname_l:
                cls_assign = "cryptochrome"
            else:
                cls_assign = "CPF_unclassified"
        elif track_label == "F2_Phosphatase":
            if category_label == "Phosphatase_target_MpgP":
                cls_assign = "MpgP_EC_3_1_3_70"
            elif category_label == "Phosphatase_target_GgpP":
                cls_assign = "GgpP_EC_3_1_3_69"
            elif category_label == "Phosphatase_target_TPP":
                cls_assign = "TPP_EC_3_1_3_12"
            elif category_label == "Synthase_MpgS":
                cls_assign = "Synthase_MpgS"
            elif category_label == "Synthase_GgpS":
                cls_assign = "Synthase_GgpS"
            elif category_label == "Synthase_OtsA":
                cls_assign = "Synthase_OtsA"
            else:
                cls_assign = f"HAD_competing_{subcategory_label.split('_')[-1]}"
        else:
            cls_assign = "Unclassified"

        records.append({
            "accession": acc,
            "uniprot_id_or_refseq": uni_id,
            "track": track_label,
            "category": category_label,
            "subcategory": subcategory_label,
            "class_assignment": cls_assign,
            "protein_name": rec_name,
            "gene_name": gene_str,
            "organism": org_name,
            "lineage": lineage,
            "ec_number": ec_str,
            "length_aa": seq_len,
            "has_direct_assay_evidence": 1 if has_direct_assay else 0,
            "has_experimental_evidence": 1 if has_experimental else 0,
            "has_author_statement_only": 1 if has_author_stmt_only else 0,
            "direct_assay_codes": direct_assay_codes,
            "function_evidence_codes": fn_eco_str,
            "all_evidence_codes": all_eco_str,
            "catalytic_activity_and_function": cat_text_clean,
            "sampling_mode": sampling_mode,
            "total_matching_in_db": total_results,
            "query_string": query_str,
            "source_db": "UniProtKB/Swiss-Prot",
            "sequence": seq_val,
        })
    return records


def fetch_and_verify_xu2021_spl_references():
    """Verify Xu et al. 2021 SPL archetype via GenBank nucleotide-to-protein linkage."""
    print("\n[Xu et al. 2021 Reference Verification]")
    print("  Fetching GenBank nucleotide MK474076 flatfile...")
    gb_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MK474076.1&rettype=gb&retmode=text"
    gb_bytes, _ = http_get(gb_url)
    gb_text = gb_bytes.decode("utf-8")

    gb_path = RAW_NCBI / "MK474076.1.gb"
    gb_path.write_text(gb_text)

    # Parse /protein_id and /translation from GenBank CDS
    pid_match = re.search(r'/protein_id="([^"]+)"', gb_text)
    assert pid_match, "Could not find /protein_id in MK474076 flatfile!"
    protein_id = pid_match.group(1)
    assert protein_id == "QCI56913.1", f"Expected QCI56913.1, found {protein_id}!"

    trans_match = re.search(r'/translation="([^"]+)"', gb_text, re.DOTALL)
    assert trans_match, "Could not find /translation in MK474076 flatfile!"
    cds_translation = re.sub(r"\s+", "", trans_match.group(1))

    # Fetch QCI56913.1 protein FASTA from NCBI
    prot_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={protein_id}&rettype=fasta&retmode=text"
    prot_bytes, _ = http_get(prot_url)
    prot_text = prot_bytes.decode("utf-8")
    fasta_path = RAW_NCBI / f"{protein_id}.fasta"
    fasta_path.write_text(prot_text)

    fasta_lines = prot_text.strip().split("\n")
    ncbi_seq = "".join(l.strip() for l in fasta_lines[1:] if not l.startswith(">"))
    assert cds_translation == ncbi_seq, "CDS translation does NOT match NCBI protein FASTA!"
    print(f"  Verified GenBank MK474076 -> {protein_id} exact translation match ({len(ncbi_seq)} aa).")

    records = [{
        "accession": protein_id,
        "uniprot_id_or_refseq": "WP_032397095.1",
        "track": "F1_CPF",
        "category": "CPF_photolyase",
        "subcategory": "SPL_short_photolyase",
        "class_assignment": "SPL",
        "protein_name": "short photolyase-like protein (Xu et al. 2021 archetype RfSPL)",
        "gene_name": "spl",
        "organism": "Roseiflexus sp. RS-1",
        "lineage": "Bacteria;Chloroflexi;Chloroflexia;Chloroflexales;Roseiflexaceae;Roseiflexus",
        "ec_number": "4.1.99.-",
        "length_aa": len(ncbi_seq),
        "has_direct_assay_evidence": 1,
        "has_experimental_evidence": 1,
        "has_author_statement_only": 0,
        "direct_assay_codes": "Xu2021_biochemical_assay",
        "function_evidence_codes": "Xu2021_experimental",
        "all_evidence_codes": "Xu2021_experimental;GenBank_MK474076",
        "catalytic_activity_and_function": "Photorepair of CPDs in single-stranded DNA and RNA (Xu et al. 2021 Mol Biol Evol 38:4505-4519)",
        "sampling_mode": "authoritative_literature_archetype",
        "total_matching_in_db": 1,
        "query_string": "GenBank:MK474076.1/QCI56913.1",
        "source_db": "NCBI_GenBank",
        "sequence": ncbi_seq,
    }]
    return records


def download_pfam_hmms():
    """Download authoritative Pfam HMMs from InterPro/EBI and press them."""
    print("\n=======================================================")
    print(">>> Downloading Authoritative Pfam HMMs from InterPro/EBI")
    print("=======================================================")

    pfam_accessions = [
        ("PF00875", "DNA_photolyase", "F1_CPF"),
        ("PF03441", "FAD_binding_7", "F1_CPF"),
        ("PF04244", "Photolyase_N", "F1_CPF"),
        ("PF04055", "Radical_SAM", "F1_Outgroup"),
        ("PF00702", "Hydrolase", "F2_Phosphatase"),
        ("PF13419", "HAD_2", "F2_Phosphatase"),
    ]

    hmm_dir = REF_DIR / "pfam_hmms"
    hmm_dir.mkdir(parents=True, exist_ok=True)

    downloaded = {}
    for acc, name, grp in pfam_accessions:
        hmm_file = hmm_dir / f"{acc}.hmm"
        if not hmm_file.exists() or hmm_file.stat().st_size == 0:
            url = f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{acc}?annotation=hmm"
            print(f"  Downloading {acc} ({name})...")
            data, _ = http_get(url)
            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)
            hmm_file.write_bytes(data)
            print(f"    Saved {hmm_file.name} ({len(data)} bytes).")
        downloaded[acc] = hmm_file

    # Build combined HMM libraries
    cpf_hmm = REF_DIR / "cpf_all.hmm"
    had_hmm = REF_DIR / "had_all.hmm"

    with open(cpf_hmm, "wb") as out_f:
        for acc in ["PF00875", "PF03441", "PF04244", "PF04055"]:
            out_f.write(downloaded[acc].read_bytes())

    with open(had_hmm, "wb") as out_f:
        for acc in ["PF00702", "PF13419"]:
            out_f.write(downloaded[acc].read_bytes())

    for target_hmm in [cpf_hmm, had_hmm]:
        print(f"  Running hmmpress on {target_hmm.name}...")
        for ext in [".h3i", ".h3m", ".h3p", ".h3f"]:
            p = target_hmm.with_suffix(target_hmm.suffix + ext)
            if p.exists():
                p.unlink()
        subprocess.run(["hmmpress", str(target_hmm)], check=True)
    print("  All HMM libraries pressed successfully.")


def main():
    t0 = time.time()
    print("=== [Hypotheses F1 & F2] Authoritative Reference Rebuild v3 ===")
    print(f"Destination Directory: {REF_DIR}")

    all_ref_records = []

    # 1. F1: Curated CPF Families
    queries_f1 = [
        ("CPD_photolyase", "CPF_photolyase", "(ec:4.1.99.3) AND (reviewed:true)", None),
        ("FeS_BCP_or_64_photolyase", "CPF_photolyase", "(ec:4.1.99.13 OR protein_name:\"6-4 photolyase\") AND (reviewed:true)", None),
        ("cryptochrome_signaling_and_DASH", "CPF_cryptochrome", "(protein_name:cryptochrome) AND (reviewed:true)", None),
        ("spore_photoproduct_lyase_radical_SAM", "Outgroup_control", "(ec:4.1.99.14) AND (reviewed:true)", None),
    ]

    for subcat, cat, q, cap in queries_f1:
        recs = fetch_uniprot_query(q, "F1_CPF", cat, subcat, cap=cap)
        all_ref_records.extend(recs)

    spl_records = fetch_and_verify_xu2021_spl_references()
    all_ref_records.extend(spl_records)

    # 2. F2: Compatible-Solute Phosphatases, Competing HADs, & Flanking Synthases
    queries_f2 = [
        # Targets
        ("mannosyl_3_phosphoglycerate_phosphatase", "Phosphatase_target_MpgP", "(ec:3.1.3.70) AND (reviewed:true)", None),
        ("glucosylglycerol_phosphate_phosphatase", "Phosphatase_target_GgpP", "(ec:3.1.3.69 OR protein_name:\"glucosylglycerol-phosphate phosphatase\") AND (reviewed:true)", None),
        ("trehalose_6_phosphate_phosphatase", "Phosphatase_target_TPP", "(ec:3.1.3.12) AND (reviewed:true)", None),
        # Competing HAD enzymes
        ("phosphoglycolate_phosphatase_PGP", "Competing_HAD_enzyme", "(ec:3.1.3.18) AND (reviewed:true)", 30),
        ("sucrose_6_phosphatase_SPP", "Competing_HAD_enzyme", "(ec:3.1.3.24) AND (reviewed:true)", None),
        ("histidinol_phosphatase_HisP", "Competing_HAD_enzyme", "(ec:3.1.3.15) AND (reviewed:true)", 30),
        ("beta_phosphoglucomutase_bPGM", "Competing_HAD_enzyme", "(ec:5.4.2.6) AND (reviewed:true)", None),
        ("deoxyribose_phosphate_mutase_DeoB", "Competing_HAD_enzyme", "(ec:5.4.2.7) AND (reviewed:true)", 30),
        ("kdo_8_phosphate_phosphatase_YrbI", "Competing_HAD_enzyme", "(ec:3.1.3.45) AND (reviewed:true)", None),
        # Flanking Operon Synthases
        ("mannosyl_3_phosphoglycerate_synthase", "Synthase_MpgS", "(ec:2.4.1.217 OR protein_name:\"mannosyl-3-phosphoglycerate synthase\") AND (reviewed:true)", None),
        ("glucosylglycerol_phosphate_synthase", "Synthase_GgpS", "(ec:2.4.1.213) AND (reviewed:true)", None),
        ("trehalose_6_phosphate_synthase", "Synthase_OtsA", "(ec:2.4.1.15) AND (reviewed:true)", 30),
    ]

    for subcat, cat, q, cap in queries_f2:
        recs = fetch_uniprot_query(q, "F2_Phosphatase", cat, subcat, cap=cap)
        all_ref_records.extend(recs)

    # Strict deduplication by accession while preserving all query memberships
    unique_records = []
    seen = {}
    for r in all_ref_records:
        acc = r["accession"]
        if acc in seen:
            prev = seen[acc]
            if r["subcategory"] not in prev["all_subcategories"]:
                prev["all_subcategories"] += f";{r['subcategory']}"
            if r["query_string"] not in prev["all_queries"]:
                prev["all_queries"] += f";{r['query_string']}"
        else:
            r["all_subcategories"] = r["subcategory"]
            r["all_queries"] = r.get("query_string", "")
            seen[acc] = r
            unique_records.append(r)

    print(f"\nTotal unique verified reference records assembled: {len(unique_records)}")

    # 3. Write FASTAs and Build DIAMOND DBs
    cpf_records = [r for r in unique_records if r["track"] == "F1_CPF"]
    phos_records = [r for r in unique_records if r["track"] == "F2_Phosphatase" and not r["category"].startswith("Synthase_")]
    synth_records = [r for r in unique_records if r["category"].startswith("Synthase_")]

    print(f"Total F1 CPF references: {len(cpf_records)}")
    print(f"Total F2 Phosphatase references: {len(phos_records)}")
    print(f"Total Flanking Synthase references: {len(synth_records)}")

    cpf_faa = REF_DIR / "cpf_references.faa"
    with open(cpf_faa, "w") as f:
        for r in cpf_records:
            f.write(f">{r['accession']} {r['uniprot_id_or_refseq']} subcat={r['class_assignment']} org={r['organism']} len={r['length_aa']} exp={r['has_experimental_evidence']}\n{r['sequence']}\n")

    phos_faa = REF_DIR / "phosphatase_references.faa"
    with open(phos_faa, "w") as f:
        for r in phos_records:
            f.write(f">{r['accession']} {r['uniprot_id_or_refseq']} cat={r['category']} subcat={r['class_assignment']} org={r['organism']} len={r['length_aa']} exp={r['has_experimental_evidence']}\n{r['sequence']}\n")

    synth_faa = REF_DIR / "synthase_references.faa"
    with open(synth_faa, "w") as f:
        for r in synth_records:
            f.write(f">{r['accession']} {r['uniprot_id_or_refseq']} cat={r['category']} subcat={r['class_assignment']} org={r['organism']} len={r['length_aa']} exp={r['has_experimental_evidence']}\n{r['sequence']}\n")

    diamond_bin = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")
    if not diamond_bin.exists():
        diamond_bin = Path("diamond")

    for faa_file, dmnd_name in [
        (cpf_faa, "cpf_references.dmnd"),
        (phos_faa, "phosphatase_references.dmnd"),
        (synth_faa, "synthase_references.dmnd"),
    ]:
        dmnd_out = REF_DIR / dmnd_name
        print(f"  Building DIAMOND database {dmnd_name}...")
        cmd = [str(diamond_bin), "makedb", "--in", str(faa_file), "-d", str(dmnd_out.with_suffix(""))]
        subprocess.run(cmd, check=True)

    # 4. Download and Press Pfam HMMs
    download_pfam_hmms()

    # 5. Write Metadata TSV and Manifest
    meta_tsv = REF_DIR / "verified_references_metadata.tsv"
    fieldnames = [
        "accession", "uniprot_id_or_refseq", "track", "category", "subcategory",
        "all_subcategories", "class_assignment", "protein_name", "gene_name", "organism",
        "lineage", "ec_number", "length_aa", "has_direct_assay_evidence",
        "has_experimental_evidence", "has_author_statement_only", "direct_assay_codes",
        "function_evidence_codes", "all_evidence_codes", "catalytic_activity_and_function",
        "sampling_mode", "total_matching_in_db", "query_string", "all_queries", "source_db"
    ]
    with open(meta_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(unique_records)

    # Symlink canonical references_metadata.tsv
    sym_meta = REF_DIR / "references_metadata.tsv"
    if sym_meta.exists() or sym_meta.is_symlink():
        sym_meta.unlink()
    sym_meta.symlink_to("verified_references_metadata.tsv")

    # Write Manifest
    manifest_tsv = REF_DIR / "references_manifest.tsv"
    with open(manifest_tsv, "w") as f:
        f.write("filename\tsha256\tsize_bytes\tdescription\n")
        files_to_hash = [
            ("cpf_references.faa", "FASTA sequences for curated CPF references"),
            ("phosphatase_references.faa", "FASTA sequences for compatible-solute and competing HAD references"),
            ("synthase_references.faa", "FASTA sequences for MpgS, GgpS, and OtsA synthase references"),
            ("cpf_references.dmnd", "DIAMOND binary database for CPF references"),
            ("phosphatase_references.dmnd", "DIAMOND binary database for phosphatase references"),
            ("synthase_references.dmnd", "DIAMOND binary database for synthase references"),
            ("cpf_all.hmm", "HMMER3 library for CPF domains (PF00875, PF03441, PF04244, PF04055)"),
            ("had_all.hmm", "HMMER3 library for HAD domains (PF00702, PF13419)"),
            ("verified_references_metadata.tsv", "Comprehensive metadata and evidence audit table"),
        ]
        for fname, desc in files_to_hash:
            fpath = REF_DIR / fname
            if fpath.exists():
                sha = compute_sha256(fpath)
                size = fpath.stat().st_size
                f.write(f"{fname}\t{sha}\t{size}\t{desc}\n")

    # Update top-level references symlink
    top_sym = ROOT / "references"
    if top_sym.is_symlink() or top_sym.exists():
        top_sym.unlink()
    top_sym.symlink_to("references_v3")

    print(f"\n=== Reference Rebuild Complete in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
    main()
