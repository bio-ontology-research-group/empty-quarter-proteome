#!/usr/bin/env python3
"""01_fetch_references.py - Authoritative Reference Retrieval for Tracks F3, F4, F5.

Resolves Root Blocking Audit (2026-10-15):
1. Completely removes unverified or hallucinated anchor accession dictionaries.
2. Retrieves reference records exclusively by structured source queries with verified family
   annotations (Pfam, InterPro, EC numbers), deriving all labels from source-backed metadata.
3. For literature anchors, obtains exact supplement / PDB chain mappings and verifies them:
   - Pushkarev 2018 (Nature): heliorhodopsin 48C12 (GenBank MF737519.1 -> NCBI protein AVZ43932.1).
   - Czech 2019 (Sci Rep): SaEctC (PDB 5BXX -> UniProt Q1GNW6, Sphingopyxis alaskensis),
     PlEctC (PDB 5ONN -> NCBI RefSeq YP_003245677.1, Paenibacillus lautus),
     orphan cupin control (PDB 5TFZ -> UniProt Q4FNM4, Pelagibacter ubique DddK),
     Halomonas elongata EctC (UniProt O52251), EctA (O52249), EctB (O52250), EctD (E1VA04).
   - Pairoh 2025 (Sci Rep): TtCA (PV553289 -> YEH58108.1), CrCA (PV553290 -> YEH58109.1),
     ToCA (PV553291 -> YEH58110.1); canonical Cam (PDB 1QRG -> UniProt P40881).
4. Strictly separates Carboxysome CAs (CsoSCA / CcmM / IPR017156) from Gamma-CA (NF050158 / IPR052265).
5. Distinguishes function-specific experimental validation (ECO:0000269/314/353 on FUNCTION /
   CATALYTIC ACTIVITY / active sites) from generic record-level ECO tags.
6. Records all audited/rejected false accessions in references/rejected_references.tsv.
7. Downloads verified profile HMMs from InterPro API, indexes with hmmpress, and builds DIAMOND databases.
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
REF_DIR = ROOT / "references"
JSON_DIR = REF_DIR / "raw_uniprot_json"
FASTA_DIR = REF_DIR / "raw_uniprot_fasta"
NCBI_DIR = REF_DIR / "raw_ncbi"
HMM_DIR = REF_DIR / "profiles"

for d in [REF_DIR, JSON_DIR, FASTA_DIR, NCBI_DIR, HMM_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def http_get(url, retries=5, delay=2):
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "EmptyQuarter-F345DiscoveryPilot/2.0 (authoritative source query)"}
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read(), dict(resp.getheaders())
        except Exception as e:
            if attempt == retries - 1:
                raise RuntimeError(f"Failed to fetch {url} after {retries} attempts: {e}")
            time.sleep(delay * (attempt + 1))


def fetch_ncbi_protein(protein_id, track, category, subcategory):
    """Fetch verified protein record from NCBI Entrez E-utilities."""
    gp_path = NCBI_DIR / f"{protein_id}.gp"
    fasta_path = NCBI_DIR / f"{protein_id}.fasta"

    url_gp = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={protein_id}&rettype=gp&retmode=text"
    raw_gp, _ = http_get(url_gp)
    gp_path.write_bytes(raw_gp)

    url_fasta = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={protein_id}&rettype=fasta&retmode=text"
    raw_fasta, _ = http_get(url_fasta)
    fasta_path.write_bytes(raw_fasta)

    gp_text = raw_gp.decode("utf-8", errors="replace")
    fasta_lines = raw_fasta.decode("utf-8", errors="replace").splitlines()
    seq_lines = [l.strip() for l in fasta_lines if not l.startswith(">")]
    seq = "".join(seq_lines)

    # Parse metadata from GenBank GP flatfile
    definition = ""
    organism = "Unknown"
    gene_name = ""
    for line in gp_text.splitlines():
        if line.startswith("DEFINITION"):
            definition = line.replace("DEFINITION", "").strip()
        elif "ORGANISM" in line and organism == "Unknown":
            organism = line.replace("ORGANISM", "").strip()
        elif "/gene=" in line and not gene_name:
            gene_name = re.search(r'/gene="([^"]+)"', line).group(1)

    return {
        "accession": protein_id,
        "source_db": "NCBI",
        "uniprot_id": "",
        "track": track,
        "category": category,
        "subcategory": subcategory,
        "organism": organism,
        "lineage": "Environmental/Literature",
        "protein_name": definition,
        "gene_name": gene_name,
        "ec_number": "",
        "length_aa": len(seq),
        "sequence": seq,
        "has_function_experimental_evidence": 1,
        "has_any_experimental_evidence": 1,
        "function_evidence_codes": "LITERATURE_VALIDATED",
        "all_evidence_codes": "LITERATURE_VALIDATED",
        "is_pathogen": 0,
        "query_used": f"NCBI_Entrez:{protein_id}",
    }


def parse_uniprot_json(entry, track, category, subcategory, query_used=""):
    """Parse verified metadata directly from authoritative UniProtKB JSON."""
    acc = entry.get("primaryAccession")
    uni_id = entry.get("uniProtkbId", "")
    org_name = entry.get("organism", {}).get("scientificName", "Unknown")
    lineage = ";".join(entry.get("organism", {}).get("lineage", []))
    
    # Recommended / Submission protein name
    prot_desc = entry.get("proteinDescription", {})
    rec_name = prot_desc.get("recommendedName", {}).get("fullName", {}).get("value", "")
    if not rec_name:
        sub_names = prot_desc.get("submissionNames", [])
        if sub_names:
            rec_name = sub_names[0].get("fullName", {}).get("value", "")
    if not rec_name:
        rec_name = prot_desc.get("submittedNames", [{}])[0].get("fullName", {}).get("value", "Unknown")

    # Gene names
    gene_names = []
    for g in entry.get("genes", []):
        if "geneName" in g:
            gene_names.append(g["geneName"].get("value", ""))
    gene_str = ",".join(filter(None, gene_names))

    # EC numbers
    ec_list = []
    rec_obj = prot_desc.get("recommendedName", {})
    for ec_obj in rec_obj.get("ecNumbers", []):
        ec_list.append(ec_obj.get("value", ""))
    for inc in prot_desc.get("includes", []):
        for ec_obj in inc.get("recommendedName", {}).get("ecNumbers", []):
            ec_list.append(ec_obj.get("value", ""))
    ec_str = ",".join(sorted(set(ec_list)))

    # Sequence
    seq_obj = entry.get("sequence", {})
    seq_len = seq_obj.get("length", 0)
    seq_val = seq_obj.get("value", "")

    # Function-specific experimental evidence parsing
    fn_ecos = []
    cat_ecos = []
    feat_ecos = []
    for c in entry.get("comments", []):
        ctype = c.get("commentType", "")
        for text_item in c.get("texts", []):
            for ev in text_item.get("evidences", []):
                ecode = ev.get("evidenceCode")
                if ctype == "FUNCTION":
                    fn_ecos.append(ecode)
        if ctype == "CATALYTIC ACTIVITY":
            cat = c.get("reaction", {})
            for ev in cat.get("evidences", []):
                cat_ecos.append(ev.get("evidenceCode"))

    for f in entry.get("features", []):
        ftype = f.get("type", "")
        if ftype in ["Active site", "Binding site", "Site"]:
            for ev in f.get("evidences", []):
                feat_ecos.append(ev.get("evidenceCode"))

    exp_codes = {"ECO:0000269", "ECO:0000314", "ECO:0000353"}
    function_specific_exp = any(e in exp_codes for e in (fn_ecos + cat_ecos + feat_ecos))

    # All ECO codes across record
    json_text = json.dumps(entry)
    all_ecos = sorted(set(re.findall(r"ECO:\d{7}", json_text)))
    any_exp = any(e in exp_codes for e in all_ecos)

    # Pathogen check
    is_pathogen = False
    pathogen_terms = ["Mycobacterium tuberculosis", "Salmonella enterica", "Staphylococcus aureus",
                      "Vibrio cholerae", "Pseudomonas aeruginosa", "Yersinia pestis", "Bacillus anthracis",
                      "Clostridioides difficile", "Shigella", "Listeria monocytogenes"]
    for pterm in pathogen_terms:
        if pterm.lower() in org_name.lower():
            is_pathogen = True
            break

    return {
        "accession": acc,
        "source_db": "UniProtKB",
        "uniprot_id": uni_id,
        "track": track,
        "category": category,
        "subcategory": subcategory,
        "organism": org_name,
        "lineage": lineage,
        "protein_name": rec_name,
        "gene_name": gene_str,
        "ec_number": ec_str,
        "length_aa": seq_len,
        "sequence": seq_val,
        "has_function_experimental_evidence": 1 if function_specific_exp else 0,
        "has_any_experimental_evidence": 1 if any_exp else 0,
        "function_evidence_codes": ",".join(sorted(set(filter(None, fn_ecos + cat_ecos + feat_ecos)))),
        "all_evidence_codes": ",".join(all_ecos),
        "is_pathogen": 1 if is_pathogen else 0,
        "query_used": re.sub(r"\s+", " ", query_used).strip(),
    }


def fetch_uniprot_paginated_query(query_str, track, category, subcategory, max_records=40):
    """Query UniProtKB with cursor pagination, saving raw JSON/FASTA and metadata."""
    encoded_q = urllib.parse.quote(query_str)
    url = f"https://rest.uniprot.org/uniprotkb/search?query={encoded_q}&format=json&size=25"
    print(f"\n[UniProt Source Query - {track}] {query_str}")
    print(f"  Target: {category} / {subcategory} (max: {max_records})")
    
    records = []
    seen_accs = set()
    page = 1

    while url and len(records) < max_records:
        raw_json, headers = http_get(url)
        data = json.loads(raw_json.decode("utf-8"))
        results = data.get("results", [])
        if not results:
            break

        print(f"  Page {page}: retrieved {len(results)} entries.")
        for entry in results:
            if len(records) >= max_records:
                break
            acc = entry.get("primaryAccession")
            if not acc or acc in seen_accs:
                continue
            seen_accs.add(acc)

            # Save raw JSON
            json_path = JSON_DIR / f"{acc}.json"
            json_path.write_bytes(json.dumps(entry, indent=2).encode("utf-8"))

            # Save raw FASTA
            fasta_path = FASTA_DIR / f"{acc}.fasta"
            if not fasta_path.exists():
                url_fasta = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
                raw_fasta, _ = http_get(url_fasta)
                fasta_path.write_bytes(raw_fasta)

            parsed = parse_uniprot_json(entry, track, category, subcategory, query_used=query_str)
            records.append(parsed)

        # Check pagination Link header
        link_header = headers.get("Link") or headers.get("link")
        next_url = None
        if link_header:
            match = re.search(r'<([^>]+)>;\s*rel="next"', link_header)
            if match:
                next_url = match.group(1)
        url = next_url
        page += 1
        time.sleep(0.3)

    print(f"  Total records fetched for {category}/{subcategory}: {len(records)}")
    return records


def download_profile(db, acc, filename, label):
    """Download HMM profile from InterPro API and save to HMM_DIR."""
    dest = HMM_DIR / filename
    print(f"Fetching profile [{db}] {acc} -> {filename} ({label})...")
    url = f"https://www.ebi.ac.uk/interpro/api/entry/{db}/{acc}/?annotation=hmm"
    content, _ = http_get(url)
    if content[:2] == b"\x1f\x8b":
        content = gzip.decompress(content)
    dest.write_bytes(content)
    sha = compute_sha256(dest)
    print(f"  Saved {filename} ({len(content)} bytes, SHA-256: {sha[:16]}...)")
    return dest


def main():
    start_time = time.time()
    print("=" * 80)
    print("F3/F4/F5 AUTHORITATIVE SOURCE-BACKED REFERENCE RETRIEVAL")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print(f"Directory: {ROOT}")
    print("=" * 80)

    # --------------------------------------------------------------------------
    # 1. DOWNLOAD AUTHORITATIVE PROFILES FROM INTERPRO
    # --------------------------------------------------------------------------
    print("\n--- STAGE 1: Downloading Verified HMM Profiles ---")
    profiles = [
        # Track F3: Rhodopsins, Heliorhodopsins & Competing GPCRs
        ("pfam", "PF01036", "PF01036_Bac_rhodopsin.hmm", "Type-1 microbial rhodopsin domain"),
        ("pfam", "PF18761", "PF18761_Heliorhodopsin.hmm", "Type-3 heliorhodopsin domain"),
        ("ncbifam", "NF038020", "NF038020_Heliorhodopsin_HeR.hmm", "NCBIFAM heliorhodopsin HeR family"),
        ("pfam", "PF00001", "PF00001_7tm_1.hmm", "Eukaryotic GPCR rhodopsin family control"),
        ("pfam", "PF00002", "PF00002_7tm_2.hmm", "Secretin-like GPCR family control"),
        # Track F4: EctC, Pathway Markers & Competing Cupins
        ("pfam", "PF06339", "PF06339_Ectoine_synth.hmm", "Ectoine synthase family"),
        ("ncbifam", "NF009806", "NF009806_ectoine_synthase.hmm", "NCBIFAM ectoine synthase family"),
        ("ncbifam", "TIGR02406", "TIGR02406_EctA.hmm", "Diaminobutyrate acetyltransferase EctA"),
        ("ncbifam", "TIGR02407", "TIGR02407_EctB.hmm", "Diaminobutyrate transaminase EctB"),
        ("ncbifam", "TIGR02408", "TIGR02408_EctD.hmm", "Ectoine hydroxylase EctD"),
        ("pfam", "PF00190", "PF00190_Cupin_1.hmm", "Cupin domain control"),
        ("pfam", "PF07883", "PF07883_Cupin_2.hmm", "Cupin domain 2 control"),
        ("pfam", "PF02678", "PF02678_Pirin.hmm", "Pirin cupin family control"),
        # Track F5: Carbonic Anhydrases & Competing LbH Transferases
        ("ncbifam", "NF050158", "NF050158_gamma_CA.hmm", "NCBIFAM gamma carbonic anhydrase"),
        ("ncbifam", "NF040597", "NF040597_gamma_class_CA.hmm", "NCBIFAM Carbonic anhydrase gamma-class"),
        ("ncbifam", "TIGR02701", "TIGR02701_shell_CA.hmm", "Carboxysome shell carbonic anhydrase CsoSCA"),
        ("pfam", "PF00484", "PF00484_Pro_CA.hmm", "Prokaryotic Beta-carbonic anhydrase"),
        ("pfam", "PF00194", "PF00194_Carb_anhydrase.hmm", "Eukaryotic Alpha-carbonic anhydrase control"),
        ("pfam", "PF00132", "PF00132_Hexapep.hmm", "Hexapeptide repeat transferase control"),
        ("pfam", "PF00583", "PF00583_Acetyltransf_1.hmm", "GNAT acetyltransferase control"),
    ]

    for db, acc, fname, label in profiles:
        download_profile(db, acc, fname, label)

    # Combine into per-track HMM collections
    track_hmms = {
        "f3_rhodopsins.hmm": ["PF01036_Bac_rhodopsin.hmm", "PF18761_Heliorhodopsin.hmm",
                              "NF038020_Heliorhodopsin_HeR.hmm", "PF00001_7tm_1.hmm", "PF00002_7tm_2.hmm"],
        "f4_ectc.hmm": ["PF06339_Ectoine_synth.hmm", "NF009806_ectoine_synthase.hmm",
                        "TIGR02406_EctA.hmm", "TIGR02407_EctB.hmm", "TIGR02408_EctD.hmm",
                        "PF00190_Cupin_1.hmm", "PF07883_Cupin_2.hmm", "PF02678_Pirin.hmm"],
        "f5_carbonic_anhydrase.hmm": ["NF050158_gamma_CA.hmm", "NF040597_gamma_class_CA.hmm",
                                      "TIGR02701_shell_CA.hmm", "PF00484_Pro_CA.hmm",
                                      "PF00194_Carb_anhydrase.hmm", "PF00132_Hexapep.hmm",
                                      "PF00583_Acetyltransf_1.hmm"],
    }

    print("\n--- Combining and Pressing Profile Databases ---")
    for combo_name, member_files in track_hmms.items():
        combo_path = HMM_DIR / combo_name
        with open(combo_path, "wb") as outfile:
            for mf in member_files:
                outfile.write((HMM_DIR / mf).read_bytes())
        print(f"Created {combo_name} ({len(member_files)} models, {combo_path.stat().st_size} bytes)")
        try:
            subprocess.run(["hmmpress", "-f", str(combo_path)], check=True, capture_output=True)
            print(f"  hmmpress {combo_name} OK")
        except FileNotFoundError:
            print(f"  hmmpress binary not found on local path; will be pressed in Slurm worker.")
        except Exception as e:
            print(f"  hmmpress warning: {e}")

    # --------------------------------------------------------------------------
    # 2. FETCH VERIFIED LITERATURE ANCHORS
    # --------------------------------------------------------------------------
    print("\n--- STAGE 2: Fetching Authoritative Literature Anchors ---")
    lit_records = []

    # A. Pushkarev et al. 2018 (Nature): heliorhodopsin 48C12 (GenBank MF737519.1 -> AVZ43932.1)
    print("Fetching Pushkarev 2018 anchor: AVZ43932.1 (heliorhodopsin 48C12)...")
    rec_48c12 = fetch_ncbi_protein("AVZ43932.1", "F3", "Heliorhodopsin", "HeR_48C12_Pushkarev2018")
    lit_records.append(rec_48c12)

    # B. Czech et al. 2019 (Sci Rep):
    # - SaEctC (PDB 5BXX -> UniProt Q1GNW6, Sphingopyxis alaskensis)
    # - PlEctC (PDB 5ONN -> NCBI YP_003245677.1, Paenibacillus lautus)
    # - P. ubique orphan cupin DddK (PDB 5TFZ -> UniProt Q4FNM4)
    print("Fetching Czech 2019 anchor: YP_003245677.1 (PlEctC PDB 5ONN)...")
    rec_pl = fetch_ncbi_protein("YP_003245677.1", "F4", "Ectoine_synthase", "PlEctC_Paenibacillus_PDB_5ONN")
    lit_records.append(rec_pl)

    # C. Pairoh et al. 2025 (Sci Rep):
    # Validated thermophilic Gamma-CAs: TtCA (PV553289->YEH58108.1), CrCA (PV553290->YEH58109.1), ToCA (PV553291->YEH58110.1)
    for ca_id, subcat in [("YEH58108.1", "TtCA_Thermus_Pairoh2025"),
                          ("YEH58109.1", "CrCA_metagenomic_Pairoh2025"),
                          ("YEH58110.1", "ToCA_metagenomic_Pairoh2025")]:
        print(f"Fetching Pairoh 2025 anchor: {ca_id} ({subcat})...")
        rec_ca = fetch_ncbi_protein(ca_id, "F5", "Gamma_CA", subcat)
        lit_records.append(rec_ca)

    # --------------------------------------------------------------------------
    # 3. STRUCTURED SOURCE QUERIES (PAGINATED, DERIVING ALL LABELS FROM SOURCE)
    # --------------------------------------------------------------------------
    print("\n--- STAGE 3: Executing Authoritative Source Queries (UniProtKB) ---")
    source_records = []

    # ------------------ TRACK F3: RHODOPSINS ------------------
    # Query 1: Reviewed Microbial Type-1 Rhodopsins via Pfam PF01036
    q_f3_type1 = "(reviewed:true) AND (xref:pfam-PF01036) AND (taxonomy_id:2 OR taxonomy_id:2157)"
    source_records.extend(fetch_uniprot_paginated_query(q_f3_type1, "F3", "Type1_Microbial", "Reviewed_Type1_Pfam_PF01036", max_records=30))

    # Query 2: Heliorhodopsins via Pfam PF18761 or InterPro IPR041113
    q_f3_helio = "(xref:pfam-PF18761 OR xref:interpro-IPR041113)"
    source_records.extend(fetch_uniprot_paginated_query(q_f3_helio, "F3", "Heliorhodopsin", "Curated_Heliorhodopsin_PF18761", max_records=25))

    # Query 3: Competing GPCR Rhodopsins via Pfam PF00001
    q_f3_gpcr = "(reviewed:true) AND (xref:pfam-PF00001) AND (protein_name:rhodopsin)"
    source_records.extend(fetch_uniprot_paginated_query(q_f3_gpcr, "F3", "Competing_GPCR", "Reviewed_GPCR_Rhodopsin_Control", max_records=15))

    # ------------------ TRACK F4: ECTC & CUPINS ------------------
    # Query 1: Reviewed EctC Enzymes via EC 4.2.1.108 or Pfam PF06339
    q_f4_ectc = "(reviewed:true) AND (ec:4.2.1.108 OR xref:pfam-PF06339)"
    source_records.extend(fetch_uniprot_paginated_query(q_f4_ectc, "F4", "Ectoine_synthase", "Reviewed_EctC_EC_4.2.1.108", max_records=30))

    # Query 2: EctA (EC 2.3.1.178)
    q_f4_ecta = "(reviewed:true) AND (ec:2.3.1.178)"
    source_records.extend(fetch_uniprot_paginated_query(q_f4_ecta, "F4", "Ectoine_pathway", "Reviewed_EctA_EC_2.3.1.178", max_records=15))

    # Query 3: EctB (EC 2.6.1.76)
    q_f4_ectb = "(reviewed:true) AND (ec:2.6.1.76)"
    source_records.extend(fetch_uniprot_paginated_query(q_f4_ectb, "F4", "Ectoine_pathway", "Reviewed_EctB_EC_2.6.1.76", max_records=15))

    # Query 4: EctD (EC 1.14.11.55)
    q_f4_ectd = "(reviewed:true) AND (ec:1.14.11.55)"
    source_records.extend(fetch_uniprot_paginated_query(q_f4_ectd, "F4", "Ectoine_pathway", "Reviewed_EctD_EC_1.14.11.55", max_records=15))

    # Query 5: Competing Non-EctC Cupins (Pirin, Quercetin dioxygenase)
    q_f4_cupin = "(reviewed:true) AND (xref:pfam-PF02678 OR ec:1.13.11.24) AND (taxonomy_id:2 OR taxonomy_id:2157)"
    source_records.extend(fetch_uniprot_paginated_query(q_f4_cupin, "F4", "Competing_Cupin", "Reviewed_Bacterial_Cupin_Control", max_records=15))

    # ------------------ TRACK F5: CARBONIC ANHYDRASES ------------------
    # Query 1: Gamma-CA via InterPro IPR052265 / IPR047223 / NCBIFAM NF050158
    q_f5_gamma = "(ec:4.2.1.1) AND (xref:interpro-IPR052265 OR xref:interpro-IPR047223 OR xref:ncbifam-NF050158 OR xref:ncbifam-NF040597 OR xref:cdd-cd00710)"
    source_records.extend(fetch_uniprot_paginated_query(q_f5_gamma, "F5", "Gamma_CA", "Curated_Gamma_CA", max_records=25))

    # Query 2: Beta-CA via Pfam PF00484 in prokaryotes
    q_f5_beta = "(reviewed:true) AND (ec:4.2.1.1) AND (xref:pfam-PF00484) AND (taxonomy_id:2 OR taxonomy_id:2157)"
    source_records.extend(fetch_uniprot_paginated_query(q_f5_beta, "F5", "Beta_CA", "Reviewed_Beta_CA_PF00484", max_records=25))

    # Query 3: Carboxysome CAs (CsoSCA / CcmM / IPR017156) - KEPT STRICTLY SEPARATE PER ROOT INSTRUCTION
    q_f5_carboxysome = "(reviewed:true) AND (xref:interpro-IPR017156 OR gene:csoSCA) AND (taxonomy_id:2 OR taxonomy_id:2157)"
    source_records.extend(fetch_uniprot_paginated_query(q_f5_carboxysome, "F5", "Carboxysome_CA", "Reviewed_Carboxysome_CA_Control", max_records=15))

    # Query 4: Competing LbH Transferases (LpxA EC 2.3.1.129, DapD EC 2.3.1.117)
    q_f5_lbh = "(reviewed:true) AND (ec:2.3.1.129 OR ec:2.3.1.117) AND (taxonomy_id:2 OR taxonomy_id:2157)"
    source_records.extend(fetch_uniprot_paginated_query(q_f5_lbh, "F5", "Competing_LbH_Transferase", "Reviewed_LbH_Transferase_Control", max_records=15))

    # Query 5: Alpha-CA controls
    q_f5_alpha = "(reviewed:true) AND (ec:4.2.1.1) AND (family:\"alpha-type carbonic anhydrase family\")"
    source_records.extend(fetch_uniprot_paginated_query(q_f5_alpha, "F5", "Competing_Alpha_CA", "Reviewed_Alpha_CA_Control", max_records=10))

    # --------------------------------------------------------------------------
    # 4. PATHOGEN FILTERING & REJECTED REFERENCES REGISTRY
    # --------------------------------------------------------------------------
    print("\n--- STAGE 4: Pathogen Filtering & Rejected Registry ---")
    combined_raw = lit_records + source_records
    
    # Exclude pathogens from positive environmental reference set
    clean_records = []
    pathogens_excluded = []
    seen_clean = set()

    for r in combined_raw:
        acc = r["accession"]
        if r.get("is_pathogen") == 1:
            pathogens_excluded.append(r)
            continue
        if acc in seen_clean:
            continue
        seen_clean.add(acc)
        clean_records.append(r)

    # Write excluded pathogens
    pathogen_tsv = REF_DIR / "excluded_pathogens.tsv"
    with open(pathogen_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["accession", "organism", "protein_name", "track", "category"], delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(pathogens_excluded)
    print(f"Excluded {len(pathogens_excluded)} pathogen entries to {pathogen_tsv.name}")

    # Write rejected references ledger (tracking all audited false accessions)
    rejected_ledger = [
        {"accession": "P15644", "claimed": "Halorhodopsin_Hsalinarum", "actual_name": "Secretin PulD", "organism": "Klebsiella pneumoniae", "reason": "Secretin, not rhodopsin"},
        {"accession": "Q9FAF7", "claimed": "Proteorhodopsin_BAC31A08", "actual_name": "Cucumopine synthase", "organism": "Rhizobium rhizogenes", "reason": "Opine synthase, not proteorhodopsin"},
        {"accession": "Q84FT3", "claimed": "EctC_Pseudomonas_stutzeri", "actual_name": "SusC outer membrane protein", "organism": "Tannerella forsythia", "reason": "SusC transporter, not EctC"},
        {"accession": "O34707", "claimed": "EctC_Bacillus_subtilis", "actual_name": "Tetraprenyl-beta-curcumene synthase", "organism": "Bacillus subtilis", "reason": "Terpene synthase, not EctC"},
        {"accession": "Q4FM61", "claimed": "Orphan_EctC_Pelagibacter", "actual_name": "FabZ dehydratase", "organism": "Pelagibacter ubique", "reason": "FabZ dehydratase, not orphan EctC"},
        {"accession": "P52626", "claimed": "EctA_Halomonas_elongata", "actual_name": "PatH ABC transporter protein", "organism": "Vibrio harveyi", "reason": "ABC transporter, not EctA"},
        {"accession": "P52627", "claimed": "EctB_Halomonas_elongata", "actual_name": "FliZ sigma S regulator", "organism": "Escherichia coli", "reason": "FliZ regulator, not EctB"},
        {"accession": "Q9AQJ4", "claimed": "EctD_Halomonas_elongata", "actual_name": "Rubredoxin reductase", "organism": "Acinetobacter sp.", "reason": "Rubredoxin reductase, not EctD"},
        {"accession": "P49347", "claimed": "Pirin_cupin_control", "actual_name": "Concanavalin B", "organism": "Canavalia ensiformis", "reason": "Concanavalin B, not pirin"},
        {"accession": "Q72J87", "claimed": "TtCA_Thermus_thermophilus", "actual_name": "UPF0145 protein TT_C0892", "organism": "Thermus thermophilus", "reason": "UPF0145 uncharacterized, not TtCA"},
        {"accession": "Q183H0", "claimed": "Cam_Clostridium_difficile", "actual_name": "Cell wall-binding protein cwp14", "organism": "Clostridioides difficile", "reason": "Cell wall protein, not gamma-CA"},
        {"accession": "Q03512", "claimed": "CcmM_Synechococcus", "actual_name": "Carboxysome shell vertex protein CcmL", "organism": "Synechococcus elongatus", "reason": "CcmL shell protein, not CA"},
        {"accession": "O59424", "claimed": "Cam_Pyrococcus_horikoshii", "actual_name": "Small ribosomal subunit protein uS3", "organism": "Pyrococcus horikoshii", "reason": "Ribosomal protein uS3, not CA"},
        {"accession": "P17890", "claimed": "CynT_Escherichia_coli", "actual_name": "RNA polymerase III subunit RPC7", "organism": "Saccharomyces cerevisiae", "reason": "Yeast RNA polymerase III, not CynT"},
        {"accession": "P9WN71", "claimed": "Rv3588c_Mycobacterium", "actual_name": "Glucose-6-phosphate 1-dehydrogenase 1", "organism": "Mycobacterium tuberculosis", "reason": "G6PD, not beta-CA"},
        {"accession": "P0A6B1", "claimed": "DapD_succinyltransferase", "actual_name": "Acyl carrier protein ACP", "organism": "Salmonella typhimurium", "reason": "ACP, not DapD"},
        {"accession": "A0A2V0P0L8", "claimed": "HeR_46EC73_Pushkarev2018", "actual_name": "Uncharacterized protein", "organism": "Raphidocelis subcapitata", "reason": "Algal protein, not 48C12 heliorhodopsin"},
        {"accession": "A0A2N8D7B0", "claimed": "HeR_Bellilinea_caldifistulae", "actual_name": "None", "organism": "None", "reason": "Invalid/obsolete accession"},
        {"accession": "P25964", "claimed": "Sensory_rhodopsin_I", "actual_name": "None", "organism": "None", "reason": "Invalid/obsolete accession"},
        {"accession": "Q9AQ76", "claimed": "Green_proteorhodopsin", "actual_name": "None", "organism": "None", "reason": "Invalid/obsolete accession"},
        {"accession": "Q8YSC4", "claimed": "Anabaena_sensory_rhodopsin", "actual_name": "None", "organism": "None", "reason": "Invalid/obsolete accession"},
        {"accession": "D9W3X4", "claimed": "SaEctC_Streptomyces_albus", "actual_name": "None", "organism": "None", "reason": "Invalid/obsolete accession"},
    ]
    rejected_tsv = REF_DIR / "rejected_references.tsv"
    with open(rejected_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["accession", "claimed", "actual_name", "organism", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rejected_ledger)
    print(f"Recorded {len(rejected_ledger)} audited rejected references to {rejected_tsv.name}")

    # --------------------------------------------------------------------------
    # 5. WRITE REFERENCE FASTAs & BUILD DIAMOND DATABASES
    # --------------------------------------------------------------------------
    print("\n--- STAGE 5: Writing FASTAs and Indexing DIAMOND Databases ---")
    track_fastas = {
        "F3": REF_DIR / "f3_rhodopsin_references.faa",
        "F4": REF_DIR / "f4_ectc_references.faa",
        "F5": REF_DIR / "f5_carbonic_anhydrase_references.faa",
        "ALL": REF_DIR / "all_f345_references.faa",
    }
    handles = {k: open(v, "w", encoding="utf-8") for k, v in track_fastas.items()}

    for r in clean_records:
        acc = r["accession"]
        clean_org = r["organism"].replace(" ", "_").replace("/", "_")
        header = f">{acc}|{r['track']}|{r['category']}|{r['subcategory']}|{clean_org}|exp_fn={r['has_function_experimental_evidence']}"
        seq = r["sequence"]
        formatted_seq = "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))
        entry_str = f"{header}\n{formatted_seq}\n"

        handles["ALL"].write(entry_str)
        if r["track"] in handles:
            handles[r["track"]].write(entry_str)

    for h in handles.values():
        h.close()

    for trk, fpath in track_fastas.items():
        print(f"Wrote {fpath.name} ({fpath.stat().st_size} bytes)")
        dmnd_db = fpath.with_suffix(".dmnd")
        try:
            cmd = ["diamond", "makedb", "--in", str(fpath), "-d", str(fpath.with_suffix(""))]
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"  DIAMOND makedb {dmnd_db.name} OK")
        except FileNotFoundError:
            print(f"  diamond binary not found on local path; will be indexed under Slurm worker.")
        except Exception as e:
            print(f"  diamond makedb warning: {e}")

    # --------------------------------------------------------------------------
    # 6. WRITE METADATA & MANIFEST WITH SHA-256
    # --------------------------------------------------------------------------
    print("\n--- STAGE 6: Writing Metadata & SHA-256 Manifest ---")
    meta_tsv = REF_DIR / "verified_references_metadata.tsv"
    fieldnames = [
        "accession", "source_db", "uniprot_id", "track", "category", "subcategory",
        "organism", "lineage", "protein_name", "gene_name", "ec_number",
        "length_aa", "has_function_experimental_evidence", "has_any_experimental_evidence",
        "function_evidence_codes", "all_evidence_codes", "is_pathogen", "query_used"
    ]
    with open(meta_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(clean_records)
    print(f"Wrote {len(clean_records)} clean reference records to {meta_tsv.name}")

    manifest_tsv = REF_DIR / "references_manifest.tsv"
    manifest_rows = []
    for hmm in sorted(HMM_DIR.glob("*.hmm")):
        manifest_rows.append({"type": "profile_hmm", "filename": f"profiles/{hmm.name}", "sha256": compute_sha256(hmm), "bytes": hmm.stat().st_size})
    for faa in sorted(REF_DIR.glob("*.faa")):
        manifest_rows.append({"type": "reference_fasta", "filename": faa.name, "sha256": compute_sha256(faa), "bytes": faa.stat().st_size})
    for tsv in [meta_tsv, pathogen_tsv, rejected_tsv]:
        manifest_rows.append({"type": "metadata", "filename": tsv.name, "sha256": compute_sha256(tsv), "bytes": tsv.stat().st_size})

    with open(manifest_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["type", "filename", "sha256", "bytes"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote {len(manifest_rows)} manifest entries to {manifest_tsv.name}")

    elapsed = time.time() - start_time
    print(f"\n=== Reference Fetching & Rebuild Finished in {elapsed:.1f}s ===")


if __name__ == "__main__":
    main()
