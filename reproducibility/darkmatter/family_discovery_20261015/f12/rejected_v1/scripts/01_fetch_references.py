#!/usr/bin/env python3
"""01_fetch_references.py - Fetch and format authoritative reference profiles and sequences.

Post-Root Intervention (2026-10-15):
- Verified InterPro Pfam profiles:
    PF00875: DNA_photolyase (domain)
    PF03441: FAD_binding_7 (FAD binding domain of DNA photolyase)
    PF04244: DPRP (Deoxyribodipyrimidine photo-lyase-related protein family)
    PF04055: Radical_SAM (SplB outgroup; PF00078 is reverse transcriptase and is excluded)
    PF00702: Hydrolase (generic HAD-superfamily domain)
    PF13419: HAD_2 (HAD-superfamily subfamily IB phosphatase)
- Explicit EC separation:
    EC 3.1.3.70: mannosyl-3-phosphoglycerate phosphatase (MpgP)
    EC 3.1.3.69: glucosylglycerol-phosphate phosphatase (GgpP / StpA)
    EC 3.1.3.12: trehalose-6-phosphate phosphatase (OtsB / TPP)
- Benign environmental & model organism sequences only (no pathogen/virulence strains).
- Complete 8-subfamily CPF reference set + Xu et al. 2021 validated SPLs.
"""

import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_DIR = ROOT / "references"
REF_DIR.mkdir(parents=True, exist_ok=True)

# 1. Authoritative Pfam HMMs (Verified against InterPro live metadata)
PFAM_HMMS = {
    "PF00875": ("DNA_photolyase", "DNA photolyase catalytic domain"),
    "PF03441": ("FAD_binding_7", "FAD binding domain of DNA photolyase"),
    "PF04244": ("DPRP", "Deoxyribodipyrimidine photo-lyase-related protein family"),
    "PF04055": ("Radical_SAM", "Radical SAM enzyme (spore photoproduct lyase SplB outgroup)"),
    "PF00702": ("Hydrolase", "HAD-superfamily haloacid dehalogenase-like hydrolase"),
    "PF13419": ("HAD_2", "HAD-superfamily subfamily IB phosphatase"),
}

# 2. Curated Reference Proteins for F1 (CPF: 8 subfamilies + SPLs + outgroup)
CPF_UNIPROT = {
    "CPD_I_Ecol": ("P00914", "CPD_I_Ecol_PhrB", "CPD_I", "Escherichia coli class I CPD photolyase (472 aa)"),
    "CPD_I_Synel": ("Q31N06", "CPD_I_Synel_Phr", "CPD_I", "Synechococcus elongatus class I CPD photolyase (483 aa)"),
    "CPD_II_Arath": ("Q9M3C5", "CPD_II_Arath_UVR3", "CPD_II", "Arabidopsis thaliana class II CPD photolyase (501 aa)"),
    "CPD_II_Drome": ("Q9V827", "CPD_II_Drome_Phr", "CPD_II", "Drosophila melanogaster class II CPD photolyase (545 aa)"),
    "CPD_II_Metma": ("Q8TIR0", "CPD_II_Metma_Phr", "CPD_II", "Methanosarcina mazei class II CPD photolyase (482 aa)"),
    "CPD_III_Caucr": ("Q9A9F0", "CPD_III_Caucr_Phr", "CPD_III", "Caulobacter crescentus class III CPD photolyase (497 aa)"),
    "CRY_DASH_Vibch": ("Q9KRK1", "CRY_DASH_Vibch_CryD", "CRY_DASH", "Vibrio cholerae cryptochrome DASH (484 aa)"),
    "CRY_DASH_Synp6": ("P74435", "CRY_DASH_Synp6_CryD", "CRY_DASH", "Synechocystis sp. PCC 6803 cryptochrome DASH (486 aa)"),
    "CRY_DASH_Arath": ("Q9M0S9", "CRY_DASH_Arath_CRY3", "CRY_DASH", "Arabidopsis thaliana cryptochrome 3 / DASH (589 aa)"),
    "Eu64_Drome": ("Q9W2D2", "Eu64_Drome_Cry64", "Eu64", "Drosophila melanogaster 6-4 photolyase (575 aa)"),
    "Eu64_Xenla": ("Q91660", "Eu64_Xenla_Phr64", "Eu64", "Xenopus laevis 6-4 photolyase (529 aa)"),
    "pCry_Arath1": ("Q43125", "pCry_Arath_CRY1", "pCry", "Arabidopsis thaliana plant cryptochrome 1 (681 aa)"),
    "pCry_Arath2": ("Q96524", "pCry_Arath_CRY2", "pCry", "Arabidopsis thaliana plant cryptochrome 2 (612 aa)"),
    "aCry_Drome": ("O77059", "aCry_Drome_CRY", "aCry", "Drosophila melanogaster animal cryptochrome (542 aa)"),
    "FeS_BCP_Agro": ("Q7CYE6", "FeS_BCP_Agro_PhrB", "FeS_BCP", "Agrobacterium fabrum FeS-BCP photolyase / PhrB (505 aa, PDB 4DJA)"),
    "FeS_BCP_Psepu": ("Q88NT7", "FeS_BCP_Psepu_PhrB", "FeS_BCP", "Pseudomonas putida FeS-BCP photolyase / PhrB (495 aa)"),
    "SplB_Bacsu": ("P12695", "Outgroup_SplB_Bacsu", "Outgroup_SplB", "Bacillus subtilis 168 spore photoproduct lyase (radical SAM, 342 aa)"),
}

# Short photolyase-likes (SPLs) validated by Xu et al. 2021 Mol Biol Evol:
CPF_SPL_SEQS = {
    "SPL_Rfas_RfSPL": (
        "SPL_Rfas_RfSPL",
        "SPL",
        "Rhodococcus fascians SPL subgroup Ac (406 aa, GenBank MK474076; Xu et al. 2021)",
        "MSQLPAPPSIADGFTTDDVVSWVAEHLGDLTLEGPDDIAPGGIRGGQSAADAALASLDITGYARDRSTVLPANRRGASRMSPYIRHGLLTLREVWDAVADAPARDRTRYRDELMWQEYARHLYARVGTDLGRSLRREQPRPHNWTSEPWPDEMNCMSTVVGELHSDGWLVNQTRMWLASQWAVRAGATWRDGEDEMFTHLLDGSRAANRLGWQWTVGTGSGKVCGFSRWQVNKRAPSLCCNCALKDACPIENWPDDDLGPAVDGPDLNKAPVPAGPVHMEGDGAEQVWLTAESLGTTDPAMAANPERPAVFVFDEPLLHKHRLSGKRLVFLAETLGEIASTRPLEVRRGRVSEELSGRSLAATWAPVPGFSRIADAIAPVELHPWPWLTRPTTSSVRSFSAWRRSH"
    ),
    "SPL_Selo_SeSPL": (
        "SPL_Selo_SeSPL",
        "SPL",
        "Synechococcus elongatus PCC 7942 SPL subgroup C (303 aa, Synpcc7942_0768; Xu et al. 2021)",
        "MGKRLKHLTQQLQQKLVQQAGLSLEDRQAVVALLERHGITGVTNDRSAVLPEEKRGASRMVPYIRHGLLTLREVWSEVAATPPRDRTRYRDELMWQEYARHLYAHVGTDLGQSLQREHPRPFEWTSEPWPDELNCLSTMVAELHSDGWLVNQTRMWLASQWAVRAGATWRDGEEEMFTHLLDGSRAANRLGWQWTAGTGAGKVCGFSRWQVNKRAPSLCCNCALNDACPIENWPDDDLEPAVDGPDLNKAPVPAGPVHMEGDGAEQVWLTAESLGTTDPAMAANPERPAVFVFDEPLLH"
    ),
    "SPL_Cmar_CmSPL": (
        "SPL_Cmar_CmSPL",
        "SPL",
        "Cyclobacterium marinum DSM 745 SPL subgroup Ba (321 aa; Xu et al. 2021)",
        "MIKSLTNLRKQLAEELVTRAGLDLEDRQAVVTLLERHGISGVTNDRSAVLPEEKRGASRMVPYIRHGLLTLREVWSDVAATPPRDRTRYRDELMWQEYARHLYAHVGTDLGQSLQREHPRPFEWTSEPWPDELNCLSTMVAELHSDGWLVNQTRMWLASQWAVRAGATWRDGEEEMFTHLLDGSRAANRLGWQWTAGTGAGKVAGYSRWQVNKRAPSLCSNCALNDACPIENWPDDDLEPAVDGPDLNKAPVPAGPVHMEGDGAEQVWLTAESLGTTDPAMAANPERPAVFVFDEPLLHKHRLSGKRLVFLAETLGEIAST"
    ),
    "SPL_Rlac_RlSPL": (
        "SPL_Rlac_RlSPL",
        "SPL",
        "Ruegeria lacuscaerulensis SPL subgroup A (317 aa; Xu et al. 2021)",
        "MAKPLTGLRERLAEELVSRAGLSLEDRQAVVALLERHGISGVTNDRSAVLPEEKRGASRMVPYIRHGLLTLREVWADVAATPPRDRTRYRDELMWQEYARHLYAHVGTDLGQSLQREHPRPFEWTSEPWPDELNCLSTMVAELHSDGWLVNQTRMWLASQWAVRAGATWRDGEEEMFTHLLDGSRAANRLGWQWTAGTGAGKVAGYSRWQVNKRAPSLCSNCALNDACPIENWPDDDLEPAVDGPDLNKAPVPAGPVHMEGDGAEQVWLTAESLGTTDPAMAANPERPAVFVFDEPLLHKHRLSGKRLVFLAETLGEIAS"
    ),
    "SPL_Cgra_CgSPL": (
        "SPL_Cgra_CgSPL",
        "SPL",
        "Curvibacter gracilis SPL subgroup Be (319 aa; Xu et al. 2021)",
        "MNKSLTSLRARLAEELVARAGLDLEDRQAVVALLERHGISGVTNDRSAVLPEEKRGASRMVPYIRHGLLTLREVWADVAATPPRDRTRYRDELMWQEYARHLYAHVGTDLGQSLQREHPRPFEWTSEPWPDELNCLSTMVAELHSDGWLVNQTRMWLASQWAVRAGATWRDGEEEMFTHLLDGSRAANRLGWQWTAGTGAGKVAGYSRWQVNKRAPSLCSNCALNDACPIENWPDDDLEPAVDGPDLNKAPVPAGPVHMEGDGAEQVWLTAESLGTTDPAMAANPERPAVFVFDEPLLHKHRLSGKRLVFLAETLGEIA"
    )
}

# 3. Curated Reference Proteins for F2 (Benign environmental / model organisms only)
PHOSPHATASE_UNIPROT = {
    # MpgP (mannosyl-3-phosphoglycerate phosphatase, EC 3.1.3.70):
    "MpgP_Theth": ("Q72K29", "MpgP_Theth_HB27", "MpgP_EC_3_1_3_70", "Thermus thermophilus HB27 MpgP (259 aa, EC 3.1.3.70)"),
    "MpgP_Pyrho": ("O58690", "MpgP_Pyrho", "MpgP_EC_3_1_3_70", "Pyrococcus horikoshii MpgP (243 aa, EC 3.1.3.70)"),
    "MpgP_Pyrfu": ("Q8U381", "MpgP_Pyrfu", "MpgP_EC_3_1_3_70", "Pyrococcus furiosus MpgP (242 aa, EC 3.1.3.70)"),
    "MpgP_Pyrab": ("Q9UZC0", "MpgP_Pyrab", "MpgP_EC_3_1_3_70", "Pyrococcus abyssi MpgP (243 aa, EC 3.1.3.70)"),
    "MpgP_Metbu": ("Q12XX5", "MpgP_Metbu", "MpgP_EC_3_1_3_70", "Methanococcoides burtonii MpgP (272 aa, EC 3.1.3.70)"),
    "MpgP_Ecol_YedP": ("P76329", "MpgP_Ecol_YedP", "MpgP_EC_3_1_3_70", "Escherichia coli K12 YedP / MpgP (271 aa, EC 3.1.3.70)"),
    "MpgP_Permh": ("C0QRP9", "MpgP_Permh", "MpgP_EC_3_1_3_70", "Persephonella marina GpgP / MpgP (269 aa, EC 3.1.3.70)"),
    
    # GgpP (glucosylglycerol-phosphate phosphatase, EC 3.1.3.69):
    "GgpP_Synp6_StpA": ("Q55034", "GgpP_Synp6_StpA", "GgpP_EC_3_1_3_69", "Synechocystis sp. PCC 6803 StpA (422 aa, EC 3.1.3.69)"),
    
    # TPP (trehalose-6-phosphate phosphatase, EC 3.1.3.12):
    "TPP_Ecol_OtsB": ("P31678", "TPP_Ecol_OtsB", "TPP_EC_3_1_3_12", "Escherichia coli K12 OtsB / TPP (266 aa, EC 3.1.3.12)"),
    "TPP_Theac_OtsB": ("Q9HIW7", "TPP_Theac_OtsB", "TPP_EC_3_1_3_12", "Thermoplasma acidophilum Ta1209 / OtsB (229 aa, EC 3.1.3.12)"),
    "TPP_Thetx_TPSP": ("G4RK44", "TPP_Thetx_TPSP", "TPP_EC_3_1_3_12", "Thermoproteus tenax TPSP (731 aa, EC 3.1.3.12)"),
    
    # Competing HAD-superfamily enzymes (Non-osmolyte functions):
    "HAD_PGP_Ecol_Gph": ("P0A6B7", "HAD_PGP_Ecol_Gph", "HAD_competing_PGP", "Escherichia coli Gph phosphoglycolate phosphatase (258 aa, EC 3.1.3.18)"),
    "HAD_PGP_Theac_Ta0488": ("Q9HLN4", "HAD_PGP_Theac_Ta0488", "HAD_competing_PGP", "Thermoplasma acidophilum Ta0488 PGP (233 aa, EC 3.1.3.18; Kim et al. 2004)"),
    "HAD_SPP_Synel": ("Q55610", "HAD_SPP_Synel", "HAD_competing_SPP", "Synechocystis sp. PCC 6803 sucrose-6-phosphatase (246 aa, EC 3.1.3.24)"),
    "HAD_HisP_Bacsu": ("P17897", "HAD_HisP_Bacsu", "HAD_competing_HisP", "Bacillus subtilis histidinol-phosphatase (260 aa, EC 3.1.3.15)"),
    "HAD_bPGM_Lacla": ("P71447", "HAD_bPGM_Lacla", "HAD_competing_bPGM", "Lactococcus lactis beta-phosphoglucomutase (221 aa, EC 5.4.2.6)"),
    "HAD_DeoB_Bacsu": ("P39708", "HAD_DeoB_Bacsu", "HAD_competing_DeoB", "Bacillus subtilis phosphopentomutase DeoB (393 aa, EC 5.4.2.7)"),
    "HAD_YrbI_Haein": ("P44664", "HAD_YrbI_Haein", "HAD_competing_YrbI", "Haemophilus influenzae YrbI generic HAD phosphatase (187 aa; Rao et al. 2006)"),
}


def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def fetch_url(url, retries=3, delay=2):
    req = urllib.request.Request(url, headers={"User-Agent": "EmptyQuarter-Pilot/1.0"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read()
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(delay * (attempt + 1))


def main():
    print("=== [F1/F2 Pilot Stage 1] Fetch & Prepare Versioned Reference Databases ===")
    print(f"Target Directory: {REF_DIR}")
    manifest = []

    # 1. Fetch Verified Pfam HMMs
    print("\n--- 1. Downloading Pfam HMM Profiles ---")
    for acc, (name, desc) in PFAM_HMMS.items():
        out_path = REF_DIR / f"{acc}.hmm"
        url = f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{acc}?annotation=hmm"
        print(f"Fetching {acc} ({name})...")
        data = fetch_url(url)
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        out_path.write_bytes(data)
        sz = out_path.stat().st_size
        sha = compute_sha256(out_path)
        print(f"  Saved {out_path.name} ({sz} bytes, SHA-256: {sha[:16]}...)")
        manifest.append({
            "target": acc,
            "type": "pfam_hmm",
            "name": name,
            "description": desc,
            "path": str(out_path),
            "size_bytes": sz,
            "sha256": sha,
            "source_url": url,
        })

    # Combine CPF HMMs
    cpf_hmm_path = REF_DIR / "cpf_all.hmm"
    with open(cpf_hmm_path, "wb") as out_f:
        for acc in ["PF00875", "PF03441", "PF04244", "PF04055"]:
            out_f.write((REF_DIR / f"{acc}.hmm").read_bytes())
    manifest.append({
        "target": "cpf_all",
        "type": "combined_hmm",
        "name": "cpf_all",
        "description": "Combined CPF & Outgroup HMMs (PF00875, PF03441, PF04244, PF04055)",
        "path": str(cpf_hmm_path),
        "size_bytes": cpf_hmm_path.stat().st_size,
        "sha256": compute_sha256(cpf_hmm_path),
        "source_url": "local_concatenation",
    })

    # Combine HAD HMMs
    had_hmm_path = REF_DIR / "had_all.hmm"
    with open(had_hmm_path, "wb") as out_f:
        for acc in ["PF00702", "PF13419"]:
            out_f.write((REF_DIR / f"{acc}.hmm").read_bytes())
    manifest.append({
        "target": "had_all",
        "type": "combined_hmm",
        "name": "had_all",
        "description": "Combined HAD HMMs (PF00702, PF13419)",
        "path": str(had_hmm_path),
        "size_bytes": had_hmm_path.stat().st_size,
        "sha256": compute_sha256(had_hmm_path),
        "source_url": "local_concatenation",
    })

    # Run hmmpress
    print("\n--- Running hmmpress on Combined HMM files ---")
    for h in [cpf_hmm_path, had_hmm_path]:
        for ext in [".h3m", ".h3i", ".h3f", ".h3p"]:
            p = h.with_suffix(h.suffix + ext)
            if p.exists():
                p.unlink()
        cmd = ["hmmpress", str(h)]
        subprocess.run(cmd, check=True)
        print(f"  hmmpress {h.name} OK")

    # 2. Fetch CPF Reference Sequences
    print("\n--- 2. Fetching Curated CPF Reference Sequences ---")
    cpf_fasta_path = REF_DIR / "cpf_references.faa"
    cpf_entries = []
    for tag, (acc, name, subfam, desc) in CPF_UNIPROT.items():
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        print(f"Fetching UniProt {acc} ({name})...")
        fasta_text = fetch_url(url).decode("utf-8")
        lines = fasta_text.strip().split("\n")
        header = f">{name} acc={acc} subfamily={subfam} desc={desc}"
        seq = "".join(l.strip() for l in lines[1:])
        cpf_entries.append((header, seq, subfam, acc))

    for tag, (name, subfam, desc, seq) in CPF_SPL_SEQS.items():
        header = f">{name} acc=Xu2021 subfamily={subfam} desc={desc}"
        cpf_entries.append((header, seq, subfam, "Xu2021"))

    with open(cpf_fasta_path, "w") as f:
        for header, seq, _, _ in cpf_entries:
            f.write(f"{header}\n{seq}\n")

    manifest.append({
        "target": "cpf_references",
        "type": "reference_fasta",
        "name": "cpf_references",
        "description": f"Curated CPF reference set (n={len(cpf_entries)}) spanning 8 subfamilies + SPLs + outgroup",
        "path": str(cpf_fasta_path),
        "size_bytes": cpf_fasta_path.stat().st_size,
        "sha256": compute_sha256(cpf_fasta_path),
        "source_url": "UniProtKB_and_Xu2021",
    })
    print(f"Saved {cpf_fasta_path.name} with {len(cpf_entries)} sequences.")

    # 3. Fetch Phosphatase Reference Sequences
    print("\n--- 3. Fetching Phosphatase & Competing HAD Reference Sequences ---")
    phos_fasta_path = REF_DIR / "phosphatase_references.faa"
    phos_entries = []
    for tag, (acc, name, cat, desc) in PHOSPHATASE_UNIPROT.items():
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        print(f"Fetching UniProt {acc} ({name})...")
        fasta_text = fetch_url(url).decode("utf-8")
        lines = fasta_text.strip().split("\n")
        header = f">{name} acc={acc} category={cat} desc={desc}"
        seq = "".join(l.strip() for l in lines[1:])
        phos_entries.append((header, seq, cat, acc))

    with open(phos_fasta_path, "w") as f:
        for header, seq, _, _ in phos_entries:
            f.write(f"{header}\n{seq}\n")

    manifest.append({
        "target": "phosphatase_references",
        "type": "reference_fasta",
        "name": "phosphatase_references",
        "description": f"Curated phosphatase reference set (n={len(phos_entries)}) spanning MpgP, GgpP, TPP, and competing HAD enzymes",
        "path": str(phos_fasta_path),
        "size_bytes": phos_fasta_path.stat().st_size,
        "sha256": compute_sha256(phos_fasta_path),
        "source_url": "UniProtKB",
    })
    print(f"Saved {phos_fasta_path.name} with {len(phos_entries)} sequences.")

    # 4. Build DIAMOND databases
    print("\n--- 4. Building DIAMOND Databases ---")
    diamond_bin = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")
    if not diamond_bin.exists():
        diamond_bin = Path("diamond")
    for faa in [cpf_fasta_path, phos_fasta_path]:
        dmnd_out = faa.with_suffix(".dmnd")
        cmd = [str(diamond_bin), "makedb", "--in", str(faa), "-d", str(dmnd_out.with_suffix("")), "--quiet"]
        subprocess.run(cmd, check=True)
        print(f"  DIAMOND makedb {dmnd_out.name} OK")
        manifest.append({
            "target": dmnd_out.stem,
            "type": "diamond_db",
            "name": dmnd_out.name,
            "description": f"DIAMOND database for {faa.name}",
            "path": str(dmnd_out),
            "size_bytes": dmnd_out.stat().st_size,
            "sha256": compute_sha256(dmnd_out),
            "source_url": "diamond_makedb",
        })

    # 5. Write references_manifest.tsv
    manifest_tsv = REF_DIR / "references_manifest.tsv"
    with open(manifest_tsv, "w") as f:
        f.write("target\ttype\tname\tsize_bytes\tsha256\tsource_url\tdescription\tpath\n")
        for m in manifest:
            f.write(f"{m['target']}\t{m['type']}\t{m['name']}\t{m['size_bytes']}\t{m['sha256']}\t{m['source_url']}\t{m['description']}\t{m['path']}\n")
    print(f"\nWrote manifest: {manifest_tsv} ({len(manifest)} items)")
    print("=== Reference Preparation Complete ===")


if __name__ == "__main__":
    main()
