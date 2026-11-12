#!/usr/bin/env python3
"""Stage 1: select proteins, attach SignalP, reuse TMbed labels by exact sequence, write FASTA files.

Scope: every member (EQ and external) of every EQ family of class dark_exported or nondark_exported
(30% identity / 80% coverage families of analysis2). Sequences/topology/domains are analysed for the EQ
members only; external members contribute only their SignalP class (route agreement).

Note: gene identifiers of the 2026-10-30 TMbed run coincide with the harmonised gene identifiers but about
a quarter of them denote different gene calls, so TMbed labels are reused ONLY for byte-identical sequences;
all other signal-peptide-positive EQ members go to work/tmbed_parts/ for a new TMbed run.

Outputs (work/): members.tsv.gz, eq_members_full.faa, eq_all_proteins.faa, tmbed_reused.tsv.gz,
tmbed_parts/part_XXX.faa (header = sequence hash), tmbed_qc_hashes.txt, prep_info.json
"""
import argparse, gzip, glob, json, re, hashlib
from pathlib import Path
import numpy as np
import pandas as pd

EQ = "EQ_Thermoplasmatota_237"
EXCLUDED = {"External_Comparators_CrossResource_Duplicate_54", "SPIRE_ProGenomes3_SpecI_Excluded_148"}
CS_RE = r"CS pos: (\d+)-(\d+)\. Pr: ([0-9.]+)"


def seqhash(s):
    return hashlib.sha256(s.encode()).hexdigest()[:24]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tmbed-dir", type=Path, required=True)
    ap.add_argument("--parts", type=int, default=32)
    ap.add_argument("--qc-n", type=int, default=1500)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    R = a.root
    info = {}
    gm = pd.read_csv(R / "metadata/harmonized_genome_manifest.tsv", sep="\t", dtype=str)
    gm = gm[~gm.cohort.isin(EXCLUDED)]
    eq_genomes = set(gm.genome_id[gm.cohort == EQ])
    assert len(eq_genomes) == 237 and len(gm) == 1957, (len(eq_genomes), len(gm))
    ft = pd.read_csv(R / "analysis2_families/results/eq_family_table_id30_c80_all_genes.tsv.gz", sep="\t")
    exported = ft[ft["class"].isin(["dark_exported", "nondark_exported"])]
    info["n_families_table"] = int(len(ft))
    info["n_exported_families"] = exported["class"].value_counts().to_dict()
    reps = set(exported.family_rep)
    fam = pd.read_csv(R / "analysis2_families/clusters/families_id30_c80.tsv", sep="\t", header=None,
                      names=["family_rep", "gene_id"], dtype=str)
    assert fam.gene_id.is_unique and len(fam) == 3617011
    fam = fam[fam.family_rep.isin(reps)].copy()
    fam["genome_id"] = fam.gene_id.str.split("::", n=1).str[0]
    assert fam.genome_id.isin(set(gm.genome_id)).all()
    fam["is_eq"] = fam.genome_id.isin(eq_genomes)
    sp = pd.read_csv(R / "signalp/harmonized_signalp6_predictions.tsv.gz", sep="\t", dtype=str,
                     usecols=["gene_id", "signalp_status", "signalp_class", "sp_any", "sp_prob", "cleavage_site_pos"])
    fam = fam.merge(sp, on="gene_id", how="left", validate="1:1")
    assert fam.signalp_class.notna().all()
    del sp
    dk = pd.read_csv(R / "mmseqs/harmonized_uniref50_per_gene.tsv.gz", sep="\t", dtype=str,
                     usecols=["gene_id", "primary_direct_uniref50_dark"])
    fam = fam.merge(dk, on="gene_id", how="left", validate="1:1")
    del dk
    fam["dark"] = fam.primary_direct_uniref50_dark == "True"
    fam = fam.drop(columns="primary_direct_uniref50_dark")
    cs = fam.cleavage_site_pos.str.extract(CS_RE)
    fam["cs_pos"] = pd.to_numeric(cs[0], errors="coerce").fillna(0).astype(int)   # last residue of the SP (1-based)
    fam["cs_prob"] = pd.to_numeric(cs[2], errors="coerce")
    man = pd.read_csv(R / "genecall/harmonized_gene_manifest.tsv.gz", sep="\t", dtype=str,
                      usecols=["gene_id", "strand", "partial_flag", "aa_length", "contig_edge_status"])
    fam = fam.merge(man, on="gene_id", how="left", validate="1:1")
    del man
    # prodigal partial flag: first char = left (contig-coordinate) end incomplete, second = right end.
    pf = fam.partial_flag.str.zfill(2)
    plus = fam.strand.astype(int) > 0
    fam["n_complete"] = np.where(plus, pf.str[0] == "0", pf.str[1] == "0")
    fam["c_complete"] = np.where(plus, pf.str[1] == "0", pf.str[0] == "0")
    fam.drop(columns=["cleavage_site_pos"]).to_csv(a.out / "members.tsv.gz", sep="\t", index=False)
    info["n_members"] = int(len(fam)); info["n_members_eq"] = int(fam.is_eq.sum())
    print("members", info, flush=True)

    eqm = fam[fam.is_eq].set_index("gene_id")
    need = set(eqm.index)
    seqs = {}
    n_eq_all = 0
    with open(R / "genecall/harmonized_proteins.faa") as fh, open(a.out / "eq_all_proteins.faa", "w") as fo:
        gid = None; keep = False; buf = []

        def flush():
            nonlocal n_eq_all
            if gid is not None and keep:
                s = "".join(buf).rstrip("*")
                fo.write(f">{gid}\n{s}\n"); n_eq_all += 1
                if gid in need:
                    seqs[gid] = s
        for line in fh:
            if line.startswith(">"):
                flush()
                gid = line[1:].split()[0]; buf = []
                keep = gid.split("::", 1)[0] in eq_genomes
            elif keep:
                buf.append(line.strip())
        flush()
    info["n_eq_proteins_all"] = n_eq_all
    assert set(seqs) == need, (len(seqs), len(need))
    with open(a.out / "eq_members_full.faa", "w") as fo:
        for g in eqm.index:
            fo.write(f">{g}\n{seqs[g]}\n")
    h2seq = {}
    for g in eqm.index:
        h2seq.setdefault(seqhash(seqs[g]), seqs[g])
    info["n_unique_eq_member_sequences"] = len(h2seq)

    reused = {}; n_tm = 0
    for f in sorted(glob.glob(str(a.tmbed_dir / "chunk_*.pred.gz"))):
        with gzip.open(f, "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline().strip().rstrip("*"); lab = fh.readline().strip()
                n_tm += 1
                k = seqhash(s)
                if k in h2seq and h2seq[k] == s:
                    reused[k] = lab[:len(s)]
    info["tmbed_entries_read"] = n_tm
    info["tmbed_reused_unique_sequences"] = len(reused)
    pd.DataFrame({"seqhash": list(reused), "label": list(reused.values())}).to_csv(a.out / "tmbed_reused.tsv.gz", sep="\t", index=False)
    sp_hash = {seqhash(seqs[g]) for g in eqm.index[eqm.signalp_class != "OTHER"]}
    todo = sorted(h for h in sp_hash if h not in reused)
    rng = np.random.default_rng(20261111)
    cov = sorted(h for h in sp_hash if h in reused)
    qc = list(rng.choice(cov, min(a.qc_n, len(cov)), replace=False)) if cov else []
    open(a.out / "tmbed_qc_hashes.txt", "w").write("\n".join(qc) + "\n")
    allh = sorted(todo + qc, key=lambda h: -len(h2seq[h]))
    info["tmbed_sp_pos_unique_sequences"] = len(sp_hash)
    info["tmbed_todo_unique_sequences"] = len(todo); info["tmbed_qc_sequences"] = len(qc)
    info["tmbed_todo_residues"] = int(sum(len(h2seq[h]) for h in allh))
    info["tmbed_todo_max_len"] = int(max(len(h2seq[h]) for h in allh)) if allh else 0
    pd_ = a.out / "tmbed_parts"; pd_.mkdir(exist_ok=True)
    outs = [open(pd_ / f"part_{i:03d}.faa", "w") for i in range(a.parts)]
    for i, h in enumerate(allh):
        outs[i % a.parts].write(f">{h}\n{h2seq[h]}\n")
    for o in outs:
        o.close()
    info["tmbed_parts"] = a.parts
    json.dump(info, open(a.out / "prep_info.json", "w"), indent=1, default=str)
    print("PREP_DONE", json.dumps(info, default=str), flush=True)


if __name__ == "__main__":
    main()
