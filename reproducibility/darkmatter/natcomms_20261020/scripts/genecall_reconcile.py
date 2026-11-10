#!/usr/bin/env python3
"""Reconcile species-representative gene calls with the assembly catalogue calls.

Revision item 3b. The 5,438 representative proteomes were re-called per MAG
(pyrodigal, run_type=Single); the assembly-wide catalogue was called per sample
assembly (pyrodigal, run_type=Metagenomic) and clustered at 90/80
(all_linclust90). Darkness (no UniRef50 hit, darkmatter_clean) is a property of
catalogue cluster representatives. Gene identifiers are therefore not assumed
interchangeable; each representative gene is joined to the catalogue

  (1) by locus: same sample contig, same strand, overlapping coordinates
      (identical / same stop, different start / same start, different stop /
      other overlap / no counterpart / contig absent), and
  (2) by exact protein sequence (MD5 of the residue string, terminal '*' removed)
      across all catalogue calls of all samples.

The catalogue darkness label of the counterpart's cluster is then propagated
and dark fractions are stratified by completeness, length, MIMAG tier and
phylum (Thermoplasmatota vs rest).

Modes (all run on Slurm workers):
  prep          parse representative proteomes -> work/rep_genes.tsv.gz, shards.tsv
  shard K       stream the catalogue calls of the samples in shard K
  merge         cluster membership + darkness join, summaries, run JSON
Set RECON_ROOT to choose the output root and RECON_TEST=1 for the small test.
"""
import gzip
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
CAT = W / "04_protein_catalog"
ASM = CAT / "assembly_proteins"
CLUSTER = CAT / "all_linclust90/all_cluster90.tsv"
DARK_ALL = CAT / "darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz"
MAG_DARK = CAT / "darkmatter_clean/05_mag_dark/mag_dark_rep_ids.txt.gz"
NONDARK_MAG = W / "06_darkmatter/darkfam_basis_clean/results/nondark_magmapped_ids.txt"
IDX = W / "06_darkmatter/heat_stress_20261016/results/genome_index.tsv"
TIERS = W / "mimag_rrna_trna_clean/results/mimag_tiers.tsv"

ROOT = Path(os.environ.get("RECON_ROOT", W / "06_darkmatter/genecall_reconcile_20261102"))
TEST = os.environ.get("RECON_TEST") == "1"
TEST_SAMPLES = ["10Dr2", "12PRr3", "11Dr3"]
NSHARD = int(os.environ.get("RECON_NSHARD", "48"))
WORK, RES, BULK = ROOT / "work", ROOT / "results", ROOT / "bulk"

LEN_BINS = [0, 100, 150, 250, 400, 10 ** 9]
LEN_LABELS = ["<100", "100-149", "150-249", "250-399", ">=400"]
HDR = re.compile(r"^(\S+) # (\d+) # (\d+) # (-?1) # .*?partial=(\d\d)")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_async(path):
    return subprocess.Popen(["sha256sum", str(path)], stdout=subprocess.PIPE, text=True)


def sha256_wait(p):
    out, _ = p.communicate()
    assert p.returncode == 0
    return out.split()[0]


def fasta(path):
    """Yield (header, sequence) from a gzipped FASTA, decompressing in a child gzip."""
    proc = subprocess.Popen(["gzip", "-dc", str(path)], stdout=subprocess.PIPE, bufsize=1 << 22)
    header, seq = None, []
    for line in proc.stdout:
        if line[:1] == b">":
            if header is not None:
                yield header, b"".join(seq).rstrip(b"*")
            header, seq = line[1:].rstrip().decode(), []
        else:
            seq.append(line.strip())
    if header is not None:
        yield header, b"".join(seq).rstrip(b"*")
    proc.stdout.close()
    assert proc.wait() == 0, f"gzip failed on {path}"


def gff_seqlens(path, wanted=None):
    lens = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                d = dict(kv.split("=", 1) for kv in line[17:].strip().split(";"))
                c = d["seqhdr"].strip('"').split()[0]
                if wanted is None or c in wanted:
                    lens[c] = int(d["seqlen"])
    return lens


# ----------------------------------------------------------------------------- prep
def prep_worker(job):
    gid, faa = job
    faa = Path(faa)
    gff = faa.with_name(faa.name[:-len(".faa.gz")] + ".gff.gz")
    rows = []
    for header, seq in fasta(faa):
        m = HDR.match(header)
        gene, start, end, strand, partial = m.groups()
        local = gene.split("::", 1)[1]
        contig = local.rsplit("_", 1)[0]
        sample = contig.rsplit("_", 1)[0]
        rows.append((gid, gene, sample, contig, int(start), int(end), int(strand), partial,
                     len(seq), int(b"*" in seq), hashlib.md5(seq).hexdigest()))
    lens = gff_seqlens(gff)
    return gid, rows, lens, dict(faa=str(faa), faa_sha256=sha256(faa), faa_bytes=faa.stat().st_size,
                                 gff=str(gff), gff_sha256=sha256(gff))


def prep():
    for d in (WORK, RES, BULK):
        d.mkdir(parents=True, exist_ok=True)
    idx = pd.read_csv(IDX, sep="\t").set_index("genome_id")
    tiers = pd.read_csv(TIERS, sep="\t").set_index("rep")
    assert len(idx) == len(tiers) == 5438 and set(idx.index) == set(tiers.index)
    if TEST:
        idx = idx[idx.source_sample.isin(TEST_SAMPLES)]
    jobs = [(g, r.proteome) for g, r in idx.sort_index().iterrows()]
    ncpu = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
    manifest, ngenes, mismatch = [], 0, 0
    src = idx.source_sample.to_dict()
    with gzip.open(WORK / "rep_genes.tsv.gz.part", "wt", compresslevel=3) as g, \
            gzip.open(WORK / "rep_contigs.tsv.gz.part", "wt", compresslevel=3) as c:
        g.write("genome\tgene\tsample\tcontig\tstart\tend\tstrand\tpartial\tlength\tinternal_stop\tmd5\n")
        c.write("genome\tcontig\tmag_len\n")
        with ProcessPoolExecutor(max_workers=ncpu) as pool:
            for i, (gid, rows, lens, man) in enumerate(pool.map(prep_worker, jobs, chunksize=4)):
                for r in rows:
                    g.write("\t".join(map(str, r)) + "\n")
                    mismatch += r[2] != src[gid]
                for contig, ln in lens.items():
                    c.write(f"{gid}\t{contig}\t{ln}\n")
                ngenes += len(rows)
                manifest.append(dict(genome=gid, n_genes=len(rows), n_contigs=len(lens), **man))
                if i % 500 == 0:
                    print("prep genomes", i, "genes", ngenes, flush=True)
    os.replace(WORK / "rep_genes.tsv.gz.part", WORK / "rep_genes.tsv.gz")
    os.replace(WORK / "rep_contigs.tsv.gz.part", WORK / "rep_contigs.tsv.gz")
    pd.DataFrame(manifest).to_csv(WORK / "proteome_manifest.tsv", sep="\t", index=False)
    # size-balanced sample shards over every catalogue sample (needed for the global sequence join)
    samples = sorted((p.parent.name, p.stat().st_size) for p in ASM.glob("*/*.faa.gz"))
    if TEST:
        samples = [s for s in samples if s[0] in TEST_SAMPLES]
    nshard = 1 if TEST else NSHARD
    load = [0] * nshard
    assign = []
    for name, size in sorted(samples, key=lambda x: -x[1]):
        k = load.index(min(load))
        load[k] += size
        assign.append((name, k, size))
    pd.DataFrame(assign, columns=["sample", "shard", "faa_bytes"]).to_csv(WORK / "shards.tsv", sep="\t", index=False)
    print(json.dumps(dict(genomes=len(jobs), genes=ngenes, gene_sample_ne_source_sample=mismatch,
                          catalogue_samples=len(samples), shards=nshard)), flush=True)


# ---------------------------------------------------------------------------- shard
def classify(r, hits, md5):
    """Best catalogue counterpart for representative gene r among same-contig catalogue calls."""
    s, e, st = r.start, r.end, r.strand
    best, best_rank, best_ov, nov = None, 9, -1, 0
    for cs, ce, cst, cpart, cid, cmd5 in hits:
        if cst != st or ce < s or cs > e:
            continue
        nov += 1
        ov = min(e, ce) - max(s, cs) + 1
        if cs == s and ce == e:
            rank = 0 if cmd5 == md5 else 1
        elif (st == 1 and ce == e) or (st == -1 and cs == s):
            rank = 2
        elif (st == 1 and cs == s) or (st == -1 and ce == e):
            rank = 3
        else:
            rank = 4
        if rank < best_rank or (rank == best_rank and ov > best_ov):
            best, best_rank, best_ov = (cs, ce, cpart, cid, cmd5), rank, ov
    names = ["identical", "same_coords_seq_differs", "same_stop_diff_start",
             "same_start_diff_stop", "overlap_other"]
    if best is None:
        return "no_counterpart", None, 0
    return names[best_rank], best, nov


def shard(k):
    shards = pd.read_csv(WORK / "shards.tsv", sep="\t")
    mine = list(shards["sample"][shards.shard == k])
    rep = pd.read_csv(WORK / "rep_genes.tsv.gz", sep="\t", dtype={"partial": str})
    repmd5 = set(bytes.fromhex(x) for x in rep.md5.unique())
    print(f"shard {k}: samples {len(mine)}; rep genes {len(rep)}; distinct rep md5 {len(repmd5)}", flush=True)
    rep = rep[rep["sample"].isin(mine)]
    maglen = pd.read_csv(WORK / "rep_contigs.tsv.gz", sep="\t").drop_duplicates("contig").set_index("contig").mag_len
    seqhit = {}  # md5 bytes -> [n catalogue copies, lexicographically smallest catalogue id]
    out_locus = gzip.open(WORK / f"shard_{k:03d}.locus.tsv.gz.part", "wt", compresslevel=3)
    out_locus.write("gene\tlocus_class\tn_overlap_same_strand\tcontig_in_catalogue\tmag_contig_len\t"
                    "asm_contig_len\tcat_id\tcat_start\tcat_end\tcat_partial\n")
    stats = []
    for sample in mine:
        faa = ASM / sample / f"{sample}.faa.gz"
        gff = ASM / sample / f"{sample}.gff.gz"
        h_faa, h_gff = sha256_async(faa), sha256_async(gff)
        r_s = rep[rep["sample"] == sample]
        need = {f"{sample}::{c}" for c in r_s.contig.unique()}
        bycontig = defaultdict(list)
        n = 0
        for header, seq in fasta(faa):
            n += 1
            m = HDR.match(header)
            cid, cs, ce, cst, cpart = m.groups()
            d = hashlib.md5(seq).digest()
            if d in repmd5:
                v = seqhit.get(d)
                if v is None:
                    seqhit[d] = [1, cid]
                else:
                    v[0] += 1
                    if cid < v[1]:
                        v[1] = cid
            contig = cid.rsplit("_", 1)[0]
            if contig in need:
                bycontig[contig].append((int(cs), int(ce), int(cst), cpart, cid, d.hex()))
        asm_len = gff_seqlens(gff, {c.split("::", 1)[1] for c in need}) if need else {}
        for r in r_s.itertuples(index=False):
            key = f"{sample}::{r.contig}"
            ml, al = maglen.get(r.contig, ""), asm_len.get(r.contig, "")
            if key not in bycontig and r.contig not in asm_len:
                out_locus.write(f"{r.gene}\tcontig_absent\t0\t0\t{ml}\t{al}\t\t\t\t\n")
                continue
            cls, best, nov = classify(r, bycontig.get(key, []), r.md5)
            if best is None:
                out_locus.write(f"{r.gene}\t{cls}\t{nov}\t1\t{ml}\t{al}\t\t\t\t\n")
            else:
                cs, ce, cpart, cid, _ = best
                out_locus.write(f"{r.gene}\t{cls}\t{nov}\t1\t{ml}\t{al}\t{cid}\t{cs}\t{ce}\t{cpart}\n")
        stats.append(dict(sample=sample, catalogue_proteins=n, rep_genes=len(r_s), rep_contigs=len(need),
                          faa=str(faa), faa_sha256=sha256_wait(h_faa), gff=str(gff), gff_sha256=sha256_wait(h_gff)))
        print(f"shard {k}: {sample} proteins {n} rep genes {len(r_s)} seqhits {len(seqhit)}", flush=True)
    out_locus.close()
    with gzip.open(WORK / f"shard_{k:03d}.seq.tsv.gz.part", "wt", compresslevel=3) as f:
        f.write("md5\tn_catalogue\texample_id\n")
        for d, (cnt, cid) in seqhit.items():
            f.write(f"{d.hex()}\t{cnt}\t{cid}\n")
    pd.DataFrame(stats).to_csv(WORK / f"shard_{k:03d}.samples.tsv", sep="\t", index=False)
    for ext in ("locus.tsv.gz", "seq.tsv.gz"):
        os.replace(WORK / f"shard_{k:03d}.{ext}.part", WORK / f"shard_{k:03d}.{ext}")
    print(f"shard {k} done", flush=True)


# ---------------------------------------------------------------------------- merge
def read_ids(path, wanted):
    hit = set()
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rb") as f:
        for line in f:
            v = line.rstrip(b"\n")
            if v in wanted:
                hit.add(v)
    return hit


def frac(n, d):
    return round(n / d, 6) if d else float("nan")


def summarise(df, keys):
    g = df.groupby(keys, observed=True, dropna=False)
    out = g.agg(n_genes=("gene", "size"),
                n_labelled=("dark", lambda x: x.notna().sum()),
                n_dark=("dark", lambda x: (x == 1).sum()),
                n_labelled_strict=("dark_strict", lambda x: x.notna().sum()),
                n_dark_strict=("dark_strict", lambda x: (x == 1).sum()),
                n_labelled_locus=("dark_locus", lambda x: x.notna().sum()),
                n_dark_locus=("dark_locus", lambda x: (x == 1).sum()),
                n_identical=("locus_class", lambda x: (x == "identical").sum())).reset_index()
    for a, b, c in [("dark_fraction", "n_dark", "n_labelled"),
                    ("dark_fraction_strict", "n_dark_strict", "n_labelled_strict"),
                    ("dark_fraction_locus", "n_dark_locus", "n_labelled_locus"),
                    ("labelled_fraction", "n_labelled", "n_genes"),
                    ("identical_fraction", "n_identical", "n_genes")]:
        out[a] = [frac(x, y) for x, y in zip(out[b], out[c])]
    return out


def merge():
    shards = pd.read_csv(WORK / "shards.tsv", sep="\t")
    nshard = shards.shard.nunique()
    missing = [k for k in range(nshard) if not (WORK / f"shard_{k:03d}.seq.tsv.gz").exists()]
    assert not missing, f"missing shards {missing}"
    h_cluster = sha256_async(CLUSTER)
    h_ids = {str(p): sha256_async(p) for p in (DARK_ALL, MAG_DARK, NONDARK_MAG, IDX, TIERS)}
    rep = pd.read_csv(WORK / "rep_genes.tsv.gz", sep="\t", dtype={"partial": str})
    locus = pd.concat([pd.read_csv(WORK / f"shard_{k:03d}.locus.tsv.gz", sep="\t",
                                   dtype={"cat_partial": str, "cat_id": str}) for k in range(nshard)])
    seq = pd.concat([pd.read_csv(WORK / f"shard_{k:03d}.seq.tsv.gz", sep="\t") for k in range(nshard)])
    seq = seq.groupby("md5").agg(n_catalogue=("n_catalogue", "sum"), example_id=("example_id", "min")).reset_index()
    samples = pd.concat([pd.read_csv(WORK / f"shard_{k:03d}.samples.tsv", sep="\t") for k in range(nshard)])
    assert locus.gene.is_unique and len(locus) == len(rep), (len(locus), len(rep))
    df = rep.merge(locus, on="gene", how="left").merge(seq, on="md5", how="left")
    df["n_catalogue"] = df.n_catalogue.fillna(0).astype(int)
    print("joined", len(df), flush=True)

    need = set(df.cat_id.dropna().str.encode("ascii")) | set(df.example_id.dropna().str.encode("ascii"))
    print("catalogue members needed", len(need), flush=True)
    m2r = {}
    n = 0
    with open(CLUSTER, "rb", buffering=1 << 24) as f:
        for line in f:
            n += 1
            t = line.find(b"\t")
            mem = line[t + 1:].rstrip(b"\n")
            if mem in need:
                m2r[mem] = line[:t]
            if n % 100_000_000 == 0:
                print("cluster lines", n, "found", len(m2r), flush=True)
    print("cluster lines", n, "members found", len(m2r), "of", len(need), flush=True)
    reps = set(m2r.values())
    dark = read_ids(DARK_ALL, reps)
    linked = read_ids(MAG_DARK, reps) | read_ids(NONDARK_MAG, reps)
    print("needed reps", len(reps), "dark", len(dark), "genome_linked", len(linked), flush=True)

    def lab(ids):
        r = [m2r.get(x.encode()) if isinstance(x, str) else None for x in ids]
        return ([x.decode() if x else None for x in r],
                [(1.0 if x in dark else 0.0) if x else float("nan") for x in r],
                [(1.0 if x in linked else 0.0) if x else float("nan") for x in r])

    df["cat_cluster_rep"], df["dark_locus"], df["linked_locus"] = lab(df.cat_id)
    df["seq_cluster_rep"], df["dark_seq"], df["linked_seq"] = lab(df.example_id)
    ident = df.locus_class == "identical"
    df["seq_class"] = "no_exact_seq"
    df.loc[df.n_catalogue > 0, "seq_class"] = "exact_seq_elsewhere"
    df.loc[ident, "seq_class"] = "exact_seq_at_locus"
    # primary label: identical locus > exact sequence anywhere > same-strand overlapping locus call
    df["label_source"] = "none"
    df.loc[df.dark_locus.notna(), "label_source"] = "locus_variant"
    df.loc[df.dark_seq.notna() & ~ident, "label_source"] = "exact_seq"
    df.loc[ident & df.dark_locus.notna(), "label_source"] = "identical"
    df["dark"] = df.dark_locus.where(df.label_source.isin(["identical", "locus_variant"]), df.dark_seq)
    df.loc[df.label_source == "none", "dark"] = float("nan")
    df["dark_strict"] = df.dark.where(df.label_source.isin(["identical", "exact_seq"]))
    df["genome_linked_cluster"] = df.linked_locus.where(df.label_source.isin(["identical", "locus_variant"]),
                                                        df.linked_seq)
    df["locus_seq_conflict"] = (df.dark_locus.notna() & df.dark_seq.notna() & (df.dark_locus != df.dark_seq)).astype(int)

    idx = pd.read_csv(IDX, sep="\t").set_index("genome_id")
    tiers = pd.read_csv(TIERS, sep="\t").set_index("rep")
    phylum = idx.gtdb_classification.str.extract(r"p__([^;]*)")[0]
    df["phylum"] = df.genome.map(phylum)
    df["phylum_group"] = (df.phylum == "Thermoplasmatota").map({True: "Thermoplasmatota", False: "rest"})
    df["mimag_tier"] = df.genome.map(tiers.mimag_tier)
    df["completeness"] = (df.partial == "00").map({True: "complete", False: "partial"})
    df["length_bin"] = pd.cut(df.length, LEN_BINS, right=False, labels=LEN_LABELS)
    df["contig_len_match"] = (pd.to_numeric(df.mag_contig_len) == pd.to_numeric(df.asm_contig_len)).astype(int)

    cols = ["genome", "gene", "contig", "start", "end", "strand", "partial", "length", "internal_stop",
            "md5", "locus_class", "n_overlap_same_strand", "contig_in_catalogue", "contig_len_match",
            "cat_id", "cat_start", "cat_end", "cat_partial", "cat_cluster_rep", "seq_class", "n_catalogue",
            "example_id", "seq_cluster_rep", "label_source", "dark", "dark_strict", "dark_locus", "dark_seq",
            "locus_seq_conflict", "genome_linked_cluster"]
    df[cols].to_csv(BULK / "rep_gene_reconciliation.tsv.gz", sep="\t", index=False,
                    compression={"method": "gzip", "compresslevel": 3})

    # mapping-class summary
    mc = df.groupby(["locus_class", "seq_class", "label_source", "completeness"], observed=True).agg(
        n_genes=("gene", "size"), n_dark=("dark", lambda x: (x == 1).sum()),
        n_conflict=("locus_seq_conflict", "sum")).reset_index()
    mc["fraction_of_genes"] = (mc.n_genes / len(df)).round(6)
    mc.to_csv(RES / "mapping_class_summary.tsv", sep="\t", index=False)

    # per-stratum summary (long): marginals, completeness x length, and the full cube
    dims = ["completeness", "length_bin", "mimag_tier", "phylum_group"]
    parts = []
    for keys in [[d] for d in dims] + [["completeness", "length_bin"],
                                       ["phylum_group", "completeness", "length_bin"],
                                       ["mimag_tier", "completeness", "length_bin"], dims]:
        s = summarise(df, keys)
        s.insert(0, "stratification", "x".join(keys))
        parts.append(s)
    tot = summarise(df.assign(all="all"), ["all"])
    tot.insert(0, "stratification", "all")
    pd.concat([tot] + parts, ignore_index=True).reindex(
        columns=["stratification"] + ["all"] + dims + [c for c in tot.columns if c not in ("stratification", "all")]
    ).drop(columns="all").to_csv(RES / "dark_fraction_strata.tsv", sep="\t", index=False)

    # per-genome table
    pg = summarise(df, ["genome"]).set_index("genome")
    cg = summarise(df[(df.completeness == "complete") & (df.length >= 100)], ["genome"]).set_index("genome")
    pg = pg.join(cg[["n_genes", "n_labelled", "n_dark", "dark_fraction"]].add_prefix("complete_ge100_"))
    cls = pd.crosstab(df.genome, df.locus_class).add_prefix("locus_")
    src = pd.crosstab(df.genome, df.label_source).add_prefix("label_")
    pg = pg.join(cls).join(src)
    pg["n_partial"] = df.groupby("genome").completeness.apply(lambda x: (x == "partial").sum())
    pg["contig_len_mismatch_genes"] = df.groupby("genome").contig_len_match.apply(lambda x: (x == 0).sum())
    pg = pg.join(idx[["gtdb_classification", "completeness", "contamination", "source_sample"]]
                 .rename(columns={"completeness": "checkm_completeness"}))
    pg["phylum_group"] = pg.index.map(lambda g: "Thermoplasmatota" if phylum[g] == "Thermoplasmatota" else "rest")
    pg = pg.join(tiers[["mimag_tier"]])
    pg.to_csv(RES / "per_genome_reconciliation.tsv", sep="\t")

    pman = pd.read_csv(WORK / "proteome_manifest.tsv", sep="\t")
    samples.to_csv(RES / "catalogue_sample_manifest.tsv", sep="\t", index=False)
    pman.to_csv(RES / "proteome_manifest.tsv", sep="\t", index=False)
    inputs = {str(CLUSTER): sha256_wait(h_cluster)}
    inputs.update({k: sha256_wait(p) for k, p in h_ids.items()})
    inputs["proteome_manifest.tsv"] = sha256(RES / "proteome_manifest.tsv")
    inputs["catalogue_sample_manifest.tsv"] = sha256(RES / "catalogue_sample_manifest.tsv")
    vc = lambda s: {str(k): int(v) for k, v in s.value_counts(dropna=False).items()}
    report = dict(job=os.environ.get("SLURM_JOB_ID"), host=platform.node(), test=TEST,
                  python=sys.version.split()[0], pandas=pd.__version__,
                  genomes=int(df.genome.nunique()), genes=len(df),
                  catalogue_samples=len(samples), catalogue_proteins=int(samples.catalogue_proteins.sum()),
                  cluster_lines=n, members_needed=len(need), members_found=len(m2r),
                  locus_class=vc(df.locus_class), seq_class=vc(df.seq_class), label_source=vc(df.label_source),
                  contig_len_mismatch_genes=int((df.contig_len_match == 0).sum()),
                  locus_seq_conflicts=int(df.locus_seq_conflict.sum()),
                  genome_linked_cluster=vc(df.genome_linked_cluster),
                  dark_fraction_primary=frac(int((df.dark == 1).sum()), int(df.dark.notna().sum())),
                  dark_fraction_strict=frac(int((df.dark_strict == 1).sum()), int(df.dark_strict.notna().sum())),
                  input_sha256=inputs,
                  rules=dict(locus="same sample contig, same strand, overlapping CDS; best = identical coords+seq > "
                                   "identical coords > same stop > same start > largest overlap",
                             sequence="MD5 of residues (terminal * removed) against every catalogue call of every sample",
                             primary_label="identical locus > exact sequence elsewhere > overlapping locus call",
                             strict_label="identical locus or exact sequence only",
                             darkness="darkmatter_clean dark_rep_ids (no UniRef50 hit) of the counterpart's 90/80 cluster rep",
                             length_bins=LEN_LABELS, completeness="partial=00 complete, else partial"),
                  scope="Representative gene calls reconciled to assembly catalogue calls; labels propagated from catalogue clusters")
    (RES / "reconcile_run.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    assert os.environ.get("SLURM_JOB_ID"), "Run on an allocated Slurm worker"
    mode = sys.argv[1]
    if mode == "prep":
        prep()
    elif mode == "shard":
        shard(int(sys.argv[2]) if len(sys.argv) > 2 else int(os.environ["SLURM_ARRAY_TASK_ID"]))
    elif mode == "merge":
        merge()
    else:
        raise SystemExit(f"unknown mode {mode}")
