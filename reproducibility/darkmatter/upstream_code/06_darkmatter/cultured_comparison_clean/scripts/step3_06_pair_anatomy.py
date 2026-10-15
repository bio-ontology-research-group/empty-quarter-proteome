#!/usr/bin/env python3
"""step3_06 - anatomy of the 23 same-species MAG / cultured pairs.

EXPLORATORY. Nothing here is registered.

For every one of the 23 MAG species representatives that also exists as a cultured
long-read genome (results/step3_paired_by_species.tsv), compare the two proteomes
protein by protein and ask three things:

  1. do the two proteomes overlap more for annotated than for dark proteins?
  2. does the MAG-unique remainder look like short-read noise (short, partial,
     contig-edge) or like real accessory content?
  3. does the culture-unique remainder look like accessory genes?

Keys and sources
  MAG side      catalogue key `<sample>::<contig>_<gene>`. The proteins of a
                representative are the pyrodigal catalogue proteins of its sample,
                M/04_protein_catalog/assembly_proteins/<sample>/<sample>.faa.gz,
                restricted to the contigs that contig_to_bin_refined.tsv.gz assigns
                to that representative. Partial flag and coordinates come from the
                pyrodigal header, contig length from the dRep winner FASTA.
  MAG darkness  protein -> 90 % linclust representative (all_cluster90.tsv, 753 M
                lines, scanned once in parallel byte ranges exactly as
                lineage_darkness_clean/scripts/v2_scan_cluster90.py does) -> dark if
                that representative is in dark_rep_ids.txt.gz. Darkness is therefore a
                cluster property, inherited by every member; that is the catalogue
                definition and it is what v2_rep_dark.tsv aggregates.
  cultured side Bakta proteins from tmp/cultured_all.faa, header `<genome>@<locus>`.
                Dark = no UniRef50 hit, read off tmp/cultured_dark_ids.txt (written by
                step3_04). Hypothetical = tmp/cultured_hypo_ids.txt. Product from the
                Bakta .tsv of the genome, in either Bakta root.

Pairing follows step3_04: every cultured genome is assigned to its best-ANI MAG
representative, so the 189 shared-species genomes collapse onto 23 representatives.
The primary comparison uses the single best-ANI cultured genome of each
representative; the union comparison pools all cultured genomes assigned to it.

Overlap rule: DIAMOND 2.0.15 blastp both directions, --very-sensitive -e 1e-5 -k 5.
Shared = reciprocal best hit (by bit score, before filtering) with >= 50 % identity
and >= 50 % coverage of both sequences. The sensitivity rule is >= 30 % identity and
>= 70 % coverage of both.

Stages, each skipped when its outputs are already on disk:
  prepare   per-pair FASTAs, MAG per-protein table, contig lengths
  cluster   the all_cluster90 scan -> MAG per-protein dark flag
  diamond   46 + 46 small blastp runs
  analyse   the tables under results/pair_anatomy/

Usage:  step3_06_pair_anatomy.py [--stage all|prepare|cluster|diamond|analyse]
                                 [--workers 32] [--force]
"""
import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from multiprocessing import Pool, Process

P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean"
E = "/ibex/scratch/projects/c2014/alelopezv/eqmgc"
W = "/ibex/scratch/projects/c2014/rob/metagenomics"
A = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/cultures/annotation"
LIN = f"{W}/06_darkmatter/lineage_darkness_clean/results"

BAKTA_ROOTS = [f"{E}/results/annotation/bakta", f"{A}/thermophilic/bakta"]
F_C2B = f"{W}/04_protein_catalog/contig_to_bin_refined.tsv.gz"
F_CLUSTER = f"{W}/04_protein_catalog/all_linclust90/all_cluster90.tsv"
F_DARK = f"{W}/04_protein_catalog/darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz"
D_ASMPROT = f"{W}/04_protein_catalog/assembly_proteins"
D_DREP = f"{W}/05_refine/_drep/out/dereplicated_genomes"
DIAMOND = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"

TMP = f"{P}/tmp/pair_anatomy"
RES = f"{P}/results/pair_anatomy"

# overlap rules: (label, min identity, min coverage of both)
RULE_MAIN = ("main", 50.0, 50.0)
RULE_SENS = ("sens", 30.0, 70.0)
EDGE_BP = 200          # "at a contig edge" = start or stop within this many bp of one
OUTFMT = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
          "evalue", "bitscore", "qcovhsp", "scovhsp"]


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def median(xs):
    s = sorted(xs)
    n = len(s)
    if not n:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def quartiles(xs):
    s = sorted(xs)
    if not s:
        return float("nan"), float("nan"), float("nan")
    return median(s[: len(s) // 2]), median(s), median(s[(len(s) + 1) // 2:])


def pct(n, d):
    return 100.0 * n / d if d else float("nan")


# ------------------------------------------------------------------ the 23 pairs
def load_pairs():
    """rep -> dict(cultured=[(genome, ani, batch, tier)] best first, best=genome)."""
    reps = []
    with open(f"{P}/results/step3_paired_by_species.tsv") as fh:
        fh.readline()
        for line in fh:
            reps.append(line.split("\t")[0])
    rep_set = set(reps)

    # step3_04 convention: each cultured genome goes to its own best-ANI representative
    best_of = {}
    info = {}
    with open(f"{P}/results/shared_species.tsv") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, im, ia = h.index("cultured_genome"), h.index("mag_rep"), h.index("ani")
        ib, iq = h.index("batch"), h.index("quality_tier")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            g, m, a = f[ic], f[im], float(f[ia])
            if g not in best_of or a > best_of[g][1]:
                best_of[g] = (m, a)
                info[g] = (f[ib], f[iq])
    pairs = {r: {"cultured": []} for r in reps}
    for g, (m, a) in best_of.items():
        if m in rep_set:
            pairs[m]["cultured"].append((g, a) + info[g])
    for r in reps:
        pairs[r]["cultured"].sort(key=lambda t: -t[1])
        pairs[r]["best"] = pairs[r]["cultured"][0][0]
        pairs[r]["ani"] = pairs[r]["cultured"][0][1]
    # MAG quality from the native-key rep table
    with open(f"{LIN}/v2_reps.tsv") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ig, ic2, ico, in50 = (h.index("genome"), h.index("completeness"),
                              h.index("contamination"), h.index("N50"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ig] in pairs:
                pairs[f[ig]].update(completeness=float(f[ic2]),
                                    contamination=float(f[ico]), n50=float(f[in50]))
    return reps, pairs


def bakta_dir(genome):
    for root in BAKTA_ROOTS:
        for d in (genome, genome + ".1"):
            if os.path.isdir(f"{root}/{d}"):
                return f"{root}/{d}", d
    return None, None


# ------------------------------------------------------------------ stage prepare
def contig_map(reps):
    """contig key '<sample>::<contig>' -> rep, for the 23 representatives only."""
    want = set(reps)
    c2r = {}
    with gzip.open(F_C2B, "rt") as fh:
        for line in fh:
            contig, _strat, b = line.rstrip("\n").split("\t")
            if b in want:
                c2r[contig] = b
    return c2r


def contig_lengths(rep):
    """contig key -> length, from the dRep winner FASTA."""
    out, name, n = {}, None, 0
    with open(f"{D_DREP}/{rep}.fa") as fh:
        for line in fh:
            if line[0] == ">":
                if name:
                    out[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name:
        out[name] = n
    return out


def extract_sample(args):
    """Pull one sample's catalogue proteins for the reps that live on its contigs."""
    sample, c2r_sample = args
    faa = f"{D_ASMPROT}/{sample}/{sample}.faa.gz"
    handles, rows = {}, []
    pr = subprocess.Popen(["zcat", faa], stdout=subprocess.PIPE, text=True,
                          bufsize=1024 * 1024)
    keep = None
    for line in pr.stdout:
        if line[0] == ">":
            head = line[1:].rstrip("\n")
            f = head.split(" # ")
            pid = f[0]
            rep = c2r_sample.get(pid[: pid.rfind("_")])
            if rep is None:
                keep = None
                continue
            if rep not in handles:
                handles[rep] = open(f"{TMP}/mag/{rep}.faa", "w")
            handles[rep].write(f">{pid}\n")
            keep = handles[rep]
            m = re.search(r"partial=(\d\d)", f[4]) if len(f) > 4 else None
            rows.append([rep, pid, pid[: pid.rfind("_")], f[1], f[2], f[3],
                         m.group(1) if m else "NA", 0])
        elif keep is not None:
            keep.write(line)
            rows[-1][7] += len(line.strip())
    pr.wait()
    for h in handles.values():
        h.close()
    return sample, rows


def stage_prepare(reps, pairs, workers):
    os.makedirs(f"{TMP}/mag", exist_ok=True)
    os.makedirs(f"{TMP}/cult", exist_ok=True)
    os.makedirs(f"{TMP}/pool", exist_ok=True)

    c2r = contig_map(reps)
    log(f"prepare: {len(c2r)} contigs on the 23 representatives")

    clen = {}
    for r in reps:
        for c, n in contig_lengths(r).items():
            clen[f"{r.split('__')[0]}::{c}"] = n
    missing = [c for c in c2r if c not in clen]
    if missing:
        sys.exit(f"ERROR: {len(missing)} binned contigs absent from the dRep FASTAs, "
                 f"e.g. {missing[:3]}")
    log(f"prepare: contig lengths for {len(clen)} contigs")

    by_sample = defaultdict(dict)
    for c, r in c2r.items():
        by_sample[c.split("::")[0]][c] = r
    with Pool(min(workers, len(by_sample))) as pool:
        allrows = []
        for sample, rows in pool.imap_unordered(extract_sample, list(by_sample.items())):
            log(f"prepare: {sample}: {len(rows)} catalogue proteins on winner contigs")
            allrows += rows
    with open(f"{TMP}/mag_proteins.tsv", "w") as out:
        out.write("rep\tprotein\tcontig\tstart\tstop\tstrand\tpartial\taa\t"
                  "contig_len\tedge_dist\n")
        for rep, pid, contig, s, e, strand, partial, aa in allrows:
            L = clen[contig]
            out.write(f"{rep}\t{pid}\t{contig}\t{s}\t{e}\t{strand}\t{partial}\t{aa}\t"
                      f"{L}\t{min(int(s) - 1, L - int(e))}\n")
    log(f"prepare: wrote {len(allrows)} MAG protein records")

    # cultured proteomes, one pass over the concatenated Bakta file
    want = {g for r in reps for g, *_ in pairs[r]["cultured"]}
    handles, keep = {}, None
    with open(f"{P}/tmp/cultured_all.faa") as fh:
        for line in fh:
            if line[0] == ">":
                pid = line[1:].split()[0]
                g = pid.split("@")[0]
                if g not in want:
                    keep = None
                    continue
                if g not in handles:
                    handles[g] = open(f"{TMP}/cult/{g}.faa", "w")
                keep = handles[g]
                keep.write(f">{pid}\n")
            elif keep is not None:
                keep.write(line)
    for h in handles.values():
        h.close()
    log(f"prepare: wrote {len(handles)} cultured proteomes of {len(want)} wanted")
    if len(handles) != len(want):
        sys.exit("ERROR: a shared-species cultured genome has no proteome")

    for r in reps:
        with open(f"{TMP}/pool/{r}.faa", "w") as out:
            for g, *_ in pairs[r]["cultured"]:
                with open(f"{TMP}/cult/{g}.faa") as fh:
                    out.write(fh.read())
    log("prepare: pooled cultured proteomes written")


# ------------------------------------------------------------------ stage cluster
def cluster_worker(k, start, end, contigs):
    """Members of the 23 reps' contigs, from one byte range of all_cluster90.tsv."""
    out = open(f"{TMP}/cl_part_{k:02d}.tsv", "wb")
    n = hit = 0
    with open(F_CLUSTER, "rb") as fh:
        fh.seek(start)
        if start > 0:
            fh.readline()                 # the previous worker owns this partial line
        pos = fh.tell()
        has = contigs.__contains__
        while pos < end:
            line = fh.readline()
            if not line:
                break
            pos += len(line)
            n += 1
            t = line.find(b"\t")
            if t < 0:
                continue
            member = line[t + 1:-1]
            if has(member[: member.rfind(b"_")]):
                hit += 1
                out.write(member + b"\t" + line[:t] + b"\n")
    out.close()
    with open(f"{TMP}/cl_stat_{k:02d}.txt", "w") as fh:
        fh.write(f"{n}\t{hit}\n")


def stage_cluster(reps, workers):
    contigs = {c.encode() for c in contig_map(reps)}
    log(f"cluster: scanning {os.path.getsize(F_CLUSTER) / 2**30:.1f} GiB with "
        f"{workers} workers for {len(contigs)} contigs")
    size = os.path.getsize(F_CLUSTER)
    step = size // workers
    procs = []
    for k in range(workers):
        s = k * step
        e = size if k == workers - 1 else (k + 1) * step
        p = Process(target=cluster_worker, args=(k, s, e, contigs))
        p.start()
        procs.append(p)
    for p in procs:
        p.join()
    n = hit = 0
    for k in range(workers):
        a, b = open(f"{TMP}/cl_stat_{k:02d}.txt").read().split()
        n, hit = n + int(a), hit + int(b)
    log(f"cluster: {n} lines scanned, {hit} members on a paired representative")

    p2c = {}
    for k in range(workers):
        with open(f"{TMP}/cl_part_{k:02d}.tsv", "rb") as fh:
            for line in fh:
                m, r = line[:-1].split(b"\t")
                p2c[m] = r
    touched = set(p2c.values())
    log(f"cluster: {len(p2c)} proteins, {len(touched)} distinct 90 % representatives")

    dark = set()
    pr = subprocess.Popen(["zcat", F_DARK], stdout=subprocess.PIPE)
    for line in pr.stdout:
        v = line[:-1]
        if v in touched:
            dark.add(v)
    pr.wait()
    log(f"cluster: {len(dark)} of those representatives are dark")

    with open(f"{TMP}/mag_protein_dark.tsv", "w") as out:
        out.write("protein\tcluster90_rep\tdark\n")
        for m, r in p2c.items():
            out.write(f"{m.decode()}\t{r.decode()}\t{1 if r in dark else 0}\n")
    for k in range(workers):
        os.remove(f"{TMP}/cl_part_{k:02d}.tsv")
        os.remove(f"{TMP}/cl_stat_{k:02d}.txt")


# ------------------------------------------------------------------ stage diamond
def run_blastp(q, d, out, threads):
    if os.path.exists(out) and os.path.getsize(out) >= 0 and os.path.exists(out + ".ok"):
        return
    dmnd = d + ".dmnd"
    if not os.path.exists(dmnd):
        subprocess.run([DIAMOND, "makedb", "--in", d, "-d", d, "--quiet"], check=True)
    subprocess.run([DIAMOND, "blastp", "-q", q, "-d", dmnd, "-o", out + ".tmp",
                    "--very-sensitive", "-e", "1e-5", "-k", "5",
                    "--outfmt", "6"] + OUTFMT +
                   ["--threads", str(threads), "--quiet", "--tmpdir", TMP], check=True)
    os.rename(out + ".tmp", out)
    open(out + ".ok", "w").close()


def stage_diamond(reps, pairs, threads):
    os.makedirs(f"{TMP}/dia", exist_ok=True)
    for i, r in enumerate(reps, 1):
        g = pairs[r]["best"]
        mag, cult, pool = f"{TMP}/mag/{r}.faa", f"{TMP}/cult/{g}.faa", f"{TMP}/pool/{r}.faa"
        run_blastp(mag, cult, f"{TMP}/dia/{r}__m2c.m8", threads)
        run_blastp(cult, mag, f"{TMP}/dia/{r}__c2m.m8", threads)
        if len(pairs[r]["cultured"]) > 1:
            run_blastp(mag, pool, f"{TMP}/dia/{r}__m2p.m8", threads)
            run_blastp(pool, mag, f"{TMP}/dia/{r}__p2m.m8", threads)
        log(f"diamond: {i}/{len(reps)} {r} vs {g}")


# ------------------------------------------------------------------ stage analyse
def best_hits(path):
    """query -> (subject, pident, qcov, scov, bitscore); ties broken by e-value then id."""
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, s = f[0], f[1]
            pid, ev, bs = float(f[2]), float(f[6]), float(f[7])
            qc, sc = float(f[8]), float(f[9])
            key = (bs, -ev, pid)
            cur = best.get(q)
            if cur is None or key > cur[0]:
                best[q] = (key, s, pid, qc, sc)
    return {q: v[1:] for q, v in best.items()}


def all_hit_subjects(path):
    out = set()
    with open(path) as fh:
        for line in fh:
            out.add(line.split("\t", 1)[0])
    return out


def rbh(fwd, rev, min_id, min_cov):
    """MAG protein -> (cultured protein, pident) for reciprocal best hits that pass."""
    out = {}
    for q, (s, pid, qc, sc) in fwd.items():
        back = rev.get(s)
        if back is None or back[0] != q:
            continue
        if pid >= min_id and qc >= min_cov and sc >= min_cov:
            out[q] = (s, pid)
    return out


def load_products(genome):
    d, name = bakta_dir(genome)
    if d is None:
        return {}
    prod = {}
    with open(f"{d}/{name}.tsv") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 7 and f[5]:
                prod[f[5]] = f[7]
    return prod


def stage_analyse(reps, pairs):
    os.makedirs(RES, exist_ok=True)

    mag = defaultdict(dict)
    with open(f"{TMP}/mag_proteins.tsv") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rec = dict(zip(h, f))
            mag[rec["rep"]][rec["protein"]] = rec
    mdark = {}
    with open(f"{TMP}/mag_protein_dark.tsv") as fh:
        fh.readline()
        for line in fh:
            p, _c, d = line.rstrip("\n").split("\t")
            mdark[p] = d == "1"
    # every catalogue protein is a member of exactly one 90 % cluster, so a protein with
    # no flag would mean the scan missed it; a silent miss would read as "annotated".
    n_prot = sum(len(v) for v in mag.values())
    n_flag = sum(1 for r in mag for p in mag[r] if p in mdark)
    log(f"analyse: dark flag for {n_flag} of {n_prot} MAG proteins")
    if n_flag < n_prot:
        sys.exit(f"ERROR: {n_prot - n_flag} MAG proteins have no 90 % cluster flag")
    # cross-check against the aggregate native-key scan
    with open(f"{LIN}/v2_rep_dark.tsv") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ir, inp, ind = h.index("rep"), h.index("n_proteins"), h.index("n_dark")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ir] in mag:
                here = len(mag[f[ir]])
                hd = sum(1 for p in mag[f[ir]] if mdark[p])
                if here != int(f[inp]) or hd != int(f[ind]):
                    log(f"WARN {f[ir]}: here {here}/{hd}, v2_rep_dark "
                        f"{f[inp]}/{f[ind]}")

    cdark = set()
    with open(f"{P}/tmp/cultured_dark_ids.txt") as fh:
        for line in fh:
            cdark.add(line.strip())
    chypo = set()
    with open(f"{P}/tmp/cultured_hypo_ids.txt") as fh:
        for line in fh:
            chypo.add(line.strip())

    cult_len = {}
    rows, prod_all = [], Counter()
    detail_pick, per_pair_products = None, []
    for r in reps:
        g = pairs[r]["best"]
        mp = mag[r]
        cl = {}
        pid_now = None
        with open(f"{TMP}/cult/{g}.faa") as fh:
            for line in fh:
                if line[0] == ">":
                    pid_now = line[1:].strip()
                    cl[pid_now] = 0
                else:
                    cl[pid_now] += len(line.strip())
        cult_len[g] = cl

        fwd = best_hits(f"{TMP}/dia/{r}__m2c.m8")
        rev = best_hits(f"{TMP}/dia/{r}__c2m.m8")
        hit_any_mag = all_hit_subjects(f"{TMP}/dia/{r}__m2c.m8")
        hit_any_cult = all_hit_subjects(f"{TMP}/dia/{r}__c2m.m8")
        pair_main = rbh(fwd, rev, RULE_MAIN[1], RULE_MAIN[2])
        pair_sens = rbh(fwd, rev, RULE_SENS[1], RULE_SENS[2])
        shared_cult = {v[0] for v in pair_main.values()}

        n_mag, n_cult = len(mp), len(cl)
        md = [p for p in mp if mdark.get(p, False)]
        ma = [p for p in mp if not mdark.get(p, False)]
        cd = [p for p in cl if p in cdark]
        ca = [p for p in cl if p not in cdark]
        md_sh = sum(1 for p in md if p in pair_main)
        ma_sh = sum(1 for p in ma if p in pair_main)
        cd_sh = sum(1 for p in cd if p in shared_cult)
        ca_sh = sum(1 for p in ca if p in shared_cult)
        full = {p for p in mp if mp[p]["partial"] == "00" and int(mp[p]["aa"]) >= 100}
        ma_full = [p for p in ma if p in full]
        md_full = [p for p in md if p in full]
        ma_full_sh = sum(1 for p in ma_full if p in pair_main)
        md_full_sh = sum(1 for p in md_full if p in pair_main)

        mu = [p for p in mp if p not in pair_main]
        cu = [p for p in cl if p not in shared_cult]
        mu_len = [int(mp[p]["aa"]) for p in mu]
        mu_part = sum(1 for p in mu if mp[p]["partial"] not in ("00", "NA"))
        mu_edge = sum(1 for p in mu if int(mp[p]["edge_dist"]) <= EDGE_BP)
        mu_dark = sum(1 for p in mu if mdark.get(p, False))
        mu_nohit = sum(1 for p in mu if p not in hit_any_mag)
        sh_len = [int(mp[p]["aa"]) for p in pair_main]
        cu_len = [cl[p] for p in cu]
        cu_dark = sum(1 for p in cu if p in cdark)
        cu_hyp = sum(1 for p in cu if p in chypo)
        cu_nohit = sum(1 for p in cu if p not in hit_any_cult)
        ids = [v[1] for v in pair_main.values()]

        prods = load_products(g)
        cu_prod = Counter(prods.get(p.split("@", 1)[1], "NA") for p in cu
                          if p not in chypo)
        for k, v in cu_prod.items():
            if k not in ("NA", "hypothetical protein"):
                prod_all[k] += v
        per_pair_products.append(
            (r, g, "; ".join(f"{k} ({v})" for k, v in cu_prod.most_common(5)
                             if k not in ("NA", "hypothetical protein"))))

        rows.append(dict(
            mag_rep=r, cultured=g, ani=pairs[r]["ani"],
            n_cultured_genomes=len(pairs[r]["cultured"]),
            mag_completeness=pairs[r].get("completeness", float("nan")),
            mag_contamination=pairs[r].get("contamination", float("nan")),
            mag_n50=pairs[r].get("n50", float("nan")),
            n_mag=n_mag, mag_dark=len(md), mag_dark_pct=pct(len(md), n_mag),
            n_cult=n_cult, cult_dark=len(cd), cult_dark_pct=pct(len(cd), n_cult),
            shared=len(pair_main), shared_sens=len(pair_sens),
            shared_pct_mag=pct(len(pair_main), n_mag),
            shared_pct_cult=pct(len(pair_main), n_cult),
            mag_annot_with_partner_pct=pct(ma_sh, len(ma)),
            mag_dark_with_partner_pct=pct(md_sh, len(md)),
            cult_annot_with_partner_pct=pct(ca_sh, len(ca)),
            cult_dark_with_partner_pct=pct(cd_sh, len(cd)),
            mag_annot=len(ma), mag_annot_shared=ma_sh,
            mag_dark_shared=md_sh, cult_annot=len(ca), cult_annot_shared=ca_sh,
            cult_dark_shared=cd_sh,
            # dark MAG proteins are systematically short and often partial calls, so the
            # annotated-vs-dark gap is repeated on complete calls of >= 100 aa only
            mag_annot_full=len(ma_full), mag_annot_full_shared=ma_full_sh,
            mag_dark_full=len(md_full), mag_dark_full_shared=md_full_sh,
            mag_unique=len(mu), mag_unique_median_aa=median(mu_len),
            mag_unique_partial=mu_part, mag_unique_edge=mu_edge,
            mag_unique_dark=mu_dark, mag_unique_nohit=mu_nohit,
            cult_unique_dark=cu_dark, cult_unique_hypothetical=cu_hyp,
            cult_unique_nohit=cu_nohit,
            mag_unique_partial_pct=pct(mu_part, len(mu)),
            mag_unique_edge_pct=pct(mu_edge, len(mu)),
            mag_unique_dark_pct=pct(mu_dark, len(mu)),
            mag_unique_nohit_pct=pct(mu_nohit, len(mu)),
            shared_median_aa=median(sh_len),
            cult_unique=len(cu), cult_unique_median_aa=median(cu_len),
            cult_unique_dark_pct=pct(cu_dark, len(cu)),
            cult_unique_hypothetical_pct=pct(cu_hyp, len(cu)),
            cult_unique_nohit_pct=pct(cu_nohit, len(cu)),
            shared_median_identity=median(ids),
            shared_ge95_pct=pct(sum(1 for i in ids if i >= 95), len(ids)),
        ))

        # union over every cultured strain of the species
        if len(pairs[r]["cultured"]) > 1:
            pf = best_hits(f"{TMP}/dia/{r}__m2p.m8")
            pr_ = best_hits(f"{TMP}/dia/{r}__p2m.m8")
            up = rbh(pf, pr_, RULE_MAIN[1], RULE_MAIN[2])
            rows[-1]["union_shared"] = len(up)
            rows[-1]["union_shared_pct_mag"] = pct(len(up), n_mag)
            rows[-1]["union_mag_annot_with_partner_pct"] = pct(
                sum(1 for p in ma if p in up), len(ma))
            rows[-1]["union_mag_dark_with_partner_pct"] = pct(
                sum(1 for p in md if p in up), len(md))
        else:
            rows[-1]["union_shared"] = len(pair_main)
            rows[-1]["union_shared_pct_mag"] = pct(len(pair_main), n_mag)
            rows[-1]["union_mag_annot_with_partner_pct"] = pct(ma_sh, len(ma))
            rows[-1]["union_mag_dark_with_partner_pct"] = pct(md_sh, len(md))

        # detail pair: highest ANI among high-quality MAGs paired to a single strain
        if (pairs[r].get("completeness", 0) >= 90 and len(pairs[r]["cultured"]) == 1
                and (detail_pick is None or pairs[r]["ani"] > detail_pick[1])):
            detail_pick = (r, pairs[r]["ani"])
        log(f"analyse: {r} vs {g}: shared {len(pair_main)}, "
            f"MAG-unique {len(mu)}, culture-unique {len(cu)}")

    cols = list(rows[0].keys())
    with open(f"{RES}/pair_anatomy_per_pair.tsv", "w") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in sorted(rows, key=lambda d: -d["ani"]):
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v)
                        for k, v in row.items()})

    # ------------------------------------------------------------ cross-pair summary
    num_cols = [c for c in cols if isinstance(rows[0][c], (int, float))
                and c not in ("ani",)]
    with open(f"{RES}/pair_anatomy_summary.tsv", "w") as out:
        out.write("statistic\tn\tq1\tmedian\tq3\tmin\tmax\n")
        for c in ["ani"] + num_cols:
            vals = [r[c] for r in rows if r[c] == r[c]]
            q1, m, q3 = quartiles(vals)
            out.write(f"{c}\t{len(vals)}\t{q1:.4f}\t{m:.4f}\t{q3:.4f}\t"
                      f"{min(vals):.4f}\t{max(vals):.4f}\n")
        # pooled over all 23 pairs, so the big proteomes are not down-weighted
        def tot(k):
            return sum(r[k] for r in rows)
        for k in ("n_mag", "mag_dark", "n_cult", "cult_dark", "shared", "shared_sens",
                  "mag_unique", "cult_unique"):
            out.write(f"pooled_{k}\t23\tNA\t{tot(k)}\tNA\tNA\tNA\n")
        # pooled rates: the question "is an annotated protein likelier to have a partner
        # than a dark one" answered once over all 87,381 MAG proteins, not per pair
        for lbl, num, den in (
                ("mag_annot_with_partner", "mag_annot_shared", "mag_annot"),
                ("mag_dark_with_partner", "mag_dark_shared", "mag_dark"),
                ("cult_annot_with_partner", "cult_annot_shared", "cult_annot"),
                ("cult_dark_with_partner", "cult_dark_shared", "cult_dark"),
                ("mag_annot_ge100aa_complete_with_partner",
                 "mag_annot_full_shared", "mag_annot_full"),
                ("mag_dark_ge100aa_complete_with_partner",
                 "mag_dark_full_shared", "mag_dark_full"),
                ("mag_unique_partial", "mag_unique_partial", "mag_unique"),
                ("mag_unique_edge", "mag_unique_edge", "mag_unique"),
                ("mag_unique_dark", "mag_unique_dark", "mag_unique"),
                ("mag_unique_nohit", "mag_unique_nohit", "mag_unique"),
                ("cult_unique_dark", "cult_unique_dark", "cult_unique"),
                ("cult_unique_hypothetical", "cult_unique_hypothetical", "cult_unique"),
                ("cult_unique_nohit", "cult_unique_nohit", "cult_unique")):
            out.write(f"pooled_{lbl}_pct\t{tot(den)}\tNA\t{pct(tot(num), tot(den)):.4f}"
                      f"\tNA\tNA\tNA\n")
        n_dark_lower = sum(1 for r in rows
                           if r["mag_dark_with_partner_pct"] < r["mag_annot_with_partner_pct"])
        out.write(f"pairs_mag_dark_partner_rate_below_annotated\t23\tNA\t"
                  f"{n_dark_lower}\tNA\tNA\tNA\n")
        n_dark_lower_c = sum(1 for r in rows
                             if r["cult_dark_with_partner_pct"] < r["cult_annot_with_partner_pct"])
        out.write(f"pairs_cult_dark_partner_rate_below_annotated\t23\tNA\t"
                  f"{n_dark_lower_c}\tNA\tNA\tNA\n")

    with open(f"{RES}/pair_anatomy_cult_unique_products.tsv", "w") as out:
        out.write("product\tn_proteins_across_23_pairs\n")
        for k, v in prod_all.most_common(40):
            out.write(f"{k}\t{v}\n")
    with open(f"{RES}/pair_anatomy_cult_unique_products_per_pair.tsv", "w") as out:
        out.write("mag_rep\tcultured\ttop5_products\n")
        for r, g, s in per_pair_products:
            out.write(f"{r}\t{g}\t{s}\n")

    # ------------------------------------------------------------------ detail pair
    r = detail_pick[0]
    g = pairs[r]["best"]
    fwd = best_hits(f"{TMP}/dia/{r}__m2c.m8")
    rev = best_hits(f"{TMP}/dia/{r}__c2m.m8")
    pair_main = rbh(fwd, rev, RULE_MAIN[1], RULE_MAIN[2])
    shared_cult = {v[0] for v in pair_main.values()}
    hit_any_cult = all_hit_subjects(f"{TMP}/dia/{r}__c2m.m8")
    prods = load_products(g)
    cl = cult_len[g]
    cand = [(cl[p], p, prods.get(p.split("@", 1)[1], "NA")) for p in cl
            if p not in shared_cult and p not in hit_any_cult
            and p not in chypo and p not in cdark
            and prods.get(p.split("@", 1)[1], "NA") not in ("NA", "hypothetical protein")]
    cand.sort(reverse=True)
    with open(f"{RES}/pair_anatomy_detail_{r}.tsv", "w") as out:
        out.write(f"# detail pair: {r} vs {g}, ANI {pairs[r]['ani']}, "
                  f"MAG completeness {pairs[r].get('completeness')}, single cultured strain\n")
        out.write("rank\tcultured_protein\taa\tproduct\n")
        for i, (L, p, pd) in enumerate(cand[:10], 1):
            out.write(f"{i}\t{p}\t{L}\t{pd}\n")

    # per-protein dump for the detail pair, so the claims can be re-checked
    mp = mag[r]
    with open(f"{RES}/pair_anatomy_detail_{r}_mag_proteins.tsv", "w") as out:
        out.write("protein\taa\tpartial\tedge_dist\tcontig_len\tdark\tstatus\t"
                  "partner\tpident\n")
        for p in sorted(mp):
            st = "shared" if p in pair_main else (
                "unique_with_hit" if p in fwd else "unique_no_hit")
            pa, pi = pair_main.get(p, ("", ""))
            out.write(f"{p}\t{mp[p]['aa']}\t{mp[p]['partial']}\t{mp[p]['edge_dist']}\t"
                      f"{mp[p]['contig_len']}\t{int(bool(mdark.get(p)))}\t{st}\t{pa}\t"
                      f"{pi if pi == '' else '%.1f' % pi}\n")
    log(f"analyse: detail pair {r} vs {g}")
    log(f"analyse: wrote {RES}")


# ------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "prepare", "cluster", "diamond", "analyse"])
    ap.add_argument("--workers", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    os.makedirs(TMP, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    reps, pairs = load_pairs()
    if len(reps) != 23:
        sys.exit(f"ERROR: expected 23 paired representatives, got {len(reps)}")
    log(f"{len(reps)} pairs; "
        f"{sum(len(pairs[r]['cultured']) for r in reps)} cultured genomes")

    if args.stage in ("all", "prepare"):
        if args.force or not os.path.exists(f"{TMP}/mag_proteins.tsv"):
            stage_prepare(reps, pairs, args.workers)
        else:
            log("prepare: outputs present, skipping")
    if args.stage in ("all", "cluster"):
        if args.force or not os.path.exists(f"{TMP}/mag_protein_dark.tsv"):
            stage_cluster(reps, args.workers)
        else:
            log("cluster: mag_protein_dark.tsv present, skipping")
    if args.stage in ("all", "diamond"):
        stage_diamond(reps, pairs, args.workers)
    if args.stage in ("all", "analyse"):
        stage_analyse(reps, pairs)
    print("DONE_STEP3_06")


if __name__ == "__main__":
    main()
