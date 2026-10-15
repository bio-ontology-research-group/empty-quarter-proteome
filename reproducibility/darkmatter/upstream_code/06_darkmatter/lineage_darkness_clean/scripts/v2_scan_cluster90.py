#!/usr/bin/env python
"""ROUND TWO, GOAL A: native-key dark fraction for all 5,438 refined species reps.

Member-level logic copied from
  W/04_protein_catalog/darkmatter_clean/scripts/darkclean_5_mag.sh
but keyed on the REFINED catalogue instead of the unrefined HQ+MQ bins, so no
rename is needed and every one of the 5,438 dRep winners is computed natively.

  contig -> winner bin   W/04_protein_catalog/contig_to_bin_refined.tsv.gz  (restricted to Wdb.csv)
  protein -> contig      protein id is  <sample>::<contig>_<gene>, contig = id up to the last "_"
  cluster membership     W/04_protein_catalog/all_linclust90/all_cluster90.tsv   rep TAB member
  darkness               cluster property, taken from the rep:
                         W/04_protein_catalog/darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz

PHASE 1 (parallel, byte ranges over the 36 GB cluster table, binary mode)
    every member whose contig belongs to a winner bin  ->  tmp/part_<k>.tsv   "cluster_rep TAB winner_index"
PHASE 2 (single process)
    unique cluster reps touched  x  dark_rep_ids  ->  dark flag  ->  per-winner totals

Usage:  v2_scan_cluster90.py [n_workers]      (default 8; re-run reuses tmp/part_*.tsv)
"""
import csv
import gzip
import os
import random
import re
import subprocess
import sys
import time
from multiprocessing import Process

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(BASE, "results")
TMP = os.path.join(BASE, "tmp")

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
F_WDB = W + "/05_refine/_drep/out/data_tables/Wdb.csv"
F_C2B = W + "/04_protein_catalog/contig_to_bin_refined.tsv.gz"
F_CLUSTER = W + "/04_protein_catalog/all_linclust90/all_cluster90.tsv"
F_DARK = W + "/04_protein_catalog/darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz"

NW = int(sys.argv[1]) if len(sys.argv) > 1 else 8


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def strip_ext(s):
    return re.sub(r"\.(fa|fasta)$", "", str(s))


# ---------------------------------------------------------------- winners
def load_winners():
    win = []
    with open(F_WDB) as fh:
        for row in csv.DictReader(fh):
            win.append(strip_ext(row["genome"]))
    return win


def build_contig_map(winners):
    """contig (bytes, 'sample::contig') -> winner index."""
    wi = {g: i for i, g in enumerate(winners)}
    c2w = {}
    with gzip.open(F_C2B, "rt") as fh:
        for line in fh:
            contig, _strat, b = line.rstrip("\n").split("\t")
            i = wi.get(b)
            if i is not None:
                c2w[contig.encode()] = i
    return c2w


# ---------------------------------------------------------------- phase 1
def worker(k, start, end, c2w, probe_idx):
    out = open(os.path.join(TMP, "part_%02d.tsv" % k), "wb")
    prb = open(os.path.join(TMP, "probe_%02d.txt" % k), "wb")
    n = 0
    hit = 0
    with open(F_CLUSTER, "rb") as fh:
        fh.seek(start)
        if start > 0:
            fh.readline()          # the previous worker owns this partial line
        pos = fh.tell()
        get = c2w.get
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
            u = member.rfind(b"_")
            i = get(member[:u])
            if i is not None:
                hit += 1
                out.write(line[:t] + b"\t" + str(i).encode() + b"\n")
                if i in probe_idx:
                    prb.write(str(i).encode() + b"\t" + member + b"\n")
    out.close()
    prb.close()
    with open(os.path.join(TMP, "stat_%02d.txt" % k), "w") as fh:
        fh.write("%d\t%d\n" % (n, hit))


def phase1(c2w, probe_idx):
    size = os.path.getsize(F_CLUSTER)
    step = size // NW
    bounds = [(k * step, size if k == NW - 1 else (k + 1) * step) for k in range(NW)]
    procs = []
    for k, (s, e) in enumerate(bounds):
        p = Process(target=worker, args=(k, s, e, c2w, probe_idx))
        p.start()
        procs.append(p)
    for p in procs:
        p.join()
    n = hit = 0
    for k in range(NW):
        a, b = open(os.path.join(TMP, "stat_%02d.txt" % k)).read().split()
        n += int(a)
        hit += int(b)
    return n, hit


# ---------------------------------------------------------------- phase 2
def phase2(winners):
    reps = []
    idxs = []
    for k in range(NW):
        with open(os.path.join(TMP, "part_%02d.tsv" % k), "rb") as fh:
            for line in fh:
                t = line.rfind(b"\t")
                reps.append(line[:t])
                idxs.append(int(line[t + 1:-1]))
    log("phase 2: %d member hits loaded" % len(reps))

    touched = set(reps)
    log("phase 2: %d distinct cluster reps touched" % len(touched))

    dark = set()
    pr = subprocess.Popen(["zcat", F_DARK], stdout=subprocess.PIPE)
    n_dark_lines = 0
    for line in pr.stdout:
        n_dark_lines += 1
        v = line[:-1]
        if v in touched:
            dark.add(v)
    pr.wait()
    log("phase 2: dark_rep_ids %d lines, %d of them touch a winner" % (n_dark_lines, len(dark)))
    del touched

    tot = [0] * len(winners)
    drk = [0] * len(winners)
    for r, i in zip(reps, idxs):
        tot[i] += 1
        if r in dark:
            drk[i] += 1
    return tot, drk, len(dark), n_dark_lines


# ---------------------------------------------------------------- main
def main():
    os.makedirs(TMP, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    rep_lines = []

    def say(s=""):
        print(s, flush=True)
        rep_lines.append(s)

    winners = load_winners()
    log("winners %d" % len(winners))

    # 20 probe reps, fixed seed, stratified over the three naming forms
    forms = {}
    for i, g in enumerate(winners):
        if "__self_SemiBin_" not in g and "__single_SemiBin_" not in g:
            forms.setdefault("no_assembly_token", []).append(i)
        elif g.endswith("_sub"):
            forms.setdefault("_sub", []).append(i)
        else:
            forms.setdefault("plain", []).append(i)
    rng = random.Random(20260910)
    probe = []
    for f, take in (("plain", 8), ("_sub", 6), ("no_assembly_token", 6)):
        probe += rng.sample(forms[f], take)
    probe_idx = set(probe)

    c2w = build_contig_map(winners)
    log("winner contigs %d" % len(c2w))

    done = all(os.path.exists(os.path.join(TMP, "part_%02d.tsv" % k)) for k in range(NW))
    if done:
        log("phase 1 outputs already present, reusing")
        n_lines = hit = -1
        try:
            n_lines = hit = 0
            for k in range(NW):
                a, b = open(os.path.join(TMP, "stat_%02d.txt" % k)).read().split()
                n_lines += int(a)
                hit += int(b)
        except OSError:
            pass
    else:
        t0 = time.time()
        n_lines, hit = phase1(c2w, probe_idx)
        log("phase 1 done in %.1f min: %d cluster-table lines, %d member hits"
            % ((time.time() - t0) / 60, n_lines, hit))
    del c2w

    tot, drk, n_dark_touch, n_dark_lines = phase2(winners)

    with open(os.path.join(RES, "v2_rep_dark.tsv"), "w") as fh:
        fh.write("rep\tn_proteins\tn_dark\tdark_fraction\n")
        for g, t, d in zip(winners, tot, drk):
            fh.write("%s\t%d\t%d\t%.6f\n" % (g, t, d, (d / t) if t else 0.0))

    say("V2 NATIVE-KEY SCAN (GOAL A)")
    say("===========================")
    say("cluster table lines scanned : %d" % n_lines)
    say("members on a winner contig  : %d" % hit)
    say("dark_rep_ids lines          : %d" % n_dark_lines)
    say("dark clusters touching a rep: %d" % n_dark_touch)
    say("winners with 0 proteins     : %d" % sum(1 for t in tot if t == 0))
    say("winners computed            : %d / %d" % (sum(1 for t in tot if t > 0), len(winners)))
    say("total proteins  %d   total dark %d   global dark fraction %.6f"
        % (sum(tot), sum(drk), sum(drk) / max(sum(tot), 1)))
    say()

    # --- GFF cross-check on the 20 probe reps
    say("GFF CROSS-CHECK, 20 probe reps (seed 20260910; 8 plain / 6 _sub / 6 no_assembly_token)")
    say("%-40s %10s %10s %8s %10s" % ("rep", "contig_n", "gff_CDS", "delta", "id_match"))
    probe_members = {}
    for k in range(NW):
        p = os.path.join(TMP, "probe_%02d.txt" % k)
        if not os.path.exists(p):
            continue
        with open(p, "rb") as fh:
            for line in fh:
                t = line.find(b"\t")
                probe_members.setdefault(int(line[:t]), set()).add(line[t + 1:-1].decode())
    n_agree = 0
    n_checked = 0
    for i in probe:
        g = winners[i]
        gff = W + "/06_darkmatter/drep_rep_proteomes/%s.gff.gz" % g
        if not os.path.exists(gff):
            say("%-40s %10d %10s %8s %10s" % (g, tot[i], "MISSING", "-", "-"))
            continue
        ids = []
        with gzip.open(gff, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) > 8 and f[2] == "CDS":
                    m = re.search(r"ID=([^;]+)", f[8])
                    if m:
                        ids.append(m.group(1))
        # gff ID is "<contig>_<gene>" and contig names are "<sample>_<number>",
        # so the catalogue-style id is "<sample>::<contig>_<gene>"
        cat = set()
        for pid in ids:
            contig = pid.rsplit("_", 1)[0]
            cat.add("%s::%s" % (contig.rsplit("_", 1)[0], pid))
        got = probe_members.get(i, set())
        inter = len(cat & got)
        n_checked += 1
        if len(ids) == tot[i]:
            n_agree += 1
        say("%-40s %10d %10d %8d %9.1f%%"
            % (g, tot[i], len(ids), len(ids) - tot[i], 100.0 * inter / max(len(cat), 1)))
    say("exact protein-count agreement: %d / %d probe reps" % (n_agree, n_checked))
    say("id_match = share of the GFF-derived catalogue-style ids that are actually")
    say("catalogue members on that bin's contigs.")

    with open(os.path.join(RES, "v2_scan_report.txt"), "w") as fh:
        fh.write("\n".join(rep_lines) + "\n")
    log("wrote results/v2_rep_dark.tsv and results/v2_scan_report.txt")


if __name__ == "__main__":
    main()
