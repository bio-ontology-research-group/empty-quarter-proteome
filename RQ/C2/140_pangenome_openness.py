#!/usr/bin/env python3
"""
RQ-C2: historical exploratory pangenome workflow for Empty Quarter MAG species.

This script reproduces an archived raw-union/dRep analysis, not the current
checksum-closed catalogue. It tests an association between historical
core/accessory labels and UniRef50 no-hit status; it cannot establish a
mechanism for taxonomy-based gene-prediction error.

Refs: Tettelin et al. 2005 (PNAS); McInerney, McNally & O'Connell 2017 (Nat Microbiol).

Pipeline (all real data, run inside an sbatch job -- NEVER on a login node):
  1. Parse dRep Cdb.csv -> species = secondary_cluster.
     Build a NON-REDUNDANT member set: one genome per (species, sample) so that
     the two binning passes (ss__ single-sample, selfsup__ self-supervised) of
     the *same* assembly cannot inject a near-duplicate genome that would
     artificially inflate the core. Prefer ss__ over selfsup__; lexicographic
     tie-break. Keep species with >= MIN_GENOMES members.
  2. Per genome: gather its proteins.  Proteins were called per-assembly
     (Prodigal, one .faa per sample); a genome's proteins = proteins whose
     contig belongs to that genome's bin.  Join key = contig id.
  3. MMseqs2 cluster the pooled member proteins of each species
     (--min-seq-id 0.5 -c 0.8 --cov-mode 1) -> gene families.
     Classify each family: core (>= CORE_FRAC of genomes), soft-core (>=0.90),
     shell (0.15..0.90), cloud/accessory (<0.15).
  4. Openness: Heaps' law fit of new-genes vs genomes-added (alpha<1 => OPEN),
     plus core/pan ratio, averaged over PERMS random genome orderings.
  5. Functional split: tag every gene family against the dark proteome
     (06_darkmatter/isdark/is_dark.tsv; join on contig_gene key) and against
     UniRef50 best-pident ("knownness"). Pool a core-vs-accessory x dark-vs-known
     2x2 across species as a historical descriptive association.
"""
import os, sys, gzip, glob, json, random, subprocess, shutil, math
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

# ----------------------------------------------------------------------------- config
M           = os.environ.get("M", "/ibex/scratch/projects/c2014/rob/metagenomics")
CDB         = os.environ.get("CDB", f"{M}/04_drep_self_vs_single/out/data_tables/Cdb.csv")
BINDIR      = os.environ.get("BINDIR", f"{M}/04_drep_self_vs_single/inputs")          # *.fa.gz member bins
FAADIR      = os.environ.get("FAADIR", f"{M}/04_protein_catalog/assembly_proteins")    # <SAMPLE>/<SAMPLE>.faa.gz
ISDARK      = os.environ.get("ISDARK", f"{M}/06_darkmatter/isdark/is_dark.tsv")        # key<TAB>pident<TAB>is_dark
OUTDIR      = os.environ.get("OUTDIR", f"{M}/06_pangenome")
WORK        = os.environ.get("WORK", f"{OUTDIR}/work")
MIN_GENOMES = int(os.environ.get("MIN_GENOMES", "10"))
CORE_FRAC   = float(os.environ.get("CORE_FRAC", "0.95"))
SOFT_FRAC   = float(os.environ.get("SOFT_FRAC", "0.90"))
SHELL_FRAC  = float(os.environ.get("SHELL_FRAC", "0.15"))
PERMS       = int(os.environ.get("PERMS", "50"))
THREADS     = int(os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("THREADS", "32")))
MIN_SEQ_ID  = os.environ.get("MIN_SEQ_ID", "0.5")
COV         = os.environ.get("COV", "0.8")
random.seed(13)

OUT  = f"{OUTDIR}/out"
os.makedirs(OUT, exist_ok=True)
os.makedirs(WORK, exist_ok=True)

def log(*a):
    print("[pangenome]", *a, flush=True)

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def sample_of(genome):
    # 'ss__F6Dr2__SemiBin_161.fa.gz' / 'selfsup__F6Dr2__SemiBin_1805.fa.gz' -> 'F6Dr2'
    return genome.split("__")[1]

def strategy_rank(genome):
    return 0 if genome.startswith("ss__") else 1   # prefer single-sample bins

# ----------------------------------------------------------------------------- 1. species/members
def load_species():
    sp_members = defaultdict(list)               # secondary_cluster -> [genome,...] (raw)
    with open(CDB) as fh:
        header = fh.readline().rstrip("\n").split(",")
        gi, si = header.index("genome"), header.index("secondary_cluster")
        for line in fh:
            f = line.rstrip("\n").split(",")
            sp_members[f[si]].append(f[gi])
    # one genome per (species, sample): prefer ss__, lexicographic tie-break
    kept = {}
    for sp, gs in sp_members.items():
        bysamp = {}
        for g in gs:
            s = sample_of(g)
            cur = bysamp.get(s)
            if cur is None or (strategy_rank(g), g) < (strategy_rank(cur), cur):
                bysamp[s] = g
        members = sorted(bysamp.values())
        if len(members) >= MIN_GENOMES:
            kept[sp] = members
    return kept

# ----------------------------------------------------------------------------- 2. contig -> genome
def bin_contigs(genome):
    p = os.path.join(BINDIR, genome)
    if not os.path.exists(p):
        p = os.path.join(BINDIR, genome if genome.endswith(".gz") else genome + ".gz")
    contigs = []
    try:
        with opener(p) as fh:
            for ln in fh:
                if ln[0] == ">":
                    contigs.append(ln[1:].split()[0])
    except FileNotFoundError:
        log("MISSING BIN", genome)
    return contigs

def build_contig_map(kept):
    jobs = sorted({g for ms in kept.values() for g in ms})
    log(f"reading contigs for {len(jobs)} member bins")
    with Pool(THREADS) as pool:
        res = pool.map(bin_contigs, jobs, chunksize=8)
    contig2genome = {}
    genome2sp = {}
    for sp, ms in kept.items():
        for g in ms:
            genome2sp[g] = sp
    for g, contigs in zip(jobs, res):
        for c in contigs:
            contig2genome[c] = g           # contigs are sample-prefixed => globally unique per kept genome
    log(f"contig->genome entries: {len(contig2genome)}")
    return contig2genome, genome2sp

# ----------------------------------------------------------------------------- 3. extract proteins per sample
def extract_sample(args):
    """Stream one sample .faa.gz once; write each protein into its genome's faa.
       Header rewritten to '>{genome}|{origid}'.  Returns (#written, set of dark-lookup keys)."""
    sample, c2g_path, genomes_for_sample = args
    import pickle
    with open(c2g_path, "rb") as fh:
        contig2genome = pickle.load(fh)
    faa = os.path.join(FAADIR, sample, f"{sample}.faa.gz")
    if not os.path.exists(faa):
        log("MISSING FAA", sample); return (0, sample)
    handles = {}
    def h(genome):
        if genome not in handles:
            spdir = os.path.join(WORK, "byspecies", SP_OF[genome])
            os.makedirs(spdir, exist_ok=True)
            handles[genome] = open(os.path.join(spdir, genome + ".faa"), "w")
        return handles[genome]
    written = 0
    cur_g = None; cur_keep = False; cur_hdr = None
    with opener(faa) as fh:
        for ln in fh:
            if ln[0] == ">":
                origid = ln[1:].split()[0]            # F6Dr2::F6Dr2_2241028_1
                after = origid.split("::", 1)[1] if "::" in origid else origid
                contig = after.rsplit("_", 1)[0]      # F6Dr2_2241028
                g = contig2genome.get(contig)
                cur_keep = g is not None
                if cur_keep:
                    cur_g = g
                    h(g).write(f">{g}|{origid}\n")
                    written += 1
            elif cur_keep:
                h(cur_g).write(ln)
    for fwh in handles.values():
        fwh.close()
    return (written, sample)

# ----------------------------------------------------------------------------- helpers for dark/pident
def load_dark(needed_keys):
    """is_dark.tsv: candidate_protein <TAB> best_uniref_pident <TAB> is_dark.
       candidate_protein = '<MAG>::<contig>_<gene>'.  Join key = part after '::'."""
    dark = {}; pid = {}
    if not os.path.exists(ISDARK):
        log("WARN no is_dark.tsv"); return dark, pid
    n = 0
    with opener(ISDARK) as fh:
        first = fh.readline()
        if not first.lower().startswith("candidate"):
            fh.seek(0)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3: continue
            key = f[0].split("::", 1)[1] if "::" in f[0] else f[0]
            if key not in needed_keys:
                continue
            try: pid[key] = float(f[1])
            except ValueError: pass
            dark[key] = (f[2] == "1")
            n += 1
    log(f"dark annotations matched to member proteins: {n}")
    return dark, pid

# ----------------------------------------------------------------------------- 4/5 per-species pangenome
def heaps(presence, genomes, perms):
    """presence: family -> set(genomes). Returns (alpha, gamma, R2, mean pan/core curves)."""
    fams = list(presence.keys())
    g2fams = defaultdict(list)
    for fam, gs in presence.items():
        for g in gs:
            g2fams[g].append(fam)
    N = len(genomes)
    pan = np.zeros(N); newg = np.zeros(N); core = np.zeros(N)
    for _ in range(perms):
        random.shuffle(genomes)
        seen = set()
        counts = defaultdict(int)
        for i, g in enumerate(genomes):
            before = len(seen)
            for fam in g2fams[g]:
                counts[fam] += 1
                seen.add(fam)
            pan[i]  += len(seen)
            newg[i] += len(seen) - before
            ncore = sum(1 for fam in seen if counts[fam] >= CORE_FRAC * (i + 1))
            core[i] += ncore
    pan /= perms; newg /= perms; core /= perms
    # Heaps' law: new genes ~ k * n^-alpha  (Tettelin); open if alpha < 1
    ns = np.arange(2, N + 1)
    y  = newg[1:]
    mask = (y > 0)
    alpha = gamma = R2 = float("nan")
    if mask.sum() >= 3:
        lx = np.log(ns[mask]); ly = np.log(y[mask])
        slope, intercept = np.polyfit(lx, ly, 1)
        alpha = float(-slope)
        pred = slope * lx + intercept
        ss_res = float(np.sum((ly - pred) ** 2)); ss_tot = float(np.sum((ly - ly.mean()) ** 2))
        R2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    # pan-genome size exponent gamma : Pan(n) ~ k n^gamma
    lx = np.log(np.arange(1, N + 1)); ly = np.log(pan)
    g_slope, _ = np.polyfit(lx, ly, 1)
    gamma = float(g_slope)
    return alpha, gamma, R2, pan.tolist(), core.tolist()

def run_species(sp, members):
    spdir = os.path.join(WORK, "byspecies", sp)
    faas = [os.path.join(spdir, g + ".faa") for g in members]
    faas = [f for f in faas if os.path.exists(f) and os.path.getsize(f) > 0]
    if len(faas) < MIN_GENOMES:
        return None
    allfaa = os.path.join(spdir, "all.faa")
    with open(allfaa, "wb") as out:
        for f in faas:
            with open(f, "rb") as fin:
                shutil.copyfileobj(fin, out)
    pref = os.path.join(spdir, "clu")
    tmp  = os.path.join(spdir, "tmp")
    cmd = ["mmseqs", "easy-cluster", allfaa, pref, tmp,
           "--min-seq-id", MIN_SEQ_ID, "-c", COV, "--cov-mode", "1",
           "-e", "1e-3", "--threads", str(THREADS), "-v", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    ctsv = pref + "_cluster.tsv"
    if r.returncode != 0 or not os.path.exists(ctsv):
        log(f"[{sp}] mmseqs FAILED rc={r.returncode}: {r.stderr[-400:]}")
        return None
    # parse clusters: rep <TAB> member ; member header = '{genome}|{origid}'
    presence = defaultdict(set)           # rep -> set(genomes)
    fam_members = defaultdict(list)       # rep -> [origid_key,...]
    with open(ctsv) as fh:
        for ln in fh:
            rep, mem = ln.rstrip("\n").split("\t")
            genome = mem.split("|", 1)[0]
            origid = mem.split("|", 1)[1]
            key = origid.split("::", 1)[1] if "::" in origid else origid   # contig_gene
            presence[rep].add(genome)
            fam_members[rep].append(key)
    # cleanup big intermediates
    for ext in ("_all_seqs.fasta", "_rep_seq.fasta"):
        try: os.remove(pref + ext)
        except OSError: pass
    shutil.rmtree(tmp, ignore_errors=True)
    try: os.remove(allfaa)
    except OSError: pass
    return presence, fam_members, len(faas), [g for g in members if os.path.exists(os.path.join(spdir, g + ".faa"))]

def classify(npres, N):
    f = npres / N
    if f >= CORE_FRAC: return "core"
    if f >= SOFT_FRAC: return "softcore"
    if f >= SHELL_FRAC: return "shell"
    return "cloud"

# ----------------------------------------------------------------------------- main
def main():
    kept = load_species()
    nsp = len(kept); nmem = sum(len(v) for v in kept.values())
    log(f"qualifying species (>= {MIN_GENOMES}, non-redundant): {nsp}  members: {nmem}")
    with open(f"{OUT}/species_members.tsv", "w") as fh:
        for sp, ms in sorted(kept.items()):
            fh.write(f"{sp}\t{len(ms)}\t{','.join(ms)}\n")

    global SP_OF
    contig2genome, SP_OF = build_contig_map(kept)

    # one pickle of contig2genome shared by extraction workers
    import pickle
    c2g_path = os.path.join(WORK, "contig2genome.pkl")
    with open(c2g_path, "wb") as fh:
        pickle.dump(contig2genome, fh, protocol=4)

    samples = sorted({sample_of(g) for ms in kept.values() for g in ms})
    log(f"streaming {len(samples)} sample faa files")
    args = [(s, c2g_path, None) for s in samples]
    with Pool(min(THREADS, len(samples))) as pool:
        for written, s in pool.imap_unordered(extract_sample, args):
            pass
    log("protein extraction done")

    # ---- per-species pangenomes
    sp_rows = []
    needed_keys = set()
    fam_records = []   # (sp, rep, npres, N, cls, keys-list)
    for i, (sp, members) in enumerate(sorted(kept.items()), 1):
        res = run_species(sp, members)
        if res is None:
            continue
        presence, fam_members, N, used = res
        npan = len(presence)
        ncore = sum(1 for r, gs in presence.items() if len(gs) >= CORE_FRAC * N)
        for rep, gs in presence.items():
            cls = classify(len(gs), N)
            keys = fam_members[rep]
            for k in keys:
                needed_keys.add(k)
            fam_records.append((sp, rep, len(gs), N, cls, keys))
        alpha, gamma, R2, pan_curve, core_curve = heaps(presence, list(used), PERMS)
        sp_rows.append(dict(species=sp, n_genomes=N, n_families=npan, n_core=ncore,
                            core_pan_ratio=round(ncore / npan, 4) if npan else 0,
                            heaps_alpha=round(alpha, 4) if alpha == alpha else None,
                            pan_exponent_gamma=round(gamma, 4) if gamma == gamma else None,
                            heaps_R2=round(R2, 4) if R2 == R2 else None,
                            open=("open" if (alpha == alpha and alpha < 1) else "closed")))
        if i % 25 == 0:
            log(f"  ...{i} species done")

    # ---- dark / pident annotation
    dark, pid = load_dark(needed_keys)

    # write species openness
    cols = ["species","n_genomes","n_families","n_core","core_pan_ratio",
            "heaps_alpha","pan_exponent_gamma","heaps_R2","open"]
    with open(f"{OUT}/species_openness.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in sp_rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    # family table + pooled enrichment (core/softcore vs shell/cloud  X  dark vs known)
    a=b=c=d=0   # a: core&dark  b: core&known  c: acc&dark  d: acc&known
    pid_core=[]; pid_acc=[]
    n_acc_total=0; n_acc_dark=0; n_core_total=0; n_core_dark=0
    with open(f"{OUT}/family_table.tsv", "w") as fh:
        fh.write("species\tfamily_rep\tn_genomes_present\tn_total\tclass\tn_annot\tdark_frac\tmean_pident\n")
        for sp, rep, npres, N, cls, keys in fam_records:
            ann = [k for k in keys if k in dark]
            ndark = sum(1 for k in ann if dark[k])
            pids = [pid[k] for k in keys if k in pid]
            dark_frac = (ndark/len(ann)) if ann else float("nan")
            mpid = (sum(pids)/len(pids)) if pids else float("nan")
            is_core = cls in ("core","softcore")
            fam_is_dark = (ndark > 0)
            fam_is_known = (len(ann)-ndark) > 0
            if ann:
                if is_core:
                    n_core_total += 1; n_core_dark += (1 if fam_is_dark else 0)
                    if fam_is_dark: a += 1
                    if fam_is_known: b += 1
                else:
                    n_acc_total += 1; n_acc_dark += (1 if fam_is_dark else 0)
                    if fam_is_dark: c += 1
                    if fam_is_known: d += 1
            if pids:
                (pid_core if is_core else pid_acc).append(mpid)
            fh.write(f"{sp}\t{rep}\t{npres}\t{N}\t{cls}\t{len(ann)}\t"
                     f"{dark_frac:.4f}\t{mpid:.2f}\n".replace("nan","NA"))

    # stats
    def fisher_or(a,b,c,d):
        a_,b_,c_,d_ = a+0.5,b+0.5,c+0.5,d+0.5
        return (a_*d_)/(b_*c_)
    summary = dict(
        min_genomes=MIN_GENOMES, core_frac=CORE_FRAC,
        mmseqs=f"--min-seq-id {MIN_SEQ_ID} -c {COV} --cov-mode 1",
        n_species=len(sp_rows),
        n_genomes_total=sum(r["n_genomes"] for r in sp_rows),
        n_open=sum(1 for r in sp_rows if r["open"]=="open"),
        n_closed=sum(1 for r in sp_rows if r["open"]=="closed"),
        median_heaps_alpha=float(np.nanmedian([r["heaps_alpha"] for r in sp_rows if r["heaps_alpha"] is not None])) if sp_rows else None,
        median_core_pan_ratio=float(np.median([r["core_pan_ratio"] for r in sp_rows])) if sp_rows else None,
        enrichment_2x2=dict(core_dark=a, core_known=b, accessory_dark=c, accessory_known=d),
        dark_odds_ratio_accessory_vs_core=round(fisher_or(c,d,a,b),3),
        accessory_dark_overlap=dict(n_accessory_annot=n_acc_total, n_accessory_dark=n_acc_dark,
                                    frac=round(n_acc_dark/n_acc_total,4) if n_acc_total else None),
        core_dark_overlap=dict(n_core_annot=n_core_total, n_core_dark=n_core_dark,
                               frac=round(n_core_dark/n_core_total,4) if n_core_total else None),
        mean_pident_core=round(float(np.mean(pid_core)),2) if pid_core else None,
        mean_pident_accessory=round(float(np.mean(pid_acc)),2) if pid_acc else None,
    )
    # optional scipy tests
    try:
        from scipy.stats import fisher_exact, mannwhitneyu
        _or, _p = fisher_exact([[c,d],[a,b]])
        summary["fisher_p_accessory_dark"] = float(_p)
        if pid_core and pid_acc:
            _u, _pu = mannwhitneyu(pid_acc, pid_core, alternative="less")
            summary["mannwhitney_p_pident_accessory_lt_core"] = float(_pu)
    except Exception as e:
        summary["scipy"] = f"unavailable: {e}"

    with open(f"{OUT}/enrichment_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    log("SUMMARY:", json.dumps(summary, indent=2))
    log("DONE. Outputs in", OUT)

if __name__ == "__main__":
    main()
