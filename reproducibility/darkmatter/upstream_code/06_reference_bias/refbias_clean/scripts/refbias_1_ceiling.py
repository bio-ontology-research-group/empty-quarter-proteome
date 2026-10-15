#!/usr/bin/env python3
"""refbias step 1 -- the no-KO CEILING on the CURRENT 5,438-representative catalogue.

The headline number of paper section 2.1: what fraction of MAG-encoded proteins carry no KEGG
Ortholog at all, and are therefore invisible to ANY KO-based 16S predictor (PICRUSt2 included)
regardless of how accurate its placement is.

The historical 56.7% was computed on the OLD 15,469 single-sample SemiBin cohort
(repo results/picrust2_bias/, June 2026) as 1 - sum(n_ko)/sum(CheckM2 Total_Coding_Sequences).
This recomputes it on the current catalogue, per protein, with two denominators:

  BAKTA   (primary)   proteins actually in the Bakta .faa -- the same gene calls eggNOG scored,
                      so numerator and denominator come from one gene-calling universe.
  CHECKM2 (secondary) CheckM2's own Total_Coding_Sequences (pyrodigal gene calls) -- reported
                      only for comparability with the historical figure, which used it.
                      The two callers disagree by ~6% per bin, so they are NOT interchangeable.

Counting, per representative:
  A = protein-coding sequences in the bin
  B = those with a KEGG_ko assignment in eggNOG
Ceiling = 1 - sum(B)/sum(A) over all 5,438 representatives (pooled, protein-weighted).

CRITICAL: `.emapper.annotations` only contains rows for queries that got a seed-ortholog hit.
Proteins with no hit are absent from the file entirely, NOT written as `-`. They are real
denominator members. Counting rows in the annotations file as the denominator would silently
delete them and understate the ceiling -- so the denominator always comes from the .faa.

_sub representatives (1,010 of 5,438) are DAS_Tool subsets of a parent bin and have no annotation
of their own. Their proteins are recovered by the SAME contig-subset construction as
repko_2_submap.py (contig_N -> original recovered from the raw bin, verified position-wise by
sequence length; refused, never guessed) -- so the ceiling inherits the rep_ko_v2 correction
rather than the parent's superset.

  python3 scripts/refbias_1_ceiling.py [n_workers]
Writes:
  results/ceiling_per_rep.tsv   per-representative counts + status
  results/ceiling_summary.tsv   pooled + per-quality ceilings, both denominators
"""
import csv, os, gzip, sys
from multiprocessing import Pool

M = "/ibex/scratch/projects/c2014/rob/metagenomics"
T = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
D = f"{M}/05_refine/_drep/out/dereplicated_genomes"
W = f"{M}/06_reference_bias/refbias_clean"
CHECKM2 = f"{M}/05_refine/_checkm2_refined/quality_report.tsv"

meta = {r["sample"]: r for r in csv.DictReader(open(f"{M}/metadata/sample_metadata.tsv"), delimiter="\t")}


def seqs(path, gz=False):
    """[(name, length)] in file order."""
    op = gzip.open if gz else open
    out, name, n = [], None, 0
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, n))
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        out.append((name, n))
    return out


def faa_locus(path):
    """Locus tags of the proteins in a Bakta .faa, in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line[1:].split()[0])
    return out


def ko_queries(emapper, keep_locus=None):
    """(queries_with_a_KO, all_queries_seen) as sets of locus tags."""
    with_ko, seen, hdr = set(), set(), None
    with open(emapper) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr = line.lstrip("#").rstrip("\n").split("\t"); continue
            if line.startswith("#") or not hdr:
                continue
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            q = row.get("query")
            if not q:
                continue
            if keep_locus is not None and q not in keep_locus:
                continue
            seen.add(q)
            ko = row.get("KEGG_ko", "-")
            if ko and ko not in ("-", ""):
                with_ko.add(q)
    return with_ko, seen


def resolve(rep):
    """(sample, mode candidates, local bin name) for a representative."""
    s = rep.split("__")[0]
    m = meta.get(s)
    if "__" not in rep or not m:
        return None
    rest = rep.split("__", 1)[1]
    if rest.startswith("single_"):   modes, local = ["single_sample"], rest[len("single_"):]
    elif rest.startswith("self_"):   modes, local = ["self_supervised"], rest[len("self_"):]
    else:                            modes, local = ["self_supervised", "single_sample"], rest
    return s, m, modes, local


def count_rep(rep):
    """-> dict of per-rep counts, or a status row explaining why not."""
    r = resolve(rep)
    if r is None:
        return {"rep": rep, "status": "NO_METADATA"}
    s, m, modes, local = r
    base = f"{T}/{m['trip']}/site{m['site']}/{s}/wgs"
    is_sub = local.endswith("_sub")

    if not is_sub:
        for mode in modes:
            faa = f"{base}/annotation/{mode}/bakta/{local}/{local}.faa"
            emap = f"{base}/annotation/{mode}/eggnog/{local}/{local}.emapper.annotations"
            if os.path.exists(faa) and os.path.exists(emap):
                prot = faa_locus(faa)
                with_ko, seen = ko_queries(emap)
                # annotations may name proteins the .faa does not (defensive: intersect)
                pset = set(prot)
                return {"rep": rep, "mode": mode, "is_sub": 0,
                        "n_prot": len(prot), "n_ko": len(with_ko & pset),
                        "n_hit": len(seen & pset), "status": "OK"}
        return {"rep": rep, "status": "ANNOTATION_NOT_FOUND"}

    # _sub: rebuild the retained-contig protein set from the parent (repko_2 logic, verbatim)
    parent = local[:-4]
    mode = modes[0] if len(modes) == 1 else "self_supervised"
    raw  = f"{base}/binning/{mode}/output_bins/{parent}.fa.gz"
    btsv = f"{base}/annotation/{mode}/bakta/{parent}/{parent}.tsv"
    bfna = f"{base}/annotation/{mode}/bakta/{parent}/{parent}.fna"
    bfaa = f"{base}/annotation/{mode}/bakta/{parent}/{parent}.faa"
    emap = f"{base}/annotation/{mode}/eggnog/{parent}/{parent}.emapper.annotations"
    if not all(os.path.exists(p) for p in (raw, btsv, bfna, bfaa, emap)):
        return {"rep": rep, "status": "PARENT_NOT_ANNOTATED"}

    R, K = seqs(raw, True), seqs(bfna)
    if len(R) != len(K) or any(R[i][1] != K[i][1] for i in range(len(R))):
        return {"rep": rep, "status": "ORDER_MISMATCH_REFUSED"}
    name_of = {K[i][0]: R[i][0] for i in range(len(R))}
    sub_contigs = {n for n, _ in seqs(f"{D}/{rep}.fa")}
    keep_ids = {c for c, orig in name_of.items() if orig in sub_contigs}

    keep_locus = set()
    with open(btsv) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            if c[0] in keep_ids:
                keep_locus.add(c[5])
    # keep_locus spans ALL feature types (tRNA/rRNA/ncRNA too) -- restrict to real proteins
    prot = set(faa_locus(bfaa)) & keep_locus
    with_ko, seen = ko_queries(emap, keep_locus)
    return {"rep": rep, "mode": mode, "is_sub": 1,
            "n_prot": len(prot), "n_ko": len(with_ko & prot),
            "n_hit": len(seen & prot), "status": "OK"}


def main():
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    reps = sorted(f[:-3] for f in os.listdir(D) if f.endswith(".fa"))
    print(f"[ceiling] representatives: {len(reps)} | workers: {nw}", flush=True)

    cm2 = {}
    with open(CHECKM2) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            cm2[row["Name"]] = row
    print(f"[ceiling] checkm2 rows: {len(cm2)}", flush=True)

    with Pool(nw) as pool:
        rows = pool.map(count_rep, reps, chunksize=8)

    os.makedirs(f"{W}/results", exist_ok=True)
    cols = ["rep", "mode", "is_sub", "quality", "n_prot", "n_ko", "n_hit", "checkm2_cds", "status"]
    ok = 0
    with open(f"{W}/results/ceiling_per_rep.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for d in rows:
            c = cm2.get(d["rep"], {})
            try:
                comp = float(c.get("Completeness", "nan")); cont = float(c.get("Contamination", "nan"))
                qual = "HQ" if (comp >= 90 and cont <= 5) else "MQ"
            except ValueError:
                qual = "NA"
            d["quality"] = qual
            d["checkm2_cds"] = c.get("Total_Coding_Sequences", "")
            if d["status"] == "OK":
                ok += 1
            o.write("\t".join(str(d.get(k, "")) for k in cols) + "\n")

    good = [d for d in rows if d["status"] == "OK"]
    def ceiling(sub, denom):
        A = sum(d["n_prot"] if denom == "bakta" else int(d["checkm2_cds"] or 0) for d in sub)
        B = sum(d["n_ko"] for d in sub)
        return A, B, (1 - B / A) if A else float("nan")

    with open(f"{W}/results/ceiling_summary.tsv", "w") as o:
        o.write("cohort\tdenominator\tn_reps\tn_proteins\tn_with_ko\tfrac_no_ko\n")
        for name, sub in [("all", good),
                          ("HQ", [d for d in good if d["quality"] == "HQ"]),
                          ("MQ", [d for d in good if d["quality"] == "MQ"]),
                          ("sub", [d for d in good if d["is_sub"] == 1]),
                          ("nonsub", [d for d in good if d["is_sub"] == 0])]:
            for denom in ("bakta", "checkm2"):
                A, B, f = ceiling(sub, denom)
                o.write(f"{name}\t{denom}\t{len(sub)}\t{A}\t{B}\t{f:.6f}\n")
                if name == "all":
                    print(f"[ceiling] {denom:8s} {A:,} proteins, {B:,} with KO -> no-KO {f*100:.2f}%", flush=True)

    bad = {}
    for d in rows:
        if d["status"] != "OK":
            bad[d["status"]] = bad.get(d["status"], 0) + 1
    print(f"[ceiling] OK {ok}/{len(reps)} | non-OK: {bad}", flush=True)
    print(f"[ceiling] -> results/ceiling_per_rep.tsv, results/ceiling_summary.tsv", flush=True)


if __name__ == "__main__":
    main()
