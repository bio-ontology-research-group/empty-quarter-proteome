#!/usr/bin/env python3
"""repko step 2 -- MAP PARENT ANNOTATIONS ONTO _sub REPRESENTATIVES.

DAS_Tool emits a `_sub` bin when a candidate bin partially overlaps an already-selected bin: it
strips the contigs claimed by the winner and keeps the remainder. A `_sub` bin is therefore a
strict SUBSET of its parent raw bin (verified: 1,010/1,010).

`gather_rep_ko.py` has no annotation for `<bin>_sub`, so it falls back to the PARENT's annotation
-- which is a superset. Those representatives currently inherit genes that refinement removed
(measured: up to 28% of a rep's KOs are spurious).

This script computes the CORRECT KO set for every `_sub` representative by keeping only the genes
that lie on contigs the sub bin actually retained. It does not re-run any tool: same Bakta gene
calls, same eggNOG version, no re-training drift.

Bakta renames contigs to `contig_1, contig_2, ...` but preserves INPUT ORDER, so the mapping
`contig_N -> original contig name` is recovered from the raw bin and VERIFIED position-by-position
by sequence length. Any bin whose lengths do not match in order is refused, not guessed.

  python3 scripts/repko_2_submap.py
Writes:
  results/sub_ko.tsv          rep <TAB> comma-separated KOs   (the corrected sets)
  results/sub_ko_report.tsv   per-rep diagnostics (contigs/features kept, KOs dropped)
"""
import csv, os, gzip, sys

M = "/ibex/scratch/projects/c2014/rob/metagenomics"
T = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
D = f"{M}/05_refine/_drep/out/dereplicated_genomes"
W = f"{M}/04_functional/repko_clean"
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


def kos_for(emapper, keep_locus=None):
    """KO set; if keep_locus is given, only genes whose locus tag is in it."""
    kos, hdr = set(), None
    with open(emapper) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr = line.lstrip("#").rstrip("\n").split("\t"); continue
            if line.startswith("#") or not hdr:
                continue
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if keep_locus is not None and row["query"] not in keep_locus:
                continue
            ko = row.get("KEGG_ko", "-")
            if ko and ko not in ("-", ""):
                kos |= {k.replace("ko:", "") for k in ko.split(",")}
    return kos


subs = sorted(r[:-3] for r in os.listdir(D) if r.endswith(".fa") and r[:-3].endswith("_sub"))
print(f"[submap] _sub representatives: {len(subs)}")

out = open(f"{W}/results/sub_ko.tsv", "w")
rep = open(f"{W}/results/sub_ko_report.tsv", "w")
rep.write("rep\tparent\tcontigs_parent\tcontigs_sub\tfeat_parent\tfeat_sub\tko_parent\tko_sub\tko_dropped\tstatus\n")
done = refused = noparent = 0
for r in subs:
    s = r.split("__")[0]; rest = r.split("__", 1)[1]
    if rest.startswith("single_"):  mode, local = "single_sample", rest[len("single_"):]
    elif rest.startswith("self_"):  mode, local = "self_supervised", rest[len("self_"):]
    else:                           mode, local = "self_supervised", rest
    parent = local[:-4]
    m = meta.get(s)
    if not m:
        refused += 1; continue
    base = f"{T}/{m['trip']}/site{m['site']}/{s}/wgs"
    raw  = f"{base}/binning/{mode}/output_bins/{parent}.fa.gz"
    btsv = f"{base}/annotation/{mode}/bakta/{parent}/{parent}.tsv"
    bfna = f"{base}/annotation/{mode}/bakta/{parent}/{parent}.fna"
    emap = f"{base}/annotation/{mode}/eggnog/{parent}/{parent}.emapper.annotations"
    if not (os.path.exists(btsv) and os.path.exists(bfna) and os.path.exists(emap)):
        noparent += 1
        rep.write(f"{r}\t{parent}\t\t\t\t\t\t\t\tPARENT_NOT_ANNOTATED\n")
        continue

    R, K = seqs(raw, True), seqs(bfna)
    # order must be preserved for contig_N -> original to be valid; verify by length
    if len(R) != len(K) or any(R[i][1] != K[i][1] for i in range(len(R))):
        refused += 1
        rep.write(f"{r}\t{parent}\t{len(R)}\t\t\t\t\t\t\tORDER_MISMATCH_REFUSED\n")
        continue
    name_of = {K[i][0]: R[i][0] for i in range(len(R))}          # contig_N -> original
    sub_contigs = {n for n, _ in seqs(f"{D}/{r}.fa")}
    keep_ids = {c for c, orig in name_of.items() if orig in sub_contigs}

    keep_locus, feat_parent = set(), 0
    with open(btsv) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            feat_parent += 1
            if c[0] in keep_ids:
                keep_locus.add(c[5])
    ko_par = kos_for(emap)
    ko_sub = kos_for(emap, keep_locus)
    out.write(f"{r}\t{','.join(sorted(ko_sub))}\n")
    rep.write(f"{r}\t{parent}\t{len(R)}\t{len(sub_contigs)}\t{feat_parent}\t{len(keep_locus)}"
              f"\t{len(ko_par)}\t{len(ko_sub)}\t{len(ko_par - ko_sub)}\tOK\n")
    done += 1
out.close(); rep.close()
print(f"[submap] corrected: {done} | parent not annotated (run step 1 first): {noparent} | refused: {refused}")
print(f"[submap] -> results/sub_ko.tsv, results/sub_ko_report.tsv")
if noparent:
    print(f"[submap] NOTE: rerun this after step 1 to pick up the {noparent} newly annotated parents.")
