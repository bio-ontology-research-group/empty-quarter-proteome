#!/usr/bin/env python3
"""repko step 3 -- REBUILD the per-representative KO table, corrected.

Rebuilds `rep_ko.tsv` for all 5,438 representatives, fixing the three defects found in the
original (`scripts/eco_probes/gather_rep_ko.py`, hohndor, 2026-06-11):

  1. NAMING     44 reps of sample 48PRr1 are named `<sample>__SemiBin_N` with no single_/self_
                token (that sample has only self-supervised binning). The original resolver splits
                on the token and failed, leaving them unannotated although the annotation exists.
                Fixed by trying BOTH modes when the token is absent.
  2. _sub BINS  1,010 reps are DAS_Tool `_sub` bins (strict subsets of a parent). The original
                fell back to the PARENT's annotation -- a superset -- attributing genes that
                refinement removed. Fixed by using the contig-subset KO sets from step 2.
  3. MISSING    124 raw bins were never annotated (80 LQ pre-refinement, 44 parents of _sub bins).
                Fixed by step 1, which writes their annotation into the per-sample tree, so the
                normal resolver now finds them.

The PI's original `04_functional/refs/rep_ko.tsv` is NEVER modified. Output goes to
`repko_clean/results/rep_ko_v2.tsv`; downstream scripts are repointed deliberately.

  python3 scripts/repko_3_build.py
"""
import csv, os, glob

M = "/ibex/scratch/projects/c2014/rob/metagenomics"
T = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
D = f"{M}/05_refine/_drep/out/dereplicated_genomes"
W = f"{M}/04_functional/repko_clean"
meta = {r["sample"]: r for r in csv.DictReader(open(f"{M}/metadata/sample_metadata.tsv"), delimiter="\t")}

# corrected KO sets for _sub reps (step 2)
sub_ko = {}
p = f"{W}/results/sub_ko.tsv"
if os.path.exists(p):
    for line in open(p):
        c = line.rstrip("\n").split("\t")
        if c[0]:
            sub_ko[c[0]] = c[1] if len(c) > 1 else ""
print(f"[build] corrected _sub KO sets available: {len(sub_ko)}")


def ann_path(rep):
    """FIXED resolver: handles the missing single_/self_ token by trying both modes."""
    s = rep.split("__")[0]
    m = meta.get(s)
    if "__" not in rep or not m:
        return None
    rest = rep.split("__", 1)[1]
    if rest.startswith("single_"):   modes, local = ["single_sample"], rest[len("single_"):]
    elif rest.startswith("self_"):   modes, local = ["self_supervised"], rest[len("self_"):]
    else:                            modes, local = ["self_supervised", "single_sample"], rest  # FIX 1
    base = f"{T}/{m['trip']}/site{m['site']}/{s}/wgs/annotation"
    for mode in modes:
        for lc in [local] + ([local[:-4]] if local.endswith("_sub") else []):
            q = f"{base}/{mode}/eggnog/{lc}/{lc}.emapper.annotations"
            if os.path.exists(q):
                return q
    return None


def parse_kos(path):
    kos, hdr = set(), None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr = line.lstrip("#").rstrip("\n").split("\t"); continue
            if line.startswith("#") or not hdr:
                continue
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            ko = row.get("KEGG_ko", "-")
            if ko and ko not in ("-", ""):
                kos |= {k.replace("ko:", "") for k in ko.split(",")}
    return kos


reps = sorted(os.path.basename(x)[:-3] for x in glob.glob(f"{D}/*.fa"))
src = {"sub_corrected": 0, "direct": 0, "missing": 0}
missing = []
with open(f"{W}/results/rep_ko_v2.tsv", "w") as o:
    for rep in reps:
        if rep in sub_ko:                       # FIX 2: contig-subset KOs
            o.write(f"{rep}\t{sub_ko[rep]}\n"); src["sub_corrected"] += 1; continue
        q = ann_path(rep)                       # FIX 1 + FIX 3
        if q:
            o.write(f"{rep}\t{','.join(sorted(parse_kos(q)))}\n"); src["direct"] += 1
        else:
            o.write(f"{rep}\t\n"); src["missing"] += 1; missing.append(rep)

print(f"[build] reps={len(reps)}")
print(f"  from contig-subset (_sub, corrected) : {src['sub_corrected']}")
print(f"  from direct annotation               : {src['direct']}")
print(f"  STILL unannotated                    : {src['missing']}")
if missing:
    open(f"{W}/results/still_unannotated.txt", "w").write("\n".join(missing) + "\n")
    print(f"    -> results/still_unannotated.txt  e.g. {missing[:3]}")

# ---- comparison against the original table (what changed, and why) ----
old = {}
for line in open(f"{M}/04_functional/refs/rep_ko.tsv"):
    c = line.rstrip("\n").split("\t")
    old[c[0]] = set(filter(None, (c[1] if len(c) > 1 else "").split(",")))
new = {}
for line in open(f"{W}/results/rep_ko_v2.tsv"):
    c = line.rstrip("\n").split("\t")
    new[c[0]] = set(filter(None, (c[1] if len(c) > 1 else "").split(",")))
gained = sum(1 for r in new if new[r] - old.get(r, set()))
lost   = sum(1 for r in new if old.get(r, set()) - new[r])
same   = sum(1 for r in new if new[r] == old.get(r, set()))
n_old_unann = sum(1 for r in old if not old[r])
n_new_unann = sum(1 for r in new if not new[r])
print(f"[build] vs original: identical {same} | gained KOs {gained} | lost (spurious removed) {lost}")
print(f"[build] unannotated reps: {n_old_unann} -> {n_new_unann}")
print(f"[build] wrote {W}/results/rep_ko_v2.tsv  (original untouched)")
