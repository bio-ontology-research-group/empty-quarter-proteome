#!/usr/bin/env python3
"""Step 5 -- the read ladder: where does the community we cannot account for get lost?

CoverM says ~51% of reads map to no representative genome. That number is CUMULATIVE: a read from
an organism that never assembled also fails to map to the catalogue. This step decomposes it.

Each sample has its own reads mapped back to its own assembly at
  <TREE>/<trip>/site<site>/<sample>/wgs/reads/<sample>.sorted.bam  (indexed)
`samtools idxstats` reports mapped reads PER CONTIG plus the unmapped total, and reads only the
index, so it is fast. Crossed with the contig->bin table this splits the reads into:

  never_assembled     unmapped in the sample's own-assembly BAM
  assembled_unbinned  on contigs that are in no bin
  binned              on contigs belonging to a bin
                      (of which in_catalogue = contigs in an HQ/MQ bin that survived refinement)

Out: results/read_ladder_by_sample.tsv , results/read_ladder_summary.tsv

Usage:  sbatch scripts/uv_5_ladder.sbatch      (or run directly for a --limit subset)
"""
import csv, os, sys, gzip, subprocess, statistics, collections

M    = "/ibex/scratch/projects/c2014/rob/metagenomics"
RES  = f"{M}/06_metabolism/uv_survival_clean/results"
TREE = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
SAMTOOLS = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/samtools"
C2B  = f"{M}/04_protein_catalog/contig_to_bin_refined.tsv.gz"

limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else None

meta = {r["sample"]: r for r in csv.DictReader(open(f"{M}/metadata/sample_metadata.tsv"), delimiter="\t")}

# ---- contig -> bin, kept per sample to bound memory --------------------------
print(f"[ladder] reading {C2B}", flush=True)
binned = collections.defaultdict(set)          # sample -> {contig}
with gzip.open(C2B, "rt") as fh:
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) < 3:
            continue
        # contig ids are written `<sample>::<contig>`; BAM reference names carry no prefix
        contig, mag = c[0], c[2]
        s = contig.split("::")[0]
        binned[s].add(contig.split("::", 1)[1] if "::" in contig else contig)
print(f"[ladder] contig->bin covers {len(binned)} samples, "
      f"{sum(len(v) for v in binned.values())} contigs", flush=True)

samples = sorted(meta)
if limit:
    samples = samples[:limit]

rows = []
for i, s in enumerate(samples, 1):
    m = meta[s]
    bam = f"{TREE}/{m['trip']}/site{m['site']}/{s}/wgs/reads/{s}.sorted.bam"
    if not os.path.exists(bam):
        print(f"[ladder] {s}: no BAM, skipped", flush=True)
        continue
    try:
        out = subprocess.run([SAMTOOLS, "idxstats", bam], capture_output=True, text=True,
                             timeout=600, check=True).stdout
    except Exception as e:
        print(f"[ladder] {s}: idxstats failed ({e})", flush=True)
        continue
    bset = binned.get(s, set())
    mapped_binned = mapped_unbinned = unmapped = 0
    # idxstats columns: reference, length, mapped reads, unmapped reads. Column 4 on a contig row
    # holds reads that are themselves unmapped but placed there because their mate mapped; they are
    # unmapped reads, so they join the `*` row in the never-assembled bucket.
    for line in out.splitlines():
        c = line.split("\t")
        if len(c) < 4:
            continue
        name, nmap, nunmap = c[0], int(c[2]), int(c[3])
        unmapped += nunmap
        if name == "*":
            continue
        if name in bset:
            mapped_binned += nmap
        else:
            mapped_unbinned += nmap
    tot = mapped_binned + mapped_unbinned + unmapped
    if tot == 0:
        continue
    rows.append(dict(sample=s, compartment=m["compartment"], total_reads=tot,
                     never_assembled=unmapped / tot * 100,
                     assembled_unbinned=mapped_unbinned / tot * 100,
                     binned=mapped_binned / tot * 100))
    if i % 25 == 0:
        print(f"[ladder] {i}/{len(samples)}", flush=True)

cols = ["sample", "compartment", "total_reads", "never_assembled", "assembled_unbinned", "binned"]
with open(f"{RES}/read_ladder_by_sample.tsv", "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

def med(key, sub=None):
    v = [r[key] for r in rows if sub is None or r["compartment"] == sub]
    return statistics.median(v) if v else float("nan")

with open(f"{RES}/read_ladder_summary.tsv", "w") as o:
    o.write("compartment\tn_samples\tnever_assembled_pct\tassembled_unbinned_pct\tbinned_pct\n")
    print(f"\n{'compartment':12s} {'n':>4s} {'never assembled':>16s} {'assembled, unbinned':>21s} {'binned':>9s}")
    for c in (None, "Surface", "Deep", "PlantRoot"):
        n = len([r for r in rows if c is None or r["compartment"] == c])
        if not n:
            continue
        o.write(f"{c or 'ALL'}\t{n}\t{med('never_assembled',c):.2f}\t"
                f"{med('assembled_unbinned',c):.2f}\t{med('binned',c):.2f}\n")
        print(f"{(c or 'ALL'):12s} {n:4d} {med('never_assembled',c):15.1f}% "
              f"{med('assembled_unbinned',c):20.1f}% {med('binned',c):8.1f}%")
print(f"\n[ladder] wrote {RES}/read_ladder_{{by_sample,summary}}.tsv")
