#!/usr/bin/env python3
# Parameterized CDS slicer (shared by dark arm [fam_3] and control arm [fam_6 prep]).
# Usage: slice_cds.py <member_fam.tsv> <out_dir> <sample>
# member_fam.tsv: member_protein_id \t family_rep. Coordinates from the per-sample
# assembly protein headers; sequence from the per-sample assembly contigs.
import sys, gzip, os
mmap_path, out_dir, s = sys.argv[1], sys.argv[2], sys.argv[3]
faa = f"/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog/assembly_proteins/{s}/{s}.faa.gz"
asm = f"/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/assemblies/{s}/{s}_contigs.fa"
out = f"{out_dir}/{s}.fna"
if os.path.exists(out) and os.path.getsize(out) > 0:
    sys.exit(0)
need = {}
for line in open(mmap_path):
    m, f = line.rstrip("\n").split("\t")
    if m.startswith(s + "::"):
        need[m] = f
if not need:
    open(out, "w").close(); sys.exit(0)
coords = {}
with gzip.open(faa, "rt") as fh:
    for line in fh:
        if not line.startswith(">"): continue
        parts = line[1:].split(" # ")
        pid = parts[0].split()[0]
        if pid in need:
            coords[pid] = (int(parts[1]), int(parts[2]), int(parts[3]))
by_contig = {}
for pid, c in coords.items():
    contig = pid.split("::")[1].rsplit("_", 1)[0]
    by_contig.setdefault(contig, []).append((pid, c))
comp = str.maketrans("ACGTacgt", "TGCAtgca")
wrote = 0
with open(out + ".part", "w") as o:
    name, seq, keep = None, [], None
    def flush():
        global wrote
        if keep and keep in by_contig:
            cs = "".join(seq)
            for pid, (st, en, strand) in by_contig[keep]:
                sub = cs[st-1:en]
                if strand == -1: sub = sub.translate(comp)[::-1]
                o.write(f">{pid} fam={need[pid]}\n{sub}\n"); wrote += 1
    for line in open(asm):
        if line.startswith(">"):
            flush(); name = line[1:].split()[0]; keep = name if name in by_contig else None; seq = []
        elif keep: seq.append(line.strip())
    flush()
os.replace(out + ".part", out)
print(f"[{s}] members {len(need)} coords {len(coords)} written {wrote} missing {len(need)-len(coords)}", flush=True)
