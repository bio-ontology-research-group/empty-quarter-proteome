#!/usr/bin/env python3
"""step3_00_prepare - build the query/target FASTAs for step three.

Outputs (all under P):
  results/step3_missing_bakta.txt      catalogue genomes with no Bakta proteome (must be empty)
  results/step3_cultured_counts.tsv    per genome: n_proteins, n_hypothetical
  tmp/cultured_all.faa                 all 722 Bakta proteomes, headers '>{genome}@{locus}'
  tmp/cultured_hypo_ids.txt            the hypothetical-protein ids in the same namespace
  tmp/cultured_chunks/chunk_NN.faa     query chunks for the UniRef50 array
  tmp/dark_chunks/chunk_NN.faa         N_DARK chunks of the 4,292,446 MAG-dark reps

Bakta roots. The proteomes live in two trees and BAKTA_ROOTS is searched in order, the
first root that holds a genome owning it. Root 0 is the eqmgc annotation tree (symlinks
into the enrichment and isolate batches, 651 genomes). Root 1 is the soil-batch
thermophilic tree, which holds the remaining 71 plus 11 further Bakta runs that are not
catalogue genomes. A directory that resolves to no catalogue id is skipped and counted;
the run fails if any catalogue genome ends up without a proteome.

Naming: catalogue.tsv calls five single-bin enrichments '<sample>' while the Bakta
directory calls them '<sample>.1'. Those five are aliases of the same genome and are
mapped back to the catalogue id here, so every downstream join is on catalogue ids.

Chunk boundaries are fixed per root, not recomputed over the concatenated file, so that
adding a root cannot shift the chunks an earlier UniRef50 array already searched. Root 0
keeps its N_CULT chunks (00 to 11, unchanged to the byte); every later root is split into
ceil(n_root / per_0) chunks of its own, numbered on from there, where per_0 is root 0's
records per chunk. Chunk files left over from a shorter previous plan are removed.
"""
import os
import sys

P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean"
E = "/ibex/scratch/projects/c2014/alelopezv/eqmgc"
W = "/ibex/scratch/projects/c2014/rob/metagenomics"
A = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/cultures/annotation"

BAKTA_ROOTS = [
    f"{E}/results/annotation/bakta",      # enrichment + isolate batches, 651 genomes
    f"{A}/thermophilic/bakta",            # soil batch, the remaining 71 genomes
]
DARK_FAA = f"{W}/06_darkmatter/darkfam_basis_clean/results/dark_reps.faa"
N_CULT = 12
N_DARK = 20

os.makedirs(f"{P}/tmp/cultured_chunks", exist_ok=True)
os.makedirs(f"{P}/tmp/dark_chunks", exist_ok=True)

# catalogue ids
cat_ids = []
with open(f"{E}/catalogue.tsv") as fh:
    header = fh.readline()
    for line in fh:
        cat_ids.append(line.split("\t")[0])
cat_set = set(cat_ids)
print(f"[info] catalogue genomes: {len(cat_ids)}")

# bakta dir -> catalogue id, root by root; a genome belongs to the first root that holds it
proteomes = []           # (root, dir, catalogue id) in write order
have = set()
for root in BAKTA_ROOTS:
    dirs = sorted(d for d in os.listdir(root) if os.path.isdir(f"{root}/{d}"))
    taken, skipped = 0, 0
    for d in dirs:
        if d in cat_set:
            g = d
        elif d.endswith(".1") and d[:-2] in cat_set:
            g = d[:-2]
        else:
            skipped += 1
            continue
        if g in have:
            continue
        have.add(g)
        proteomes.append((root, d, g))
        taken += 1
    print(f"[info] {root}: {len(dirs)} dirs, {taken} catalogue genomes taken, "
          f"{skipped} not in the catalogue")
missing = [g for g in cat_ids if g not in have]
print(f"[info] catalogue genomes covered {len(have)}; missing {len(missing)}")

with open(f"{P}/results/step3_missing_bakta.txt", "w") as out:
    for g in missing:
        out.write(g + "\n")
if missing:
    sys.exit(f"ERROR: {len(missing)} catalogue genomes have no Bakta proteome "
             f"(see results/step3_missing_bakta.txt)")

# concatenate the proteomes, tagging every protein with its genome
counts = {}
n_prot = 0
root_prot = {root: 0 for root in BAKTA_ROOTS}
with open(f"{P}/tmp/cultured_all.faa", "w") as out, \
     open(f"{P}/tmp/cultured_hypo_ids.txt", "w") as hyo:
    for root, d, g in proteomes:
        npro = 0
        with open(f"{root}/{d}/{d}.faa") as fh:
            for line in fh:
                if line[0] == ">":
                    locus = line[1:].split()[0]
                    out.write(f">{g}@{locus}\n")
                    npro += 1
                else:
                    out.write(line)
        nhyp = 0
        hp = f"{root}/{d}/{d}.hypotheticals.faa"
        if os.path.exists(hp):
            with open(hp) as fh:
                for line in fh:
                    if line[0] == ">":
                        hyo.write(f"{g}@{line[1:].split()[0]}\n")
                        nhyp += 1
        counts[g] = (npro, nhyp)
        n_prot += npro
        root_prot[root] += npro
print(f"[info] cultured genomes: {len(counts)}; cultured proteins: {n_prot}")

with open(f"{P}/results/step3_cultured_counts.tsv", "w") as out:
    out.write("genome\tn_proteins\tn_hypothetical\n")
    for g in cat_ids:
        if g in counts:
            out.write(f"{g}\t{counts[g][0]}\t{counts[g][1]}\n")


def split_sizes(total, n):
    """n chunks of ceil(total/n) records, the last one taking what is left."""
    per = (total + n - 1) // n
    sizes, left = [], total
    for _ in range(n):
        sizes.append(min(per, left))
        left -= sizes[-1]
    return [s for s in sizes if s > 0]


def write_chunks(src, outdir, sizes):
    """Write src into len(sizes) chunks with exactly those record counts."""
    idx, seen = 0, 0
    out = open(f"{outdir}/chunk_{idx:02d}.faa", "w")
    with open(src) as fh:
        for line in fh:
            if line[0] == ">":
                if seen == sizes[idx] and idx + 1 < len(sizes):
                    out.close()
                    idx += 1
                    seen = 0
                    out = open(f"{outdir}/chunk_{idx:02d}.faa", "w")
                seen += 1
            out.write(line)
    out.close()
    for stale in sorted(f for f in os.listdir(outdir)
                        if f.startswith("chunk_") and f.endswith(".faa")):
        if int(stale[len("chunk_"):-len(".faa")]) >= len(sizes):
            os.remove(f"{outdir}/{stale}")
            print(f"[info] removed stale {outdir}/{stale}")
    print(f"[info] {src}: {sum(sizes)} records -> {len(sizes)} chunks, sizes {sizes}")


def count_records(path):
    with open(path) as fh:
        return sum(1 for line in fh if line[0] == ">")


# root 0 keeps its N_CULT chunks; later roots are split at the same records per chunk
cult_sizes = split_sizes(root_prot[BAKTA_ROOTS[0]], N_CULT)
per0 = cult_sizes[0]
for root in BAKTA_ROOTS[1:]:
    n = root_prot[root]
    if n:
        cult_sizes += split_sizes(n, (n + per0 - 1) // per0)
write_chunks(f"{P}/tmp/cultured_all.faa", f"{P}/tmp/cultured_chunks", cult_sizes)
write_chunks(DARK_FAA, f"{P}/tmp/dark_chunks", split_sizes(count_records(DARK_FAA), N_DARK))
print("DONE_STEP3_00")
