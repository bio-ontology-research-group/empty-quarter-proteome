#!/usr/bin/env python3
"""v2_thermo_tree — phylogenomic placement of the 237 desert Thermoplasmatota.

Purpose: the archaeal section calls this group a radiation. That word needs a tree. This builds one
from the ar53 marker alignment GTDB-Tk already produced, which contains our genomes and the GTDB
reference set aligned together, so no new alignment is needed.

Taxa: our 237 Thermoplasmatota representatives + all 1,647 GTDB r232 Thermoplasmatota species
representatives + an outgroup of 30 Halobacteriota and Thermoproteota representatives.
Context worth knowing: in GTDB r232 the class SW-10-69-26 holds 21 genomes worldwide and E2 holds
197, and Penumbrarchaeia, the class in Maeke 2025, holds 37.

Tree: FastTreeMP, LG model with gamma, approximate likelihood ratio supports. This is an
approximate-ML tree and should be described as such; it answers whether the desert genomes form
coherent groups, not fine branching order.

Inputs (read-only): the GTDB-Tk ar53 MSA and the r232 archaeal taxonomy.
Outputs: results/thermo_tree.{fasta,nwk,taxa.tsv,log}
"""
import os, gzip, csv, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
W = "/ibex/scratch/projects/c2014/rob/metagenomics"
MSA = (f"{W}/rebuild_20260729/coassembly_novelty/legacy_taxonomic_screen_r232/"
       "gtdbtk_2.7.2_r232/run/align/gtdbtk.ar53.msa.fasta.gz")  # r232 run (2026-09-21); was the r226 run's MSA
TAX = f"{W}/rebuild_20260729/references/gtdbtk_r232/taxonomy/ar53_taxonomy_r232_reps.tsv"
OURS = os.path.join(RES, "v2_reps.tsv")
FASTTREE = "/home/velazqal/miniconda3/envs/phabox2/bin/FastTreeMP"
N_OUTGROUP = 30


def main():
    ours = set()
    with open(OURS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r["phylum"] or "").replace("p__", "") == "Thermoplasmatota":
                ours.add(r["genome"])
    print(f"our Thermoplasmatota representatives: {len(ours)}", flush=True)

    ref_thermo, outgroup, lineage = set(), [], {}
    with open(TAX) as fh:
        for line in fh:
            acc, tax = line.rstrip("\n").split("\t")[:2]
            lineage[acc] = tax
            if "p__Thermoplasmatota" in tax:
                ref_thermo.add(acc)
            elif ("p__Halobacteriota" in tax or "p__Thermoproteota" in tax) and len(outgroup) < N_OUTGROUP:
                outgroup.append(acc)
    keep = ours | ref_thermo | set(outgroup)
    print(f"GTDB Thermoplasmatota reps {len(ref_thermo)}, outgroup {len(outgroup)}, "
          f"target taxa {len(keep)}", flush=True)

    out_fa = os.path.join(RES, "thermo_tree.fasta")
    written, name = 0, None
    with gzip.open(MSA, "rt") as fh, open(out_fa, "w") as o:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if name in keep:
                    o.write(f">{name}\n")
                    written += 1
                    name = name
                else:
                    name = None
            elif name:
                o.write(line)
    print(f"sequences written: {written} of {len(keep)} requested", flush=True)

    with open(os.path.join(RES, "thermo_tree_taxa.tsv"), "w") as o:
        o.write("taxon\tset\tgtdb_lineage\n")
        for g in sorted(keep):
            s = ("desert" if g in ours else
                 "gtdb_thermoplasmatota" if g in ref_thermo else "outgroup")
            o.write(f"{g}\t{s}\t{lineage.get(g, '')}\n")

    # The tree itself is a Slurm job, never a login-node run: submit
    # `sbatch scripts/v2_thermo_tree.sbatch` once this alignment exists.
    print(f"alignment ready: {out_fa}\nnow run: sbatch scripts/v2_thermo_tree.sbatch", flush=True)


if __name__ == "__main__":
    main()
