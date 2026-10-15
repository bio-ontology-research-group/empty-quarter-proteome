#!/bin/bash
# REBUILD step 3 — member CDS extraction (shared input for CoverM leg + dN/dS leg).
# For every member variant of every ADMITTED family (~540k proteins), cut its nucleotide gene
# sequence out of the sample's assembly contigs, using the prodigal coordinates carried in the
# protein FASTA headers (start/end/strand). Output: results/darkfam_members.fna + member→family map.
#
# Submit: sbatch scripts/fam_3_cds.sh
#
#SBATCH --job-name=ale_fam_cds
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam3_%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=batch

set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
AP=/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog/assembly_proteins
ASM=/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/assemblies
R=$B/results

# member -> family map for admitted families only
export LC_ALL=C
awk -F'\t' 'NR==FNR{adm[$1];next} $1 in adm{print $2"\t"$1}' \
    "$R/FAMILIES_FINAL.txt" "$R/fam_cluster.tsv" > "$R/member_fam.tsv"
echo "[$(date -Iseconds)] admitted-family members: $(wc -l < "$R/member_fam.tsv")"

mkdir -p "$R/cds_parts"
cut -f1 "$R/member_fam.tsv" | awk -F'::' '{print $1}' | sort -u > "$R/_samples_needed.txt"
echo "[$(date -Iseconds)] samples with members: $(wc -l < "$R/_samples_needed.txt")"

slice_one(){  # sample — python: harvest coords from faa.gz headers, slice contigs.fa
python3 - "$1" <<'PY'
import sys, gzip, os
s = sys.argv[1]
B = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results"
faa = f"/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog/assembly_proteins/{s}/{s}.faa.gz"
asm = f"/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/assemblies/{s}/{s}_contigs.fa"
out = f"{B}/cds_parts/{s}.fna"
if os.path.exists(out) and os.path.getsize(out) > 0: sys.exit(0)
need = {}
for line in open(f"{B}/member_fam.tsv"):
    m, f = line.rstrip("\n").split("\t")
    if m.startswith(s + "::"): need[m] = f
if not need: open(out, "w").close(); sys.exit(0)
coords = {}
with gzip.open(faa, "rt") as fh:
    for line in fh:
        if not line.startswith(">"): continue
        parts = line[1:].split(" # ")
        pid = parts[0].split()[0]
        if pid in need:
            coords[pid] = (int(parts[1]), int(parts[2]), int(parts[3]))
missing = len(need) - len(coords)
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
print(f"[{s}] members {len(need)} coords {len(coords)} written {wrote} missing-coords {missing}", flush=True)
PY
}
export -f slice_one
xargs -P "$SLURM_CPUS_PER_TASK" -I{} bash -c 'slice_one "$@"' _ {} < "$R/_samples_needed.txt"

cat "$R"/cds_parts/*.fna > "$R/darkfam_members.fna"
n=$(grep -c '^>' "$R/darkfam_members.fna")
echo "[$(date -Iseconds)] total member CDS: $n (expect ~$(wc -l < "$R/member_fam.tsv"))"
echo "[$(date -Iseconds)] step 3 done"
