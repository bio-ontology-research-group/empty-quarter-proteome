#!/bin/bash
# REBUILD step 1 — dark-family universe from the paper's own dark call (FINAL design 2026-08-13):
# extract the 4,292,446 MAG-mapped dark cluster90 reps (darkclean_5; §1.2's shipped object) from
# all_reps90.faa.gz, then mmseqs easy-cluster 50% id / 80% cov (same convention as the June
# orthogroups, 74_mmseqs_orthogroups.sh). Families = clusters; recurrence layers added in step 2.
# Supersedes the strict re-search chain (basis_1..3, never run) — darkness needs no new search here.
#
# Submit: sbatch scripts/fam_1_cluster.sh
#
#SBATCH --job-name=ale_fam_cluster
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam1_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=batch

set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
CAT=/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog
IDS=$CAT/darkmatter_clean/05_mag_dark/mag_dark_rep_ids.txt.gz
FAA=$CAT/all_linclust90/all_reps90.faa.gz
R=$B/results
mkdir -p "$R/fam_tmp"

DARK=$R/dark_reps.faa
if [ ! -s "$DARK" ]; then
    echo "[$(date -Iseconds)] extracting 4,292,446 dark reps from all_reps90.faa.gz (one pass)"
    zcat "$FAA" | awk 'NR==FNR{ids[$1];next} /^>/{h=substr($1,2); p=(h in ids)} p' \
        <(zcat "$IDS") - > "$DARK.part"
    mv "$DARK.part" "$DARK"
fi
n=$(grep -c '^>' "$DARK")
echo "[$(date -Iseconds)] extracted: $n (expect 4292446)"
[ "$n" -eq 4292446 ] || { echo "ERROR: extraction count mismatch"; exit 1; }

source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate mmseqs2
mmseqs version || true
echo "[$(date -Iseconds)] mmseqs easy-cluster 50% id / 80% cov"
mmseqs easy-cluster "$DARK" "$R/fam" "$R/fam_tmp" \
    --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads "$SLURM_CPUS_PER_TASK"

echo "[$(date -Iseconds)] family size distribution (members = dark cluster90 reps)"
cut -f1 "$R/fam_cluster.tsv" | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c \
    > "$R/family_size_hist.txt"
awk '{tot+=$1; if($2>=3) ge3+=$1} END{print "clusters total:", tot, "  with >=3 members:", ge3}' \
    "$R/family_size_hist.txt"
echo "[$(date -Iseconds)] step 1 done"
