#!/bin/bash
#SBATCH --job-name=darkfam_taxa
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology/taxa/taxa_%j.out
# Resolve member->MAG->GTDB taxon for every multi-member dark family.
# orthogroups_cluster.tsv is a 2-col mmseqs file (rep<TAB>member); both ids are
# <mag>::<protein>, so the member's MAG is the substring before '::'. We stream the
# 20.4M-line file once with a 23,382-key hash, then join MAGs to GTDB-Tk and emit a
# compact per-family taxonomy summary (broad/narrow, dominant phylum/genus, top lineage).
set -euo pipefail
BASE=/ibex/scratch/projects/c2014/rob/metagenomics
DM=$BASE/06_darkmatter
ECO=$DM/darkfam_ecology
OUT=$ECO/taxa
mkdir -p "$OUT"
ORTHO=$DM/orthogroups/orthogroups_cluster.tsv
META=$ECO/catalog/darkfam_meta.tsv
# NB: the orthogroup members are dRep representative proteomes named in the
# self_/single_SemiBin namespace -> use the REFINED GTDB summary (same namespace,
# 5,438 reps), NOT 02_taxonomy/gtdbtk_out (renamed full 15,469-MAG set -> only 45% match).
GBAC=$BASE/05_refine/_gtdbtk_refined/gtdbtk.bac120.summary.tsv
GAR=$BASE/05_refine/_gtdbtk_refined/gtdbtk.ar53.summary.tsv

echo "[$(date)] extracting members for $(($(wc -l < "$META")-1)) families from $ORTHO"
awk -F'\t' 'FNR==NR{ if(FNR>1) want[$1]=1; next }
            ($1 in want){ p=index($2,"::"); print $1"\t"substr($2,1,p-1) }' \
    "$META" "$ORTHO" > "$OUT/fam_mag.tsv"
echo "[$(date)] fam_mag rows: $(wc -l < "$OUT/fam_mag.tsv")"

python3 "$OUT/99_summarize_taxa.py" "$OUT/fam_mag.tsv" "$GBAC" "$GAR" "$OUT/fam_taxa_summary.tsv"
echo "[$(date)] DONE -> $OUT/fam_taxa_summary.tsv"
touch "$OUT/TAXA_DONE"
