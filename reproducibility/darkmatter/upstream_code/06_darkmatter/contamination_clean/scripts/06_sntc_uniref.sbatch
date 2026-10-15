#!/bin/bash
# 06 - what organism is the negative control's 11 kb of assembled sequence? Step 05 can only say
# whether SNTC contigs are in the catalogue; it cannot name anything outside it. This calls genes on
# the 11 contigs and searches them against UniRef50, whose headers carry a Tax= field, so the blank
# gets a name rather than only an absence.
#
#   sbatch scripts/06_sntc_uniref.sbatch
#
# UniRef50 mmseqs2 DB and search parameters are the paper's own, from
# 06_darkmatter/cultured_comparison_clean/scripts/step3_03_cultured_vs_uniref50.sbatch, which in turn
# takes them verbatim from darkclean_3b_uniref.sh, the rule that defines MAG darkness:
# -s 7.5 -e 1e-3 --max-seqs 10. The query here is ~a dozen proteins, not 2.7 M, so it runs unchunked.
# theader is added to the output format because the organism name is the point of this step.
#SBATCH --job-name=ale_sntc_uniref
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/sntc_uniref_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/sntc_uniref_%j.err
set -euo pipefail

M=/ibex/scratch/projects/c2014/rob/metagenomics
P=$M/06_darkmatter/contamination_clean
ASM=/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/assemblies
UNIREF_DB=$M/04_structure/uniref50/uniref50
WORK=$P/tmp/sntc_uniref
mkdir -p "$WORK/tmp" "$P/results/identify"

[ -f "${UNIREF_DB}.dbtype" ] || { echo "ERROR: UniRef50 DB not found at $UNIREF_DB" >&2; exit 2; }

source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh

# ---------------------------------------------------------------- genes on 11 kb of contigs
# -p meta: the contigs are 561 to 3,171 bp, far too short for prodigal to train on.
conda activate gtdbtk
prodigal -v 2>&1 | head -2
prodigal -p meta -i "$ASM/SNTC/SNTC_contigs.fa" \
    -a "$WORK/SNTC.faa" -d "$WORK/SNTC.ffn" -f gff -o "$WORK/SNTC.gff" -q
echo "proteins called: $(command grep -c '^>' "$WORK/SNTC.faa")"

# ---------------------------------------------------------------- UniRef50
conda activate mmseqs2
mmseqs version
rm -rf "$WORK/query.mmseqsdb"* "$WORK/result.mmseqsdb"*
mmseqs createdb "$WORK/SNTC.faa" "$WORK/query.mmseqsdb"
mmseqs search "$WORK/query.mmseqsdb" "$UNIREF_DB" "$WORK/result.mmseqsdb" "$WORK/tmp" \
    -s 7.5 -e 1e-3 --max-seqs 10 --threads "${SLURM_CPUS_PER_TASK:-16}" --split-memory-limit 200G
mmseqs convertalis "$WORK/query.mmseqsdb" "$UNIREF_DB" "$WORK/result.mmseqsdb" \
    "$P/results/identify/sntc_uniref50_hits.tsv" \
    --format-output 'query,target,fident,alnlen,evalue,bits,theader' \
    --threads "${SLURM_CPUS_PER_TASK:-16}"
rm -rf "$WORK/tmp"
echo "hit rows: $(wc -l < "$P/results/identify/sntc_uniref50_hits.tsv")"
echo "=== DONE ==="
