#!/bin/bash
# 05 - what are the two control libraries? Three independent identifications of the control
#      assemblies, run in one job:
#
#   A  GTDB-Tk 2.7.2 classify_wf against the project's own R232 reference. Same software version,
#      same reference and the same flags as the run that produced 05_refine/_gtdbtk_r232, so the
#      answer is directly comparable to the catalogue's own taxonomy.
#   B  skani whole-genome ANI against the 5,438 species representatives. Answers "is this organism
#      already in the catalogue" independently of GTDB.
#   C  minimap2 of the control contigs against the catalogue's own concatenated reference. This is
#      the only one of the three that can say anything about SNTC, whose 11 contigs total 11 kb and
#      are far below the length where whole-genome ANI or marker placement mean anything.
#
#   sbatch scripts/05_identify_controls.sbatch
#
# Provenance of the reused invocations:
#   A  rebuild_20260729/coassembly_novelty/legacy_taxonomic_screen_r232/gtdbtk_2.7.2_r232/run/gtdbtk.json
#      command_line: classify_wf --genome_dir ... --extension fa --cpus 64 --pplacer_cpus 32
#      --keep_intermediates --write_single_copy_genes. CPU counts scaled to this job; the two
#      intermediate flags are dropped because nothing downstream reads their output for 2 genomes.
#   B  06_darkmatter/cultured_comparison_clean/scripts/01_skani.sbatch, skani 0.3.1, -s 75, whole
#      files via --ql (no --qi, which in 0.3.1 means per-contig queries).
#   C  the index and the contig-to-genome map built by rebuild_20260729/scripts/prepare_current_coverm.py.
#SBATCH --job-name=ale_ctrl_id
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/identify_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/identify_%j.err
set -uo pipefail

M=/ibex/scratch/projects/c2014/rob/metagenomics
P=$M/06_darkmatter/contamination_clean
ASM=/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/assemblies
REPS=$M/05_refine/_drep/out/dereplicated_genomes
PREP=$M/rebuild_20260729/picrust2_current_20260803/coverm_current/prepared
SK=/ibex/sw/rl9c/gtdb-tk/2.5.2/rl9.4_conda3/Miniconda3/envs/gtdbtk_v2.5.2/bin/skani
MM2=/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/minimap2
T=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$P/tmp/control_genomes" "$P/results/identify"

# ---------------------------------------------------------------- stage inputs
cp -f "$ASM/SNTC/SNTC_contigs.fa"         "$P/tmp/control_genomes/SNTC.fa"
cp -f "$ASM/SPositive/SPositive_contigs.fa" "$P/tmp/control_genomes/SPositive.fa"
ls -l "$P/tmp/control_genomes"
find "$REPS" -maxdepth 1 -name '*.fa' | sort > "$P/tmp/mag_paths.txt"
n=$(wc -l < "$P/tmp/mag_paths.txt")
[ "$n" -eq 5438 ] || { echo "FAIL: expected 5438 representative .fa, got $n" >&2; exit 1; }
printf '%s\n' "$P/tmp/control_genomes/SNTC.fa" "$P/tmp/control_genomes/SPositive.fa" \
  > "$P/tmp/control_paths.txt"

# ---------------------------------------------------------------- A: GTDB-Tk 2.7.2 / R232
echo "=== A: GTDB-Tk classify_wf, R232 ==="; date
module load gtdb-tk/2.7.2
export GTDBTK_DATA_PATH=$M/rebuild_20260729/references/gtdbtk_r232
gtdbtk --version
rm -rf "$P/results/identify/gtdbtk" "$P/tmp/gtdbtk_tmp"
mkdir -p "$P/tmp/gtdbtk_tmp"
gtdbtk classify_wf \
    --genome_dir "$P/tmp/control_genomes" \
    --out_dir "$P/results/identify/gtdbtk" \
    --extension fa \
    --cpus "$T" \
    --pplacer_cpus 16 \
    --tmpdir "$P/tmp/gtdbtk_tmp"
echo "gtdbtk exit $?"
for f in "$P/results/identify/gtdbtk"/gtdbtk.*.summary.tsv; do
    [ -s "$f" ] && { echo "--- $f"; cut -f1,2,6,7,15 "$f"; }
done
[ -s "$P/results/identify/gtdbtk/gtdbtk.bac120.filtered.tsv" ] && \
    { echo "--- filtered out (too few markers):"; cat "$P/results/identify/gtdbtk/gtdbtk.bac120.filtered.tsv"; }
date

# ---------------------------------------------------------------- B: skani vs the 5,438
echo "=== B: skani ANI, controls vs the 5,438 representatives ==="; date
$SK --version
$SK dist --ql "$P/tmp/control_paths.txt" --rl "$P/tmp/mag_paths.txt" \
    -o "$P/results/identify/skani_controls_vs_reps.tsv" -t "$T" -s 75
echo "skani exit $?; rows: $(wc -l < "$P/results/identify/skani_controls_vs_reps.tsv")"
head -20 "$P/results/identify/skani_controls_vs_reps.tsv"
date

# ---------------------------------------------------------------- C: contigs vs the catalogue
# Same reps.mmi (minimap2 -x sr -I 64G over the concatenated 5,438) that CoverM used in step 04, so
# a contig hit and a read hit are resolved against exactly the same sequence set. The sr preset is
# fixed by the prebuilt index and cannot be changed here; it is a short-read preset applied to
# contigs, which is why this stage is read as presence/absence of catalogue sequence in a contig,
# not as an alignment-quality statistic.
echo "=== C: minimap2, control contigs vs the concatenated catalogue ==="; date
for S in SNTC SPositive; do
    $MM2 -x sr -t "$T" --secondary=no -c \
        "$PREP/reps.mmi" "$P/tmp/control_genomes/$S.fa" \
        > "$P/results/identify/${S}_contigs_vs_reps.paf" 2> "$P/logs/mm2_${S}.log"
    echo "$S: $(wc -l < "$P/results/identify/${S}_contigs_vs_reps.paf") PAF rows"
done
date
echo "=== DONE ==="
