#!/bin/bash
# v2_thermo_relatives - ANI and AAI of the 237 desert Thermoplasmatota against GTDB r232 and MGnify.
#   step 0 (login, network): python scripts/v2_thermo_relatives_fetch.py
#   step 1: J1=$(sbatch --parsable scripts/v2_thermo_relatives.sbatch prep)
#           skani 0.3.1 triangle over all genomes (same binary and -s 75 as
#           cultured_comparison_clean/scripts/01_skani.sbatch), Pyrodigal proteins, DIAMOND db
#   step 2: J2=$(sbatch --parsable --dependency=afterok:$J1 --array=0-19 --cpus-per-task=16 \
#                  --mem=48G scripts/v2_thermo_relatives.sbatch blast)
#   step 3: sbatch --dependency=afterok:$J2 --cpus-per-task=4 --mem=96G \
#                  scripts/v2_thermo_relatives.sbatch aai
#   step 4 (login, seconds): python scripts/v2_thermo_relatives.py
#SBATCH --job-name=ale_thermo_rel
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean/logs/thermo_relatives_%x_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean/logs/thermo_relatives_%x_%A_%a.err
set -euo pipefail

STAGE=${1:?stage required: prep | blast | aai}
P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean
T=$P/tmp/thermo_relatives
PY=/home/velazqal/miniconda3/envs/darkfunc/bin/python          # pyrodigal 3.7.1, numpy
SK=/ibex/sw/rl9c/gtdb-tk/2.5.2/rl9.4_conda3/Miniconda3/envs/gtdbtk_v2.5.2/bin/skani
N=${SLURM_CPUS_PER_TASK:-32}
test -s $T/genome_paths.tsv || { echo "run v2_thermo_relatives_fetch.py first" >&2; exit 1; }
date; echo "stage $STAGE, cpus $N"

case $STAGE in
prep)
    echo "skani: $($SK --version 2>&1 | head -1)"
    tail -n +2 $T/genome_paths.tsv | cut -f4 > $T/skani_paths.txt
    echo "genomes: $(wc -l < $T/skani_paths.txt)"
    $SK triangle -l $T/skani_paths.txt --sparse -t $N -s 75 > $T/skani_triangle.tsv
    echo "skani rows: $(($(wc -l < $T/skani_triangle.tsv) - 1))"
    gzip -c $T/skani_triangle.tsv > $P/results/v2_thermo_relatives_skani.tsv.gz
    $PY $P/scripts/v2_thermo_relatives_aai.py proteins $N
    module load diamond/2.1.16
    diamond --version
    diamond makedb --in $T/faa/all.faa -d $T/faa/all --threads $N
    ;;
blast)
    module load diamond/2.1.16
    C=$(printf "%02d" ${SLURM_ARRAY_TASK_ID:?array task id required})
    diamond blastp -q $T/faa/chunk_$C.faa -d $T/faa/all.dmnd --very-sensitive \
        -e 1e-5 --id 30 --query-cover 70 --subject-cover 70 --max-target-seqs 0 \
        --outfmt 6 qseqid sseqid pident bitscore --threads $N --tmpdir $T --quiet \
      | $PY $P/scripts/v2_thermo_relatives_aai.py reduce $C
    ;;
aai)
    $PY $P/scripts/v2_thermo_relatives_aai.py aai
    ;;
*)
    echo "unknown stage $STAGE" >&2; exit 1 ;;
esac
date; echo "stage $STAGE done"
