#!/bin/bash
# 09 - control of the control. Step 04 reports 0 of 119,683,914 SPositive reads mapped to the
# catalogue, while step 05 finds 77 contig alignments of at least 150 bp and at least 90 % identity
# between the SPositive assembly and the catalogue, several of them to the E. coli representative's
# own contigs. Those two facts are not obviously compatible, so the zero is checked directly before
# it is written down.
#
#   sbatch scripts/09_sanity_spositive_zero.sbatch
#
# Three mappings of 1,000,000 R1 reads each, same index, minimap2 called directly so that CoverM's
# accounting is not in the path:
#   SNTC        positive control for the method. Step 04 puts it at 8.56 %, so a comparable figure
#               here means the index, the preset and the counting all work.
#   SPositive   the number in question.
#   SPositive contigs, reads simulated end to end from the 77 good alignment blocks would be the
#               stronger test; instead the blocks themselves are re-aligned at read length by
#               cutting the SPositive assembly into 98 bp fragments, the library's own read length.
#               If 98 bp of SPositive sequence cannot match the catalogue even when taken straight
#               from the assembly, the zero is a property of short-read identity, not a failure.
#SBATCH --job-name=ale_sanity_zero
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=04:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/sanity_zero_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/logs/sanity_zero_%j.err
set -euo pipefail

M=/ibex/scratch/projects/c2014/rob/metagenomics
P=$M/06_darkmatter/contamination_clean
PREP=$M/rebuild_20260729/picrust2_current_20260803/coverm_current/prepared
READS=/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads
MM2=/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/minimap2
SAM=/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/samtools
T=${SLURM_CPUS_PER_TASK:-16}
WORK=$P/tmp/sanity; mkdir -p "$WORK"

# 1,000,000 reads = 4,000,000 fastq lines
for S in SNTC SPositive; do
    zcat "$READS/${S}_R1.fastq.gz" | head -4000000 > "$WORK/${S}_1M.fq" || true
    echo "$S subsample reads: $(( $(wc -l < "$WORK/${S}_1M.fq") / 4 ))"
done

# SPositive assembly cut into 98 bp non-overlapping fragments, the library's own read length
awk '/^>/{next}{seq=seq $0} END{
        n=0
        for (i = 1; i + 97 <= length(seq); i += 98) { n++; print ">frag_" n; print substr(seq, i, 98) }
     }' "$P/tmp/control_genomes/SPositive.fa" > "$WORK/SPositive_frag98.fa"
echo "SPositive 98 bp fragments: $(command grep -c '^>' "$WORK/SPositive_frag98.fa")"

for Q in "$WORK/SNTC_1M.fq" "$WORK/SPositive_1M.fq" "$WORK/SPositive_frag98.fa"; do
    echo "=== query: $(basename "$Q") ==="; date
    $MM2 -a -x sr -t "$T" "$PREP/reps.mmi" "$Q" 2>> "$P/logs/sanity_mm2.log" \
        | $SAM view -c -F 0x904 - > "$WORK/$(basename "$Q").mapped"
    echo "primary mapped records: $(cat "$WORK/$(basename "$Q").mapped")"
done
date
echo "=== DONE ==="
