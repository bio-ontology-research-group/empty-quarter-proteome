#!/bin/bash
# mimag_rrna_trna_clean/scripts/00_prepare_inputs.sh
#
# Builds the two small inputs the array job reads:
#   results/genome_list.txt  — the 5,438 representative names, one per line, sorted
#                              (fixed order = the chunking order for the array job)
#   results/domain_map.tsv   — rep <TAB> domain (Bacteria|Archaea), from the GTDB-Tk r232
#                              bac120/ar53 summaries that already carry all 5,438 reps
#
# Read-only against 05_refine/; writes only inside mimag_rrna_trna_clean/. Idempotent
# (regenerates both files each time it is run; cheap, login-node only).
set -euo pipefail

REFINE=/ibex/scratch/projects/c2014/rob/metagenomics/05_refine
BASE=/ibex/scratch/projects/c2014/rob/metagenomics/mimag_rrna_trna_clean
FASTA_DIR="$REFINE/_drep/out/dereplicated_genomes"
BAC="$REFINE/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
AR="$REFINE/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"

mkdir -p "$BASE/results"

# genome list: strip .fa, sort for a stable, reproducible chunk order
ls "$FASTA_DIR" | sed 's/\.fa$//' | sort > "$BASE/results/genome_list.txt"

n=$(wc -l < "$BASE/results/genome_list.txt")
echo "genome_list.txt: $n representatives"

# domain map: user_genome -> Bacteria / Archaea, from the classification column's d__ prefix
{
  awk -F'\t' 'NR>1{split($2,a,";"); dom=a[1]; sub("d__","",dom); print $1"\t"dom}' "$BAC"
  awk -F'\t' 'NR>1{split($2,a,";"); dom=a[1]; sub("d__","",dom); print $1"\t"dom}' "$AR"
} | sort -k1,1 > "$BASE/results/domain_map.tsv"

m=$(wc -l < "$BASE/results/domain_map.tsv")
echo "domain_map.tsv: $m rows"

# sanity: every genome in genome_list.txt must have a domain row, and counts must match 5,438
missing=$(comm -23 "$BASE/results/genome_list.txt" <(cut -f1 "$BASE/results/domain_map.tsv") | wc -l)
echo "genomes with no domain match: $missing"
if [ "$n" -ne "$m" ] || [ "$missing" -ne 0 ]; then
  echo "FATAL: genome_list ($n) and domain_map ($m) do not match 1:1 (missing=$missing)" >&2
  exit 1
fi
echo "OK: $n genomes, $m domain rows, 0 missing"
