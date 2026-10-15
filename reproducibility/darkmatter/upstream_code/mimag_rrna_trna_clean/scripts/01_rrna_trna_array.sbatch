#!/bin/bash -l
#SBATCH --job-name=mimag_rrna_trna
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-108
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/mimag_rrna_trna_clean/logs/rrna_trna_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/mimag_rrna_trna_clean/logs/rrna_trna_%A_%a.err
#
# mimag_rrna_trna_clean/scripts/01_rrna_trna_array.sbatch
#
# rRNA (barrnap 0.9) and tRNA (tRNAscan-SE 2.0.12) detection on the 5,438 species
# representative genomes, chunked ~50 genomes per array task (109 tasks x 50 = 5,450 >= 5,438).
# Kingdom/search-mode is set per genome from its GTDB-Tk r232 domain (domain_map.tsv) — every
# one of the 5,438 has a domain call, so no both-kingdoms fallback is needed here.
#
# Idempotent: a genome is skipped if its two output files (rRNA GFF, tRNA table) already exist
# and are non-empty. Safe to resubmit for failed/incomplete tasks.
#
# Inputs (read-only):
#   FASTA_DIR   = 05_refine/_drep/out/dereplicated_genomes/<rep>.fa   (5,438 representative FASTAs)
#   GENOME_LIST = results/genome_list.txt   (5,438 rep names, sorted; fixed chunk order)
#   DOMAIN_MAP  = results/domain_map.tsv    (rep -> Bacteria|Archaea)
#
# Outputs (small, per-genome only — no genome copies):
#   results/rrna/<rep>.gff            barrnap GFF3
#   results/trna/<rep>.trnascan.txt   tRNAscan-SE tabular output

set -uo pipefail

BASE=/ibex/scratch/projects/c2014/rob/metagenomics/mimag_rrna_trna_clean
FASTA_DIR=/ibex/scratch/projects/c2014/rob/metagenomics/05_refine/_drep/out/dereplicated_genomes
GENOME_LIST="$BASE/results/genome_list.txt"
DOMAIN_MAP="$BASE/results/domain_map.tsv"
RRNA_OUT="$BASE/results/rrna"
TRNA_OUT="$BASE/results/trna"
CHUNK_SIZE=50

mkdir -p "$RRNA_OUT" "$TRNA_OUT"

TOTAL=$(wc -l < "$GENOME_LIST")
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE + 1 ))
END=$(( START + CHUNK_SIZE - 1 ))
if [ "$END" -gt "$TOTAL" ]; then END=$TOTAL; fi
if [ "$START" -gt "$TOTAL" ]; then
  echo "Task $SLURM_ARRAY_TASK_ID: start ($START) beyond total ($TOTAL), nothing to do"
  exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: genomes $START-$END of $TOTAL"

mapfile -t CHUNK < <(sed -n "${START},${END}p" "$GENOME_LIST")

# barrnap (perl/hmmer/bedtools) and tRNAscan-SE (its own conda perl) collide if both modules
# are loaded in the same shell (mismatched perl library handshake). Each tool therefore runs
# in its own subshell that does a fresh module purge + load; subshell env never leaks to the
# parent loop or to the other tool (verified interactively before submitting this array).
echo "barrnap: $(module purge >/dev/null 2>&1; module load barrnap/0.9 >/dev/null 2>&1; echo "$(command -v barrnap) ($(barrnap --version 2>&1))")"
echo "tRNAscan-SE: $(module purge >/dev/null 2>&1; module load trnascan-se/2.0.12 >/dev/null 2>&1; command -v tRNAscan-SE)"

fail=0
for rep in "${CHUNK[@]}"; do
  [ -z "$rep" ] && continue
  fasta="$FASTA_DIR/${rep}.fa"
  gff_out="$RRNA_OUT/${rep}.gff"
  trna_out="$TRNA_OUT/${rep}.trnascan.txt"

  if [ -s "$gff_out" ] && [ -s "$trna_out" ]; then
    echo "[skip] $rep — outputs already present"
    continue
  fi
  if [ ! -s "$fasta" ]; then
    echo "[ERROR] $rep — missing FASTA $fasta" >&2
    fail=1
    continue
  fi

  domain=$(awk -F'\t' -v r="$rep" '$1==r{print $2}' "$DOMAIN_MAP")
  case "$domain" in
    Bacteria) bkingdom=bac; tflag=-B ;;
    Archaea)  bkingdom=arc; tflag=-A ;;
    *)
      echo "[ERROR] $rep — no/unknown domain ('$domain') in $DOMAIN_MAP" >&2
      fail=1
      continue
      ;;
  esac

  # --- rRNA: barrnap (own subshell: fresh purge + load, isolated from tRNAscan-SE's env) ---
  if [ ! -s "$gff_out" ]; then
    (
      module purge >/dev/null 2>&1
      module load barrnap/0.9 >/dev/null 2>&1
      barrnap --kingdom "$bkingdom" --quiet --threads "${SLURM_CPUS_PER_TASK:-4}" "$fasta"
    ) > "$gff_out.tmp" 2> "$BASE/logs/${rep}.barrnap.log"
    if [ $? -eq 0 ] && [ -s "$gff_out.tmp" ]; then
      mv "$gff_out.tmp" "$gff_out"
    else
      echo "[ERROR] $rep — barrnap failed (domain=$domain)" >&2
      rm -f "$gff_out.tmp"
      fail=1
    fi
  fi

  # --- tRNA: tRNAscan-SE (own subshell: fresh purge + load) ---
  if [ ! -s "$trna_out" ]; then
    (
      module purge >/dev/null 2>&1
      module load trnascan-se/2.0.12 >/dev/null 2>&1
      tRNAscan-SE "$tflag" -Q -q --thread "${SLURM_CPUS_PER_TASK:-4}" \
        -o "$trna_out.tmp" "$fasta"
    ) > "$BASE/logs/${rep}.trnascan.log" 2>&1
    rc=$?
    if [ -s "$trna_out.tmp" ]; then
      # tRNAscan-SE writes a header-only file when zero tRNAs are found — that is a valid result
      mv "$trna_out.tmp" "$trna_out"
    else
      echo "[ERROR] $rep — tRNAscan-SE failed (domain=$domain, rc=$rc)" >&2
      rm -f "$trna_out.tmp"
      fail=1
    fi
  fi

  echo "[done] $rep domain=$domain"
done

echo "Task $SLURM_ARRAY_TASK_ID complete, fail=$fail"
exit $fail
