#!/bin/bash
# §3.1 rarefaction control — the leg that licenses any "content beats depth" statement.
#
# WHY: DIAMOND E-values scale with database size, so a fixed --evalue 1e-3 is more permissive against a
# small catalogue than a large one — which flatters exactly the small biome-matched desert catalogues.
# Measured on our own data: the five desert DBs searched SEPARATELY find 63,413 dark families, the same
# sequences concatenated into one 147.5M DB find 60,869 — 2,544 families (2.40% of the universe) hit
# only because the databases were split. So the raw cross-catalogue ordering is DESCRIPTIVE only, and
# the claim has to ride on this control.
#
# HOW: one search against full MGnify90 at --evalue 1 -k 500 --max-hsps 1, keeping every E_full, then
# ANALYTIC subsampling in raref_curve.py: a hit survives a p-subsample iff E_full <= 1e-3/p, and
# P(no hit) = (1-p)^h(p). NO SUBSAMPLED DATABASES ARE BUILT — that was the design trap, since building
# them would reintroduce the very size-dependent threshold we are trying to correct for.
# Both query sets go through (dark + matched annotated control), so the control curve shows what the
# size effect looks like for proteins that plainly have homologs.
#
# Reading it: at p = n_target/717738164 the curve gives the size-ONLY expectation for that depth.
# The gap between a real catalogue's observed absence and the curve at its matched p is the part that
# size cannot explain.
#   sbatch spec_4_rarefaction.sbatch
#SBATCH --job-name=ale_spec_raref
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_raref_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_raref_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=48:00:00
#SBATCH --partition=largemem
#SBATCH --account=pi-hohndor
set -uo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
S=$WD/06_darkmatter/specificity_clean
B=$WD/06_darkmatter/darkfam_basis_clean/results
OUT=$S/results
QUERY=$OUT/query_all.faa
DB=$WD/databases/mgnify90/mgy_clusters
WORK=/ibex/user/velazqal/eqmags_31/raref     # ~3 GB hit table: own quota, not the 99%-full volume
M8=$WORK/raref_mgnify90_e1_k500.m8
mkdir -p "$WORK" "$S/logs"
source $WD/../miniconda3/etc/profile.d/conda.sh; conda activate eggnog
T=${SLURM_CPUS_PER_TASK:-32}
log(){ echo "[$(date -Iseconds)] $*"; }
log "diamond: $(diamond --version 2>&1 | head -1)"

A=$(df -B1G --output=avail "$WORK" | tail -1 | tr -d ' ')
log "free on $WORK: ${A} GB (hit table needs ~5 GB)"
[ "$A" -lt 50 ] && { log "ABORT: <50 GB free"; exit 1; }

# ---------------- one permissive search, every E_full retained ----------------
if [ -s "$M8" ]; then
  log "hit table exists ($(du -h "$M8" | cut -f1)) — skipping search"
else
  log "diamond blastp vs mgnify90 (--sensitive --evalue 1 -k 500 --max-hsps 1)"
  # --max-hsps 1 matters: h(p) counts distinct SUBJECT SEQUENCES, because subjects are what get
  # retained with probability p. Multiple HSPs per subject would double-count.
  /usr/bin/time -v diamond blastp -q "$QUERY" -d "$DB" --sensitive --evalue 1 -k 500 --max-hsps 1 \
     -p "$T" -f 6 qseqid evalue -o "$M8" 2> "$WORK/raref_search.log" \
     || { log "search FAILED — see $WORK/raref_search.log"; exit 1; }
fi
log "hit rows: $(wc -l < "$M8")"

# ---------------- analytic e-value-corrected rarefaction ----------------
CURVE=$OUT/rarefaction_curve.tsv
python3 "$S/scripts/raref_curve.py" "$M8" "$B/FAMILIES_FINAL.txt" "$B/CONTROL_FAMILIES.txt" "$CURVE" \
  || { log "curve computation FAILED"; exit 1; }
log "curve -> $CURVE"
column -t "$CURVE" 2>/dev/null || cat "$CURVE"
log "DONE."
