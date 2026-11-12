#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=eq_loc_search
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis3_localisation/logs/%x-%j.out
# Stage 2a: HMMER (Pfam, dbCAN on EQ member proteins; TIGRFAM/Pfam sorting-signal + sortase subset on all EQ
# proteins) and DIAMOND (MEROPS pepunit). Full-length proteins; mature-region filtering is done in 03_combine.py.
# Temp files live in node RAM (/dev/shm) and are removed at exit; Lustre gets a handful of packed tables.
set -euo pipefail
ROOT=/scratch/project/k10018/eq-protein-external-20261110
OUT=$ROOT/analysis3_localisation
WORK=$OUT/work
PY=/project/k10018/sw/envs/mapping/bin/python
IMG=/project/k10018/sw/images/eq-annot
PFAM=/project/k10018/sw/db/pfam-ibex/pfam/Pfam-A.hmm
TIGR=/project/k10018/sw/db/interproscan-5.61-93.0-ibex/interproscan-5.61-93.0/data/tigrfam/15.0/TIGRFAMs_HMM.LIB
LDB=/project/k10018/sw/db/eq-localisation-20261111
module load singularity >/dev/null 2>&1
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C
sif(){ local n=$1; shift; singularity exec -B /project,/scratch,/dev/shm "$IMG/$n.sif" "$@"; }
T=/dev/shm/eqloc_${SLURM_JOB_ID}; rm -rf "$T"; mkdir -p "$T"; trap 'rm -rf "$T"' EXIT
NCHUNK=90
mkdir -p $WORK

[ -s $WORK/prep_info.json ]
NSEQ=$(grep -c '>' $WORK/eq_members_full.faa)
NALL=$(grep -c '>' $WORK/eq_all_proteins.faa)
echo "member seqs $NSEQ ; all EQ proteins $NALL"

# --- sorting-signal / sortase / pilin / S-layer model subset (TIGRFAM 15.0 + Pfam) -------------------
grep -A2 '^NAME' $TIGR | awk '/^NAME/{n=$2} /^DESC/{print n"\t"$0}' \
  | grep -i -E 'sortase|CTERM|PGF|PEF|sorting|flagellin|pilin|S-layer|S_layer|class III signal|LPXTG|prepilin' \
  | grep -v -E '_Cterm:|_cterm:|_Cterm_|Cterm_rpt|biosynthesis' > $T/tigr_keys.tsv || true
cp $T/tigr_keys.tsv $WORK/tigrfam_subset_keys.tsv
cut -f1 $T/tigr_keys.tsv | sort -u > $T/tigr_keys.txt
grep -E '^NAME' $PFAM | awk '{print $2}' \
  | grep -E '^(Class_IIIsignal|Exosortase_EpsH|Flagellin_arch-type|PEF-CTERM_ARCH|PEP-CTERM|PGF-CTERM|PID_CTERM|IPTL-CTERM|VPDSG-CTERM|VPEID-CTERM|VPxxxP-CTERM_ARCH|Pilin_N|S-layer|S_layer_C|S_layer_N|Sortase|Arch_flagellin|Archaeosortase.*|PGF-pre-PGF|PGF_pre_PGF|DUF.*PGF.*)$' > $T/pfam_keys.txt || true
sif hmmer hmmfetch -f $TIGR $T/tigr_keys.txt > $T/sort.hmm
sif hmmer hmmfetch -f $PFAM $T/pfam_keys.txt >> $T/sort.hmm
grep -E '^(NAME|ACC|DESC|LENG|GA|TC|NC)' $T/sort.hmm > $WORK/sorting_models_info.txt
echo "sorting subset models: $(grep -c '^NAME' $T/sort.hmm)"
( sif hmmer hmmsearch --cpu 16 -Z $NALL --noali -o /dev/null --domtblout $T/sort.domtbl $T/sort.hmm $WORK/eq_all_proteins.faa \
  && grep -v '^#' $T/sort.domtbl | gzip -c > $WORK/sorting_models_vs_eq_all.domtbl.gz.tmp \
  && mv $WORK/sorting_models_vs_eq_all.domtbl.gz.tmp $WORK/sorting_models_vs_eq_all.domtbl.gz && echo SORT_DONE ) &
PSORT=$!

# --- MEROPS pepunit (DIAMOND): see sbatch_02b_merops.sh. In job 15955356 the MEROPS step of this script failed
#     (DIAMOND could not parse the raw pepunit.lib), after the HMMER outputs had been written; the HMMER results
#     of 15955356 are the ones used, and the MEROPS search was rerun on a normalised FASTA as job 15955359.

# --- Pfam + dbCAN on member sequences, chunked ---------------------------------------------------------
$PY - "$WORK/eq_members_full.faa" "$T" $NCHUNK <<'EOF'
import sys
src, t, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
outs = [open(f"{t}/chunk_{i:03d}.faa", "w") for i in range(n)]
i = -1
for line in open(src):
    if line.startswith(">"):
        i += 1
    outs[i % n].write(line)
for o in outs:
    o.close()
EOF
run_chunk(){
  c=$1
  singularity exec -B /project,/scratch,/dev/shm $IMG/hmmer.sif hmmsearch --cpu 1 -Z $NSEQ --noali -o /dev/null \
      --domtblout $T/pfam_$c.domtbl $PFAM $T/chunk_$c.faa
  singularity exec -B /project,/scratch,/dev/shm $IMG/hmmer.sif hmmsearch --cpu 1 -Z $NSEQ --noali -o /dev/null \
      --domtblout $T/dbcan_$c.domtbl $LDB/dbCAN-HMMdb-V14.txt $T/chunk_$c.faa
  echo ok > $T/done_$c
}
export -f run_chunk; export T IMG PFAM LDB NSEQ
seq -f '%03g' 0 $((NCHUNK-1)) | xargs -P $NCHUNK -I{} bash -c 'run_chunk {}'
[ "$(ls $T/done_* | wc -l)" -eq $NCHUNK ]
cat $T/pfam_*.domtbl | grep -v '^#' | gzip -c > $WORK/pfam_vs_members.domtbl.gz.tmp && mv $WORK/pfam_vs_members.domtbl.gz.tmp $WORK/pfam_vs_members.domtbl.gz
cat $T/dbcan_*.domtbl | grep -v '^#' | gzip -c > $WORK/dbcan_vs_members.domtbl.gz.tmp && mv $WORK/dbcan_vs_members.domtbl.gz.tmp $WORK/dbcan_vs_members.domtbl.gz
# per-model cutoffs and descriptions for post-hoc gathering-threshold filtering
awk '/^NAME/{n=$2} /^ACC/{a=$2} /^DESC/{d=$0; sub(/^DESC  /,"",d)} /^LENG/{l=$2} /^GA/{g1=$2; g2=$3; sub(/;/,"",g2); print n"\t"a"\t"l"\t"g1"\t"g2"\t"d}' $PFAM > $WORK/pfam_models.tsv
awk '/^NAME/{n=$2; a="-"; d="-"; g1="NA"; g2="NA"; t1="NA"; t2="NA"} /^ACC/{a=$2} /^DESC/{d=$0; sub(/^DESC  /,"",d)} /^LENG/{l=$2} /^GA/{g1=$2; g2=$3; sub(/;/,"",g2)} /^TC/{t1=$2; t2=$3; sub(/;/,"",t2)} /^HMM /{print n"\t"a"\t"l"\t"g1"\t"g2"\t"t1"\t"t2"\t"d}' $T/sort.hmm > $WORK/sorting_models.tsv
wait $PSORT
ls -la $WORK
sif hmmer hmmsearch -h | sed -n 2p
sif diamond diamond version
echo SEARCH_DONE
