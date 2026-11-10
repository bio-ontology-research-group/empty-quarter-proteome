#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --job-name=eq_sci_stg2
#SBATCH --output=/scratch/project/k10018/eq-science-20261108/logs/%x_%j.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261108/logs/%x_%j.err

set -euo pipefail
export EQ_SCIENCE_OUT=/scratch/project/k10018/eq-science-20261108
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping
module load singularity >/dev/null 2>&1 || true
export SINGULARITY_CACHEDIR=/scratch/project/k10018/cache/singularity
export SINGULARITYENV_LC_ALL=C SINGULARITYENV_LANG=C

T=/dev/shm/eq_stage2_${SLURM_JOB_ID}
MM_TMP="$EQ_SCIENCE_OUT/results/mmseqs_tmp_${SLURM_JOB_ID}"
mkdir -p "$T" "$MM_TMP"
trap 'rm -rf "$T" "$MM_TMP"' EXIT

M8_OUT="$EQ_SCIENCE_OUT/results/thermo538_vs_uniref50_s75_e1e3.m8.gz"

if [ ! -s "$M8_OUT" ]; then
  echo "[Stage 2] Extracting and prefixing 538 Thermoplasmatota proteomes to $T/thermo538.faa at $(date -Is)..."
  python3 - "$T/thermo538.faa" << 'EOF'
import sys, tarfile
from pathlib import Path
out_faa = sys.argv[1]
tar_p = "/project/k10018/work/eq-annot-20261030/results/c7_thermo_evol/thermo_proteomes.tar.gz"
n_seqs = 0
with tarfile.open(tar_p, "r:gz") as tf, open(out_faa, "w") as out:
    for m in tf:
        if not m.isfile() or not m.name.endswith(".faa"):
            continue
        gid = Path(m.name).stem
        f = tf.extractfile(m)
        if f is None:
            continue
        for raw in f:
            line = raw.decode("utf-8", errors="replace")
            if line.startswith(">"):
                tok = line[1:].strip().split()[0]
                out.write(f">{gid}|{tok}\n")
                n_seqs += 1
            else:
                out.write(line)
print(f"Wrote {n_seqs} prefixed sequences to {out_faa}", flush=True)
EOF

  echo "[Stage 2] Running MMseqs2 easy-search (-s 7.5 -e 1e-3 --max-accept 10, no identity/cov gate) vs UniRef50 2026_01 at $(date -Is)..."
  singularity exec -B /project,/scratch,/dev/shm /project/k10018/sw/images/eq-annot/mmseqs2.sif \
    mmseqs easy-search "$T/thermo538.faa" /project/k10018/sw/db/eq-c2/uniref50/uniref50 \
    "$T/thermo538_uniref50.m8" "$MM_TMP" \
    -s 7.5 -e 1e-3 --max-accept 10 --threads "${SLURM_CPUS_ON_NODE:-192}" --split-memory-limit 250G \
    --format-output query,target,pident,alnlen,evalue,bits >/dev/null
  gzip -c "$T/thermo538_uniref50.m8" > "${M8_OUT}.tmp" && mv "${M8_OUT}.tmp" "$M8_OUT"
  rm -rf "$MM_TMP"
fi

echo "[Stage 2] Running stage2_uniref50_harmonised.py at $(date -Is)..."
python3 "$EQ_SCIENCE_OUT/code/stage2_uniref50_harmonised.py" "$M8_OUT"
echo "[Stage 2] Done at $(date -Is)."
