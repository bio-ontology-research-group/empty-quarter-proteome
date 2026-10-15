#!/bin/bash
# 05b_fold — ESMFold the to-fold family reps (one shard per array task) on IBEX A100.
# Adapted from the proven scripts/95_esmfold_fold.py + pilot_af2/03b (per-protein
# checkpointing: skip if the PDB already exists, so requeues are cheap). Env + weights
# provisioned by scripts/94_esmfold_env.sh (conda env `esmfold`, HF cache on scratch).
# Submit after 05a with the shard count it prints:  sbatch --array=1-<NSH> scripts/05b_fold.sh
#SBATCH --job-name=sc5b_fold
#SBATCH --account=pi-hohndor
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5b_fold_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5b_fold_%A_%a.err
#SBATCH --partition=batch
#SBATCH --gres=gpu:a100:1
#SBATCH --exclude=gpu201-23-r
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
V=$M/04_structure/00_structure_clean/05_validate
SH=$(printf "%s/shards/shard_%04d.faa" "$V" $((SLURM_ARRAY_TASK_ID-1)))
[ -s "$SH" ] || { echo "missing shard $SH"; exit 2; }
O=$V/models/fold; mkdir -p "$O"
TSV=$V/models/plddt_$(basename "$SH" .faa).tsv

export HF_HOME=$M/06_darkmatter/esmfold/hf
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate esmfold

python3 - "$SH" "$O" "$TSV" <<'PY'
import sys, os, torch
from transformers import AutoTokenizer, EsmForProteinFolding
FAA, PDB_DIR, PLDDT_TSV = sys.argv[1:4]

def read_faa(p):
    name, seq = None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name: yield name, "".join(seq)
                name = line[1:].strip().split()[0]; seq = []
            else:
                seq.append(line.strip())
    if name: yield name, "".join(seq)

todo = [(n, s.replace("*", "")[:1022]) for n, s in read_faa(FAA)]
todo = [(n, s) for n, s in todo
        if s and not os.path.exists(os.path.join(PDB_DIR, n + ".pdb"))]
print(f"to fold this run: {len(todo)}", flush=True)
if todo:
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.cuda().eval()
    model.esm = model.esm.half()
    torch.backends.cuda.matmul.allow_tf32 = True
    model.trunk.set_chunk_size(64)

mode = "a" if os.path.exists(PLDDT_TSV) else "w"
with open(PLDDT_TSV, mode) as out:
    if mode == "w":
        out.write("protein\tlength\tmean_plddt\n")
    for i, (name, s) in enumerate(todo):
        try:
            with torch.no_grad():
                ids = tok([s], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
                o = model(ids)
            pdb = model.output_to_pdb(o)[0]
            plddt = float(o["plddt"][0, :, 1].mean().item())
            with open(os.path.join(PDB_DIR, name + ".pdb"), "w") as fh:
                fh.write(pdb)
            out.write(f"{name}\t{len(s)}\t{plddt:.4f}\n"); out.flush()
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            out.write(f"{name}\t{len(s)}\tOOM\n"); out.flush()
        if (i + 1) % 100 == 0:
            print(f"{i+1}/{len(todo)}", flush=True)
PY
echo "DONE shard $(basename "$SH")  pdbs=$(grep -c '^>' "$SH") -> $(ls "$O" | wc -l) total in fold/"
