#!/bin/bash
# fam_6b — per-family FNA split + family_index for the dN/dS worker (mirrors June's layout,
# repointed into darkfam_basis_clean/results/selection per review fix). ARM=dark|ctrl.
#   sbatch --export=ARM=dark scripts/fam_6b_prep.sh
#
#SBATCH --job-name=ale_dnds_prep
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam6b_%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --partition=batch
set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
R=$B/results
ARM=${ARM:?set ARM=dark|ctrl}
if [ "$ARM" = dark ]; then FNA=$R/darkfam_members.fna; else FNA=$R/ctrl_members.fna; fi
[ -s "$FNA" ] || { echo "ERROR: $FNA missing"; exit 1; }

python3 - "$R" "$ARM" "$FNA" <<'PY'
import sys, os, hashlib
R, ARM, FNA = sys.argv[1:4]
base = f"{R}/selection/members/{ARM}"
os.makedirs(base, exist_ok=True)
for i in range(256): os.makedirs(f"{base}/{i:02x}", exist_ok=True)
fams = {}
cur_fam, cur_fh = None, None
n = 0
with open(FNA) as fh:
    for line in fh:
        if line.startswith(">"):
            fam = line.rsplit("fam=", 1)[1].strip()
            if fam != cur_fam:
                if cur_fh: cur_fh.close()
                shard = hashlib.md5(fam.encode()).hexdigest()[:2]
                path = f"{base}/{shard}/{fam}.fna"
                cur_fh = open(path, "a"); cur_fam = fam
                fams.setdefault(fam, [0, path])
            fams[fam][0] += 1
            cur_fh.write(line.split(" ")[0] + "\n"); n += 1
        else:
            cur_fh.write(line)
if cur_fh: cur_fh.close()
idx = f"{R}/selection/family_index_{ARM}.tsv"
with open(idx + ".part", "w") as out:
    out.write("family_id\tset\tx\tn_members\tfna\n")
    for fam, (cnt, path) in sorted(fams.items()):
        out.write(f"{fam}\t{ARM}\t-\t{cnt}\t{path}\n")
os.replace(idx + ".part", idx)
print(f"[{ARM}] members {n}, families {len(fams)} -> {idx}", flush=True)
PY
echo "[$(date -Iseconds)] fam_6b $ARM done"
