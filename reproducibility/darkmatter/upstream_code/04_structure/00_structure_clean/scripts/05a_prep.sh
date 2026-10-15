#!/bin/bash
# 05a_prep — pick one fold-representative per ADMITTED family (>=3 members, from 04),
# preferring members that already have a PI ESMFold model (06_esmfold_novelfold/pdb,
# filename = <protein_id>.pdb) so we only GPU-fold what does not exist. Extract AA
# sequences for the to-fold reps from the merged candidate DB (main data = AA) and
# shard them for the 05b GPU array. Assembles the reuse-model symlink dir.
# Submit after 04: sbatch scripts/05a_prep.sh
#SBATCH --job-name=sc5a_prep
#SBATCH --account=pi-hohndor
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5a_prep_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5a_prep_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
SC=$M/04_structure/00_structure_clean; W=$SC/04_cluster; V=$SC/05_validate
FS=$M/tools/foldseek/bin/foldseek
PIMODELS=$M/06_esmfold_novelfold/pdb
PER_SHARD=500   # small shards + 4h walltime = backfillable at low fairshare priority (2026-08-16)

[ -s "$W/families_ge3_reps.txt" ] || { echo "ABORT: run 04 first"; exit 2; }
rm -rf "$V"; mkdir -p "$V/models/reuse" "$V/shards"

# name<TAB>header map for the merged DB (no lookup after concatdbs; headers carry names)
"$FS" prefixid "$W/merged_h" "$V/merged_h.tsv" --tsv >/dev/null

ls "$PIMODELS" | sed 's/\.pdb$//' > "$V/pi_model_ids.txt"

python3 - "$W" "$V" <<'PY'
import sys, os
W, V = sys.argv[1], sys.argv[2]
admitted = set(l.strip() for l in open(f"{W}/families_ge3_reps.txt"))
pi = set(l.strip() for l in open(f"{V}/pi_model_ids.txt"))
# clusters.tsv: rep_name \t member_name
fams = {}
for line in open(f"{W}/clusters.tsv"):
    rep, mem = line.rstrip("\n").split("\t")[:2]
    if rep in admitted:
        fams.setdefault(rep, []).append(mem)
# key map from merged headers: key \t name(first token of header)
name2key = {}
for line in open(f"{V}/merged_h.tsv"):
    p = line.rstrip("\n").split("\t")
    name2key[p[1].split()[0] if len(p) > 1 else ""] = p[0]
out = open(f"{V}/fam_reps.tsv", "w")
out.write("family\tchosen_rep\tn_members\tsource\n")
foldkeys, reuse = [], []
miss_key = 0
for rep, members in sorted(fams.items()):
    withmodel = sorted(m for m in members if m in pi)
    if withmodel:
        chosen, src = withmodel[0], "REUSE"
        reuse.append(chosen)
    else:
        chosen, src = rep, "FOLD"
        k = name2key.get(chosen)
        if k is None:
            miss_key += 1
            continue
        foldkeys.append(k)
    out.write(f"{rep}\t{chosen}\t{len(members)}\t{src}\n")
out.close()
open(f"{V}/fold_rep_keys.txt", "w").write("".join(k + "\n" for k in foldkeys))
open(f"{V}/reuse_rep_ids.txt", "w").write("".join(r + "\n" for r in reuse))
print(f"admitted_families\t{len(fams)}")
print(f"reuse_existing_model\t{len(reuse)}")
print(f"to_fold\t{len(foldkeys)}")
print(f"missing_key_ERROR\t{miss_key}")
assert miss_key == 0, "chosen FOLD reps missing from merged DB"
assert len(reuse) + len(foldkeys) == len(fams)
PY

# reuse symlinks
while read -r id; do ln -s "$PIMODELS/$id.pdb" "$V/models/reuse/$id.pdb"; done < "$V/reuse_rep_ids.txt"
NRE=$(ls "$V/models/reuse" | wc -l); echo "reuse links: $NRE"
[ "$NRE" -eq "$(wc -l < "$V/reuse_rep_ids.txt")" ] || { echo "ABORT: reuse link mismatch"; exit 2; }

# AA fasta for to-fold reps, then shard by record count
"$FS" createsubdb "$V/fold_rep_keys.txt" "$W/merged"   "$V/repDB"   >/dev/null
"$FS" createsubdb "$V/fold_rep_keys.txt" "$W/merged_h" "$V/repDB_h" >/dev/null
"$FS" convert2fasta "$V/repDB" "$V/fold_targets.faa" >/dev/null
NFA=$(grep -c '^>' "$V/fold_targets.faa")
[ "$NFA" -eq "$(wc -l < "$V/fold_rep_keys.txt")" ] || { echo "ABORT: fasta count != keys ($NFA)"; exit 2; }
awk -v per=$PER_SHARD -v dir="$V/shards" '/^>/{if(n%per==0){f=sprintf("%s/shard_%04d.faa",dir,n/per)} n++} {print > f}' "$V/fold_targets.faa"
NSH=$(ls "$V/shards" | wc -l)
{
  echo -e "admitted_families\t$(($(wc -l < "$V/fam_reps.tsv")-1))"
  echo -e "reuse_existing_model\t$NRE"
  echo -e "to_fold\t$NFA"
  echo -e "shards_of_${PER_SHARD}\t$NSH"
} | tee "$V/SUMMARY_05a.txt"
echo "next: sbatch --array=1-$NSH scripts/05b_fold.sh"
echo "DONE_05a"
