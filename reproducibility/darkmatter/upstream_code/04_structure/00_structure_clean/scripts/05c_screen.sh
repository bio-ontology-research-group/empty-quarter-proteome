#!/bin/bash
# ⚠ SUPERSEDED 2026-08-18 by 05c_afdb.sh (split AFDB leg) + 05d_integrate.sh — this single job TIMED OUT
# at 24h (job 50603810) mid AFDB split 2/2. Kept for provenance; the integration Python now lives
# once in 05_integrate.py. Do not resubmit this script.
# 05c_screen — full-coordinate FoldSeek TM search of ALL 17,221 admitted-family rep models
# (our 05b folds + reused PI models) PLUS the PI's 3,684 blinded known-fold control models as
# co-queries, vs the frozen current references (PDB 2026-07-30 + AFDB v6). Params verbatim from
# rebuild_20260729/scripts/167_postfold (foldseek/10-941cd33, --alignment-type 1 -e 10 -s 9.5
# --max-seqs 50, both TM normalisations in the output). Controls searched under the IDENTICAL
# invocation give our own measured false-dual-novel floor. Then integrates: gates (length>=100,
# mean pLDDT>=0.70), dual qTM<0.5 call, control floor with Wilson CI, alignment-TM sensitivity.
# Submit after the 05b array completes: sbatch scripts/05c_screen.sh
#SBATCH --job-name=sc5c_screen
#SBATCH --account=pi-hohndor
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=48
#SBATCH --mem=768G
#SBATCH --time=1-00:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5c_screen_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5c_screen_%j.err
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
V=$M/04_structure/00_structure_clean/05_validate
RUN=$M/rebuild_20260729
PDB=$RUN/structural_calibration/references/pdb_20260730
AFDB=$RUN/structural_calibration/references/afdb_uniprot_v6_uniprot_2025_03
CALDIR=$RUN/structural_calibration
MANIFEST=$CALDIR/current_reference_control_manifest_20260802.tsv
PIREPS=$CALDIR/current_structural_novelty_20260730/per_candidate_current_structures.tsv
FORMAT='query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,qcov,tcov,qlen,tlen,alnlen'
test -f "$PDB/REFERENCE_COMPLETE"; test -f "$AFDB/REFERENCE_COMPLETE"
test -f "$MANIFEST"; test -f "$PIREPS"

# --- gate: every 05b shard produced its pLDDT table; fold PDB count == foldable rows ---
NSH=$(ls "$V/shards" | wc -l)
NTSV=$(ls "$V"/models/plddt_shard_*.tsv 2>/dev/null | wc -l)
[ "$NTSV" -eq "$NSH" ] || { echo "ABORT: $NTSV/$NSH pLDDT shard tables (05b incomplete)"; exit 2; }
NOK=$(awk -F'\t' 'FNR>1 && $3!="OOM"' "$V"/models/plddt_shard_*.tsv | wc -l)
NPDB=$(ls "$V/models/fold" | wc -l)
NOOM=$(awk -F'\t' 'FNR>1 && $3=="OOM"' "$V"/models/plddt_shard_*.tsv | wc -l)
echo "folded ok: $NOK  pdbs: $NPDB  OOM-skipped: $NOOM"
[ "$NPDB" -eq "$NOK" ] || { echo "ABORT: pdb count != ok rows"; exit 2; }

# --- assemble query dir: fold + reuse + KF controls ---
Q=$V/query_pdb; rm -rf "$Q"; mkdir -p "$Q"
for f in "$V"/models/fold/*.pdb "$V"/models/reuse/*.pdb; do ln -s "$f" "$Q/$(basename "$f")"; done
NKF=0
while IFS=$'\t' read -r blind rel; do
  [ "$blind" = "blind_id" ] && continue
  ln -s "$CALDIR/$rel" "$Q/$blind.pdb"; NKF=$((NKF+1))
done < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="blind_id")b=i;if($i=="fold_status")s=i;if($i=="model_pdb_relative_path")p=i}} NR>1&&$s=="ok"{print $b"\t"$p}' "$MANIFEST")
NQ=$(ls "$Q" | wc -l)
echo "queries: $NQ (folds $NPDB + reuse $(ls "$V/models/reuse" | wc -l) + controls $NKF)"

# --- search both refs, params verbatim from 167_postfold ---
# Working dir on /ibex/user (1.5T personal quota): the c2014 filesystem is ~99% full and job
# 50603505 died on transient ENOSPC mid-search (2026-08-17). Only compact final outputs go back.
module load foldseek/10-941cd33
foldseek version | tee "$V/foldseek_version.txt"
rm -rf /ibex/user/velazqal/sc5c_screen_*
SW=/ibex/user/velazqal/sc5c_screen_${SLURM_JOB_ID}; mkdir -p "$SW/pdb_tmp" "$SW/afdb_tmp"
foldseek easy-search "$Q" "$PDB/foldseek_pdb" "$SW/pdb_hits.raw.tsv" "$SW/pdb_tmp" \
  --alignment-type 1 -e 10 -s 9.5 --max-seqs 50 \
  --threads "${SLURM_CPUS_PER_TASK:-48}" --format-output "$FORMAT"
foldseek easy-search "$Q" "$AFDB/database/afdb" "$SW/afdb_hits.raw.tsv" "$SW/afdb_tmp" \
  --alignment-type 1 -e 10 -s 9.5 --max-seqs 50 --db-load-mode 2 --split-memory-limit 650G \
  --threads "${SLURM_CPUS_PER_TASK:-48}" --format-output "$FORMAT"
rm -rf "$V/screen"; mkdir -p "$V/screen"
mv "$SW/pdb_hits.raw.tsv" "$SW/afdb_hits.raw.tsv" "$V/screen/"
rm -rf "$SW"

# --- integrate: gates + dual-novel call + control floor + alignment-TM sensitivity ---
python3 - "$V" "$PIREPS" "$MANIFEST" <<'PY'
import sys, json, csv, math, os
V, PIREPS, MANIFEST = sys.argv[1:4]

def norm_plddt(x):
    return x / 100.0 if x > 1.5 else x

# rep metadata: chosen_rep -> (family, n_members, source)
meta = {}
with open(f"{V}/fam_reps.tsv") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        meta[row["chosen_rep"]] = (row["family"], int(row["n_members"]), row["source"])

# pLDDT + length per query
info = {}   # id -> (length, plddt)
for fn in sorted(os.listdir(f"{V}/models")):
    if fn.startswith("plddt_shard_") and fn.endswith(".tsv"):
        with open(f"{V}/models/{fn}") as fh:
            next(fh)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[2] == "OOM": continue
                info[p[0]] = (int(p[1]), norm_plddt(float(p[2])))
with open(PIREPS) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["rep"] in meta and row["rep"] not in info:
            info[row["rep"]] = (int(row["length"]), norm_plddt(float(row["mean_plddt"])))
controls = {}
with open(MANIFEST) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["fold_status"] == "ok":
            controls[row["blind_id"]] = (int(row["fold_length"]), norm_plddt(float(row["fold_mean_plddt"])))

def best_scores(path):
    best = {}   # query -> [qtm, ttm, alntm]
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
            alntm, qtm, ttm = float(p[2]), float(p[3]), float(p[4])
            b = best.setdefault(q, [0.0, 0.0, 0.0])
            b[0] = max(b[0], qtm); b[1] = max(b[1], ttm); b[2] = max(b[2], alntm)
    return best

pdb = best_scores(f"{V}/screen/pdb_hits.raw.tsv")
afdb = best_scores(f"{V}/screen/afdb_hits.raw.tsv")
Z = [0.0, 0.0, 0.0]

def wilson(k, n, z=1.959964):
    if n == 0: return (0.0, 0.0)
    p = k / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d; h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return (c - h, c + h)

def classify(ids, source_info, out_tsv):
    rows, stats = [], dict(searched=0, conf_q=0, dual_qtm=0, conf_dual_qtm=0, dual_alntm=0, conf_dual_alntm=0)
    for q in sorted(ids):
        L, pl = source_info[q]
        pq, aq = pdb.get(q, Z), afdb.get(q, Z)
        conf = (L >= 100 and pl >= 0.70)
        dual_q = (pq[0] < 0.5 and aq[0] < 0.5)
        dual_a = (pq[2] < 0.5 and aq[2] < 0.5)
        stats["searched"] += 1
        stats["conf_q"] += conf; stats["dual_qtm"] += dual_q
        stats["conf_dual_qtm"] += (conf and dual_q)
        stats["dual_alntm"] += dual_a; stats["conf_dual_alntm"] += (conf and dual_a)
        rows.append([q, L, f"{pl:.4f}", f"{pq[0]:.4f}", f"{pq[1]:.4f}", f"{pq[2]:.4f}",
                     f"{aq[0]:.4f}", f"{aq[1]:.4f}", f"{aq[2]:.4f}", int(conf), int(dual_q), int(conf and dual_q)])
    with open(out_tsv, "w") as fh:
        fh.write("id\tlength\tmean_plddt\tbest_pdb_qtm\tbest_pdb_ttm\tbest_pdb_alntm\tbest_afdb_qtm\tbest_afdb_ttm\tbest_afdb_alntm\tconfidence_qualified\tdual_novel_qtm\tconf_qualified_dual_novel\n")
        for r_ in rows: fh.write("\t".join(map(str, r_)) + "\n")
    return stats

rep_ids = [q for q in meta if q in info]
missing = [q for q in meta if q not in info]
rs = classify(rep_ids, info, f"{V}/per_rep_structures.tsv")
cs = classify(controls, controls, f"{V}/per_control_structures.tsv")
flo, fhi = wilson(cs["conf_dual_qtm"], cs["conf_q"])
summary = {
    "admitted_families": len(meta),
    "reps_searched": rs["searched"], "reps_missing_model": len(missing),
    "reps_confidence_qualified_Q": rs["conf_q"],
    "reps_dual_novel_qtm_ungated": rs["dual_qtm"],
    "reps_conf_qualified_dual_novel_K": rs["conf_dual_qtm"],
    "K_over_Q_pct": round(100.0 * rs["conf_dual_qtm"] / rs["conf_q"], 2) if rs["conf_q"] else None,
    "sensitivity_alntm_conf_qualified_dual_novel": rs["conf_dual_alntm"],
    "controls_searched": cs["searched"], "controls_confidence_qualified": cs["conf_q"],
    "control_false_dual_novel_conf_qualified": cs["conf_dual_qtm"],
    "control_floor_pct": round(100.0 * cs["conf_dual_qtm"] / cs["conf_q"], 2) if cs["conf_q"] else None,
    "control_floor_wilson95_pct": [round(100*flo, 2), round(100*fhi, 2)],
    "control_floor_alntm": cs["conf_dual_alntm"],
    "gates": {"length_min": 100, "mean_plddt_min": 0.70, "tm_threshold": 0.5, "primary_norm": "qtm"},
}
json.dump(summary, open(f"{V}/summary_05c.json", "w"), indent=1)
with open(f"{V}/SUMMARY_05c.txt", "w") as fh:
    for k, v in summary.items(): fh.write(f"{k}\t{v}\n")
print(json.dumps(summary, indent=1))
if missing: print(f"WARNING: {len(missing)} reps had no model/plddt (first: {missing[:3]})")
PY
echo "DONE_05c"
