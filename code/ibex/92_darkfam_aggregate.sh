#!/bin/bash
# Phase A.3 — aggregate per-sample CoverM into the dark-family x sample matrices.
# Each catalogue gene IS one dark family rep, so rows = families directly.
#   sbatch --dependency=afterok:<91arrayid> scripts/ibex/92_darkfam_aggregate.sh
#SBATCH --job-name=hohndor_darkfam_agg
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_agg_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_agg_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$WD/06_darkmatter/darkfam_ecology
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eco_probes
COVF_MIN=0.5   # covered_fraction threshold for presence

python3 - "$OUT" "$COVF_MIN" <<'PY'
import sys, os, csv
OUT, COVF_MIN = sys.argv[1], float(sys.argv[2])
man = os.path.join(OUT, "samples_darkfam.tsv")
samples = [r["sample"] for r in csv.DictReader(open(man), delimiter="\t")]

# per-sample: family -> (tpm, covf); keep families detected (covf>0) anywhere
tpm = {}   # family -> {sample: tpm}
prs = {}   # family -> {sample: 1}
done = []
for s in samples:
    cov = os.path.join(OUT, "cov", s, "cov.tsv")
    if not (os.path.exists(cov) and os.path.exists(os.path.join(OUT, "cov", s, "DONE"))):
        continue
    done.append(s)
    with open(cov) as fh:
        next(fh, None)  # header
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            fam = c[0]
            if fam in ("Contig", ""):
                continue
            try:
                t = float(c[1]); cf = float(c[2])
            except ValueError:
                continue
            if t > 0:
                tpm.setdefault(fam, {})[s] = t
            if cf >= COVF_MIN:
                prs.setdefault(fam, {})[s] = 1

fams = sorted(set(tpm) | set(prs))
print(f"[agg] samples processed: {len(done)}/{len(samples)} ; families detected: {len(fams)}", flush=True)

def write_matrix(path, data, fmt):
    with open(path, "w") as o:
        o.write("family_id\t" + "\t".join(done) + "\n")
        for f in fams:
            row = data.get(f, {})
            o.write(f + "\t" + "\t".join(fmt(row.get(s, 0)) for s in done) + "\n")

write_matrix(os.path.join(OUT, "darkfam_abundance.tsv"), tpm, lambda v: f"{v:.4f}")
write_matrix(os.path.join(OUT, "darkfam_presence.tsv"),  prs, lambda v: str(int(v)))

# summary: prevalence (n samples present at covf>=thr) per family + sample richness
prev = {f: len(prs.get(f, {})) for f in fams}
rich = {s: sum(1 for f in fams if prs.get(f, {}).get(s)) for s in done}
with open(os.path.join(OUT, "darkfam_summary.txt"), "w") as o:
    def w(x): o.write(x + "\n"); print(x, flush=True)
    w("=== Phase A: dark-family x sample matrix ===")
    w(f"samples processed: {len(done)} / {len(samples)}")
    w(f"dark families detected (>=1 sample, TPM>0): {len(fams)}")
    pv = sorted(prev.values())
    if pv:
        import statistics as st
        w(f"prevalence (samples present, covf>={COVF_MIN}): median={st.median(pv)} "
          f"max={max(pv)} singletons(1 sample)={sum(1 for v in pv if v==1)}")
        core = sum(1 for v in pv if v >= 0.5*len(done))
        w(f"core families (present in >=50% of samples): {core}")
    if rich:
        rv = sorted(rich.values())
        import statistics as st
        w(f"per-sample dark-family richness: median={st.median(rv)} min={min(rv)} max={max(rv)}")
PY

echo "=== summary ==="; cat "$OUT/darkfam_summary.txt"
ls -lh "$OUT"/darkfam_abundance.tsv "$OUT"/darkfam_presence.tsv
date -Iseconds > "$OUT/AGG_DONE"
echo "[$(date -Iseconds)] aggregate done -> $OUT"
