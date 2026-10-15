#!/bin/bash
# Phase C.3 — aggregate per-family dN/dS and run the selection tests:
#  (1) are DARK families under purifying selection (dN/dS << 1) -> real, functional genes?
#  (2) DARK vs size-matched NON-DARK control (the null): Mann-Whitney on per-family median dN/dS.
#  (3) if a Phase-B per-family label file is present (structured/endemic), test structured vs cosmopolitan.
#   sbatch --dependency=afterok:<94array> scripts/ibex/96_darkfam_dnds_aggregate.sh
#SBATCH --job-name=hohndor_darkfam_dnds_agg
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_dnds_agg_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_dnds_agg_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$WD/06_darkmatter/darkfam_ecology/selection
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eco_probes

# concatenate parts -> dnds.tsv
head -1 "$(ls "$OUT"/dnds_parts/part_*.tsv | head -1)" > "$OUT/dnds.tsv"
for f in "$OUT"/dnds_parts/part_*.tsv; do tail -n +2 "$f"; done >> "$OUT/dnds.tsv"

python3 - "$OUT" <<'PY'
import os, sys, csv, statistics as st
from scipy.stats import mannwhitneyu
OUT = sys.argv[1]
rows = list(csv.DictReader(open(os.path.join(OUT, "dnds.tsv")), delimiter="\t"))
def vals(rs): return [float(r["median_dNdS"]) for r in rs if r["median_dNdS"]]
dark = [r for r in rows if r["set"] == "dark"]
ctrl = [r for r in rows if r["set"] == "control"]
vd, vc = vals(dark), vals(ctrl)

out = []
def w(x): out.append(x); print(x)
w("# Phase C — selection on dark protein families (NG86 dN/dS)")
w(f"families with dN/dS: dark={len(vd)}  control(size-matched)={len(vc)}\n")
def desc(v, name):
    if not v: w(f"{name}: none"); return
    frac1 = sum(1 for x in v if x < 1)/len(v)
    w(f"{name}: median dN/dS={st.median(v):.3f}  mean={st.mean(v):.3f}  "
      f"frac<1 (purifying)={frac1:.3f}  n={len(v)}")
desc(vd, "DARK families")
desc(vc, "CONTROL (non-dark, size-matched)")
if vd and vc:
    U, p = mannwhitneyu(vd, vc, alternative="two-sided")
    w(f"\nDARK vs CONTROL median dN/dS: Mann-Whitney p={p:.2e} "
      f"(dark {'<' if st.median(vd)<st.median(vc) else '>='} control)")
    w("Interpretation: dN/dS << 1 in dark families => under purifying selection => real, functional")
    w("genes (not spurious ORFs). If dark <= control, dark matter is at least as constrained as typical genes.")

# (3) optional structured/endemic labels from Phase B
lab = os.path.join(OUT, "..", "phase_b_family_labels.tsv")
if os.path.exists(lab):
    L = {r["family_id"]: r for r in csv.DictReader(open(lab), delimiter="\t")}
    dd = {r["family_id"]: float(r["median_dNdS"]) for r in dark if r["median_dNdS"]}
    structured = [dd[f] for f in dd if L.get(f, {}).get("structured") == "1"]
    cosmo = [dd[f] for f in dd if L.get(f, {}).get("structured") == "0"]
    if structured and cosmo:
        U, p = mannwhitneyu(structured, cosmo, alternative="two-sided")
        w(f"\nSTRUCTURED vs COSMOPOLITAN dark families: median dN/dS "
          f"{st.median(structured):.3f} vs {st.median(cosmo):.3f}  MWU p={p:.2e}  "
          f"(n={len(structured)}/{len(cosmo)})")
        w("=> if structured families are MORE constrained, the env-responsive dark matter is ADAPTIVE.")
else:
    w("\n[note] Phase-B per-family labels (phase_b_family_labels.tsv) not present yet -> "
      "structured-vs-cosmopolitan test deferred.")
open(os.path.join(OUT, "SELECTION_SUMMARY.md"), "w").write("\n".join(out) + "\n")
PY
echo "[$(date -Iseconds)] aggregate done -> $OUT/SELECTION_SUMMARY.md"
date -Iseconds > "$OUT/DNDS_DONE"
