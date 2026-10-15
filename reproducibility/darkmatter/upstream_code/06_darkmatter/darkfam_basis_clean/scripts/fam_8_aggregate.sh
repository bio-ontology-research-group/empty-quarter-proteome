#!/bin/bash
# Aggregate fam_4 per-sample CoverM (member level) -> family x sample presence/abundance matrices
# + leg-3 headline stats. Rerunnable in minutes when late samples land.
#SBATCH --job-name=ale_fam_agg
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam8_%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
python3 - <<'PY'
import os, glob
B="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean"
R=f"{B}/results"
mfam = {}
for line in open(f"{R}/member_fam_final.tsv"):
    m, f = line.rstrip("\n").split("\t"); mfam[m] = f
fams = sorted(set(mfam.values())); fidx = {f: i for i, f in enumerate(fams)}
files = sorted(glob.glob(f"{R}/coverm/*.tsv"))
print(f"[agg] {len(fams)} families, {len(files)} samples", flush=True)
import array
tpm_rows, pres_rows, samples = [], [], []
for fp in files:
    s = os.path.basename(fp)[:-4]
    tpm = array.array('d', [0.0]) * 0
    tpm = [0.0]*len(fams); pres = [0]*len(fams)
    with open(fp) as fh:
        head = fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            f = mfam.get(c[0])
            if f is None: continue
            i = fidx[f]
            t = float(c[1]); cf = float(c[2])
            tpm[i] += t
            if cf >= 0.5: pres[i] = 1
    samples.append(s); tpm_rows.append(tpm); pres_rows.append(pres)
    print(f"[agg] {s} done", flush=True)
with open(f"{R}/fam_abundance.tsv.part","w") as ot, open(f"{R}/fam_presence.tsv.part","w") as op:
    ot.write("family_rep\t" + "\t".join(samples) + "\n")
    op.write("family_rep\t" + "\t".join(samples) + "\n")
    for f, i in fidx.items():
        ot.write(f + "\t" + "\t".join(f"{tpm_rows[k][i]:.4f}" for k in range(len(samples))) + "\n")
        op.write(f + "\t" + "\t".join(str(pres_rows[k][i]) for k in range(len(samples))) + "\n")
os.replace(f"{R}/fam_abundance.tsv.part", f"{R}/fam_abundance.tsv")
os.replace(f"{R}/fam_presence.tsv.part", f"{R}/fam_presence.tsv")
import statistics as st
prev = [sum(pres_rows[k][i] for k in range(len(samples))) for i in range(len(fams))]
peak = [max(tpm_rows[k][i] for k in range(len(samples))) for i in range(len(fams))]
single = sum(1 for p in prev if p <= 1)
print(f"[agg] LEG3: median prevalence {st.median(prev)}/{len(samples)} samples; single-sample {single} ({100*single/len(fams):.1f}%); median peak TPM {st.median(peak):.1f}")
PY
