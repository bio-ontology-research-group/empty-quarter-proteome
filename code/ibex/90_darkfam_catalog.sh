#!/bin/bash
# Phase A.1 (dark-matter ecology) — build the nucleotide gene catalogue for read recruitment.
# Unit = orthogroups whose REPRESENTATIVE is sequence-dark (no UniRef50 / <30%) AND that have >=3
# members (multi-member dark families). For each such rep we extract its pyrodigal CDS nucleotide
# sequence from the dereplicated genome FASTA (the .faa/.gff are amino-acid; Bakta .ffn IDs don't
# match the pyrodigal orthogroup IDs, so we slice the contig directly). Also writes the 551-sample
# read manifest used by the CoverM array.
#   sbatch scripts/ibex/90_darkfam_catalog.sh
#SBATCH --job-name=hohndor_darkfam_catalog
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_catalog_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_catalog_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
DM=$WD/06_darkmatter
OG=$DM/orthogroups
GEN=$WD/05_refine/_drep/out/dereplicated_genomes
DRP=$DM/drep_rep_proteomes
DARK=$DM/reallydark/isdark_ids.txt
READS=/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads
OUT=$DM/darkfam_ecology
CAT=$OUT/catalog
mkdir -p "$CAT"
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eco_probes
MIN_MEMBERS=3

echo "[$(date -Iseconds)] building dark-family nucleotide catalogue"
python3 - "$OG" "$GEN" "$DRP" "$DARK" "$CAT" "$MIN_MEMBERS" <<'PY'
import sys, os, gzip
OG, GEN, DRP, DARK, CAT, MINM = sys.argv[1:7]
MINM = int(MINM)

def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

# 1) dark protein ids (rep ids that are sequence-novel)
dark = set(l.strip() for l in open(DARK) if l.strip())
print(f"[catalog] dark proteins: {len(dark)}", flush=True)

# 2) og_recurrence.tsv : rep \t n_members \t n_samples  -> keep dark reps with >= MINM members
wanted = {}   # rep_full -> (n_members, n_samples)
rec = os.path.join(OG, "og_recurrence.tsv")
for line in open(rec):
    c = line.rstrip("\n").split("\t")
    if len(c) < 3:
        continue
    rep = c[0]
    if rep not in dark:
        continue
    try:
        nm, ns = int(c[1]), int(c[2])
    except ValueError:
        continue
    if nm >= MINM:
        wanted[rep] = (nm, ns)
print(f"[catalog] dark multi-member (>= {MINM}) families: {len(wanted)}", flush=True)

# group wanted reps by mag; locus_id = part after '::'
by_mag = {}
for rep in wanted:
    mag, locus = rep.split("::", 1)
    by_mag.setdefault(mag, {})[locus] = rep

# meta table
mt = open(os.path.join(CAT, "darkfam_meta.tsv"), "w")
mt.write("family_id\tmag\tlocus\tn_members\tn_samples_catalog\n")
for rep, (nm, ns) in wanted.items():
    mag, locus = rep.split("::", 1)
    mt.write(f"{rep}\t{mag}\t{locus}\t{nm}\t{ns}\n")
mt.close()

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s): return s.translate(COMP)[::-1]

def read_fasta(path):
    name, seq = None, []
    for line in opn(path):
        if line.startswith(">"):
            if name is not None: yield name, "".join(seq)
            name = line[1:].split()[0]; seq = []
        else:
            seq.append(line.strip())
    if name is not None: yield name, "".join(seq)

out = open(os.path.join(CAT, "darkfam_catalog.fna"), "w")
n_fam = n_done = n_missing_mag = n_missing_cds = 0
total_bp = 0
for mag, loci in by_mag.items():
    n_fam += len(loci)
    gpath = os.path.join(GEN, mag + ".fa")
    gff = os.path.join(DRP, mag + ".gff.gz")
    if not (os.path.exists(gpath) and os.path.exists(gff)):
        n_missing_mag += len(loci); continue
    contigs = {n: s for n, s in read_fasta(gpath)}
    # parse CDS coords for the wanted loci of this mag
    coords = {}  # locus -> (contig,start,end,strand)
    for line in opn(gff):
        if line.startswith("#") or "\tCDS\t" not in line:
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9: continue
        attrs = c[8]
        i = attrs.find("ID=")
        if i < 0: continue
        locus = attrs[i+3:].split(";")[0]
        if locus in loci:
            coords[locus] = (c[0], int(c[3]), int(c[4]), c[6])
    for locus, rep in loci.items():
        cc = coords.get(locus)
        if not cc:
            n_missing_cds += 1; continue
        contig, s, e, strand = cc
        seq = contigs.get(contig)
        if not seq:
            n_missing_cds += 1; continue
        nt = seq[s-1:e]
        if strand == "-":
            nt = revcomp(nt)
        if not nt:
            n_missing_cds += 1; continue
        out.write(f">{rep}\n{nt}\n")
        n_done += 1; total_bp += len(nt)
out.close()
print(f"[catalog] families={n_fam} extracted={n_done} missing_mag={n_missing_mag} "
      f"missing_cds={n_missing_cds} total_bp={total_bp}", flush=True)
open(os.path.join(CAT, "CATALOG_STATS.txt"), "w").write(
    f"families={n_fam}\nextracted={n_done}\nmissing_mag={n_missing_mag}\n"
    f"missing_cds={n_missing_cds}\ntotal_bp={total_bp}\nmin_members={MINM}\n")
PY

echo "[$(date -Iseconds)] building 551-sample read manifest"
python3 - "$READS" "$OUT" <<'PY'
import sys, os, glob, re
READS, OUT = sys.argv[1], sys.argv[2]
rows = []
for r1 in sorted(glob.glob(os.path.join(READS, "*_R1.fastq.gz"))):
    s = os.path.basename(r1)[:-len("_R1.fastq.gz")]
    # exclude controls / non-soil per CLAUDE.md naming convention
    if s.startswith("C-") or s in ("qc","SNTC","SPositive") or s.startswith("P34") \
       or s.startswith("SNTC") or s.startswith("SPositive"):
        continue
    r2 = os.path.join(READS, f"{s}_R2.fastq.gz")
    if not os.path.exists(r2):
        continue
    se = os.path.join(READS, f"{s}_SE.fastq.gz")
    rows.append((s, r1, r2, se if os.path.exists(se) else ""))
with open(os.path.join(OUT, "samples_darkfam.tsv"), "w") as o:
    o.write("sample\tR1\tR2\tSE\n")
    for r in rows:
        o.write("\t".join(r) + "\n")
print(f"[manifest] samples: {len(rows)}", flush=True)
PY

echo "=== catalog stats ==="; cat "$CAT/CATALOG_STATS.txt"
echo "=== n samples ==="; tail -n +2 "$OUT/samples_darkfam.tsv" | wc -l
ls -lh "$CAT/darkfam_catalog.fna"
echo "[$(date -Iseconds)] catalogue + manifest ready -> $OUT"
date -Iseconds > "$OUT/CATALOG_DONE"
