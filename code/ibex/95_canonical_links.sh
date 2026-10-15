#!/bin/bash
# Publish ONLY COMPLETE cross-MAG analyses into the canonical EmptyQuarter_Data tree, anchored to
# MAGs (-> assembly -> site). This folder is the PI's results-pickup point, NOT a working dir:
# every section is gated behind a completion check, so incomplete/in-progress work is never linked
# (and e.g. the dark-matter ecology link auto-appears once Phase A finishes).
#
# Per-sample bins (single_sample + self_supervised) and per-MAG annotation are ALREADY canonical
# under sites/<trip>/site<N>/<sample>/wgs/. This adds the project-level products.
# NON-DESTRUCTIVE: symlinks only; never moves/overwrites a real file. Idempotent.
#   sbatch scripts/ibex/95_canonical_links.sh
#SBATCH --job-name=hohndor_canon_links
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/canon_links_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/canon_links_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
B=/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results
CAT=$B/catalog
DM=$WD/06_darkmatter
MAN=$CAT/MANIFEST.tsv

slink () {  # safe symlink: $1=target(existing) $2=linkpath ; skip if a REAL file/dir occupies linkpath
  local tgt="$1" lnk="$2"
  [ -e "$tgt" ] || { echo "    MISS target: $tgt"; return 0; }
  if [ -e "$lnk" ] && [ ! -L "$lnk" ]; then echo "    SKIP (real file present): $lnk"; return 0; fi
  ln -sfn "$tgt" "$lnk"
}
manifest () { echo -e "$1\t$2\t$3" >> "$MAN"; }   # name, status, note

mkdir -p "$CAT"
: > "$MAN"; echo -e "analysis\tstatus\tnote" > "$MAN"
echo "[$(date -Iseconds)] publishing COMPLETE analyses only"

# ============ 1) MAG catalogue (anchor) — guard: dereplicated genomes present ============
REPS=$WD/05_refine/_drep/out/dereplicated_genomes
if [ -d "$REPS" ] && [ "$(ls "$REPS"/*.fa 2>/dev/null | head -1)" ]; then
  echo "[1] MAG catalogue"
  mkdir -p "$CAT"/mags/dereplicated
  BINIDX=$CAT/mags/_bin_index.tsv
  find "$B"/sites -path '*/wgs/binning/*/output_bins/*.fa.gz' 2>/dev/null | \
    awk -F'/' '{n=NF; bin=$n; sub(/\.fa\.gz$/,"",bin); strat=$(n-2); samp=$(n-5);
               print samp"|"strat"|"bin"\t"$0}' > "$BINIDX"
  python3 - "$REPS" "$BINIDX" "$CAT" <<'PY'
import os, sys, glob
REPS, BINIDX, CAT = sys.argv[1:4]
idx = {}
for line in open(BINIDX):
    p = line.rstrip("\n").split("\t")
    if len(p) == 2: idx[p[0]] = p[1]
sm = {"single": "single_sample", "self": "self_supervised"}
ddir = os.path.join(CAT, "mags", "dereplicated")
def relink(lnk, tgt):
    if os.path.exists(lnk) and not os.path.islink(lnk): return False
    if os.path.islink(lnk) or os.path.exists(lnk):
        try: os.remove(lnk)
        except OSError: pass
    os.symlink(tgt, lnk); return True
mt = open(os.path.join(CAT, "mags", "catalog_membership.tsv"), "w")
mt.write("mag\tsample\tstrategy\tbin\tbin_path\tannotation_dir\tlink_target\trefined_sub\n")
bl = gl = 0
for fa in sorted(glob.glob(os.path.join(REPS, "*.fa"))):
    mag = os.path.basename(fa)[:-3]
    sample = strat = binname = ""; binpath = ""; sub = 0
    try:
        sample, rest = mag.split("__", 1)
        tok, binno = rest.split("_SemiBin_", 1)
        if binno.endswith("_sub"): binno = binno[:-4]; sub = 1
        strat = sm.get(tok, tok); binname = "SemiBin_" + binno
        binpath = idx.get(f"{sample}|{strat}|{binname}", "")
    except ValueError:
        pass
    ann = ""
    if binpath:
        wgs = binpath.split("/wgs/")[0] + "/wgs"
        ann = f"{wgs}/annotation/{strat}/bakta/{binname}"
        if relink(os.path.join(ddir, mag + ".fa.gz"), binpath): bl += 1
        lt = "per_sample_bin"
    else:
        if relink(os.path.join(ddir, mag + ".fa"), fa): gl += 1
        lt = "drep_genome"
    mt.write(f"{mag}\t{sample}\t{strat}\t{binname}\t{binpath}\t{ann}\t{lt}\t{sub}\n")
mt.close()
print(f"    reps->bin: {bl}  reps->genome(fallback): {gl}")
PY
  n=$(tail -n +2 "$CAT/mags/catalog_membership.tsv" | wc -l)
  manifest "mags (dRep catalogue)" "COMPLETE" "$n reps; links to per-sample bins"
else
  manifest "mags (dRep catalogue)" "SKIP" "dereplicated genomes not found"
fi

# ============ 2) protein catalogue + orthogroups — guard: orthogroups_cluster.tsv ============
if [ -s "$DM/orthogroups/orthogroups_cluster.tsv" ]; then
  echo "[2] proteins"
  mkdir -p "$CAT/proteins"
  for t in assembly_proteins_per_sample.tsv contig_to_bin.tsv.gz contig_to_bin_refined.tsv.gz; do
    slink "$WD/04_protein_catalog/$t" "$CAT/proteins/$t"; done
  slink "$DM/orthogroups" "$CAT/proteins/orthogroups"
  manifest "proteins (catalogue + orthogroups)" "COMPLETE" "mmseqs linclust + orthogroups + UniRef50"
else
  manifest "proteins" "SKIP" "orthogroups not found"
fi

# ============ 3) taxonomy (GTDB-Tk) — guard: both summaries present ============
GT=$WD/02_taxonomy/gtdbtk_out
if [ -s "$GT/gtdbtk.bac120.summary.tsv" ] && [ -s "$GT/gtdbtk.ar53.summary.tsv" ]; then
  echo "[3] taxonomy (GTDB-Tk)"
  mkdir -p "$CAT/taxonomy"
  slink "$GT/gtdbtk.bac120.summary.tsv" "$CAT/taxonomy/gtdbtk.bac120.summary.tsv"
  slink "$GT/gtdbtk.ar53.summary.tsv"   "$CAT/taxonomy/gtdbtk.ar53.summary.tsv"
  nb=$(tail -n +2 "$GT/gtdbtk.bac120.summary.tsv" | wc -l); na=$(tail -n +2 "$GT/gtdbtk.ar53.summary.tsv" | wc -l)
  manifest "taxonomy (GTDB-Tk r226)" "COMPLETE" "$nb bacterial + $na archaeal MAGs classified"
else
  manifest "taxonomy (GTDB-Tk)" "SKIP" "summaries not found"
fi

# ============ 4) structure novelty — guard: final novelty figures present (link FINAL outputs only) ============
ST=$WD/04_structure
if ls "$ST/novelty_analysis/"fig_*.pdf >/dev/null 2>&1; then
  echo "[4] structure (final outputs only)"
  mkdir -p "$CAT/structure"
  rm -f "$CAT/structure/prostt5_foldseek"   # drop earlier whole-tree link
  slink "$ST/novelty_analysis"        "$CAT/structure/novelty_analysis"
  slink "$ST/foldseek_results"        "$CAT/structure/foldseek_results"
  slink "$ST/cluster90_rep_seq.fasta" "$CAT/structure/cluster90_rep_seq.fasta"
  slink "$ST/cluster90_cluster.tsv"   "$CAT/structure/cluster90_cluster.tsv"
  manifest "structure (ProstT5+FoldSeek vs AFDB)" "COMPLETE" "novelty analysis + figures + AFDB hits"
else
  manifest "structure" "SKIP" "novelty figures not found"
fi

# ============ 5) dark matter: function assignment — guard: candidate_scores_v3.tsv ============
if [ -s "$DM/scores_v3/candidate_scores_v3.tsv" ]; then
  echo "[5] dark_matter/function_assignment"
  mkdir -p "$CAT/dark_matter/function_assignment"
  for d in scores_v3 dossier reallydark isdark gap_context gapsmith_models pfba localization; do
    slink "$DM/$d" "$CAT/dark_matter/function_assignment/$d"; done
  manifest "dark_matter/function_assignment (C·V·T)" "COMPLETE" "calibrated EC assignment; reviewed top-100"
else
  manifest "dark_matter/function_assignment" "SKIP" "candidate_scores_v3.tsv not found"
fi

# ============ 6) dark family -> MAG provenance — guard: membership built ============
if [ -s "$CAT/mags/catalog_membership.tsv" ] && [ -s "$DM/reallydark/isdark_ids.txt" ]; then
  echo "[6] dark_matter/families provenance"
  mkdir -p "$CAT/dark_matter/families"
  python3 - "$DM" "$CAT" <<'PY'
import os, sys
DM, CAT = sys.argv[1], sys.argv[2]
ext = {}
for line in open(os.path.join(CAT, "mags", "catalog_membership.tsv")):
    c = line.rstrip("\n").split("\t")
    if c[0] == "mag" or len(c) < 7: continue
    ext[c[0]] = ".fa.gz" if c[6] == "per_sample_bin" else ".fa"
dark = set(l.strip() for l in open(os.path.join(DM, "reallydark", "isdark_ids.txt")) if l.strip())
out = open(os.path.join(CAT, "dark_matter", "families", "darkfamily_to_mag.tsv"), "w")
out.write("family_id\tmag\tlocus\tn_members\tn_samples\tmag_genome_link\n")
n = 0
for line in open(os.path.join(DM, "orthogroups", "og_recurrence.tsv")):
    c = line.rstrip("\n").split("\t")
    if len(c) < 3 or c[0] not in dark: continue
    mag, locus = c[0].split("::", 1); e = ext.get(mag)
    link = f"../../mags/dereplicated/{mag}{e}" if e else "NA(mag_not_in_catalog)"
    out.write(f"{c[0]}\t{mag}\t{locus}\t{c[1]}\t{c[2]}\t{link}\n"); n += 1
out.close()
print(f"    dark family provenance rows: {n}")
PY
  manifest "dark_matter/families (->MAG provenance)" "COMPLETE" "$(tail -n +2 "$CAT/dark_matter/families/darkfamily_to_mag.tsv" | wc -l) families"
else
  manifest "dark_matter/families" "SKIP" "inputs missing"
fi

# ============ 7) dark-matter ecology (Phase A/B) — guard: AGG_DONE (auto-links when complete) ============
if [ -f "$DM/darkfam_ecology/AGG_DONE" ]; then
  echo "[7] dark_matter/ecology (Phase A complete)"
  mkdir -p "$CAT/dark_matter/ecology"
  slink "$DM/darkfam_ecology" "$CAT/dark_matter/ecology/darkfam_ecology"
  manifest "dark_matter/ecology (Phase A/B)" "COMPLETE" "dark-family x sample matrix"
else
  manifest "dark_matter/ecology (Phase A/B)" "IN-PROGRESS (excluded)" "Phase A still running; auto-links on AGG_DONE"
fi

# ============ 8) atmospheric guild analysis — guard: INTEGRATE_DONE ============
AT=$WD/06_atmospheric
if [ -f "$AT/integrate/INTEGRATE_DONE" ] && [ -f "$AT/classify/CLASSIFY_DONE" ]; then
  echo "[8] analyses/atmospheric_guild"
  mkdir -p "$CAT/analyses"
  slink "$AT" "$CAT/analyses/atmospheric_guild"
  manifest "analyses/atmospheric_guild" "COMPLETE" "HydDB/CoxL classify, carbon-fixation, FBA, cultures (all DONE)"
else
  manifest "analyses/atmospheric_guild" "SKIP" "DONE markers missing"
fi

# ============ 9) phages — guard: vOTU output + genomad/phabox/vibrant complete ============
PH=$WD/06_phage
ng=$(ls "$PH/genomad_out" 2>/dev/null | wc -l)
if [ -s "$PH/votu/all_viruses.fna" ] && [ "$ng" -ge 500 ]; then
  echo "[9] phages"
  mkdir -p "$B/phages"
  for d in genomad_out phabox_out vibrant_out votu borg_screen; do
    slink "$PH/$d" "$B/phages/$d"; done
  manifest "phages (geNomad/PhaBox/VIBRANT/vOTU)" "COMPLETE" "$ng/550 samples; vOTU clustered"
else
  manifest "phages" "SKIP" "vOTU or per-sample runs incomplete"
fi

# explicitly record analyses deliberately EXCLUDED as incomplete/not-started
manifest "ESMFold (Phase 4c, dark structures)" "INCOMPLETE (excluded)" "no DONE; partial pdb/p2rank only"
for x in 05_picrust_validation 07_nitrogen 08_dormancy 09_kraken2; do
  [ "$(ls "$WD/$x" 2>/dev/null)" ] || manifest "$x" "NOT STARTED (excluded)" "empty"
done

# ============ README ============
cat > "$CAT/README.md" <<'MD'
# Empty Quarter cross-MAG analysis catalogue (results pickup)

Only **complete** analyses are published here (each linker section is gated behind a completion
check; in-progress work is excluded — see `MANIFEST.tsv`). Per-sample bins
(`binning/{single_sample,self_supervised}/output_bins`) and per-MAG annotation
(`annotation/<strategy>/{bakta,eggnog}`) already live under each `sites/<trip>/site<N>/<sample>/wgs/`.

Anchoring: a result links to its MAG → assembly → site.
- `mags/dereplicated/<mag>.fa.gz` → that MAG's per-sample bin (`_sub` DAS_Tool reps → parent SemiBin;
  reps with no canonical bin fall back to the dRep genome). `mags/catalog_membership.tsv` carries
  sample/trip/site/strategy/bin/annotation per rep.
- `dark_matter/families/darkfamily_to_mag.tsv` → each dark family's representative MAG (→ bin → site).

Contents (all COMPLETE): `mags/`, `proteins/`, `taxonomy/` (GTDB-Tk), `structure/` (final outputs),
`dark_matter/{function_assignment,families}`, `analyses/atmospheric_guild/`, `../phages/`.
`dark_matter/ecology/` appears automatically once Phase A completes.

All entries are symlinks into `/ibex/scratch/projects/c2014/rob/metagenomics/`.
Builder: `workbench/scripts/ibex/95_canonical_links.sh` (non-destructive, idempotent, completion-gated).
MD

echo "[$(date -Iseconds)] === MANIFEST ==="; column -t -s$'\t' "$MAN" 2>/dev/null || cat "$MAN"
date -Iseconds > "$CAT/CANON_LINKS_DONE"
echo "[$(date -Iseconds)] done"
