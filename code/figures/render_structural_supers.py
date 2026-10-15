#!/usr/bin/env python3
"""PyMOL (open-source) driver — renders the Section 1.3 fold-gallery assets.
Run:  source ~/miniconda3/etc/profile.d/conda.sh && conda activate pymol-oss
      pymol -cq code/figures/render_structural_supers.py

Tooling history (2026-08-27): both IBEX pymol modules (2.4.1, 3.1.6) are unlicensed Schrodinger
builds that watermark every render; the chimerax/1.7 module cannot render headless (no OSMesa, and
no GL under Xvfb on login nodes). Open-source PyMOL (conda-forge, env `pymol-oss`) ray-traces on
pure CPU with no watermark — the working path.

Renders 12 exemplar candidate fold-family models (conf_qualified_dual_novel==1, pLDDT>=0.85,
quartile-mirrored lengths incl. one >500 aa showpiece) as cartoons coloured by per-residue pLDDT
(B-factor 0-1: orange 0.5 -> yellow -> marine 1.0), and two superpositions of a candidate onto its
best local AFDB match (candidate marine, match transparent wheat): F34PRr2 (alnTM 0.96, qTM 0.47)
and V45Dr1 (alnTM 0.92, qTM 0.42). AFDB v6 targets fetched from EBI (accessions from
screen/afdb_hits.raw.tsv), stored beside the renders. Consumed by make_structural_calibration_figure.py.
FOLD models live in models/fold/, REUSE models in models/reuse/ (fam_reps.tsv source column).
"""
import os, sys
from pymol import cmd

OUT = "/ibex/project/c2014/alelopezv/repos/empty-quarter-mags/figures/renders_structural"
MODELS = "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/00_structure_clean/05_validate/models/fold"

SUPER = [("F1Dr1::F1Dr1_5662978_6", "7q83", "D", "super_F1Dr1"),
         ("V45Dr1::V45Dr1_703031_2", "AF-A0AAW7WCR9", None, "super_V45Dr1"),
         ("22PRr3::22PRr3_1915643_1", "AF-A0A533VC05", None, "super_22PRr3"),
         ("CONTROL:query_pdb/KF000280.pdb", "1crw", "G", "super_control_KF000280")]

def model_path(pid):
    for sub in ("", "../reuse"):
        p = os.path.normpath(os.path.join(MODELS, sub, pid + ".pdb"))
        if os.path.exists(p):
            return p
    sys.exit(f"missing model: {pid}")

os.makedirs(OUT, exist_ok=True)
cmd.set("ray_opaque_background", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.bg_color("white")
# paper palette C (see code/figures/palette.py): teal ramp for pLDDT, teal candidate, grey known
cmd.set_color("plow", [0.910, 0.949, 0.945])
cmd.set_color("pmid", [0.388, 0.686, 0.678])
cmd.set_color("phigh", [0.059, 0.365, 0.373])
cmd.set_color("cteal", [0.122, 0.541, 0.549])
cmd.set_color("refviolet", [0.541, 0.435, 0.722])

# fixed A-per-pixel across the four superpositions: every render uses the same viewport width
# (largest candidate CA extent + 8 A), so equal-size cells show true relative scale; scale.json
# records the viewport for the compose script's 20 A scale bar.
import math, json
def ca_extent(path):
    c = [(float(l[30:38]), float(l[38:46]), float(l[46:54]))
         for l in open(path) if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    c = c[::max(1, len(c)//80)]
    return max(math.dist(p, q) for i, p in enumerate(c) for q in c[i+1:])
def super_model_path(pid):
    return os.path.normpath(os.path.join(MODELS, "..", "..", pid.split(":", 1)[1]))         if pid.startswith("CONTROL:") else model_path(pid)
_exts = {pid: ca_extent(super_model_path(pid)) for pid, _, _, _ in SUPER}
VIEW_A = max(_exts.values()) + 8.0
json.dump({"viewport_A": VIEW_A}, open(os.path.join(OUT, "scale.json"), "w"))

for pid, acc, chain, stem in SUPER:
    cmd.delete("all")
    cmd.load(super_model_path(pid), "cand")
    cmd.load(os.path.join(OUT, acc + ".pdb"), "ref")
    if chain:
        cmd.remove(f"ref and not chain {chain}")
    cmd.remove("ref and not polymer")
    cmd.hide("everything"); cmd.show("cartoon")
    cmd.color("cteal", "cand")
    cmd.color("refviolet", "ref")
    cmd.set("cartoon_transparency", 0.2, "ref")
    # CE alignment locks onto the structurally matched region (pymol's default super may fit a
    # different region than the FoldSeek match); reference residues far from the candidate are
    # trimmed so only the matching part of the known structure is drawn (stated in the caption).
    cmd.cealign("cand", "ref")
    cmd.remove("ref and not (ref within 8 of cand)")
    cmd.orient("cand")
    cmd.zoom("cand", buffer=(VIEW_A - _exts[pid]) / 2.0, complete=1)
    cmd.ray(1100, 1100)
    cmd.png(os.path.join(OUT, stem + ".png"), dpi=200)

print("done ->", OUT)
