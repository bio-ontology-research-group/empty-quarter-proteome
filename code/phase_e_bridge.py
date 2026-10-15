#!/usr/bin/env python3
"""Phase E — decoding<->ecology bridge. Of the decodable dark slice (C·V·T best assignment per family),
which desert-adaptive functional categories are present, and do they associate with environmental
structure / the aridity gradient? Assignments are hypotheses (calibrated P ~0.22-0.43), so this is a
SUPPORTING line, framed as enrichment/association, not certainty.
  (1) functional composition of the decoded dark slice (desert-adaptive categories)
  (2) category vs environmentally-structured families (Fisher exact)
  (3) category abundance vs aridity (do stress/osmolyte functions track aridity?)
"""
import os, csv, re
import numpy as np
from scipy.stats import fisher_exact, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")

# desert-adaptive functional categories (keyword match on rxn_name + pathways)
CATS = {
    "osmolyte":        r"trehalose|glucosylglycerol|ectoine|mannosylglycerate|betaine|glycerophosphodiester|sucrose-6|proline.*synth|glutamate.*synth",
    "oxidative_stress":r"catalase|peroxidase|peroxiredoxin|superoxide|thioredoxin|glutathione|glutaredoxin|rubrerythrin|dismutase",
    "desiccation_dormancy": r"dipicolinate|sporulation|spore|stringent|ppgpp|guanosine.*tetraphosphate|trehalose-phosphate",
    "dna_repair_uv":   r"photolyase|excinuclease|recombinase|uvr|endonuclease.*iv|glycosylase|photoreactiv",
    "c1_oligotrophy":  r"formate|formyl|methanol|methylamine|carbon monoxide|formaldehyde|methenyltetrahydro",
    "metal_handling":  r"vanadium|molybdo|molybdenum|tungsten|metallochaperone|siderophore|ferric|cobalamin",
    "aromatic_scav":   r"dioxygenase|gentisate|benzoate|catechol|protocatechuate|aromatic|phenol|salicylate",
}
CATC = {k: re.compile(v, re.I) for k, v in CATS.items()}

def categorize(rxn, path):
    s = (rxn or "") + " " + (path or "")
    return [k for k, rx in CATC.items() if rx.search(s)]

fec = {r["family_id"]: r for r in csv.DictReader(open(os.path.join(ECO, "family_ec.tsv")), delimiter="\t")}
lab_path = os.path.join(ECO, "selection", "phase_b_family_labels.tsv")
if not os.path.exists(lab_path): lab_path = os.path.join(ECO, "phase_b_family_labels.tsv")
struct = {r["family_id"]: r["structured"] == "1" for r in csv.DictReader(open(lab_path), delimiter="\t")}

# assign categories per family
fam_cats = {}
for f, r in fec.items():
    if r.get("rxn_name"):
        cs = categorize(r["rxn_name"], r.get("pathways"))
        if cs: fam_cats[f] = cs

out = []
def w(x): out.append(x); print(x)
w("# Phase E — decoding<->ecology bridge")
w(f"families with a decoded function: {sum(1 for r in fec.values() if r.get('rxn_name'))}")
w(f"families matching a desert-adaptive category: {len(fam_cats)}\n")

# (1) composition
from collections import Counter
cc = Counter(c for cs in fam_cats.values() for c in cs)
w("## (1) desert-adaptive functional composition of the decoded dark slice")
for k, n in cc.most_common():
    w(f"  {k:20s} {n:6d} families")

# (2) category enrichment among environmentally-structured families (Fisher)
w("\n## (2) is each category over-represented among environmentally-STRUCTURED dark families?")
families = [f for f in fec if f in struct]
nS = sum(1 for f in families if struct[f]); nC = len(families) - nS
w(f"   (structured={nS}, cosmopolitan={nC})")
for k in CATS:
    inS = sum(1 for f in families if k in fam_cats.get(f, []) and struct[f])
    inC = sum(1 for f in families if k in fam_cats.get(f, []) and not struct[f])
    if inS + inC < 10: continue
    odds, p = fisher_exact([[inS, nS - inS], [inC, nC - inC]], alternative="two-sided")
    tag = "ENRICHED in structured" if odds > 1 else "depleted in structured"
    w(f"  {k:20s} struct={inS:5d} cosmo={inC:5d}  OR={odds:.2f}  p={p:.2e}  {tag if p<0.05 else ''}")

# (3) category abundance vs aridity
env = {r["sample"]: r for r in csv.DictReader(open(os.path.join(ECO, "env_matrix.tsv")), delimiter="\t")}
def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None
with open(os.path.join(ECO, "darkfam_abundance.tsv")) as fh:
    samples = fh.readline().rstrip("\n").split("\t")[1:]
    famrow = {}
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if c[0] in fam_cats: famrow[c[0]] = np.array([float(x) if x else 0.0 for x in c[1:]], float)
arid = np.array([fnum(env.get(s, {}).get("aridity_index")) for s in samples], float)
am = ~np.isnan(arid)
w("\n## (3) does category total abundance track the aridity gradient? (Spearman vs aridity index; low index = more arid)")
for k in CATS:
    fams = [f for f in fam_cats if k in fam_cats[f] and f in famrow]
    if len(fams) < 5: continue
    tot = np.sum([famrow[f] for f in fams], axis=0)
    if np.ptp(tot[am]) == 0: continue
    r, p = spearmanr(tot[am], arid[am])
    trend = "increases toward arid" if r < 0 else "decreases toward arid"
    w(f"  {k:20s} n_fam={len(fams):5d}  Spearman r={r:+.3f} p={p:.2e}  {trend if p<0.05 else ''}")

open(os.path.join(ECO, "DECODE_ECOLOGY_BRIDGE.md"), "w").write("\n".join(out) + "\n")
print("\nwritten -> results/darkmatter_ecology/DECODE_ECOLOGY_BRIDGE.md")
