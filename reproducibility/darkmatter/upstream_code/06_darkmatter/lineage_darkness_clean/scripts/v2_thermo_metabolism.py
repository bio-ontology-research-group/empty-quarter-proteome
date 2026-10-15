#!/usr/bin/env python3
"""v2_thermo_metabolism - what the 237 Thermoplasmatota representatives CAN encode.

The ecology script (`v2_thermo_ecology.py`) established clean ABSENCES. This one reads the
positive content of the same genomes: which KEGG modules they encode, at what fraction of
their steps, and how that compares with archaeal and completeness-matched bacterial controls.

Module definitions and the presence rule are reused verbatim from the differential-function
analysis (`06_metabolism/differential_function_clean/scripts/00_build_features.py` and its
`results/feature_definitions.json`): 399 modules in DRAM's KEGG step form, a module is
PRESENT when at least 2/3 of its top-level steps are satisfied, a step being satisfied if
any of its alternative KOs is annotated. The raw fraction of steps encoded is reported
alongside the call, and for the single-step respiratory complexes (M00144, M00153, M00155,
M00157, M00159 ... one step with up to 17 alternative subunits, so the 2/3 rule reduces to
"any one subunit") the fraction of the module's KOs encoded is reported as well, and strict
subunit-count features are defined separately under `ko_set`.

Groups. thermo_237 and its two classes (SW-10-69-26 185, E2 51; the single Thermoplasmata rep
is kept in thermo_237 but not given its own row). Comparators are reported twice: over the
whole pool and over a 1:1 completeness-matched subset (greedy nearest neighbour without
replacement, caliper 1.0 completeness points), because the 237 sit at a HIGHER median
completeness (85.42) than the archaeal pools, so unmatched absences would flatter them.
`archaea_other` follows the brief and means all 591 archaea minus the 237, so it contains the
269 Halobacteriota; `archaea_nonhalo` is the 85 that remain.

Inputs (read-only)
  results/v2_reps.tsv                      quality, size, lineage for the 5,438 reps
  results/v2_thermo_ecology_by_rep.tsv     GTDB class per rep (same 5,438)
  W/04_functional/repko_clean/results/rep_ko_v2.tsv        rep -> comma-joined KO list
  W/06_metabolism/differential_function_clean/results/feature_definitions.json
  W/databases/dram/kofam_ko_list.tsv       KO names, used for the EC-class counts

Outputs
  results/v2_thermo_metabolism_summary.tsv   long: feature x group
  results/v2_thermo_metabolism_by_rep.tsv    wide: 5,438 reps x focus features
  results/v2_thermo_metabolism_report.txt    contrasts, standouts, auxotrophy, caveats
"""
import os
import sys
import json
import csv
import math
import statistics as st
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
W = "/ibex/scratch/projects/c2014/rob/metagenomics/"

REPKO = W + "04_functional/repko_clean/results/rep_ko_v2.tsv"
FEATDEF = W + "06_metabolism/differential_function_clean/results/feature_definitions.json"
KOLIST = W + "databases/dram/kofam_ko_list.tsv"
REPS = os.path.join(RES, "v2_reps.tsv")
ECOL = os.path.join(RES, "v2_thermo_ecology_by_rep.tsv")

MODULE_FRACTION = 2.0 / 3.0     # same rule as 00_build_features.py
CALIPER = 1.0                   # completeness points, 1:1 matching

# --------------------------------------------------------------------------- KO feature sets
# Every KO identity below was read off databases/dram/kofam_ko_list.tsv, not from memory.
# rule: ("any", n) -> present when at least n of the KOs are annotated
#       ("all",)   -> present when every KO is annotated
KO_FEATURES = {
    # --- respiration, strict subunit counts (the modules are one step of many alternatives)
    "complexI_nuo_half":      ("respiration", ("any", 7),
                               ["K00330", "K00331", "K00332", "K00333", "K00334", "K00335",
                                "K00336", "K00337", "K00338", "K00339", "K00340", "K00341",
                                "K00342", "K00343"]),
    "complexII_sdh_2sub":     ("respiration", ("any", 2),
                               ["K00239", "K00240", "K00241", "K00242",
                                "K00234", "K00235", "K00236", "K00237"]),
    "complexIII_bc1_2sub":    ("respiration", ("any", 2),
                               ["K00410", "K00411", "K00412", "K00413", "K03886", "K03887",
                                "K03888", "K03889", "K03890", "K03891"]),
    "complexIV_coxA_coxB":    ("respiration", ("all",), ["K02274", "K02275"]),
    "complexIV_coxA":         ("respiration", ("any", 1), ["K02274"]),
    "bd_oxidase":             ("respiration", ("any", 1), ["K00425", "K00426"]),
    "cbb3_oxidase":           ("respiration", ("any", 1), ["K00404", "K00405", "K00406", "K00407"]),
    "ATPase_F_half":          ("respiration", ("any", 4),
                               ["K02108", "K02109", "K02110", "K02111", "K02112", "K02113",
                                "K02114", "K02115"]),
    "ATPase_VA_half":         ("respiration", ("any", 4),
                               ["K02107", "K02117", "K02118", "K02119", "K02120", "K02121",
                                "K02122", "K02123", "K02124"]),
    # --- central carbon, single diagnostic enzymes
    "pyruvate_ferredoxin_OR": ("central_carbon", ("any", 1), ["K00169", "K00170", "K00171",
                                                              "K00172", "K03737"]),
    "pyruvate_dehydrogenase": ("central_carbon", ("any", 1), ["K00161", "K00162", "K00163"]),
    "oxoglutarate_ferredoxin_OR": ("central_carbon", ("any", 1), ["K00174", "K00175"]),
    "oxoglutarate_dehydrogenase": ("central_carbon", ("any", 1), ["K00164", "K00658"]),
    "acetyl_CoA_synthetase":  ("central_carbon", ("any", 1), ["K01895"]),
    "pta_ackA":               ("central_carbon", ("all",), ["K00625", "K00925"]),
    "pyruvate_kinase":        ("central_carbon", ("any", 1), ["K00873"]),
    "pep_carboxykinase":      ("central_carbon", ("any", 1), ["K01610", "K01596"]),
    "rieske_bc_K03886":       ("respiration", ("any", 1), ["K03886"]),
    # branched-chain amino acids: the carbon-skeleton route (ilvBN/ilvC/ilvD) against the
    # terminal transamination step, which is what tells an auxotroph from a prototroph
    "ilv_carbon_skeleton":    ("aa_biosynthesis", ("any", 2), ["K01652", "K01653", "K11258",
                                                               "K00053", "K01687"]),
    "bcaa_aminotransferase":  ("aa_biosynthesis", ("any", 1), ["K00826"]),
    # --- nitrogen
    "urease_ureABC":          ("nitrogen", ("any", 2), ["K01428", "K01429", "K01430", "K14048"]),
    "ammonium_transporter_amt": ("nitrogen", ("any", 1), ["K03320"]),
    "glutamine_synthetase":   ("nitrogen", ("any", 1), ["K01915"]),
    "glutamate_synthase_GOGAT": ("nitrogen", ("any", 1), ["K00265", "K00266", "K00284"]),
    "glutamate_dehydrogenase": ("nitrogen", ("any", 1), ["K00260", "K00261", "K00262"]),
    "nitrate_reductase_nar_nap": ("nitrogen", ("any", 1), ["K00370", "K00371", "K00374",
                                                           "K02567", "K02568"]),
    "nitrate_reductase_assim": ("nitrogen", ("any", 1), ["K00360", "K00367", "K00372", "K10534"]),
    "nitrite_reductase_any":  ("nitrogen", ("any", 1), ["K00362", "K00363", "K00366", "K00368",
                                                        "K03385", "K15864", "K17877"]),
    "nitric_oxide_reductase": ("nitrogen", ("any", 1), ["K04561", "K02305", "K15877"]),
    "nitrous_oxide_reductase": ("nitrogen", ("any", 1), ["K00376"]),
    # --- sulfur
    "sulfate_adenylyltransferase": ("sulfur", ("any", 1), ["K00955", "K00956", "K00957",
                                                           "K00958", "K13811"]),
    "paps_reductase":         ("sulfur", ("any", 1), ["K00390"]),
    "sulfite_reductase":      ("sulfur", ("any", 1), ["K00380", "K00381", "K00392"]),
    "sqr_sulfide_oxidation":  ("sulfur", ("any", 1), ["K17218"]),
    "sox_thiosulfate":        ("sulfur", ("any", 2), ["K17222", "K17223", "K17224", "K17225",
                                                      "K17226", "K17227"]),
    "sulfate_transport_cysP": ("sulfur", ("any", 1), ["K02048"]),
    # --- osmolytes
    "trehalose_otsAB":        ("osmolyte", ("all",), ["K00697", "K01087"]),
    "trehalose_treYZ":        ("osmolyte", ("all",), ["K06044", "K01236"]),
    "trehalose_treS":         ("osmolyte", ("any", 1), ["K05343"]),
    "betaine_synthesis_betAB": ("osmolyte", ("any", 1), ["K00108", "K00130", "K11440", "K17755"]),
    "betaine_transport_proVWX": ("osmolyte", ("any", 1), ["K02000", "K02001", "K02002", "K05020"]),
    "ectoine_ectABC":         ("osmolyte", ("any", 2), ["K00836", "K06718", "K06720"]),
    "ectoine_ectC":           ("osmolyte", ("any", 1), ["K06720"]),
    # --- stress
    "photolyase":             ("stress", ("any", 1), ["K01669"]),
    "catalase":               ("stress", ("any", 1), ["K03781", "K03782", "K07217"]),
    "superoxide_dismutase":   ("stress", ("any", 1), ["K04564", "K04565", "K16627", "K00518"]),
    "kaiC":                   ("stress", ("any", 1), ["K08482"]),
    "dps":                    ("stress", ("any", 1), ["K04047"]),
    "rhodopsin":              ("stress", ("any", 1), ["K04641", "K04642", "K04643"]),
    # --- motility, adhesion, defence
    "bacterial_flagellum":    ("motility", ("any", 2), ["K02406", "K02400", "K02410", "K02416",
                                                        "K02556", "K02557", "K02387", "K02392"]),
    "archaellum_fla":         ("motility", ("any", 2), ["K07324", "K07325", "K07327", "K07328",
                                                        "K07329", "K07330", "K07331", "K07332",
                                                        "K07333"]),
    "chemotaxis_cheAY":       ("motility", ("any", 1), ["K03407", "K03413"]),
    "type_IV_pilus":          ("motility", ("any", 2), ["K02650", "K02652", "K02653", "K02655",
                                                        "K02656", "K02661", "K02662", "K02663",
                                                        "K02666", "K02669"]),
    "crispr_cas1_cas2":       ("defence", ("any", 1), ["K15342", "K09951"]),
    "crispr_cas3_typeI":      ("defence", ("any", 1), ["K07012", "K19046"]),
    "crispr_cas9_typeII":     ("defence", ("any", 1), ["K09952"]),
    "crispr_any":             ("defence", ("any", 1), ["K15342", "K09951", "K07012", "K09952",
                                                       "K07464", "K09002", "K19046"]),
    # --- uptake: what the cell can pull in, which is the positive read on carbon source
    "peptide_nickel_ABC":     ("uptake", ("any", 2), ["K02033", "K02034", "K02035"]),
    "oligopeptide_POT":       ("uptake", ("any", 1), ["K03305"]),
    "branched_chain_aa_ABC":  ("uptake", ("any", 2), ["K01995", "K01996", "K01997", "K01998",
                                                      "K01999"]),
    "sugar_ABC_msm":          ("uptake", ("any", 2), ["K10117", "K10118", "K10119", "K02025",
                                                      "K02026", "K02027"]),
}

# --------------------------------------------------------------------------- module groupings
FOCUS_MODULES = {
    "central_carbon": ["M00001", "M00002", "M00003", "M00004", "M00006", "M00007", "M00580",
                       "M00008", "M00308", "M00309", "M00633", "M00009", "M00010", "M00011",
                       "M00012", "M00307", "M00579", "M00173", "M00374", "M00375", "M00377",
                       "M00165"],
    "respiration": ["M00144", "M00148", "M00149", "M00150", "M00151", "M00152", "M00153",
                    "M00154", "M00155", "M00156", "M00157", "M00159", "M00416", "M00417"],
    "carbon_source": ["M00014", "M00061", "M00081", "M00550", "M00552", "M00631", "M00632",
                      "M00854", "M00855", "M00549", "M00554", "M00086", "M00087", "M00098",
                      "M00082", "M00083", "M00032", "M00036", "M00035", "M00044", "M00045",
                      "M00038", "M00046", "M00546", "M00027", "M00013", "M00373", "M00741"],
    "nitrogen": ["M00028", "M00029", "M00175", "M00528", "M00529", "M00530", "M00531",
                 "M00615", "M00804"],
    "sulfur": ["M00176", "M00595", "M00596"],
    "osmolyte": ["M00033", "M00555", "M00565"],
    "aa_biosynthesis": ["M00015", "M00016", "M00017", "M00018", "M00019", "M00020", "M00021",
                        "M00022", "M00023", "M00024", "M00025", "M00026", "M00028", "M00030",
                        "M00031", "M00040", "M00338", "M00432", "M00433", "M00525", "M00526",
                        "M00527", "M00535", "M00570", "M00609", "M00763", "M00844", "M00845"],
    "cofactor": ["M00095", "M00096", "M00115", "M00116", "M00117", "M00118", "M00119", "M00120",
                 "M00121", "M00122", "M00123", "M00124", "M00125", "M00126", "M00127", "M00140",
                 "M00358", "M00378", "M00573", "M00577", "M00836", "M00840", "M00841", "M00846",
                 "M00847", "M00849"],
    "nucleotide": ["M00048", "M00049", "M00050", "M00051", "M00052", "M00053", "M00005"],
    "lipid_envelope": ["M00060", "M00866", "M00093", "M00364", "M00365", "M00093"],
}

# amino acid -> the KEGG modules whose presence would make the cell a prototroph for it.
# Ala, Asp, Asn, Glu, Gly and Gln have no KEGG biosynthesis module and are therefore
# untestable this way; that is stated in the report rather than scored as an auxotrophy.
AA_MODULES = {
    "Arg": ["M00844", "M00845"],
    "Cys": ["M00021", "M00338", "M00609"],
    "His": ["M00026"],
    "Ile": ["M00019", "M00535", "M00570"],
    "Leu": ["M00432"],
    "Lys": ["M00016", "M00030", "M00031", "M00433", "M00525", "M00526", "M00527"],
    "Met": ["M00017"],
    "Phe": ["M00024"],
    "Pro": ["M00015"],
    "Ser": ["M00020"],
    "Thr": ["M00018"],
    "Trp": ["M00023"],
    "Tyr": ["M00025", "M00040"],
    "Val": ["M00019"],
    "chorismate": ["M00022"],
    "ornithine": ["M00028", "M00763"],
}

COFACTOR_MODULES = {
    "NAD": ["M00115"], "CoA": ["M00120"], "pantothenate": ["M00119"], "riboflavin": ["M00125"],
    "thiamine": ["M00127"], "pyridoxal": ["M00124"], "biotin": ["M00123", "M00573", "M00577"],
    "folate": ["M00126", "M00840", "M00841"], "cobalamin": ["M00122"],
    "heme": ["M00121", "M00847", "M00868"], "siroheme": ["M00846"], "menaquinone": ["M00116"],
    "ubiquinone": ["M00117"], "glutathione": ["M00118"], "F420": ["M00378"],
    "coenzyme_M": ["M00358"], "isoprenoid": ["M00095", "M00096", "M00849"],
}


# --------------------------------------------------------------------------- loaders
def load_ko_sets():
    """rep_ko_v2.tsv has no header: genome TAB comma-joined KOs."""
    out = {}
    with open(REPKO) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if not p[0]:
                continue
            out[p[0]] = set(k for k in (p[1].split(",") if len(p) > 1 and p[1] else []) if k)
    return out


def load_reps():
    reps = pd.read_csv(REPS, sep="\t",
                       usecols=["genome", "completeness", "contamination", "length",
                                "phylum", "domain", "dark_fraction"])
    eco = pd.read_csv(ECOL, sep="\t", usecols=["genome", "class"])
    reps = reps.merge(eco, on="genome", how="left")
    reps["domain"] = reps["domain"].astype(str).str.replace("d__", "", regex=False)
    reps["phylum"] = reps["phylum"].astype(str).str.replace("p__", "", regex=False)
    return reps


def load_ko_names():
    kl = pd.read_csv(KOLIST, sep="\t", usecols=["knum", "definition"])
    return dict(zip(kl["knum"], kl["definition"].astype(str)))


# --------------------------------------------------------------------------- module presence
def module_tables(ko_sets, step_groups):
    """Per genome: satisfied steps and encoded KOs per module, via a KO -> (module, step) index."""
    idx = defaultdict(list)
    n_steps, n_kos = {}, {}
    for mid, steps in step_groups.items():
        n_steps[mid] = len(steps)
        n_kos[mid] = len(set(k for s in steps for k in s))
        for si, step in enumerate(steps):
            for ko in step:
                idx[ko].append((mid, si))

    genomes = sorted(ko_sets)
    mids = sorted(step_groups)
    mpos = {m: i for i, m in enumerate(mids)}
    sat = np.zeros((len(genomes), len(mids)), dtype=np.int16)
    enc = np.zeros((len(genomes), len(mids)), dtype=np.int16)
    for gi, g in enumerate(genomes):
        hit_steps = defaultdict(set)
        hit_kos = defaultdict(set)
        for ko in ko_sets[g]:
            for mid, si in idx.get(ko, ()):
                hit_steps[mid].add(si)
                hit_kos[mid].add(ko)
        for mid, s in hit_steps.items():
            sat[gi, mpos[mid]] = len(s)
        for mid, s in hit_kos.items():
            enc[gi, mpos[mid]] = len(s)

    tot_steps = np.array([n_steps[m] for m in mids], dtype=float)
    tot_kos = np.array([n_kos[m] for m in mids], dtype=float)
    steps_frac = sat / tot_steps
    ko_frac = enc / tot_kos
    present = steps_frac >= MODULE_FRACTION - 1e-9
    return genomes, mids, steps_frac, ko_frac, present


# --------------------------------------------------------------------------- matching
def match_on_completeness(target, pool, caliper=CALIPER):
    """Greedy 1:1 nearest-neighbour on completeness without replacement."""
    cand = []
    for tg, tc in target:
        for pg, pc in pool:
            d = abs(tc - pc)
            if d <= caliper:
                cand.append((d, tg, pg))
    cand.sort()
    used, done, pairs = set(), set(), []
    for d, tg, pg in cand:
        if tg in done or pg in used:
            continue
        done.add(tg)
        used.add(pg)
        pairs.append((tg, pg))
    return pairs


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


# --------------------------------------------------------------------------- main
def main():
    os.makedirs(RES, exist_ok=True)
    log = []

    def say(s=""):
        print(s)
        log.append(s)

    reps = load_reps()
    ko_sets = load_ko_sets()
    assert set(reps["genome"]) == set(ko_sets), "rep set and KO table disagree"
    spec = json.load(open(FEATDEF))
    step_groups = spec["module"]["step_groups"]
    mod_names = spec["module"]["names"]
    ko_names = load_ko_names()

    genomes, mids, steps_frac, ko_frac, present = module_tables(ko_sets, step_groups)
    gpos = {g: i for i, g in enumerate(genomes)}
    mpos = {m: i for i, m in enumerate(mids)}
    say(f"modules {len(mids)}  genomes {len(genomes)}  "
        f"median modules present per genome "
        f"{np.median(present.sum(axis=1)):.0f}")

    # ---------------- KO-derived per-genome counts
    gh_kos = {k for k, d in ko_names.items() if "EC:3.2.1" in d}
    pep_kos = {k for k, d in ko_names.items() if "EC:3.4." in d}
    tr_kos = {k for k, d in ko_names.items()
              if "transport system" in d.lower() or "abc transporter" in d.lower()}
    carb_words = ("sugar", "monosaccharide", "oligosaccharide", "maltose", "glucose", "xylose",
                  "arabinose", "rhamnose", "fructose", "trehalose", "cellobiose", "ribose",
                  "glycerol", "mannose", "galactose", "fucose", "chitobiose", "lactose")
    pep_words = ("amino acid", "peptide", "oligopeptide", "dipeptide", "branched-chain",
                 "polar amino", "glutamine transport", "arginine transport", "histidine transport")
    carb_tr = {k for k in tr_kos if any(w in ko_names[k].lower() for w in carb_words)}
    peptr = {k for k in tr_kos if any(w in ko_names[k].lower() for w in pep_words)}
    say(f"KO classes: glycoside hydrolase (EC:3.2.1) {len(gh_kos)}, peptidase (EC:3.4.) "
        f"{len(pep_kos)}, transport-system {len(tr_kos)} (carbohydrate {len(carb_tr)}, "
        f"peptide/amino-acid {len(peptr)})")

    counts = {}
    for g in genomes:
        ks = ko_sets[g]
        counts[g] = {
            "ko_richness": len(ks),
            "n_modules_present": int(present[gpos[g]].sum()),
            "gh_ko_count": len(ks & gh_kos),
            "peptidase_ko_count": len(ks & pep_kos),
            "transporter_ko_count": len(ks & tr_kos),
            "carb_transporter_ko_count": len(ks & carb_tr),
            "pep_aa_transporter_ko_count": len(ks & peptr),
        }

    # ---------------- KO-set features
    kof_present = {}
    for name, (_cat, rule, kos) in KO_FEATURES.items():
        col = np.zeros(len(genomes), dtype=bool)
        frac = np.zeros(len(genomes))
        for gi, g in enumerate(genomes):
            hit = len(ko_sets[g].intersection(kos))
            frac[gi] = hit / len(kos)
            col[gi] = (hit == len(kos)) if rule[0] == "all" else (hit >= rule[1])
        kof_present[name] = (col, frac)

    # ---------------- amino acid / cofactor prototrophy
    aa_present, cof_present = {}, {}
    for aa, mods in AA_MODULES.items():
        mm = [m for m in mods if m in mpos]
        aa_present[aa] = present[:, [mpos[m] for m in mm]].any(axis=1)
    for cf, mods in COFACTOR_MODULES.items():
        mm = [m for m in mods if m in mpos]
        cof_present[cf] = present[:, [mpos[m] for m in mm]].any(axis=1)
    aa_score = np.vstack([aa_present[a] for a in AA_MODULES]).sum(axis=0)
    cof_score = np.vstack([cof_present[c] for c in COFACTOR_MODULES]).sum(axis=0)

    # ---------------- groups
    comp = dict(zip(reps["genome"], reps["completeness"]))
    is_th = reps["phylum"] == "Thermoplasmatota"
    thermo = list(reps.loc[is_th, "genome"])
    halo = list(reps.loc[reps["phylum"] == "Halobacteriota", "genome"])
    arch_other = list(reps.loc[(reps["domain"] == "Archaea") & ~is_th, "genome"])
    arch_nonhalo = list(reps.loc[(reps["domain"] == "Archaea") & ~is_th
                                 & (reps["phylum"] != "Halobacteriota"), "genome"])
    bact = list(reps.loc[reps["domain"] == "Bacteria", "genome"])
    cls = dict(zip(reps["genome"], reps["class"]))

    tgt = [(g, comp[g]) for g in thermo]
    groups = {
        "thermo_237": thermo,
        "thermo_SW-10-69-26": [g for g in thermo if cls.get(g) == "SW-10-69-26"],
        "thermo_E2": [g for g in thermo if cls.get(g) == "E2"],
        "halo_all": halo,
        "archaea_other_all": arch_other,
        "archaea_nonhalo_all": arch_nonhalo,
        "bacteria_all": bact,
    }
    matched_pairs = {}
    for label, pool in (("halo", halo), ("archaea_other", arch_other), ("bacteria", bact)):
        pairs = match_on_completeness(tgt, [(g, comp[g]) for g in pool])
        matched_pairs[label] = pairs
        groups[label + "_matched"] = [p for _, p in pairs]
        say(f"matched {label}: {len(pairs)} of 237 pairs, caliper {CALIPER} completeness points; "
            f"median completeness thermo {st.median([comp[t] for t, _ in pairs]):.2f} vs control "
            f"{st.median([comp[p] for _, p in pairs]):.2f}")
    # thermo subset actually used in each matched comparison
    for label in matched_pairs:
        groups["thermo_for_" + label] = [t for t, _ in matched_pairs[label]]

    for name, gl in groups.items():
        if not name.startswith("thermo_for_"):
            say(f"group {name:24s} n={len(gl):5d}  median completeness "
                f"{st.median([comp[g] for g in gl]):.2f}")

    # ---------------- long summary
    rows = []
    mcat = {}
    for cat, mm in FOCUS_MODULES.items():
        for m in mm:
            mcat.setdefault(m, cat)

    def add_rows(fid, ftype, fname, cat, vec_present, vec_sfrac, vec_kfrac):
        for gname, gl in groups.items():
            ii = [gpos[g] for g in gl]
            if not ii:
                continue
            p = vec_present[ii]
            rows.append({
                "feature": fid, "feature_type": ftype, "feature_name": fname, "category": cat,
                "group": gname, "n": len(ii), "n_present": int(p.sum()),
                "share": round(float(p.mean()), 6),
                "median_steps_frac": (round(float(np.median(vec_sfrac[ii])), 6)
                                      if vec_sfrac is not None else ""),
                "median_ko_frac": (round(float(np.median(vec_kfrac[ii])), 6)
                                   if vec_kfrac is not None else ""),
                "median_value": "", "median_value_per100comp": "",
            })

    for m in mids:
        j = mpos[m]
        add_rows(m, "module", mod_names.get(m, ""), mcat.get(m, "other"),
                 present[:, j], steps_frac[:, j], ko_frac[:, j])
    for name, (col, frac) in kof_present.items():
        add_rows(name, "ko_set", "+".join(KO_FEATURES[name][2]), KO_FEATURES[name][0],
                 col, None, frac)
    for aa, vec in aa_present.items():
        add_rows("AA_" + aa, "aa_prototrophy", "+".join(AA_MODULES[aa]), "aa_biosynthesis",
                 vec, None, None)
    for cf, vec in cof_present.items():
        add_rows("COF_" + cf, "cofactor", "+".join(COFACTOR_MODULES[cf]), "cofactor",
                 vec, None, None)

    count_cols = ["ko_richness", "n_modules_present", "gh_ko_count", "peptidase_ko_count",
                  "transporter_ko_count", "carb_transporter_ko_count",
                  "pep_aa_transporter_ko_count"]
    for cname in count_cols + ["aa_prototrophy_score", "cofactor_score"]:
        for gname, gl in groups.items():
            if cname == "aa_prototrophy_score":
                vals = [float(aa_score[gpos[g]]) for g in gl]
            elif cname == "cofactor_score":
                vals = [float(cof_score[gpos[g]]) for g in gl]
            else:
                vals = [float(counts[g][cname]) for g in gl]
            per100 = [v / (comp[g] / 100.0) for v, g in zip(vals, gl)]
            rows.append({
                "feature": cname, "feature_type": "count", "feature_name": cname,
                "category": "aggregate", "group": gname, "n": len(gl), "n_present": "",
                "share": "", "median_steps_frac": "", "median_ko_frac": "",
                "median_value": round(st.median(vals), 4),
                "median_value_per100comp": round(st.median(per100), 4),
            })

    summ = pd.DataFrame(rows)
    out_s = os.path.join(RES, "v2_thermo_metabolism_summary.tsv")
    summ.to_csv(out_s, sep="\t", index=False)
    say(f"wrote {out_s}  rows {len(summ)}")

    # ---------------- by-rep table
    focus = [m for cat in FOCUS_MODULES for m in FOCUS_MODULES[cat]]
    focus = sorted(set(m for m in focus if m in mpos))
    br = pd.DataFrame({"genome": genomes})
    br["domain"] = br["genome"].map(dict(zip(reps["genome"], reps["domain"])))
    br["phylum"] = br["genome"].map(dict(zip(reps["genome"], reps["phylum"])))
    br["class"] = br["genome"].map(cls)
    br["completeness"] = br["genome"].map(comp)
    br["contamination"] = br["genome"].map(dict(zip(reps["genome"], reps["contamination"])))
    br["length"] = br["genome"].map(dict(zip(reps["genome"], reps["length"])))
    br["dark_fraction"] = br["genome"].map(dict(zip(reps["genome"], reps["dark_fraction"])))
    gassign = {}
    for gname in ["bacteria_all", "archaea_nonhalo_all", "halo_all", "thermo_237"]:
        for g in groups[gname]:
            gassign[g] = gname.replace("_all", "")
    br["group"] = br["genome"].map(gassign)
    for gname in ["halo_matched", "archaea_other_matched", "bacteria_matched"]:
        s = set(groups[gname])
        br["in_" + gname] = br["genome"].map(lambda g: int(g in s))
    for c in count_cols:
        br[c] = [counts[g][c] for g in genomes]
    br["aa_prototrophy_score"] = aa_score
    br["cofactor_score"] = cof_score
    for aa in AA_MODULES:
        br["AA_" + aa] = aa_present[aa].astype(int)
    for cf in COFACTOR_MODULES:
        br["COF_" + cf] = cof_present[cf].astype(int)
    for name, (col, frac) in kof_present.items():
        br["KO_" + name] = col.astype(int)
    for m in focus:
        br["M_" + m] = present[:, mpos[m]].astype(int)
        br["Mfrac_" + m] = np.round(steps_frac[:, mpos[m]], 4)
    out_b = os.path.join(RES, "v2_thermo_metabolism_by_rep.tsv")
    br.to_csv(out_b, sep="\t", index=False)
    say(f"wrote {out_b}  rows {len(br)}  cols {br.shape[1]}")

    # ---------------- contrasts: thermo vs each matched control
    say()
    say("=== contrasts, Thermoplasmatota against completeness-matched controls "
        "(Fisher exact, BH q over all features x 3 controls) ===")
    feats = [(m, "module", mod_names.get(m, "")) for m in mids]
    feats += [(k, "ko_set", "") for k in KO_FEATURES]
    feats += [("AA_" + a, "aa_prototrophy", "") for a in AA_MODULES]
    feats += [("COF_" + c, "cofactor", "") for c in COFACTOR_MODULES]

    def vec_of(fid, ftype):
        if ftype == "module":
            return present[:, mpos[fid]]
        if ftype == "ko_set":
            return kof_present[fid][0]
        if ftype == "aa_prototrophy":
            return aa_present[fid[3:]]
        return cof_present[fid[4:]]

    contrasts = []
    for label in ["halo", "archaea_other", "bacteria"]:
        tlist = groups["thermo_for_" + label]
        clist = groups[label + "_matched"]
        ti = [gpos[g] for g in tlist]
        ci = [gpos[g] for g in clist]
        for fid, ftype, fname in feats:
            v = vec_of(fid, ftype)
            a, b = int(v[ti].sum()), len(ti) - int(v[ti].sum())
            c, d = int(v[ci].sum()), len(ci) - int(v[ci].sum())
            if a + c == 0 or b + d == 0:
                p = 1.0
            else:
                p = fisher_exact([[a, b], [c, d]])[1]
            contrasts.append({"control": label, "feature": fid, "feature_type": ftype,
                              "feature_name": fname, "n_thermo": len(ti), "n_control": len(ci),
                              "share_thermo": a / len(ti), "share_control": c / len(ci),
                              "delta": a / len(ti) - c / len(ci), "p": p})
    cdf = pd.DataFrame(contrasts)
    cdf["q"] = bh(cdf["p"].values)
    cdf = cdf.sort_values("delta")

    for label in ["halo", "archaea_other", "bacteria"]:
        sub = cdf[(cdf["control"] == label) & (cdf["q"] < 0.01)]
        for direction, s in (("LOWER in the 237", sub.nsmallest(12, "delta")),
                             ("HIGHER in the 237", sub.nlargest(12, "delta"))):
            say()
            say(f"--- vs {label}_matched, {direction}")
            say(f"{'feature':<12}{'thermo':>8}{'control':>9}{'delta':>8}{'q':>11}  name")
            for _, r in s.iterrows():
                say(f"{r['feature']:<12}{r['share_thermo']:>8.3f}{r['share_control']:>9.3f}"
                    f"{r['delta']:>8.3f}{r['q']:>11.2e}  {r['feature_name'][:58]}")

    # standouts: consistent sign against all three matched controls
    piv = cdf.pivot_table(index="feature", columns="control", values="delta")
    piv = piv.dropna()
    up = piv[(piv > 0.10).all(axis=1)].copy()
    dn = piv[(piv < -0.10).all(axis=1)].copy()
    up["min_delta"] = up.min(axis=1)
    dn["max_delta"] = dn.max(axis=1)
    say()
    say("=== consistent standouts: |delta| > 0.10 against ALL THREE matched controls ===")
    say(f"HIGHER in the 237: {len(up)} features")
    for f, r in up.sort_values("min_delta", ascending=False).head(20).iterrows():
        nm = mod_names.get(f, f)
        say(f"  {f:<12} halo {r['halo']:+.3f}  arch {r['archaea_other']:+.3f}  "
            f"bact {r['bacteria']:+.3f}  {nm[:60]}")
    say(f"LOWER in the 237: {len(dn)} features")
    for f, r in dn.sort_values("max_delta").head(25).iterrows():
        nm = mod_names.get(f, f)
        say(f"  {f:<12} halo {r['halo']:+.3f}  arch {r['archaea_other']:+.3f}  "
            f"bact {r['bacteria']:+.3f}  {nm[:60]}")

    # ---------------- the two classes
    say()
    say("=== SW-10-69-26 (185) vs E2 (51): features differing by more than 0.25 in share ===")
    si = [gpos[g] for g in groups["thermo_SW-10-69-26"]]
    ei = [gpos[g] for g in groups["thermo_E2"]]
    # the control that matters here: the two classes sit at the SAME completeness but NOT at the
    # same genome size, so the gaps below are reduction, not recovery
    ln = dict(zip(reps["genome"], reps["length"]))
    csz = {g: ln[g] / (comp[g] / 100.0) for g in ln}
    swg, e2g = groups["thermo_SW-10-69-26"], groups["thermo_E2"]
    say(f"completeness  SW {st.median([comp[g] for g in swg]):.2f}  E2 "
        f"{st.median([comp[g] for g in e2g]):.2f}  MWU p="
        f"{mannwhitneyu([comp[g] for g in swg], [comp[g] for g in e2g]).pvalue:.2e}")
    say(f"corrected size  SW {st.median([csz[g] for g in swg]) / 1e6:.2f} Mb  E2 "
        f"{st.median([csz[g] for g in e2g]) / 1e6:.2f} Mb  MWU p="
        f"{mannwhitneyu([csz[g] for g in swg], [csz[g] for g in e2g]).pvalue:.2e}")
    e2pairs = match_on_completeness([(g, comp[g]) for g in e2g], [(g, comp[g]) for g in swg])
    say(f"1:1 completeness-matched SW control for E2: {len(e2pairs)} of {len(e2g)} pairs; "
        f"median KO richness E2 "
        f"{st.median([counts[t]['ko_richness'] for t, _ in e2pairs]):.0f} vs SW "
        f"{st.median([counts[p]['ko_richness'] for _, p in e2pairs]):.0f}")
    crows = []
    for fid, ftype, fname in feats:
        v = vec_of(fid, ftype)
        ps, pe = v[si].mean(), v[ei].mean()
        if abs(ps - pe) >= 0.25:
            a, b = int(v[si].sum()), len(si) - int(v[si].sum())
            c, d = int(v[ei].sum()), len(ei) - int(v[ei].sum())
            p = fisher_exact([[a, b], [c, d]])[1]
            crows.append((fid, fname, ps, pe, ps - pe, p))
    crows.sort(key=lambda r: -abs(r[4]))
    say(f"{'feature':<12}{'SW':>7}{'E2':>7}{'delta':>8}{'p':>11}  name")
    for fid, fname, ps, pe, d, p in crows[:30]:
        say(f"{fid:<12}{ps:>7.3f}{pe:>7.3f}{d:>+8.3f}{p:>11.2e}  {fname[:58]}")

    # ---------------- aggregates per group
    say()
    say("=== aggregate counts per group (median; per100comp = divided by completeness/100) ===")
    say(f"{'group':<24}{'n':>5}{'comp':>7}" + "".join(f"{c.replace('_count','').replace('_ko','')[:11]:>12}"
                                                      for c in count_cols)
        + f"{'AAproto/16':>12}{'COF/17':>8}")
    for gname, gl in groups.items():
        if gname.startswith("thermo_for_"):
            continue
        vals = [st.median([counts[g][c] for g in gl]) for c in count_cols]
        say(f"{gname:<24}{len(gl):>5}{st.median([comp[g] for g in gl]):>7.2f}"
            + "".join(f"{v:>12.0f}" for v in vals)
            + f"{st.median([aa_score[gpos[g]] for g in gl]):>12.0f}"
            + f"{st.median([cof_score[gpos[g]] for g in gl]):>8.0f}")

    # ---------------- amino acid / cofactor detail
    for title, table, keys in (("amino-acid prototrophy (module-encoded)", aa_present, AA_MODULES),
                               ("cofactor biosynthesis", cof_present, COFACTOR_MODULES)):
        say()
        say(f"=== {title}: share of genomes per group ===")
        gl_order = ["thermo_237", "thermo_SW-10-69-26", "thermo_E2", "halo_matched",
                    "archaea_other_matched", "bacteria_matched", "bacteria_all"]
        say(f"{'feature':<14}" + "".join(f"{g.replace('_matched','_m')[:13]:>14}" for g in gl_order))
        for k in keys:
            v = table[k]
            say(f"{k:<14}" + "".join(
                f"{v[[gpos[g] for g in groups[go]]].mean():>14.3f}" for go in gl_order))

    # ---------------- which hydrolases and transporters, by name
    say()
    say("=== most prevalent enzymes of each class in the 237, with matched-control shares ===")
    classes = [("glycoside hydrolase (EC:3.2.1)", gh_kos), ("peptidase (EC:3.4.)", pep_kos),
               ("carbohydrate transport system", carb_tr),
               ("peptide / amino-acid transport system", peptr)]
    ti = [gpos[g] for g in groups["thermo_237"]]
    for title, kos in classes:
        prev = []
        for k in kos:
            n = sum(1 for g in groups["thermo_237"] if k in ko_sets[g])
            if n:
                prev.append((n / len(ti), k))
        prev.sort(reverse=True)
        say()
        say(f"--- {title}: {len(prev)} of {len(kos)} seen in the 237")
        say(f"{'KO':<9}{'thermo':>8}{'SW':>7}{'E2':>7}{'halo_m':>8}{'arch_m':>8}{'bact_m':>8}  name")
        for share, k in prev[:12]:
            cells = []
            for go in ["thermo_SW-10-69-26", "thermo_E2", "halo_matched",
                       "archaea_other_matched", "bacteria_matched"]:
                gl = groups[go]
                cells.append(sum(1 for g in gl if k in ko_sets[g]) / len(gl))
            say(f"{k:<9}{share:>8.3f}" + "".join(f"{c:>7.3f}" if i < 2 else f"{c:>8.3f}"
                                                 for i, c in enumerate(cells))
                + f"  {ko_names.get(k, '')[:52]}")

    out_r = os.path.join(RES, "v2_thermo_metabolism_report.txt")
    with open(out_r, "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("wrote", out_r)


if __name__ == "__main__":
    main()
