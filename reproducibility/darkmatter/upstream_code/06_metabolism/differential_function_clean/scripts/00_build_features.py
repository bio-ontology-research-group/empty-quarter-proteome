#!/usr/bin/env python3
"""Build the two feature universes and write them to results/feature_definitions.json.

MODULES
-------
Source: DRAM's KEGG module step form
  /ibex/scratch/projects/c2014/rob/metagenomics/databases/dram/module_step_form.20260623.tsv
One row per (module, KO, path). The `path` column's FIRST component is the top-level
STEP index; rows sharing a step are ALTERNATIVE enzymes for that step (M00001 step 0 has
six alternative hexokinases). 399 modules, 276 of which have more KOs than steps.

Because of those alternatives, the completeness rule is applied over STEPS, not over the
raw KO list: a module is present in a genome when at least 2/3 of its steps are satisfied,
a step being satisfied if ANY of its alternative KOs is annotated. Taken literally over the
raw KO list the 2/3 rule would be unsatisfiable for the alternative-rich modules (the
bacterial ribosome module is one step with 55 alternatives, so it would demand 37 of 55),
which is why the step form is used instead.

The JSON therefore stores, per module, the list of STEP GROUPS. `10_differential.py`
consumes a flat KO list plus a fraction, so the step structure is pre-resolved here into
a compact per-step representative encoding: each module is emitted as its step groups and
the presence test is done on groups. See `step_groups` in the JSON.

KOs
---
Every KO annotated in rep_ko_v2.tsv, later restricted by 10_differential.py to those
carried by 5-95 % of the 5,438 species representatives. Names come from the kofamscan
ko_list shipped with DRAM.
"""
import os
import sys
import json
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import W, RES, load_ko_sets

MODULE_STEP_FORM = W + "databases/dram/module_step_form.20260623.tsv"
KO_LIST = W + "databases/dram/kofam_ko_list.tsv"


def main():
    ms = pd.read_csv(MODULE_STEP_FORM, sep="\t", usecols=["ko", "module", "module_name", "path"])
    ms["step"] = ms["path"].astype(str).str.split(",").str[0]

    module_steps, module_names = {}, {}
    for mid, grp in ms.groupby("module"):
        steps = [sorted(set(g["ko"].dropna())) for _, g in grp.groupby("step")]
        steps = [s for s in steps if s]
        if not steps:
            continue
        module_steps[mid] = steps
        module_names[mid] = str(grp["module_name"].iloc[0])

    ko_sets = load_ko_sets()
    all_kos = sorted(set().union(*ko_sets.values()))

    kl = pd.read_csv(KO_LIST, sep="\t", usecols=["knum", "definition"])
    ko_names = dict(zip(kl["knum"], kl["definition"].astype(str)))

    spec = {
        "module": {
            "step_groups": module_steps,
            "definitions": {m: sorted(set(k for s in st for k in s))
                            for m, st in module_steps.items()},
            "names": module_names,
            "source": MODULE_STEP_FORM,
            "rule": "present when >= 2/3 of top-level steps are satisfied; "
                    "a step is satisfied if any of its alternative KOs is annotated",
        },
        "ko": {
            "definitions": {k: [k] for k in all_kos},
            "names": {k: ko_names.get(k, "") for k in all_kos},
            "source": "rep_ko_v2.tsv annotations; names from " + KO_LIST,
            "rule": "present when the KO is annotated; universe restricted downstream "
                    "to KOs carried by 5-95 % of species representatives",
        },
    }

    os.makedirs(RES, exist_ok=True)
    out = os.path.join(RES, "feature_definitions.json")
    with open(out, "w") as fh:
        json.dump(spec, fh)

    nsteps = [len(v) for v in module_steps.values()]
    nkos = [len(v) for v in spec["module"]["definitions"].values()]
    print(f"modules: {len(module_steps)}  steps/module median {pd.Series(nsteps).median():.0f} "
          f"(max {max(nsteps)})  KOs/module median {pd.Series(nkos).median():.0f}")
    print(f"KOs annotated across {len(ko_sets)} genomes: {len(all_kos)}  "
          f"named: {sum(1 for k in all_kos if ko_names.get(k))}")
    print("wrote", out)


if __name__ == "__main__":
    main()
