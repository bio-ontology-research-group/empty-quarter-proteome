#!/usr/bin/env python3
"""Print the RQ-A3 confidence/length sensitivity summary.

The large structures and FoldSeek outputs remain on IBEX. Inputs used to
produce the summary are:

* ``best_afdb_tm.tsv`` and ``best_pdb_tm.tsv`` in this directory; and
* ``$M/06_esmfold_novelfold/pdb/plddt_shard*.tsv`` on IBEX.

Dual-novel means best TM < 0.5 against both AFDB and the curated PDB; a missing
hit counts as novel. ``finalize_tm.py`` recomputes the summary from those
inputs. This script only prints the checked key/value result.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = os.path.join(HERE, "fold_validation_summary.tsv")


def load(path):
    data = {}
    with open(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2 and not row[0].startswith("#"):
                data[row[0]] = row[1]
    return data


def main():
    data = load(SUMMARY)
    print(f"pLDDT shard records                         : {int(data['n_plddt_records']):,}")
    print(f"  scored models (status=ok)                 : {int(data['n_plddt_ok']):,}")
    print(f"  excluded as too_long                      : {int(data['n_plddt_too_long']):,}")
    print("--- dual-novel sensitivity: best TM<0.5 vs AFDB AND PDB ---")
    for threshold in ("0.00", "0.40", "0.50", "0.60", "0.70", "0.80"):
        stem = f"dual_novel_plddt_ge_{threshold}"
        count = int(data[f"{stem}_n"])
        denominator = int(data[f"{stem}_denominator"])
        rate = float(data[f"{stem}_pct"])
        print(f"  pLDDT>={threshold}: {count:6,} / {denominator:6,} ({rate:4.1f}%)")
    print("--- overall medians: dual-novel / structural match ---")
    print(f"  pLDDT                                    : "
          f"{data['dual_novel_median_plddt']} / {data['matched_median_plddt']}")
    print(f"  length (aa)                              : "
          f"{data['dual_novel_median_length_aa']} / {data['matched_median_length_aa']}")
    print("--- pLDDT>=0.60 novelty rate by length ---")
    for label, stem in (
        ("<100", "lt100"),
        ("100-150", "100_150"),
        ("150-250", "150_250"),
        ("250-400", "250_400"),
        (">=400", "ge400"),
    ):
        n_models = int(data[f"plddt_ge_0.60_length_{stem}_n"])
        rate = float(data[f"plddt_ge_0.60_length_{stem}_dual_novel_pct"])
        print(f"  {label:>7} aa: n={n_models:5,}; dual-novel={rate:4.1f}%")
    print("  dual-novel with length>=100 aa            : "
          f"{int(data['plddt_ge_0.60_length_ge100_dual_novel_n']):,}")
    print("Matched known-fold false-novelty control    : "
          f"{data['known_fold_matched_false_novelty_control']}")


if __name__ == "__main__":
    main()
