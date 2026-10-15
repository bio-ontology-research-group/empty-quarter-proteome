# Curated MAG-paper result artifacts

This directory contains compact inputs and outputs used by the public reproducibility harness. Files come from
different, explicitly labelled analysis universes; directory membership alone does not make a statistic an active
manuscript claim. The active values and extraction tiers are registered in `../reproducibility/claims.tsv` and
checked by `../reproducibility/verify.py`.

## Catalogue scopes

- Current taxonomy and recurrent-family analyses use the frozen 14,426-refined-bin / 5,438-ANI95-representative
  catalogue and GTDB-Tk 2.7.2/R232 artifacts.
- Several functional tables are a historical 15,469 single-sample-MAG snapshot generated in June 2026 from
  Bakta/EggNOG, AMRFinderPlus and antiSMASH workbench outputs. They remain usable only where the manuscript and
  claim ledger explicitly label that denominator.
- The 753,305,533-protein assembly catalogue is a separate discovery arm and is not a denominator for the tables
  listed here.

The detailed universe crosswalk is in `../docs/binning-and-annotation-scope.md`.

## Principal compact tables

| File or directory | Content and boundary |
|---|---|
| `amr_class_summary.tsv`, `amr_long.tsv` | historical single-sample-MAG intrinsic-resistome census; gene presence is genomic capacity, not phenotype |
| `bgc_product_summary.tsv`, `bgc_per_mag.tsv` | historical antiSMASH `--minimal` region/product census; not a conservative BGC-novelty estimate |
| `metabolism_marker_long.tsv`, `metabolism_*_by_phylum.tsv` | historical EggNOG marker census; atmospheric forms are adjudicated separately under `atmospheric/` |
| `darkmatter_ecology/` | recurrent sequence-dark-family ecology and bounded historical diagnostics; read each file's status/caveat |
| `global_specificity/` | no-hit results relative to named catalogues; never evidence of universal absence or endemism |
| `picrust2_bias/` | historical PICRUSt2/MAG comparison with denominator sensitivity retained explicitly |
| `smetana/` | 16-sample pilot, not biome-scale network inference |

## Selected registered historical functional counts

The harness currently checks, among other values, atmospheric-hydrogenase capacity in 3,386 MAGs, form-I CODH
capacity in 3,828 MAGs, a 1,204-MAG dual-marker set, and antiSMASH product counts such as terpene regions in
9,984 MAGs. These are gene/region-presence summaries in the labelled historical catalogue; they do not directly
measure expression, flux, ecological activity, resistance phenotype or chemical production.

To reproduce all promoted figures, claims, manifests and the public manuscript:

```bash
bash reproducibility/run_all.sh
```
