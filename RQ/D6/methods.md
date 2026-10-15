# RQ-D6 — Methods

**Question.** How prevalent is rhodopsin-based light-energy harvesting across the 15,469 Empty
Quarter soil MAGs, and does it track a light niche (sun-exposed Surface vs shielded Deep)?

## Annotation source
Per-MAG functional annotations are the **existing EggNOG-mapper outputs** (`.emapper.annotations`,
one per MAG) generated during the main pipeline — **no new search was run**. File paths are derived
from `06_mobilome/mag_index.tsv` (col1 = MAG id, col2 = Bakta `.faa`); the sibling EggNOG file is the
same path with `/bakta/`→`/eggnog/` and `.faa`→`.emapper.annotations`. All 15,469 MAGs had an EggNOG
file (0 missing). Relevant columns: 8 = Description, 9 = Preferred_name, 12 = KEGG_ko, 21 = PFAMs.

## Census (IBEX, Slurm batch — not login node)
`scripts/ibex/D6_rho_census.sh`, submitted with `sbatch --account=pi-hohndor --partition=batch`
(job **47843530**; 4 CPU, 8 GB, ~6 min wall). For each MAG it scans the EggNOG table and flags a
protein as a **microbial type-1 rhodopsin (retinylidene opsin)** if:
- PFAMs (col 21) contains a real rhodopsin Pfam — **Bac_rhodopsin (PF01036)** or Heliorhodopsin; **or**
- KEGG_ko (col 12) ∈ {K04641, K04642, K04643, K04644}; **or**
- Description / Preferred_name contains "rhodopsin",

**excluding** ASRT / "…rhodopsin **transducer**" proteins — the Anabaena Sensory Rhodopsin Transducer
is a soluble cytoplasmic signalling protein (Pfam `ASRT`), **not an opsin**, and only matches on the
word "rhodopsin". Output: `mag · protein · ko · pname · pfams · description`
(`06_metabolism/rhodopsin/rhodopsin_hits.tsv`). 398 candidate proteins in 341 MAGs were pulled.

## Local QC + analysis (`run.py`)
1. **Drop eukaryotic false positives.** 7 hits had Pfam `7tm_1` / `7tm_1,Orexin_rec2` / `adh_short`
   ("7-transmembrane receptor (rhodopsin family)" — animal/plant **GPCR** opsin-family, plus one
   retinol dehydrogenase RDH12). All came from **4 PlantRoot MAGs** and are eukaryotic contamination,
   not microbial opsins. Retained set = proteins with Pfam **Bac_rhodopsin** (or Heliorhodopsin):
   **391 proteins in 337 MAGs**.
2. **Prevalence** = distinct rhodopsin-bearing MAGs / 15,469.
3. **Sub-typing** is assigned **only** from an informative KEGG ortholog or an explicit
   proteo/xantho/helio name. EggNOG's generic description "Bacteriorhodopsin-like protein" /
   "COG5524 Bacteriorhodopsin" is the descriptor for the **entire** type-1 superfamily and does not
   resolve subtype, so such proteins are labelled "type-1 rhodopsin (subtype unresolved)".
4. **Light-niche test.** MAG→sample via the id prefix before `__` (e.g. `10Dr2__SemiBin_0`→`10Dr2`);
   sample→compartment from `results/darkmatter_ecology/env_matrix.tsv`. The denominator is **all**
   MAGs in each compartment (rhodopsin is rare, so per-compartment totals are required, not just
   rhodopsin-bearing MAGs). Surface vs Deep enrichment tested with a **Fisher exact test**.
   Benchmarked against the light-driven photolyase pattern from
   `results/uv_survival/compartment.tsv` (Surface 22.4% vs Deep 13.1%).

## Inputs / outputs
- In: `rhodopsin_hits.tsv`, `all_mags.txt` (15,469 ids, both pulled from IBEX),
  `results/darkmatter_ecology/env_matrix.tsv`.
- Out: `rhodopsin_census.tsv`, `fig_rhodopsin.{pdf,png}`.

## Caveats
- **Homology-based annotation only** (EggNOG/Pfam); no HMM-confirmed retinal-binding-lysine check or
  phylogenetic placement, so subtype calls (proton-pump proteorhodopsin vs sensory) are coarse and
  some divergent soil opsins may be missed. A targeted Pfam/HMMER + retinal-residue validation is the
  natural follow-up. Prevalence is robust (Bac_rhodopsin PF01036 is a specific, well-curated family).
