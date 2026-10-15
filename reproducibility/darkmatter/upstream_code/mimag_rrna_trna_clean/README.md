# mimag_rrna_trna_clean: rRNA + tRNA detection for MIMAG quality tiers on the 5,438 species reps

Owner: velazqal.

## STATUS (2026-09-24): COMPLETE and REGISTERED (CATALOGUE_NUMBERS cohorts table; R1 + Methods sentences in main.tex)

Array 52468548 (109 tasks) all COMPLETED, ~28 min per task; `02_summarise.py` run on the login node.
Full MIMAG tiers on the 5,438 reps: **high 126**, medium 5,312, low 0. Of the 2,259 completeness/contamination
HQ reps, 2,133 drop to medium. All three rRNAs in 157; 16S in 495; 23S in 333; 5S in 2,537; ≥18 tRNA isotypes in
3,629 (median 19). Paper keeps completeness/contamination tiers as operational and discloses the MIMAG count.
Outputs: `results/mimag_tiers.tsv`, `results/rrna_trna_per_rep.tsv`; per-genome barrnap GFF and tRNAscan tables in
`results/rrna/`, `results/trna/`.

**Purpose.** The 5,438-representative catalogue is currently tiered as HQ/MQ/LQ on completeness
and contamination alone (2,259 HQ + 3,179 MQ, §1.1). MIMAG's actual high-quality definition also
requires all three rRNA genes (5S, 16S, 23S) and >=18 tRNAs covering the standard amino acids.
This folder runs that rRNA/tRNA detection and reports how the tiers change under the full rule.

**Inputs (read-only, never modified):**
- Representative FASTAs: `05_refine/_drep/out/dereplicated_genomes/*.fa` (5,438 files, confirmed
  count matches the catalogue)
- CheckM2 completeness/contamination: `05_refine/_checkm2_refined/quality_report.tsv` (14,426
  refined bins; all 5,438 reps present; confirmed the `>=90/<=5` vs `>=50/<=10` split reproduces
  the catalogue's 2,259 HQ / 3,179 MQ exactly)
- GTDB domain per rep: `05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv` (4,847, Bacteria) +
  `gtdbtk.ar53.summary.tsv` (591, Archaea) = 5,438; domain read from the `classification` column's
  `d__` prefix. Every one of the 5,438 reps has a domain call, so no both-kingdoms fallback run
  was needed.

**Tools:**
- barrnap **0.9** (`module load barrnap/0.9`) — rRNA, run with `--kingdom bac` or `--kingdom arc`
  per the rep's GTDB domain above.
- tRNAscan-SE **2.0.12** (`module load trnascan-se/2.0.12`) — a module exists on IBEX, so this was
  used directly (no conda env, no aragorn fallback needed); run with `-B` (bacterial) or `-A`
  (archaeal) search mode per domain.
- ⚠ barrnap and tRNAscan-SE **cannot** be module-loaded together in one shell — tRNAscan-SE's own
  conda perl collides with barrnap's module perl (`ListUtil.c: loadable library and perl binaries
  are mismatched`). The array script runs each tool in its own subshell with a fresh
  `module purge; module load ...`; verified interactively that this does not leak between calls.
- Summariser: `/home/velazqal/miniconda3/envs/python-utils/bin/python` (3.9), stdlib only (csv, re).

**Scripts (run in order):**
1. `scripts/00_prepare_inputs.sh` — login node, seconds. Builds `results/genome_list.txt` (5,438
   rep names, sorted, fixed chunk order) and `results/domain_map.tsv` (rep -> Bacteria/Archaea).
   Already run; asserts 5,438 = 5,438 = 0 missing.
2. `scripts/01_rrna_trna_array.sbatch` — Slurm array, 109 tasks x ~50 genomes, 4 CPUs / 8G / 2h
   each, account `pi-hohndor`, partition `batch`. Idempotent: skips a genome whose GFF + tRNA
   output already both exist and are non-empty. Smoke-tested on 2 genomes (1 rep completeness
   79.18/contam 3.82 domain Bacteria; 1 rep completeness 91.42/contam 0.19 domain Bacteria) before
   submission — both tools ran correctly under the subshell fix.
3. `scripts/02_summarise.py` — **run later, after the array finishes.** Not part of the submitted
   job. Exact command:
   ```
   /home/velazqal/miniconda3/envs/python-utils/bin/python \
     /ibex/scratch/projects/c2014/rob/metagenomics/mimag_rrna_trna_clean/scripts/02_summarise.py
   ```

**Job IDs:** array job **52468548** (109 tasks, `52468548_[0-108]`), submitted 2026-09-24,
account `pi-hohndor` / partition `batch`, 4 CPUs / 8G / 2h per task. Not waited on; check with
`squeue -j 52468548` or `sacct -j 52468548`.

**Outputs:**
- `results/rrna/<rep>.gff` — barrnap GFF3 (rRNA hits), one per genome
- `results/trna/<rep>.trnascan.txt` — tRNAscan-SE tabular output, one per genome
- `results/rrna_trna_per_rep.tsv` — rep, domain, has_5S, has_16S, has_23S, n_tRNA_unique_aa,
  n_tRNA_total (from the summariser, not yet run)
- `results/mimag_tiers.tsv` — joins CheckM2 completeness/contamination with the rRNA/tRNA calls;
  `old_tier` (completeness/contamination only, matches the catalogue's 2,259/3,179 split) vs
  `mimag_tier` (full MIMAG rule: completeness >=90, contamination <=5, all three rRNAs, >=18
  tRNAs = high; completeness >=50, contamination <=10 = medium)

No genome copies are kept anywhere in this folder — only per-genome GFF/text outputs (a few KB
each) and the two small summary tables.
