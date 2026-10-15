# GEMINI.md — Empty Quarter MAGs Project Memory

## Project Principles & Compute Rules
- Always run python code through `uv`.
- Heavy computation, bulk parsing, database extraction, and search belong strictly on guarded IBEX Slurm workers (account `pi-hohndor`, partition `batch`/`debug`) or ORIX GPU nodes; NEVER on login nodes.
- Guard scripts with `require_slurm_worker.sh` to enforce host execution within allocated Slurm jobs.
- Never delete or remove source files without granular integrity validation and explicit user confirmation.
- Never guess or invent ontology or family identifiers (e.g. Pfam, ChEBI); always validate programmatically against source databases (InterPro API, etc.).

## Dark Matter / Family Discovery (2026-10-15)
- F3 Genomic Context Test artifacts:
  - Local: `reproducibility/darkmatter/family_discovery_20261015/f3_context/`
  - Remote: `/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f3_context/`
  - Sensory HMM database: `f3_sensory_context.hmm` (17 Pfam signaling models + 2 rhodopsin controls) with official Gathering Thresholds (`--cut_ga`).
  - Key finding: Type-1 microbial rhodopsins carry dedicated MCP/HAMP transducers (30.0% vs 2.1% in controls, $p = 3.63 \times 10^{-11}$), matching halobacterial Sensory Rhodopsins I/II; Heliorhodopsins (Type-3) lack transducers (0/65) and show no two-component enrichment over background (12.3% vs 21.5%, $p = 0.24$), justifying downgrading the sensory hypothesis for EQ heliorhodopsins.
