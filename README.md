# Rub’ al-Khali proteome: reproducibility package

This paper asks how unknown protein families are distributed among desert microbial lineages and shared with other environments. The manuscript leads with the disproportionate unknown-family contribution of widespread Thermoplasmatota, their subsurface distribution, cross-soil sequence sharing and heliorhodopsin diversity.

## Reproduce the published analysis snapshot

Requirements: Python 3 with `pip install -r requirements.txt`, curl, latexmk and a TeX distribution including natbib and authblk.

```bash
git clone https://github.com/bio-ontology-research-group/empty-quarter-proteome.git
cd empty-quarter-proteome
bash reproducibility/fetch_large_files.sh
bash reproducibility/run_all.sh
```

Downloads use public GitHub release assets and verify their SHA-256 checksums. No IBEX credentials or other repository checkout is needed for these commands. The harness regenerates figures, verifies frozen numerical evidence and builds `manuscript/main.pdf`. It reports 62 COMPUTE checks, 288 FROZEN checks, six TABLE sources and 25 historical provenance entries. TABLE and historical entries have different evidential roles from locally recomputed results. The historical registry and supporting analyses are retained for traceability; the active dark-proteome checks are also available through `bash reproducibility/darkmatter/run_current.sh`.

## Contents and scope

- `manuscript/`: LaTeX, full author list, bibliography and figures.
- `reproducibility/`: entry points, checksum manifests, claim registry, retained analysis scripts and compact data.
- `code/`, `results/`, `RQ/`: analysis sources, derived evidence and historical supporting analyses.
- `reviews/`: review reports, responses and verification evidence.

Raw assembly, binning, large sequence searches and model inference are documented for IBEX/ORIX. The frozen-output harness is distinct from rerunning those upstream jobs. Missing raw-file paths, permissions, run-merge provenance, ENA reviewer access, MAG deposition and archival DOI remain documented prerequisites in the manuscript. Public derived inputs do not resolve those separate issues.

## Manuscript and companion paper

Manuscript editing: https://www.overleaf.com/project/6a40af16a879c0f96d4c1761, local `~/Documents/papers/empty-quarter-proteome/`. This package is checked out at `~/Public/software/empty-quarter-proteome/`.

The Life on Air paper has its own independent package: https://github.com/bio-ontology-research-group/empty-quarter-lifeonair. See `MIGRATION.md` for source provenance.
