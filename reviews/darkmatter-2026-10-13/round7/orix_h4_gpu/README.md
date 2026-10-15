# H4 ORIX GPU comparison — 2026-10-14

Requested by Robert. GPU Foldseek searches move to ORIX; CPU-only preparation stays on IBEX workers. Existing CPU results remain frozen.

- IBEX preparation job: **53760480**, debug, 8 CPUs, 32 GB. Guarded worker execution. Downloads release 10-941cd33 GPU binary, verifies version and every frozen PDB component, verifies 5,851 whole-model and 7,773 crop hashes, builds query databases and a padded target. Archives release URL and SHA256, reference hashes and query manifest.
- ORIX array: **315233**, tasks 0 (whole) and 1 (crop), freecycle-h200, QoS freecycle, one explicitly typed nvidia_h200 GPU and 16 CPUs per task, 96 GB, 12 hours, requeue enabled. Submitted held pending transfer; no GPU allocated while held.
- Local relay: `stage_and_release.py`, process record `relay_process.json`; live status `transfer_status.json`, log `stage_and_release.log`. Waits for successful preparation, rsyncs self-contained files through disk cache `/home/leechuck/.cache/darkmatter-h4-orix-20261014/payload`, then releases the array. Errors leave it held. Full checksums are verified on ORIX workers before searching. Cache retained until verification/collection.
- Remote results: `/mnt/data/u/hohndor/eq-current-structures/h4-round7-gpu-20261014/results/{wc,crop}`. Partial attempts have separate restart directories; completed results atomically promoted with output hashes. GPU utilization recorded.

## Interpretation

Same Foldseek commit 941cd33ff0771cd2e3f144e3293e22a2b87e9fda, coordinate inputs and frozen PDB release as CPU run. Alignment type 1, E<=10, max-seqs 1000, identical output columns. GPU uses ungapped prefilter mode 1, which differs from CPU k-mer prefilter; `-s` is not applied in GPU mode. Results are a GPU-retrieval sensitivity analysis and must be compared with CPU results before promotion. A GPU/host change does not fix the known alnTM normalization defect. The strict qTM/ttm checks and an explicit alnTM anomaly disposition still apply. Adaptive cap checks and downstream five-tier classifications remain required before review.

## Validation so far

Python compilation and shell syntax checks passed. ORIX sbatch test-only accepted the GPU configuration, then held submission returned 315233. Worker guard previously tested on actual IBEX login and compute hosts. Input preparation and GPU output validity are not yet complete. Flash read-only review attempt timed out without findings (empty response despite receipt status SUCCESS); no completed Flash review is claimed. Root inspected the scripts; scientific validation and Argon review remain pending.

Sources: https://github.com/steineggerlab/foldseek#faster-search-with-gpu-acceleration and release https://github.com/steineggerlab/foldseek/releases/tag/10-941cd33. ORIX live Slurm partition settings verified; website was inaccessible through browser tool.
