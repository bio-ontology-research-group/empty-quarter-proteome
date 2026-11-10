#!/usr/bin/env bash
# Regenerate the four reframed main manuscript figures and per-panel TSVs from
# compact bundled inputs, run the figure/panel verification tests, and build
# the manuscript PDF.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

python3 "$HERE/plot_proteome_reframed.py" \
  --out "$REPO/manuscript/figures" \
  --panel-dir "$HERE/panel_data" \
  --comp-res "$HERE/inputs/compute_results" \
  --repro-dm "$REPO/reproducibility/darkmatter" \
  --arch-res "$REPO/reproducibility/integration_20261108/results"

cp -f "$REPO/manuscript/figures/"*_reframed.* "$REPO/figures/"

pytest -q "$HERE/test_proteome_reframed_figures.py"
if [[ -d "/home/leechuck/Documents/papers/empty-quarter-data-current" && -d "/home/leechuck/Documents/papers/empty-quarter-ecology-current" ]]; then
  pytest -q "$HERE/test_three_paper_consistency.py"
fi
bash "$REPO/reproducibility/build_manuscript.sh"
