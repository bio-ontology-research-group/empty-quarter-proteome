#!/usr/bin/env bash
# Reproduce the current compact evidence package. This does not rerun upstream searches.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
python3 reproducibility/darkmatter/verify_review_inputs.py
python3 reproducibility/darkmatter/check_current_tables.py
python3 reproducibility/darkmatter/round1_sensitivity.py
python3 reproducibility/darkmatter/round2_checks.py
python3 reproducibility/darkmatter/match_selection.py
python3 reproducibility/darkmatter/selection_estimands.py
python3 reproducibility/darkmatter/rebuild_environment.py
python3 reproducibility/darkmatter/overlap_sensitivity.py
python3 reproducibility/darkmatter/produce_culture_diagnostics.py
python3 reproducibility/darkmatter/current_revision_checks.py
python3 reproducibility/darkmatter/round4_checks.py
python3 reproducibility/darkmatter/round5_extensions/checks.py
python3 reproducibility/darkmatter/make_review_figures.py
python3 reproducibility/darkmatter/check_round7_registered.py
python3 reproducibility/darkmatter/round7_figures/make_structure.py
bash reproducibility/build_manuscript.sh
