#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

cd "${REPO_ROOT}"

echo "[1/5] Verifying shared_foundation frozen outputs..."
python3 reproducibility/mechanism_20261109/shared_foundation/verify_shared_foundation.py

echo "[2/5] Verifying exp1_metabolic frozen outputs and unit tests..."
python3 reproducibility/mechanism_20261109/exp1_metabolic/verify_frozen_outputs.py
python3 -m pytest reproducibility/mechanism_20261109/exp1_metabolic/interpretation/test_descriptive_interpretation.py reproducibility/mechanism_20261109/exp1_metabolic/covariates/test_covariate_metadata_supplement.py -q

echo "[3/5] Verifying exp3_viruses frozen outputs and unit tests..."
python3 reproducibility/mechanism_20261109/exp3_viruses/verify_frozen_outputs.py
python3 -m pytest reproducibility/mechanism_20261109/exp3_viruses/tests/test_exp3_viruses.py -q

echo "[4/5] Verifying exp4_surface frozen outputs and unit tests..."
python3 reproducibility/mechanism_20261109/exp4_surface/verify_frozen_outputs.py
python3 -m pytest reproducibility/mechanism_20261109/exp4_surface/tests/test_exp4_surface.py -q

echo "[5/5] Verifying independent expression eligibility..."
python3 reproducibility/expression_20261110/verify_expression_20261110.py

echo "ALL MECHANISTIC VALIDATION CHECKS PASSED."
