#!/usr/bin/env bash
# Master reproducibility driver: regenerate every figure, recompute & verify every paper number,
# refresh the claims manifest, build the manuscript, and write REPORT.md. Exit non-zero if any
# gate, COMPUTE/FROZEN claim, TABLE source, figure check, or manuscript build fails.
#
#   bash reproducibility/run_all.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPORT="$HERE/REPORT.md"

echo "# Reproducibility report" > "$REPORT"
echo >> "$REPORT"

echo "===================  1/4  figures  ==================="
bash "$HERE/regenerate_figures.sh" | tee "$HERE/.fig.log"
figrc=${PIPESTATUS[0]}
{ echo "## Figures"; echo '```'; cat "$HERE/.fig.log"; echo '```'; echo; } >> "$REPORT"

echo; echo "===================  2/4  numbers  ==================="
python3 "$HERE/verify.py" | tee "$HERE/.verify.log"
verrc=${PIPESTATUS[0]}
{ echo "## Numbers (recompute vs paper)"; echo '```'; cat "$HERE/.verify.log"; echo '```'; echo; } >> "$REPORT"

if [[ -f "$HERE/tests/test_manuscript_scoping.py" ]]; then
  echo "== running manuscript scoping & include tests =="
  python3 -m pytest -q "$HERE/tests/test_manuscript_scoping.py" || verrc=1
fi
if [[ -f "$HERE/reframed_20261108/test_proteome_reframed_figures.py" ]]; then
  echo "== running reframed figure verification tests =="
  python3 -m pytest -q "$HERE/reframed_20261108/test_proteome_reframed_figures.py" || verrc=1
fi
if [[ -f "$HERE/integration_20261110/surface_figure/test_surface_figure.py" ]]; then
  echo "== running surface modules figure verification tests =="
  python3 -m pytest -q "$HERE/integration_20261110/surface_figure/test_surface_figure.py" || verrc=1
fi

echo; echo "===================  3/4  manifest  ==================="
python3 "$HERE/verify.py" --manifest > "$HERE/claims.tsv"
echo "wrote $(wc -l < "$HERE/claims.tsv") lines to claims.tsv"

echo; echo "===================  4/4  manuscript  ==================="
bash "$HERE/build_manuscript.sh" | tee "$HERE/.build.log"
buildrc=${PIPESTATUS[0]}
{ echo "## Manuscript build"; echo '```'; cat "$HERE/.build.log"; echo '```'; echo; } >> "$REPORT"

rm -f "$HERE/.fig.log" "$HERE/.verify.log" "$HERE/.build.log"
echo; echo "REPORT.md written. figures rc=$figrc, verify rc=$verrc, manuscript rc=$buildrc"
[[ $figrc -eq 0 && $verrc -eq 0 && $buildrc -eq 0 ]] && echo "REPRODUCIBILITY: OK" || echo "REPRODUCIBILITY: ISSUES (see above)"
exit $(( figrc || verrc || buildrc ))
