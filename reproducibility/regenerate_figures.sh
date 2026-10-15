#!/usr/bin/env bash
# Regenerate derived result prose and every manuscript figure from curated tables
# under ../results/ plus checksum-pinned compact inputs under this directory;
# separately rebuild the bounded RQ/E2 products. Run from anywhere.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
cd "$REPO"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1785801600}"
export TZ="${TZ:-UTC}"

echo "== regenerating derived result prose and figures =="
rc=0
run() { echo "  -> python3 code/$1"; python3 "code/$1" >/dev/null 2>"$HERE/.figerr" || { echo "     FAILED: $1"; sed 's/^/       /' "$HERE/.figerr" | tail -4; rc=1; }; }
run make_figures.py
run phase_d_endemism.py
run phase_d_global.py
run make_darkmatter_figures.py
run make_global_specificity_figure.py
run make_picrust2_bias_figure.py
echo "  -> python3 RQ/E2/run.py"
python3 RQ/E2/run.py >/dev/null 2>"$HERE/.figerr" || {
  echo "     FAILED: RQ/E2/run.py"
  sed 's/^/       /' "$HERE/.figerr" | tail -4
  rc=1
}
rm -f "$HERE/.figerr"

echo "== checking expected figure PDFs exist =="
EXPECT=(
  fig_amr_class fig_bgc_product fig_lifeonair_by_phylum
  fig_darkmatter_dnds fig_darkmatter_pcoa fig_darkmatter_permanova
  fig_darkmatter_endemism fig_darkmatter_global_endemism fig_darkmatter_endemism_sites
  fig_picrust2_bias fig_darkmatter_global_specificity
)
miss=0
for f in "${EXPECT[@]}"; do
  if [[ -s "figures/$f.pdf" ]]; then echo "  OK   $f.pdf"; else echo "  MISS $f.pdf"; miss=1; fi
done
[[ $miss -eq 0 ]] && echo "All ${#EXPECT[@]} figures present." || { echo "Some figures missing."; rc=1; }

echo "== checking B4 exploratory-screen products exist =="
E2_EXPECT=(
  RQ/E2/fig_posel_function.pdf
  RQ/E2/table_positive_families_summary.tsv
  RQ/E2/table_structural_tiers.tsv
  RQ/E2/table_functional_categories.tsv
  RQ/E2/table_boltz_hypotheses.tsv
)
for path in "${E2_EXPECT[@]}"; do
  if [[ -s "$path" ]]; then echo "  OK   $path"; else echo "  MISS $path"; miss=1; rc=1; fi
done

if [[ $miss -eq 0 ]]; then
  echo "== synchronising generated PDFs into manuscript/figures/ =="
  mkdir -p manuscript/figures
  for f in "${EXPECT[@]}"; do
    cp -f "figures/$f.pdf" "manuscript/figures/$f.pdf"
    if cmp -s "figures/$f.pdf" "manuscript/figures/$f.pdf"; then
      echo "  SYNC $f.pdf"
    else
      echo "  FAIL $f.pdf"
      rc=1
    fi
  done
  if [[ "${EQ_SYNC_OVERLEAF:-0}" == "1" && -d paper/figures ]]; then
    echo "== synchronising generated PDFs into optional Overleaf submodule =="
    for f in "${EXPECT[@]}"; do
      cp -f "figures/$f.pdf" "paper/figures/$f.pdf"
    done
  fi
fi
python3 "$HERE/darkmatter/round7_figures/make_structure.py" || rc=1
python3 "$HERE/darkmatter/family_followup_20261015/plot_family.py" || rc=1
exit $rc
