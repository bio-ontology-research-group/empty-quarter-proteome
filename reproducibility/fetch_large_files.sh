#!/usr/bin/env bash
# Fetch only this paper's bulk inputs from IBEX, then verify their original checksums.
set -euo pipefail
IBEX_USER="${IBEX_USER:-}"
IBEX_HOST="${IBEX_HOST:-ilogin.ibex.kaust.edu.sa}"
SOURCE="/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/mags-paper-large-files"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -z "$IBEX_USER" ]]; then
  echo 'Set IBEX_USER to your KAUST account; IBEX project c2014 access is required.' >&2
  exit 2
fi
LIST="$(mktemp)"
trap 'rm -f "$LIST"' EXIT
awk '{print $2}' "$ROOT/reproducibility/large_files.sha256" > "$LIST"
rsync -av --files-from="$LIST" "${IBEX_USER}@${IBEX_HOST}:${SOURCE}/" "$ROOT/"
(cd "$ROOT" && sha256sum -c reproducibility/large_files.sha256)
echo 'All institutional inputs verified.'
