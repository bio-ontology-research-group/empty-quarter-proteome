#!/usr/bin/env bash
# Build the manuscript with its BibTeX bibliography and fail on unresolved citations/references.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
MANUSCRIPT="$ROOT/manuscript"
PAPER="$ROOT/paper"

if [[ ! -s "$MANUSCRIPT/main.tex" || ! -s "$MANUSCRIPT/references.bib" ]]; then
  echo "public manuscript snapshot is missing" >&2
  exit 2
fi
# The public snapshot is authoritative; an optional private submodule may lag it.
# Set EQ_CHECK_PRIVATE_COPY=1 when explicitly validating synchronisation.
if [[ "${EQ_CHECK_PRIVATE_COPY:-0}" == "1" && -s "$PAPER/main.tex" ]]; then
  cmp -s "$PAPER/main.tex" "$MANUSCRIPT/main.tex" || {
    echo "paper/main.tex differs from manuscript/main.tex; reconcile the private copy with authoritative manuscript/main.tex" >&2
    exit 2
  }
  cmp -s "$PAPER/references.bib" "$MANUSCRIPT/references.bib" || {
    echo "paper/references.bib differs from manuscript/references.bib; reconcile with authoritative manuscript/references.bib" >&2
    exit 2
  }
fi
cd "$MANUSCRIPT"

latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex >/dev/null
test -s main.pdf

if grep -Eq \
    'Citation .* undefined|Reference .* undefined|There were undefined references|There were undefined citations' \
    main.log; then
    echo "manuscript contains unresolved citations or references" >&2
    exit 1
fi

pages="$(pdfinfo main.pdf | awk '/^Pages:/ {print $2}')"
echo "MANUSCRIPT: OK (${pages} pages; BibTeX resolved; public snapshot)"
