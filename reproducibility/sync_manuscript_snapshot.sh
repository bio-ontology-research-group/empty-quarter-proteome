#!/usr/bin/env bash
# Explicit refresh from the dedicated manuscript checkout; existing source edits should be reconciled first.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE="${1:-$HOME/Documents/papers/empty-quarter-proteome}"
for file in main.tex authors.tex references.bib; do
  test -s "$SOURCE/$file"
done
python3 - "$SOURCE" "$ROOT/manuscript" <<'PYCODE'
import re, shutil, sys
from pathlib import Path
source, target = map(Path, sys.argv[1:])
text = (source/'main.tex').read_text()
figures = re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', text)
for relative in figures:
    if not (source/relative).is_file():
        raise SystemExit('Missing figure: ' + relative)
for relative in ['main.tex', 'authors.tex', 'references.bib'] + figures:
    destination = target/relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source/relative, destination)
print('Manuscript source, authors, bibliography and referenced figures synchronised.')
PYCODE
