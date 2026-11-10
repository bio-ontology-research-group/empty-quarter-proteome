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
pending, sources, figures = ['main.tex'], set(), set()
while pending:
    relative = pending.pop()
    if relative in sources:
        continue
    path = source/relative
    if not path.resolve().is_relative_to(source.resolve()) or not path.is_file():
        raise SystemExit('Missing or external manuscript source: ' + relative)
    sources.add(relative)
    text = path.read_text()
    figures.update(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', text))
    for included in re.findall(r'\\(?:input|include)\{([^}]+)\}', text):
        pending.append(included if included.endswith('.tex') else included + '.tex')
for relative in figures:
    if not (source/relative).resolve().is_relative_to(source.resolve()) or not (source/relative).is_file():
        raise SystemExit('Missing or external figure: ' + relative)
for relative in sorted(sources | figures | {'authors.tex', 'references.bib'}):
    destination = target/relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source/relative, destination)
print('Manuscript source, authors, bibliography and referenced figures synchronised.')
PYCODE
