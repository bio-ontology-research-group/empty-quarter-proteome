#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
r=Path('reviews/darkmatter-2026-10-13/round7')
roots=[Path('resources/literature/darkmatter-2026-10-13'),Path('manuscript'),r,Path('reviews/darkmatter-2026-10-13/round6/snapshots/h1_final'),Path('reproducibility/darkmatter/round7_final'),Path('reproducibility/darkmatter/round7_analyses'),Path('reproducibility/darkmatter/round7_figures')]
files=set()
for root in roots:
 for p in root.rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts and p.suffix not in {'.aux','.fls','.fdb_latexmk','.blg'} and 'ibex_mirror' not in p.name:files.add(str(p))
for root in ['reproducibility','reproducibility/darkmatter']:
 for p in Path(root).glob('*'):
  if p.is_file():files.add(str(p))
files.update(['CLAUDE.md','notes/CORRECTIONS.md',str(r/'manuscript/verify_ibex_mirror.sbatch'),str(r/'manuscript/prepare_ibex_mirror.py')])
lines=[]
for f in sorted(files):
 h=hashlib.sha256()
 with open(f,'rb') as inp:
  for b in iter(lambda:inp.read(1024*1024),b''):h.update(b)
 lines.append(h.hexdigest()+'  '+f)
p=r/'manuscript/ibex_mirror.sha256';p.write_text('\n'.join(lines)+'\n')
files.add(str(p))
(r/'manuscript/ibex_mirror_files.txt').write_text('\n'.join(sorted(files))+'\n')
print(len(lines),'files pinned for IBEX review mirror')
