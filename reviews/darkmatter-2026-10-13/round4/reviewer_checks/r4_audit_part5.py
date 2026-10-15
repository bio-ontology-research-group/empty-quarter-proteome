from pathlib import Path

ROOT = Path('/home/leechuck/Documents/papers/empty-quarter-mags')
DM = ROOT / 'reproducibility/darkmatter'
SNAP = DM / 'source_snapshot'

print('=== 1. JSON SUMMARIES UNDER reproducibility/darkmatter ===')
for p in sorted(DM.rglob('*.json')):
  rel = p.relative_to(DM)
  if (
      'environments/' in str(rel)
      or 'sha256' in str(rel)
      or 'gff_checks' in str(rel)
      or 'source_snapshot' in str(rel)
  ):
    continue
  print(f'\n--- {rel} ---')
  print(p.read_text()[:4000])

print('\n=== 2. METABOLISM / BGC / AMR SUMMARIES IN source_snapshot ===')
for p in sorted((SNAP / '06_metabolism').rglob('*summary*')) + sorted(
    (SNAP / '06_metabolism').rglob('*report*')
):
  rel = p.relative_to(SNAP)
  if p.is_file() and p.stat().st_size < 25000:
    print(f'\n--- {rel} ---')
    print(p.read_text()[:4000])
