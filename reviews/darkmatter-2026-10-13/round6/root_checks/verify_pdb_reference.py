import pathlib, hashlib, json, time
root=pathlib.Path('/ibex/scratch/projects/c2014/rob/metagenomics/rebuild_20260729/structural_calibration/references/pdb_20260730')
manifest=root/'foldseek_pdb.sha256'
rows=[]
for line in manifest.read_text().splitlines():
    expected, original=line.split(maxsplit=1)
    p=root/pathlib.Path(original.lstrip('*')).name
    if p == manifest:
        continue  # Historical checksum list includes itself; pin its current digest separately.
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    rows.append(dict(path=str(p),bytes=p.stat().st_size,expected=expected,sha256=h.hexdigest(),match=h.hexdigest()==expected))
result=dict(manifest_path=str(manifest),manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),checked_unix=time.time(),all_match=all(r['match'] for r in rows),files=rows)
pathlib.Path('pdb_reference_verification.json').write_text(json.dumps(result,indent=2)+'\n')
assert result['all_match']
print(len(rows),'reference components verified')
