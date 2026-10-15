#!/usr/bin/env python3
"""Validate active H1-H4 frozen claims without invoking historical manuscript guards.

This is compact-table extraction and arithmetic, not a rerun of the Slurm analyses.
Raw audits and their execution receipts are retained separately in the review package.
"""
import csv, hashlib, importlib.util, math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('paper_verifier',ROOT/'reproducibility/verify.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
claims=[c for c in v.CLAIMS if c[1] in {'H1','H2','H3','H4'}]
assert claims, 'No active extension claims registered'
assert len({c[0] for c in claims})==len(claims), 'Duplicate claim identifiers'
manifest={r['artifact']:r for r in csv.DictReader(open(ROOT/'reproducibility/current_manuscript_artifacts.tsv'),delimiter='\t')}
checked=set();out=[]
for rel,row in manifest.items():
    if rel.startswith('darkmatter/round7_final/'):
        assert hashlib.sha256((ROOT/'reproducibility'/rel).read_bytes()).hexdigest()==row['sha256'],rel
        checked.add(rel)
for cid,section,description,expected,tier,source,fn,tol in claims:
    assert tier=='FROZEN' and fn is not None,(cid,tier)
    p=ROOT/source;rel=str(p.relative_to(ROOT/'reproducibility'))
    assert rel in manifest,(cid,'unpinned source',rel)
    if rel not in checked:
        assert hashlib.sha256(p.read_bytes()).hexdigest()==manifest[rel]['sha256'],rel
        checked.add(rel)
    actual=float(fn());assert math.isfinite(actual),(cid,actual)
    assert abs(actual-float(expected))<=tol,(cid,actual,expected,tol)
    out.append({'claim_id':cid,'section':section,'description':description,'paper_value':expected,'extracted_value':actual,'tolerance':tol,'tier':tier,'source':source,'validation':'PASS: pinned compact extraction; upstream execution separately audited'})
with open(ROOT/'reproducibility/darkmatter/round7_registered_claims.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
print(f'PASS: {len(out)} H1-H4 FROZEN claims from {len(checked)} checksum-pinned source files. No upstream rerun claimed.')
