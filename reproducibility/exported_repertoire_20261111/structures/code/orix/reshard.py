#!/usr/bin/env python3
"""Login-node helper (file handling only): regroup not-yet-claimed shards by alignment depth.
Alignments with <= 1500 rows go to SET_shallow (run with a reduced MSA padding, --max-msa 512:1024 for monomers
and 508:1024 for complexes: nothing is subsampled, only the zero padding of the MSA tensors shrinks);
deeper alignments go to SET_deep and run with the ColabFold defaults (512:5120 / 508:2048).
usage: reshard.py SET ITEMS_PER_SHARD [skip_shard ...]"""
import sys, json, shutil
from pathlib import Path
D = Path("/mnt/data/u/hohndor/eq-exported-dark-20261111")
SET, per = sys.argv[1], int(sys.argv[2]); skip = set(sys.argv[3:])
items = {"shallow": [], "deep": []}
for sd in sorted((D / "inputs" / SET).glob("shard_*")):
    if sd.name in skip or (D / "locks" / f"{SET}_{sd.name}").exists(): continue
    for f in sorted(sd.glob("*.a3m")):
        txt = f.read_text(); lines = txt.splitlines()
        rows = sum(1 for l in lines if l.startswith(">"))
        L = sum(int(x) for x in lines[0][1:].split("\t")[0].split(",")) if lines[0].startswith("#") else len(lines[1])
        items["shallow" if rows <= 1500 else "deep"].append((L, f.name, f, rows))
    (D / "locks" / f"{SET}_{sd.name}").mkdir(parents=True)          # retire the old shard
    (D / "locks" / f"{SET}_{sd.name}" / "owner").write_text("resharded\n")
out = {}
for cls, its in items.items():
    its.sort(key=lambda x: -x[0])
    n = max(1, -(-len(its) // per)); load = [0.0] * n; sh = [[] for _ in range(n)]
    for it in its:
        k = load.index(min(load)); load[k] += it[0] ** 2 + 40000; sh[k].append(it)
    for k in range(n):
        d = D / "inputs" / f"{SET}_{cls}" / f"shard_{k:02d}"; d.mkdir(parents=True, exist_ok=True)
        for L, name, f, rows in sh[k]: shutil.copy2(f, d / name)
    out[cls] = {"n": len(its), "shards": [len(x) for x in sh], "max_rows": max([x[3] for x in its], default=0)}
(D / "inputs" / f"reshard_{SET}.json").write_text(json.dumps(out, indent=1))
print(SET, json.dumps(out))
