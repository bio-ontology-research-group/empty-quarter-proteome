#!/usr/bin/env python3
"""Login-node helper (file handling only): move the long, not yet started monomer queries of the running
shards into separate shards (SET=monomer_tail) for additional GPUs. A .done.txt marker in the original
output directory makes the running colabfold_batch skip them (it checks the marker before each query).
usage: offload_tail.py MINLEN NSHARDS"""
import sys, json
from pathlib import Path
D = Path("/mnt/data/u/hohndor/eq-exported-dark-20261111")
minlen, n = int(sys.argv[1]), int(sys.argv[2])
items = []
for sd in sorted((D / "inputs/monomer").glob("shard_*")):
    od = D / "predictions/monomer" / sd.name
    for f in sorted(sd.glob("*.a3m")):
        L = len(f.read_text().splitlines()[1])
        name = f.stem
        if L >= minlen and not (od / f"{name}.done.txt").exists() and not list(od.glob(f"{name}_*")):
            items.append((L, name, f, od))
items.sort(reverse=True, key=lambda x: x[0])
load = [0.0] * n; sh = [[] for _ in range(n)]
for it in items:
    k = load.index(min(load)); load[k] += it[0] ** 2 + 40000; sh[k].append(it)
moved = []
for k in range(n):
    d = D / "inputs/monomer_tail" / f"shard_{k:02d}"; d.mkdir(parents=True, exist_ok=True)
    for L, name, f, od in sh[k]:
        (d / f.name).write_text(f.read_text())
        (od / f"{name}.done.txt").write_text("offloaded to monomer_tail\n")
        moved.append({"id": name, "len": L, "from": f.parent.name, "to": d.name})
(D / "inputs/monomer_tail/offload.json").write_text(json.dumps(moved, indent=1))
print(len(moved), "queries offloaded;", [len(x) for x in sh])
