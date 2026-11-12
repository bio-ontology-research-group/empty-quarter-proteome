#!/usr/bin/env python3
"""Login-node helper (file handling only): unpack an a3m tar and split it into length-balanced shards.
usage: prep_shards.py TAR OUTDIR NSHARDS MAXLEN"""
import sys, tarfile, os, json
from pathlib import Path
tar, out, n, maxlen = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
items = []
with tarfile.open(tar) as tf:
    for m in tf.getmembers():
        if not m.isfile(): continue
        txt = tf.extractfile(m).read().decode()
        lines = txt.splitlines()
        L = sum(int(x) for x in lines[0][1:].split("\t")[0].split(",")) if lines[0].startswith("#") else len(lines[1])
        items.append((L, os.path.basename(m.name), txt))
items.sort(reverse=True)
load = [0.0] * n
shard = [[] for _ in range(n)]
long_ = []
for L, name, txt in items:
    if L > maxlen:
        long_.append((L, name, txt)); continue
    k = load.index(min(load)); load[k] += L ** 2 + 40000; shard[k].append((L, name, txt))
for k in range(n):
    d = out / f"shard_{k:02d}"; d.mkdir(parents=True, exist_ok=True)
    for L, name, txt in shard[k]: (d / name).write_text(txt)
d = out / "too_long"; d.mkdir(parents=True, exist_ok=True)
for L, name, txt in long_: (d / name).write_text(txt)
man = {"tar": str(tar), "n_items": len(items), "n_shards": n, "max_len": maxlen,
       "per_shard": [len(s) for s in shard], "too_long": [(L, nm) for L, nm, _ in long_]}
(out / "shards.json").write_text(json.dumps(man, indent=1))
print(json.dumps(man))
