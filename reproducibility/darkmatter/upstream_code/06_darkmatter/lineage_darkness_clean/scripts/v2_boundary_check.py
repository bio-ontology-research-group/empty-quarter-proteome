#!/usr/bin/env python
"""ROUND TWO, GOAL A audit: no cluster-table line was lost at a worker boundary.

v2_scan_cluster90.py splits all_cluster90.tsv into 8 byte ranges. Worker k stops
once it has passed its end offset; worker k+1 seeks to that offset and throws away
one readline() because it assumes the byte lands mid-line. If an offset happens to
land exactly on a line start, that line belongs to nobody. This checks the 7
internal boundaries and reports any line that fell through.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import v2_scan_cluster90 as v2

NW = int(sys.argv[1]) if len(sys.argv) > 1 else 8

size = os.path.getsize(v2.F_CLUSTER)
step = size // NW
lost = []
with open(v2.F_CLUSTER, "rb") as fh:
    for k in range(1, NW):
        off = k * step
        fh.seek(off - 1)
        prev = fh.read(1)
        if prev == b"\n":                 # the offset is a line start: nobody read it
            line = fh.readline()
            lost.append((k, off, line.rstrip(b"\n").decode()))

print("boundaries checked: %d" % (NW - 1))
print("lines that fell through: %d" % len(lost))
if lost:
    win = v2.load_winners()
    c2w = v2.build_contig_map(win)
    for k, off, line in lost:
        rep, member = line.split("\t")
        contig = member[:member.rfind("_")].encode()
        i = c2w.get(contig)
        print("  boundary %d @ %d  rep=%s member=%s  winner=%s"
              % (k, off, rep, member, win[i] if i is not None else "no"))
else:
    print("every boundary offset lands mid-line, so all %d bytes were covered exactly once" % size)
