#!/usr/bin/env python3
"""Inventory gene-call completeness across the species representatives.

These are retained representative proteomes. Their identifiers and gene calls
are not assumed interchangeable with the assembly catalogue's original calls.
"""
import gzip
import hashlib
import json
import os
import platform
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "results"
BULK = HERE / "bulk"


def fasta(path):
    header, seq = None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq).rstrip("*")
                header, seq = line[1:].rstrip(), []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq).rstrip("*")


def worker(job):
    gid, filename, high = job
    p = Path(filename)
    counters = dict(genome=gid, proteins=0, complete=0, length_ge100=0,
                    complete_ge100=0, internal_stops=0, residues=0)
    selected = []
    for header, seq in fasta(p):
        complete = re.search(r"partial=(\d\d)", header).group(1) == "00"
        counters["proteins"] += 1
        counters["complete"] += complete
        counters["length_ge100"] += len(seq) >= 100
        counters["complete_ge100"] += complete and len(seq) >= 100 and "*" not in seq
        counters["internal_stops"] += "*" in seq
        counters["residues"] += len(seq)
        if high and complete and len(seq) >= 100 and "*" not in seq:
            selected.append(f">{header.split()[0]}\n{seq}\n")
    counters["source_bytes"] = p.stat().st_size
    with open(p, "rb") as f:
        digest = hashlib.file_digest(f, "sha256").hexdigest()
    counters["source_sha256"] = digest
    counters["source_path"] = str(p)
    return counters, "".join(selected)


if __name__ == "__main__":
    assert os.environ.get("SLURM_JOB_ID"), "Run on an allocated Slurm worker"
    OUT.mkdir(exist_ok=True)
    BULK.mkdir(exist_ok=True)
    idx_path = W / "06_darkmatter/heat_stress_20261016/results/genome_index.tsv"
    tiers_path = W / "mimag_rrna_trna_clean/results/mimag_tiers.tsv"
    idx = pd.read_csv(idx_path, sep="\t").set_index("genome_id")
    tiers = pd.read_csv(tiers_path, sep="\t").set_index("rep")
    assert len(idx) == len(tiers) == 5438
    assert set(idx.index) == set(tiers.index)
    high = set(tiers.index[tiers.mimag_tier == "high"])
    assert len(high) == 126
    jobs = [(gid, row.proteome, gid in high) for gid, row in idx.sort_index().iterrows()]
    rows = []
    with open(BULK / "high_quality_complete_ge100.faa", "w") as fasta_out:
        with ProcessPoolExecutor(max_workers=1) as pool:
            for i, (row, seq) in enumerate(pool.map(worker, jobs)):
                rows.append(row)
                fasta_out.write(seq)
                if i % 500 == 0:
                    print("genomes", i, flush=True)
    table = pd.DataFrame(rows).set_index("genome")
    table = table.join(idx[["gtdb_classification", "completeness", "contamination"]])
    table = table.join(tiers[["mimag_tier"]])
    table.to_csv(OUT / "representative_gene_inventory.tsv", sep="\t")
    report = dict(job=os.environ["SLURM_JOB_ID"], host=platform.node(),
                  genomes=len(table),
                  totals={c: int(table[c].sum()) for c in
                          ["proteins", "complete", "length_ge100", "complete_ge100",
                           "internal_stops", "residues"]},
                  high_quality_complete_ge100=int(table.loc[list(high), "complete_ge100"].sum()),
                  scope="Retained representative calls; no darkness assignment inferred from gene identifiers")
    (OUT / "inventory_run.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
