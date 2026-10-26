#!/usr/bin/env python3
"""NCBI taxon (id + name) for every bin in wgs_bins.tsv from its GTDB r232 lineage, checked against ENA's
taxonomy service (only "submittable" taxa can be used for a sample). Fills ncbi_taxid / ncbi_name; idempotent.

  srun -c 2 --mem 16G -t 02:00:00 python3 15_resolve_bin_taxa.py        # metadata parsing: not on a login node

Resolution order per lineage (the most specific NCBI name ENA accepts for an organism-level environmental sample,
ENA taxonomy FAQ: "Genus sp." / "<Rank> bacterium|archaeon" are the organism-level environmental names):
  1. GTDB species that is a real name (not sp######): "<Genus> <epithet>" via ENA scientific-name lookup; if ENA
     does not know it, the majority NCBI species taxid of the GTDB species cluster (GTDB r232 bac120/ar53 metadata)
  2. genus: "<Genus> sp."                                 3. family / order / class / phylum: "<Rank> bacterium|archaeon"
  4. "uncultured bacterium" (77133) / "uncultured archaeon" (115547)
GTDB polyphyly suffixes (_A, _B ...) are stripped before lookup. Placeholder ranks (UBA..., JA..., numeric, with
"-") are skipped. Lookups are cached in ena_wgs/bins/taxa_cache.tsv (name, taxid, submittable, rank).
Output summary: bins per resolution level. No email address is ever sent.
"""
from __future__ import annotations

import csv
import gzip
import json
import re
import sys
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from ena_common import HERE, WORK

TABLE = HERE / "wgs_bins.tsv"
CACHE = WORK / "bins" / "taxa_cache.tsv"
META = {"d__Bacteria": WORK / "bins" / "reference" / "bac120_metadata_r232.tsv.gz",
        "d__Archaea": WORK / "bins" / "reference" / "ar53_metadata_r232.tsv.gz"}
REST = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"
FALLBACK = {"d__Bacteria": ("77133", "uncultured bacterium"), "d__Archaea": ("115547", "uncultured archaeon")}
PLACEHOLDER = re.compile(r"^(UBA|JA|GCA|GCF|RBG|CG|SW-|DTU|2-|Fen-|UBA|QH|SCGC|sp\d)|[-\d]")


def clean(name: str) -> str:
    """Strip GTDB suffixes: 'Bacillus_A' -> 'Bacillus', 'Bacillus_A flava' -> 'Bacillus flava'."""
    return re.sub(r"_[A-Z]+(?=\s|$)", "", name)


def is_placeholder(rank_name: str) -> bool:
    return not rank_name or bool(PLACEHOLDER.search(rank_name))


class Ena:
    def __init__(self) -> None:
        self.cache: dict[str, tuple[str, str, str]] = {}
        if CACHE.is_file():
            for r in csv.DictReader(open(CACHE), delimiter="\t"):
                self.cache[r["name"]] = (r["taxid"], r["submittable"], r["rank"])
        self.fh = open(CACHE, "a")
        if CACHE.stat().st_size == 0:
            self.fh.write("name\ttaxid\tsubmittable\trank\n")

    def lookup(self, name: str) -> str | None:
        """taxid if ENA knows the exact scientific name and it is submittable, else None."""
        if name not in self.cache:
            url = REST + urllib.parse.quote(name)
            try:
                with urllib.request.urlopen(url, timeout=60) as resp:
                    body = resp.read().decode()
                hits = json.loads(body) if body.strip() else []
            except Exception:
                hits = []
            hit = next((h for h in hits if h.get("scientificName") == name), None)
            self.cache[name] = (hit["taxId"], hit["submittable"], hit.get("rank", "")) if hit else ("", "false", "")
            self.fh.write("\t".join((name,) + self.cache[name]) + "\n")
            self.fh.flush()
        taxid, submittable, _ = self.cache[name]
        return taxid if taxid and submittable == "true" else None


def species_from_metadata(domain: str) -> dict[str, tuple[str, str]]:
    """GTDB species -> majority (ncbi_species_taxid, ncbi_organism_name) over the cluster's genomes, species-level only."""
    votes: dict[str, Counter] = defaultdict(Counter)
    names: dict[str, str] = {}
    with gzip.open(META[domain], "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            sp = r["gtdb_taxonomy"].split(";")[-1]
            tid, org = r["ncbi_species_taxid"], r["ncbi_organism_name"]
            if tid and tid != "none" and " sp." not in org and not re.search(r"\b(bacterium|archaeon)\b", org):
                votes[sp][tid] += 1
                names.setdefault(tid, " ".join(org.split()[:2]))
    return {sp: (c.most_common(1)[0][0], names[c.most_common(1)[0][0]]) for sp, c in votes.items()}


def main() -> None:
    rows = list(csv.DictReader(open(TABLE), delimiter="\t"))
    cols = list(rows[0].keys())
    ena = Ena()
    meta: dict[str, dict] = {}
    level = Counter()
    resolved: dict[str, tuple[str, str, str]] = {}
    for r in rows:
        lin = r["gtdb_lineage"]
        if not lin:
            level["no lineage yet"] += 1
            continue
        if lin in resolved:
            r["ncbi_taxid"], r["ncbi_name"], lv = resolved[lin]
            level[lv] += 1
            continue
        ranks = dict(p.split("__", 1) for p in lin.split(";"))
        domain = "d__" + ranks["d"]
        word = "bacterium" if domain == "d__Bacteria" else "archaeon"
        taxid = name = None
        lv = ""
        sp = ranks.get("s", "")
        if sp and not re.search(r"\bsp\d", sp):
            cand = clean(sp)
            taxid = ena.lookup(cand)
            if taxid:
                name, lv = cand, "species (ENA name)"
            else:
                if domain not in meta:
                    meta[domain] = species_from_metadata(domain)
                hit = meta[domain].get("s__" + sp)
                if hit and ena.lookup(hit[1]) == hit[0]:
                    taxid, name, lv = hit[0], hit[1], "species (GTDB metadata majority)"
        if not taxid and not is_placeholder(ranks.get("g", "")):
            cand = f"{clean(ranks['g'])} sp."
            taxid = ena.lookup(cand)
            if taxid:
                name, lv = cand, "genus sp."
        if not taxid:
            for rk in ("f", "o", "c", "p"):
                if not is_placeholder(ranks.get(rk, "")):
                    cand = f"{clean(ranks[rk])} {word}"
                    taxid = ena.lookup(cand)
                    if taxid:
                        name, lv = cand, f"{rk}__ {word}"
                        break
        if not taxid:
            taxid, name = FALLBACK[domain]
            lv = "uncultured fallback"
        resolved[lin] = (taxid, name, lv)
        r["ncbi_taxid"], r["ncbi_name"], _ = resolved[lin]
        level[lv] += 1
    with open(TABLE, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} bins, {len(resolved)} distinct lineages, {len(ena.cache)} cached ENA lookups")
    for k, v in level.most_common():
        print(f"  {v:6d}  {k}")


if __name__ == "__main__":
    main()
