#!/usr/bin/env python3
"""NCBI taxon (id + name) for every bin in wgs_bins.tsv from its GTDB r232 lineage, checked against ENA's
taxonomy service (only "submittable" taxa can be used for a sample). Fills ncbi_taxid / ncbi_name; idempotent.

  srun -c 2 --mem 16G -t 02:00:00 python3 15_resolve_bin_taxa.py        # metadata parsing: not on a login node

Resolution order per lineage (the most specific NCBI name ENA accepts for an organism-level environmental sample,
ENA taxonomy FAQ: "Genus sp." / "<Rank> bacterium|archaeon" are the organism-level environmental names). Every
candidate is used only if ENA's taxonomy service returns it as submittable:
  1. GTDB species that is a real name (not sp######): "<Genus> <epithet>"; if ENA does not know it, the majority
     NCBI species of the GTDB species cluster (GTDB r232 bac120/ar53 metadata, ncbi_species_taxid)
  2. genus: "<Genus> sp."
  3. family, then order, class, phylum, each tried two ways: (a) "<GTDB rank name> bacterium|archaeon";
     (b) the NCBI organism names that genomes of that GTDB rank carry in the r232 metadata, most frequent first,
     restricted to environmental-style names (… bacterium / … archaeon / uncultured …). (b) is what resolves
     GTDB placeholder lineages whose NCBI name differs (e.g. GTDB Patescibacteriota → "Patescibacteria group
     bacterium"); the same mapping RAK-PR's gtdb_to_ncbi_organism.py drew from the r226 metadata, extended to all ranks
  4. "uncultured bacterium" (77133) / "uncultured archaeon" (115547)
GTDB polyphyly suffixes (_A, _B ...) are stripped before lookup. Placeholder rank names (UBA..., JA..., numeric,
with "-") are skipped in (a) but still used as keys in (b). Lookups are cached in ena_wgs/bins/taxa_cache.tsv.
Columns written: ncbi_taxid, ncbi_name, ncbi_level (which step resolved it). No email address is ever sent.
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


ENV_NAME = re.compile(r"^(uncultured |Candidatus )?[A-Z][A-Za-z0-9 .-]* (bacterium|archaeon)$|^uncultured ")


def load_metadata(domain: str) -> tuple[dict[str, tuple[str, str]], dict[str, Counter]]:
    """From the GTDB r232 metadata of one domain:
    species: GTDB species -> majority (ncbi_species_taxid, name) over the cluster's genomes, species-level names only;
    rank_names: GTDB rank name (p__/c__/o__/f__ value) -> Counter of environmental-style NCBI organism names carried
    by the genomes placed in that GTDB rank (the data-driven GTDB->NCBI synonym map for placeholder ranks)."""
    votes: dict[str, Counter] = defaultdict(Counter)
    names: dict[str, str] = {}
    rank_names: dict[str, Counter] = defaultdict(Counter)
    with gzip.open(META[domain], "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ranks = r["gtdb_taxonomy"].split(";")
            sp, org, tid = ranks[-1], r["ncbi_organism_name"], r["ncbi_species_taxid"]
            if tid and tid != "none" and " sp." not in org and not re.search(r"\b(bacterium|archaeon)\b", org):
                votes[sp][tid] += 1
                names.setdefault(tid, " ".join(org.split()[:2]))
            if ENV_NAME.match(org):
                for rk in ranks[1:5]:          # p__ c__ o__ f__
                    rank_names[rk][org] += 1
    species = {sp: (c.most_common(1)[0][0], names[c.most_common(1)[0][0]]) for sp, c in votes.items()}
    return species, rank_names


def main() -> None:
    rows = list(csv.DictReader(open(TABLE), delimiter="\t"))
    cols = list(rows[0].keys())
    if "ncbi_level" not in cols:
        cols.insert(cols.index("ncbi_name") + 1, "ncbi_level")
    ena = Ena()
    meta: dict[str, tuple[dict, dict]] = {}
    level = Counter()
    resolved: dict[str, tuple[str, str, str]] = {}
    for r in rows:
        lin = r["gtdb_lineage"]
        if not lin:
            r["ncbi_taxid"] = r["ncbi_name"] = r["ncbi_level"] = ""
            level["no lineage yet"] += 1
            continue
        if lin in resolved:
            r["ncbi_taxid"], r["ncbi_name"], r["ncbi_level"] = resolved[lin]
            level[r["ncbi_level"]] += 1
            continue
        ranks = dict(p.split("__", 1) for p in lin.split(";"))
        domain = "d__" + ranks["d"]
        word = "bacterium" if domain == "d__Bacteria" else "archaeon"
        if domain not in meta:
            meta[domain] = load_metadata(domain)
        species_map, rank_names = meta[domain]
        taxid = name = None
        lv = ""
        sp = ranks.get("s", "")
        if sp and not re.search(r"\bsp\d", sp):
            cand = clean(sp)
            taxid = ena.lookup(cand)
            if taxid:
                name, lv = cand, "1 species (ENA name)"
            else:
                hit = species_map.get("s__" + sp)
                if hit and ena.lookup(hit[1]) == hit[0]:
                    taxid, name, lv = hit[0], hit[1], "1 species (GTDB metadata majority)"
        if not taxid and not is_placeholder(ranks.get("g", "")):
            cand = f"{clean(ranks['g'])} sp."
            taxid = ena.lookup(cand)
            if taxid:
                name, lv = cand, "2 genus sp."
        if not taxid:
            for rk in ("f", "o", "c", "p"):
                rname = ranks.get(rk, "")
                cands = [] if is_placeholder(rname) else [f"{clean(rname)} {word}"]
                cands += [n for n, _ in rank_names.get(f"{rk}__{rname}", Counter()).most_common(5)]
                for cand in cands:
                    taxid = ena.lookup(cand)
                    if taxid:
                        name = cand
                        lv = f"3 {rk}__ {'GTDB name' if cand == cands[0] and not is_placeholder(rname) else 'NCBI name of its genomes'}"
                        break
                if taxid:
                    break
        if not taxid:
            taxid, name = FALLBACK[domain]
            lv = "4 uncultured fallback"
        resolved[lin] = (taxid, name, lv)
        r["ncbi_taxid"], r["ncbi_name"], r["ncbi_level"] = resolved[lin]
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
