#!/usr/bin/env python3
"""Derived BioSamples for the bins (ENA requires one per binned metagenome and one per MAG), by XML.

  python3 16_register_bin_samples.py --tier binned --service test --dryrun          # write XML only
  python3 16_register_bin_samples.py --tier binned --service test --only 10Dr2__self_SemiBin_190 --submit
  python3 16_register_bin_samples.py --tier binned --service prod --submit          # IRREVERSIBLE (her word)
  python3 16_register_bin_samples.py --tier mag    --service prod --submit          # after the binned tier

Tier binned: checklist ERC000050 (ENA binned metagenome), one sample per row of wgs_bins.tsv, alias = bin name,
  "sample derived from" = the soil BioSample, description "This sample represents a metagenomic bin from the
  metagenomic sample ERSxxx. GTDB release 232: <lineage>." (ENA's prescribed phrase + the taxon).
Tier mag: checklist ERC000047 (GSC MIMAGS), one sample per representative (is_representative = 1), alias =
  <bin>_MAG, "sample derived from" = the bin's binned sample accession (so the binned tier must be registered
  first), description "This sample represents a MAG derived from the metagenomic sample ERSxxx through binned
  sample ERSyyy. GTDB release 232: <lineage>."
Field values (decisions: notes/ENA_DEPOSIT.md, Step 3): constants below; collection date, coordinates, ENVO
contexts, country, depth and elevation copied from the soil BioSample record (fetched once, cached in
ena_wgs/bins/cache/soil/); quality, rRNA/tRNA, taxon from wgs_bins.tsv. Samples go in SAMPLE_SETs of --batch
(default 200) per submission with HOLD = the study hold; accessions are folded into wgs_bins.tsv
(binned_sample_accession / mag_sample_accession). Resume: aliases OK in the tier's log are skipped.
"""
from __future__ import annotations

import argparse
import csv
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from ena_common import HERE, WORK, fetch_record, password, print_receipt, submission_xml, submit

TABLE = HERE / "wgs_bins.tsv"
OUT = WORK / "bins"
HOLD = "2026-11-25"
CHECKLIST = {"binned": "ERC000050", "mag": "ERC000047"}
PROJECT_NAME = "Rub' al-Khali Desert Microbiome"
METAGENOMIC_SOURCE = "soil metagenome"
ISOLATION_SOURCE = "Rub' al-Khali desert soil"
BINNING_SOFTWARE = "SemiBin2 2.2.1; DAS Tool 1.1.7"
BINNING_PARAMETERS = ("SemiBin2 single_easy_bin on coverage and k-mer composition, pretrained global model "
                      "(--environment global) and self-supervised (--self-supervised -m 1000) runs on each assembly; "
                      "DAS Tool refinement, score threshold 0.5")
ASSEMBLY_SOFTWARE = "MEGAHIT;1.2.9;--presets meta-large --min-contig-len 500 --min-count 2 --prune-level 1 --low-local-ratio 0.3"
COMPLETENESS_SOFTWARE = "CheckM2 1.1.0"
TAXONOMIC_MARKER = "multi-marker approach (GTDB-Tk 2.7.2, 120 bacterial / 53 archaeal marker genes, GTDB release 232)"
RRNA_SOFTWARE = "barrnap 0.9"
TRNA_SOFTWARE = "tRNAscan-SE 2.0.12"
COPY_FROM_SOIL = ("collection date", "geographic location (latitude)", "geographic location (longitude)",
                  "geographic location (country and/or sea)", "broad-scale environmental context",
                  "local environmental context", "environmental medium")


def tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def soil_attrs(acc: str, pw: str) -> dict[str, tuple[str, str | None]]:
    """tag -> (value, units) of the soil BioSample, from a cached copy of the live record."""
    cache = OUT / "cache" / "soil" / f"{acc}.xml"
    if cache.is_file():
        rec = ET.parse(cache).getroot()
    else:
        rec = fetch_record("samples", acc, pw)
        if rec is None:
            sys.exit(f"ERROR: could not fetch soil sample {acc}")
        cache.parent.mkdir(parents=True, exist_ok=True)
        ET.ElementTree(rec).write(cache, encoding="UTF-8", xml_declaration=True)
    sample = rec.find(".//SAMPLE") if rec.tag != "SAMPLE" else rec
    return {a.findtext("TAG"): (a.findtext("VALUE") or "", a.findtext("UNITS")) for a in sample.iter("SAMPLE_ATTRIBUTE")}


def add(attrs: ET.Element, tag: str, value: str, units: str | None = None) -> None:
    if value == "":
        sys.exit(f"ERROR: empty value for {tag!r}")
    el = ET.SubElement(attrs, "SAMPLE_ATTRIBUTE")
    ET.SubElement(el, "TAG").text = tag
    ET.SubElement(el, "VALUE").text = value
    if units:
        ET.SubElement(el, "UNITS").text = units


def build(row: dict, tier: str, soil: dict[str, tuple[str, str | None]]) -> ET.Element:
    b = row["bin"]
    alias = b if tier == "binned" else f"{b}_MAG"
    if not row["ncbi_taxid"]:
        sys.exit(f"ERROR: {b} has no NCBI taxon (run 15_resolve_bin_taxa.py)")
    s = ET.Element("SAMPLE", {"alias": alias, "center_name": "KAUST"})
    ET.SubElement(s, "TITLE").text = alias
    name = ET.SubElement(s, "SAMPLE_NAME")
    ET.SubElement(name, "TAXON_ID").text = row["ncbi_taxid"]
    ET.SubElement(name, "SCIENTIFIC_NAME").text = row["ncbi_name"]
    lineage = f" GTDB release 232: {row['gtdb_lineage']}." if row["gtdb_lineage"] else ""
    if tier == "binned":
        derived = row["biosample_accession"]
        desc = f"This sample represents a metagenomic bin from the metagenomic sample {derived}.{lineage}"
    else:
        derived = row["binned_sample_accession"]
        if not derived:
            sys.exit(f"ERROR: {b} has no binned sample accession yet (register the binned tier first)")
        desc = (f"This sample represents a MAG derived from the metagenomic sample {row['biosample_accession']} "
                f"through binned sample {derived}.{lineage}")
    ET.SubElement(s, "DESCRIPTION").text = desc
    attrs = ET.SubElement(s, "SAMPLE_ATTRIBUTES")
    add(attrs, "project name", PROJECT_NAME)
    add(attrs, "metagenomic source", METAGENOMIC_SOURCE)
    add(attrs, "sample derived from", derived)
    add(attrs, "isolation_source", ISOLATION_SOURCE)
    for tag in COPY_FROM_SOIL:
        if tag in soil and soil[tag][0]:
            add(attrs, tag, soil[tag][0], soil[tag][1])
        else:
            sys.exit(f"ERROR: soil sample {row['biosample_accession']} lacks {tag!r}")
    add(attrs, "binning software", BINNING_SOFTWARE)
    add(attrs, "binning parameters", BINNING_PARAMETERS)
    add(attrs, "assembly software", ASSEMBLY_SOFTWARE)
    add(attrs, "completeness score", row["completeness"], "%")
    add(attrs, "completeness software", COMPLETENESS_SOFTWARE)
    add(attrs, "contamination score", row["contamination"], "%")
    add(attrs, "taxonomic identity marker", TAXONOMIC_MARKER)
    if row["assembly_quality"]:
        add(attrs, "assembly quality", row["assembly_quality"])
        add(attrs, "16s recovered", "Yes" if row["has_16S"] == "1" else "No")
        add(attrs, "16S recovery software", RRNA_SOFTWARE)
        add(attrs, "number of standard tRNAs extracted", row["n_trna_total"])
        add(attrs, "tRNA extraction software", TRNA_SOFTWARE)
    elif tier == "mag":
        sys.exit(f"ERROR: {b}: assembly quality is mandatory for a MAG sample (rRNA/tRNA not available)")
    add(attrs, "ENA-CHECKLIST", CHECKLIST[tier])
    return s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", required=True, choices=["binned", "mag"])
    ap.add_argument("--service", required=True, choices=["test", "prod"])
    ap.add_argument("--batch", type=int, default=200)
    ap.add_argument("--only", help="comma-separated bin names (e.g. a test-service rehearsal set)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--dryrun", action="store_true", help="write the XML, send nothing")
    g.add_argument("--submit", action="store_true", help="POST to the chosen service")
    a = ap.parse_args()

    pw = password()
    rows = tsv(TABLE)
    cols = list(rows[0].keys())
    acc_col = "binned_sample_accession" if a.tier == "binned" else "mag_sample_accession"
    # MAG tier = species representatives with contamination < 10 (checklist: "If contamination >= 10% then please
    # submit as a binned assembly"); the one representative at exactly 10 stays binned only
    only = set(a.only.split(",")) if a.only else None
    todo = [r for r in rows if (a.tier == "binned" or (r["is_representative"] == "1" and float(r["contamination"]) < 10))
            and not r[acc_col] and (only is None or r["bin"] in only)]
    log = OUT / "logs" / f"register_{a.tier}.{a.service}.tsv"
    log.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if log.is_file():
        done = {l.split("\t")[0] for l in log.read_text().splitlines()[1:] if "\tOK\t" in l}
    else:
        log.write_text("alias\tstatus\taccession\tdetail\n")
    todo = [r for r in todo if (r["bin"] if a.tier == "binned" else f"{r['bin']}_MAG") not in done]
    print(f"{a.tier}: {len(todo)} samples to register on {a.service} (batch {a.batch})")

    soil_cache: dict[str, dict] = {}
    xml_dir = OUT / "xml" / a.tier
    xml_dir.mkdir(parents=True, exist_ok=True)
    sub_xml = OUT / "xml" / f"submission.{a.tier}.xml"
    sub_xml.write_text(submission_xml("ADD", HOLD))
    by_alias = {}
    for i in range(0, len(todo), a.batch):
        batch = todo[i:i + a.batch]
        sset = ET.Element("SAMPLE_SET")
        for r in batch:
            acc = r["biosample_accession"]
            if acc not in soil_cache:
                soil_cache[acc] = soil_attrs(acc, pw)
            s = build(r, a.tier, soil_cache[acc])
            sset.append(s)
            by_alias[s.get("alias")] = r
        ET.indent(sset)
        path = xml_dir / f"batch_{i // a.batch:04d}.xml"
        ET.ElementTree(sset).write(path, encoding="UTF-8", xml_declaration=True)
        if a.dryrun:
            print(f"[dryrun] {path} ({len(batch)} samples)")
            continue
        ok, info = submit({"SUBMISSION": sub_xml, "SAMPLE": path}, pw, a.service,
                          OUT / "receipts" / a.tier / f"{a.service}.{path.stem}.receipt.xml")
        print_receipt(path.stem, ok, info)
        with open(log, "a") as fh:
            got = {alias: acc for tag, acc, alias in info["accessions"] if tag == "SAMPLE" and acc}
            for r in batch:
                alias = r["bin"] if a.tier == "binned" else f"{r['bin']}_MAG"
                if ok and alias in got:
                    r[acc_col] = got[alias]
                    fh.write(f"{alias}\tOK\t{got[alias]}\t{a.service}\n")
                else:
                    fh.write(f"{alias}\tFAILED\t\t{'; '.join(t for _, t in info['messages'])[:300]}\n")
        if a.service == "prod":
            with open(TABLE, "w") as fh:
                w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
                w.writeheader()
                w.writerows(rows)
    if a.dryrun:
        print(f"XML in {xml_dir}; first sample:\n" + ET.tostring(sset[0], encoding="unicode")[:1500] if todo else "nothing to do")


if __name__ == "__main__":
    main()
