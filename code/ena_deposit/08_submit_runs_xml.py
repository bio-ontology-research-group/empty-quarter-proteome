#!/usr/bin/env python3
"""Register experiments + runs by XML (drop-box) for reads already uploaded to the root of the Webin upload area.

Why: Webin-CLI ignores pre-uploaded files and re-uploads each pair into webin-cli/reads/<run>/ (2026-09-06:
6.6 TB would have been transferred twice). The XML route references the uploaded files by name + md5.
Same metadata as the manifests (wgs_runs.tsv); aliases follow Webin-CLI's convention (webin-reads-<run_alias>)
so the two routes cannot create the same run twice. Before anything is sent, the account's live experiment and
run lists are fetched: runs already on ENA are skipped; an experiment without a run gets the run only.

  python 08_submit_runs_xml.py --service test|prod --mode dryrun|submit [--only RUN_ALIAS] [--limit N]

Outputs under WORK: xml/runs/<run_alias>.{experiment,run}.xml, receipts/runs/<run_alias>.<service>.receipt.xml,
logs/runs_xml.<service>.tsv (alias, status, ERX, ERR, detail; OK aliases skipped on rerun).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import xml.etree.ElementTree as ET

from ena_common import HERE, WORK, curl, password, print_receipt, submission_xml, submit

STUDY = "PRJEB125153"
HOLD = "2026-11-01"
CENTER = "KAUST"
RUNS = HERE / "wgs_runs.tsv"
MD5 = WORK / "logs" / "md5sums.tsv"
REPORT = "https://www.ebi.ac.uk/ena/submit/report/{kind}?format=json&max-results=10000"


def live(kind: str, pw: str) -> dict[str, dict]:
    """alias -> report for the account's experiments or runs under the study."""
    rc, body = curl([REPORT.format(kind=kind)], pw)
    if rc != 0 or not body.strip().startswith("["):
        sys.exit(f"ERROR: could not list {kind} (report API)")
    return {r["report"]["alias"]: r["report"] for r in json.loads(body) if r["report"].get("studyId") == "ERP204894"}


def experiment_xml(r: dict) -> ET.Element:
    alias = f"webin-reads-{r['run_alias']}"
    e = ET.Element("EXPERIMENT", {"alias": alias, "center_name": CENTER})
    ET.SubElement(e, "TITLE").text = r["run_alias"]
    ET.SubElement(e, "STUDY_REF", {"accession": STUDY})
    d = ET.SubElement(e, "DESIGN")
    desc = f"Shotgun metagenome of soil sample {r['sample_id']}"
    if r["sample_id"] == "48Dr2":
        desc += "; fastp-processed reads (raw R1 file lost)"
    ET.SubElement(d, "DESIGN_DESCRIPTION").text = desc
    ET.SubElement(d, "SAMPLE_DESCRIPTOR", {"accession": r["biosample_accession"]})
    ld = ET.SubElement(d, "LIBRARY_DESCRIPTOR")
    ET.SubElement(ld, "LIBRARY_NAME").text = r["run_alias"]
    ET.SubElement(ld, "LIBRARY_STRATEGY").text = "WGS"
    ET.SubElement(ld, "LIBRARY_SOURCE").text = "METAGENOMIC"
    ET.SubElement(ld, "LIBRARY_SELECTION").text = "RANDOM"
    layout = ET.SubElement(ld, "LIBRARY_LAYOUT")
    paired = {"NOMINAL_LENGTH": str(round(float(r["insert_size"])))} if r.get("insert_size") else {}
    ET.SubElement(layout, "PAIRED", paired)
    ET.SubElement(ld, "LIBRARY_CONSTRUCTION_PROTOCOL").text = r["library_kit"]
    p = ET.SubElement(e, "PLATFORM")
    ET.SubElement(ET.SubElement(p, "ILLUMINA"), "INSTRUMENT_MODEL").text = r["instrument_model"]
    return e


def run_xml(r: dict, md5: dict[str, str], exp_ref: dict[str, str]) -> ET.Element:
    ru = ET.Element("RUN", {"alias": f"webin-reads-{r['run_alias']}", "center_name": CENTER})
    ET.SubElement(ru, "TITLE").text = r["run_alias"]
    ET.SubElement(ru, "EXPERIMENT_REF", exp_ref)
    files = ET.SubElement(ET.SubElement(ru, "DATA_BLOCK"), "FILES")
    for name in (r["upload_r1"], r["upload_r2"]):
        if name not in md5:
            sys.exit(f"ERROR: no md5 for {name} in {MD5}")
        ET.SubElement(files, "FILE", {"filename": name, "filetype": "fastq", "checksum_method": "MD5", "checksum": md5[name]})
    return ru


def write_set(tag: str, el: ET.Element, path) -> None:
    s = ET.Element(tag)
    s.append(el)
    ET.indent(s)
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(s).write(path, encoding="UTF-8", xml_declaration=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--service", required=True, choices=["test", "prod"])
    ap.add_argument("--mode", required=True, choices=["dryrun", "submit"])
    ap.add_argument("--only")
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    pw = password()

    with open(RUNS) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(MD5) as fh:
        md5 = {r["upload_name"]: r["md5"] for r in csv.DictReader(fh, delimiter="\t")}
    exps, runs = live("experiments", pw), live("runs", pw)
    print(f"on ENA under {STUDY}: {len(exps)} experiments, {len(runs)} runs")

    log = WORK / "logs" / f"runs_xml.{args.service}.tsv"
    done: set[str] = set()
    if log.exists():
        done = {l.split("\t")[0] for l in open(log) if "\tOK\t" in l}
    else:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text("alias\tstatus\terx\terr\tdetail\n")

    remaining = [r for r in rows if f"webin-reads-{r['run_alias']}" not in runs]
    print(f"runs without a run record on ENA: {len(remaining)} of {len(rows)}")
    targets = remaining
    if args.only:
        targets = [r for r in targets if r["run_alias"] == args.only]
    if args.limit:
        targets = targets[: args.limit]
    print(f"this invocation handles {len(targets)} run(s)")

    sxml = WORK / "xml" / "runs.submission.xml"
    sxml.write_text(submission_xml("ADD", hold_until=HOLD))
    n_ok = n_fail = 0
    for r in targets:
        alias = f"webin-reads-{r['run_alias']}"
        if r["run_alias"] in done and args.mode == "submit":
            continue
        parts = {"SUBMISSION": sxml}
        xdir = WORK / "xml" / "runs"
        if alias in exps:                       # experiment exists (Webin-CLI created it, run failed): run only
            exp_ref = {"accession": exps[alias]["id"]}
        else:
            exp_ref = {"refname": alias}
            write_set("EXPERIMENT_SET", experiment_xml(r), xdir / f"{r['run_alias']}.experiment.xml")
            parts["EXPERIMENT"] = xdir / f"{r['run_alias']}.experiment.xml"
        write_set("RUN_SET", run_xml(r, md5, exp_ref), xdir / f"{r['run_alias']}.run.xml")
        parts["RUN"] = xdir / f"{r['run_alias']}.run.xml"
        if args.mode == "dryrun":
            print(f"[dry] {r['run_alias']}: {'run only (experiment ' + exps[alias]['id'] + ')' if alias in exps else 'experiment + run'}")
            continue
        ok, info = submit(parts, pw, args.service, WORK / "receipts" / "runs" / f"{r['run_alias']}.{args.service}.receipt.xml")
        erx = next((a for t, a, _ in info["accessions"] if t == "EXPERIMENT"), exps.get(alias, {}).get("id", ""))
        err = next((a for t, a, _ in info["accessions"] if t == "RUN"), "")
        with open(log, "a") as fh:
            fh.write(f"{r['run_alias']}\t{'OK' if ok else 'FAILED'}\t{erx}\t{err}\t{' | '.join(t for _, t in info['messages'])[:200]}\n")
        print_receipt(r["run_alias"], ok, info)
        n_ok += ok
        n_fail += not ok
    if args.mode == "dryrun":
        print(f"DRY RUN - XML in {WORK / 'xml' / 'runs'}, nothing sent.")
    else:
        print(f"{args.service}: {n_ok} OK, {n_fail} failed; log {log}")


if __name__ == "__main__":
    main()
