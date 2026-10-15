#!/usr/bin/env python3
"""Soil BioSamples: `register` the missing ones (ADD), `cleanup` the registered ones (fetch-then-patch MODIFY).
Field values come from the ecology-data curated field table and per-trip sample sheets; the rules (title,
description, project name, checklist, kept/dropped fields, corrections) are recorded in notes/ENA_DEPOSIT.md.

  python 03_register_biosamples.py register|cleanup --service test|prod --mode dryrun|submit [--only ID] [--limit N] [--refresh]

Outputs under WORK: cache/samples/, xml/{register,cleanup}/, receipts/, logs/{register,cleanup}.<service>.tsv
(resume log: OK aliases are skipped), logs/cleanup_diff.tsv. MODIFY of production samples cannot be rehearsed
on the test service; run one sample first. A new cleanup rule needs --refresh and a renamed resume log.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from ena_common import HERE, WORK, curl, fetch_record, password, print_receipt, service_url, submission_xml, submit

ECO = Path("/ibex/scratch/projects/c2014/alelopezv/repos/empty-quarter-ecology-data/metadata")
ENV_TABLE = ECO / "environmental" / "environmental_measurements_curated.tsv"
ALIASES = ECO / "samples" / "site_aliases.tsv"
SHEETS = {"1": "trip1-2023.tsv", "2": "trip2-2023.tsv", "3": "trip3-2024.tsv", "4": "trip4-2024.tsv", "5": "trip5-2025.tsv"}
# Field-sheet `feature` (site x trip) -> local environmental context. The seven terms are the ones already
# used on the registered records (labels verified in OLS 2026-09-03). Named-location features (road,
# campground, water well, hot spring) have no ENVO mapping here and fall back to sibling/template.
FEATURE_ENVO = {"dune": "ENVO:00000170", "sand dune": "ENVO:00000170", "saline pan": "ENVO:00000279",
                "desert oasis": "ENVO:00000156", "gravel": "ENVO:01000018", "dune slack": "ENVO:00000308",
                "oil spill": "ENVO:00002061", "oilspill": "ENVO:00002061", "aeolian lake": "ENVO:01001105",
                "aerolian lake": "ENVO:01001105"}
BIOSAMPLE_MAP = HERE / "wgs_biosample_map.tsv"
CACHE = WORK / "cache" / "samples"

TEMPLATE_ACC = "ERS28393005"      # 10Dr2: full ERC000022 soil record, source of the constant attributes
CHECKLIST = "ERC000022"
PROJECT_NAME = "Rub' al-Khali Desert Microbiome"
DROP_NEW = {"target gene", "pcr primers", "nucleic acid amplification", "nucleic acid extraction"}
DROP_CLEAN = {"nucleic acid amplification"}
FIELD_TAGS = ("collection date", "geographic location (latitude)", "geographic location (longitude)",
              "temperature", "pressure", "humidity")
PREFIX_TRIP = {"": ("1", "Trip 1 (2023)"), "T": ("2", "Trip 2 (2023)"), "F": ("3", "Trip 3 (2024)"),
               "S": ("4", "Trip 4 (2024)"), "V": ("5", "Trip 5 (2025)")}
DEPTH_MM = {"D": "40", "PR": "40", "S": "0"}
COMPARTMENT = {"D": "deep", "PR": "root-adjacent", "S": "surface"}
ID_RE = re.compile(r"^([FSTV]?)(\d+)(D|S|PR)r(\d)$")
# Sites whose REGISTERED coordinates are kept as they are (decision 2026-09-03): the field table places
# site 52 ~19 km away in trips 1 and 3 (20.8515 N, 53.7579 E) vs trips 4-5 and every registered record
# (20.8278 N, 53.5783 E). Open question for the ecology-data owner; not corrected here.
COORD_HOLD_SITES = {"52"}


# ----------------------------------------------------------------------------- field table
def parse_coords(s: str) -> tuple[str, str]:
    m = re.match(r"^\s*([\d.]+)\s*N\s*,\s*([\d.]+)\s*E\s*$", s or "")
    if not m:
        raise ValueError(f"unparsable coordinates {s!r}")
    return f"{float(m.group(1)):.4f}", f"{float(m.group(2)):.4f}"


def iso_date(s: str) -> str:
    return datetime.strptime(s.strip(), "%d/%m/%Y").strftime("%Y-%m-%d")


def load_field_table() -> dict[tuple[str, str], dict]:
    """(expedition label, site id) -> curated row. Sites 61-64 are trip-1 named locations (alias table)."""
    with open(ENV_TABLE) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(ALIASES) as fh:
        aliases = {r["source_site_id"]: r for r in csv.DictReader(fh, delimiter="\t") if r["status"] == "confirmed"}
    table = {(r["expedition"], r["site"]): r for r in rows if r["record_role"] == "primary_transect_site"}
    for sid, a in aliases.items():
        for r in rows:
            if r["expedition"] == "Trip 1 (2023)" and r["record_role"] == "named_special_record":
                try:
                    lat, lon = parse_coords(r["coordinates"])
                except ValueError:
                    continue
                if abs(float(lat) - float(a["source_latitude"])) < 1e-4 and abs(float(lon) - float(a["source_longitude"])) < 1e-4:
                    table[("Trip 1 (2023)", sid)] = dict(r, alias_label=a["canonical_site_label"])
    return table


def load_sheet_features() -> dict[tuple[str, str], str]:
    """(trip prefix, site id) -> ENVO id of the field-sheet `feature` for that site and expedition."""
    out: dict[tuple[str, str], str] = {}
    prefix_of = {v[0]: k for k, v in PREFIX_TRIP.items()}
    for trip_no, fname in SHEETS.items():
        with open(ECO / "samplesheets" / fname) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                site = (r.get("site") or "").strip()
                envo = FEATURE_ENVO.get((r.get("feature") or "").strip().lower())
                if site and envo:
                    out[(prefix_of[trip_no], site)] = envo
    return out


def field_row(table: dict, sid: str) -> dict | None:
    m = ID_RE.match(sid)
    return table.get((PREFIX_TRIP[m.group(1)][1], m.group(2))) if m else None


def field_values(field: dict) -> dict[str, tuple[str, str | None]]:
    lat, lon = parse_coords(field["coordinates"])
    return {
        "collection date": (iso_date(field["date"]), None),
        "geographic location (latitude)": (lat, "DD"),
        "geographic location (longitude)": (lon, "DD"),
        "temperature": (field.get("temperature_c", "").strip(), "ºC"),
        "pressure": (field.get("pressure_mbar", "").strip(), "mbar"),
        "humidity": (field.get("relative_humidity_pct", "").strip(), "%"),
    }


SITE_LABEL = {"62": "site road (location 1)", "64": "site campground"}   # named trip-1 locations, house style


def description(sid: str, field: dict) -> str:
    prefix, site, comp, _ = ID_RE.match(sid).groups()
    site_label = SITE_LABEL.get(site) or f"site {site}"
    return f"Rub' al-Khali desert soil, {site_label}, {COMPARTMENT[comp]} compartment, expedition {PREFIX_TRIP[prefix][0]}"


# ----------------------------------------------------------------------------- live records
def account_soil_samples(pw: str, refresh: bool = False) -> list[dict]:
    """All registered samples of the account whose alias is a soil id (report API, cached as JSON)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cache = WORK / "cache" / "account_samples.json"
    if refresh or not cache.exists():
        rc, body = curl(["https://www.ebi.ac.uk/ena/submit/report/samples?format=json&max-results=10000"], pw)
        if rc != 0 or not body.strip().startswith("["):
            sys.exit("ERROR: could not list the account's samples (report API)")
        cache.write_text(body)
    recs = [r["report"] for r in json.loads(cache.read_text())]
    return [r for r in recs if ID_RE.match(r["alias"] or "")]


def cached_record(acc: str, pw: str, refresh: bool = False) -> ET.Element | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{acc}.xml"
    if refresh or not path.exists():
        root = fetch_record("samples", acc, pw)
        if root is None:
            return None
        path.write_bytes(ET.tostring(root, encoding="utf-8"))
    root = ET.parse(path).getroot()
    return root.find("SAMPLE") if root.tag != "SAMPLE" else root


def attrs_of(sample: ET.Element) -> dict[str, str]:
    return {a.findtext("TAG"): (a.findtext("VALUE") or "") for a in sample.iter("SAMPLE_ATTRIBUTE")}


def site_context(records: dict[str, ET.Element]) -> dict[tuple[str, str], str]:
    """(prefix, site) -> local environmental context from registered ERC000022 records; falls back to
    the same site in any expedition (the term is a site property; it never varies by compartment)."""
    by_trip: dict[tuple[str, str], Counter] = defaultdict(Counter)
    by_site: dict[str, Counter] = defaultdict(Counter)
    for alias, rec in records.items():
        m = ID_RE.match(alias)
        v = attrs_of(rec).get("local environmental context")
        if m and v:
            by_trip[(m.group(1), m.group(2))][v] += 1
            by_site[m.group(2)][v] += 1
    out = {k: c.most_common(1)[0][0] for k, c in by_trip.items()}
    for alias in records:
        m = ID_RE.match(alias)
        key = (m.group(1), m.group(2))
        if key not in out and by_site.get(m.group(2)):
            out[key] = by_site[m.group(2)].most_common(1)[0][0]
    return out


# ----------------------------------------------------------------------------- builders
def add_attr(attrs: ET.Element, tag: str, value: str, unit: str | None = None) -> None:
    el = ET.SubElement(attrs, "SAMPLE_ATTRIBUTE")
    ET.SubElement(el, "TAG").text = tag
    ET.SubElement(el, "VALUE").text = value
    if unit:
        ET.SubElement(el, "UNITS").text = unit


def build_new(sid: str, template: ET.Element, field: dict, local_ctx: str | None,
              keep_from: ET.Element | None = None) -> ET.Element:
    """A full ERC000022 record for `sid`: constants from the template, field values from the table.
    `keep_from` (an existing record) contributes its extraction kit, if it has one."""
    prefix, site, comp, rep = ID_RE.match(sid).groups()
    sample = ET.Element("SAMPLE", {"alias": sid})
    ET.SubElement(sample, "TITLE").text = sid
    name = ET.SubElement(sample, "SAMPLE_NAME")
    ET.SubElement(name, "TAXON_ID").text = template.findtext("SAMPLE_NAME/TAXON_ID")
    ET.SubElement(name, "SCIENTIFIC_NAME").text = template.findtext("SAMPLE_NAME/SCIENTIFIC_NAME")
    ET.SubElement(sample, "DESCRIPTION").text = description(sid, field)
    attrs = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")
    values = field_values(field)
    kept_kit = attrs_of(keep_from).get("nucleic acid extraction") if keep_from is not None else None
    if keep_from is not None and site in COORD_HOLD_SITES:      # keep registered coordinates, 4 decimals
        old = attrs_of(keep_from)
        for tag in ("geographic location (latitude)", "geographic location (longitude)"):
            if old.get(tag):
                values[tag] = (f"{float(old[tag]):.4f}", "DD")
    for a in template.iter("SAMPLE_ATTRIBUTE"):
        tag = a.findtext("TAG")
        if tag in DROP_NEW and not (tag == "nucleic acid extraction" and kept_kit):
            continue
        if tag in values:
            val, unit = values[tag]
            if val:
                add_attr(attrs, tag, val, unit)
        elif tag == "depth":
            add_attr(attrs, tag, DEPTH_MM[comp], "mm")
        elif tag == "local environmental context":
            add_attr(attrs, tag, local_ctx or a.findtext("VALUE"))
        elif tag == "project name":
            add_attr(attrs, tag, PROJECT_NAME)
        elif tag == "nucleic acid extraction":
            add_attr(attrs, tag, kept_kit)
        elif tag == "ENA-CHECKLIST":
            add_attr(attrs, tag, CHECKLIST)
        else:
            add_attr(attrs, tag, a.findtext("VALUE"), a.findtext("UNITS"))
    return sample


def patch_existing(rec: ET.Element, sid: str, template: ET.Element, field: dict,
                   local_ctx: str | None) -> tuple[ET.Element, list[tuple[str, str, str]]]:
    """Cleaned copy of a registered record + list of (field, old, new)."""
    old = attrs_of(rec)
    changes: list[tuple[str, str, str]] = []
    if not old.get("ENA-CHECKLIST", "").startswith(CHECKLIST):
        new = build_new(sid, template, field, local_ctx, keep_from=rec)
        new.set("accession", rec.get("accession"))
        if rec.get("center_name"):
            new.set("center_name", rec.get("center_name"))
        newa = attrs_of(new)
        changes.append(("ENA-CHECKLIST", old.get("ENA-CHECKLIST", ""), CHECKLIST))
        changes.append(("TITLE", rec.findtext("TITLE") or "", sid))
        changes.append(("DESCRIPTION", rec.findtext("DESCRIPTION") or "", new.findtext("DESCRIPTION")))
        for tag in sorted(set(old) | set(newa)):
            if tag in ("ENA-CHECKLIST",):
                continue
            if old.get(tag, "") != newa.get(tag, ""):
                changes.append((tag, old.get(tag, ""), newa.get(tag, "")))
        return new, changes

    new = copy.deepcopy(rec)
    new.set("alias", sid)
    if new.findtext("TITLE") != sid:
        changes.append(("TITLE", new.findtext("TITLE") or "", sid))
        new.find("TITLE").text = sid
    desc = description(sid, field)
    d = new.find("DESCRIPTION")
    if d is None:
        d = ET.SubElement(new, "DESCRIPTION")
        new.remove(d)
        new.insert(list(new).index(new.find("SAMPLE_NAME")) + 1, d)
    if (d.text or "") != desc:
        changes.append(("DESCRIPTION", d.text or "", desc))
        d.text = desc
    attrs = new.find("SAMPLE_ATTRIBUTES")
    values = field_values(field)
    for a in list(attrs):
        tag = a.findtext("TAG")
        if tag in DROP_CLEAN:
            changes.append((tag, a.findtext("VALUE") or "", ""))
            attrs.remove(a)
        elif tag == "project name" and a.findtext("VALUE") != PROJECT_NAME:
            changes.append((tag, a.findtext("VALUE") or "", PROJECT_NAME))
            a.find("VALUE").text = PROJECT_NAME
        elif tag == "local environmental context" and local_ctx and a.findtext("VALUE") != local_ctx:
            changes.append((tag, a.findtext("VALUE") or "", local_ctx))       # field sheet wins
            a.find("VALUE").text = local_ctx
        elif tag in ("collection date", "geographic location (latitude)", "geographic location (longitude)"):
            if tag != "collection date" and ID_RE.match(sid).group(2) in COORD_HOLD_SITES:
                continue
            want = values[tag][0]
            have = a.findtext("VALUE") or ""
            same = have == want or (tag != "collection date" and have and abs(float(have) - float(want)) < 2e-4)
            if not same:
                changes.append((tag, have, want))
                a.find("VALUE").text = want
    return new, changes


# ----------------------------------------------------------------------------- driver
def write_xml(sample: ET.Element, path: Path) -> None:
    sset = ET.Element("SAMPLE_SET")
    sset.append(sample)
    ET.indent(sset)
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(sset).write(path, encoding="UTF-8", xml_declaration=True)


def resume_log(path: Path) -> set[str]:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("alias\tstatus\taccession\tdetail\n")
        return set()
    with open(path) as fh:
        return {l.split("\t")[0] for l in fh if "\tOK\t" in l}


def do_submit(sid: str, xml_path: Path, action: str, service: str, pw: str, log: Path, mode_dir: str) -> tuple[bool, str]:
    sxml = WORK / "xml" / f"{mode_dir}.submission.xml"
    sxml.write_text(submission_xml(action))
    receipt = WORK / "receipts" / mode_dir / f"{sid}.{service}.receipt.xml"
    ok, info = submit({"SUBMISSION": sxml, "SAMPLE": xml_path}, pw, service, receipt)
    acc = next((a for t, a, _ in info["accessions"] if t == "SAMPLE" and a), "")
    with open(log, "a") as fh:
        fh.write(f"{sid}\t{'OK' if ok else 'FAILED'}\t{acc}\t{' | '.join(t for _, t in info['messages'])[:200]}\n")
    print_receipt(sid, ok, info)
    return ok, acc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", choices=["register", "cleanup"])
    ap.add_argument("--service", required=True, choices=["test", "prod"])
    ap.add_argument("--mode", required=True, choices=["dryrun", "submit"])
    ap.add_argument("--only", help="a single sample id")
    ap.add_argument("--limit", type=int, help="first N targets only")
    ap.add_argument("--refresh", action="store_true", help="re-fetch cached live records")
    args = ap.parse_args()

    pw = password()
    table = load_field_table()
    template = cached_record(TEMPLATE_ACC, pw, args.refresh)
    if template is None:
        sys.exit(f"ERROR: could not fetch template {TEMPLATE_ACC}")

    # every registered soil sample of the account (live records, cached on disk); the account list itself
    # is re-read for cleanup so samples registered earlier the same day are included
    account = account_soil_samples(pw, refresh=(args.task == "cleanup" or args.refresh))
    records: dict[str, ET.Element] = {}
    print(f"registered soil samples in Webin-65194: {len(account)}; loading live records ...")
    for r in account:
        rec = cached_record(r["id"], pw, args.refresh)
        if rec is not None:
            records[r["alias"]] = rec
    sheet = load_sheet_features()
    sib_ctx = site_context(records)
    ctx = {k: sheet.get(k, v) for k, v in sib_ctx.items()}     # field sheet first, siblings as fallback
    ctx.update(sheet)

    with open(BIOSAMPLE_MAP) as fh:
        map_rows = list(csv.DictReader(fh, delimiter="\t"))

    if args.task == "register":
        targets = [r["sample_id"] for r in map_rows if r["status"].startswith("NO_BIOSAMPLE")]
        already = {a for a in targets if a in records}
        if already:
            sys.exit(f"ERROR: {len(already)} 'missing' samples are now registered: {sorted(already)[:5]} - update the map")
    else:
        targets = sorted(records)
    if args.only:
        targets = [t for t in targets if t == args.only]
    if args.limit:
        targets = targets[: args.limit]

    log = WORK / "logs" / f"{args.task}.{args.service}.tsv"
    done = resume_log(log) if args.mode == "submit" else set()
    xml_dir = WORK / "xml" / args.task
    n_ok = n_fail = n_skip = 0

    if args.task == "register":
        print(f"{'sample':<9} {'exp':<4} {'date':<11} {'lat':<9} {'lon':<9} {'depth':<5} {'temp':<5} {'press':<5} {'hum':<5} local_ctx      description")
        for sid in targets:
            field = field_row(table, sid)
            if field is None:
                print(f"{sid:<9} NEEDS FIELD DATA")
                n_fail += 1
                continue
            m = ID_RE.match(sid)
            sample = build_new(sid, template, field, ctx.get((m.group(1), m.group(2))))
            a = attrs_of(sample)
            xml_path = xml_dir / f"{sid}.sample.xml"
            write_xml(sample, xml_path)
            print(f"{sid:<9} {PREFIX_TRIP[m.group(1)][0]:<4} {a['collection date']:<11} {a['geographic location (latitude)']:<9} "
                  f"{a['geographic location (longitude)']:<9} {a['depth']:<5} {a.get('temperature', '-'):<5} "
                  f"{a.get('pressure', '-'):<5} {a.get('humidity', '-'):<5} {a['local environmental context']:<14} {sample.findtext('DESCRIPTION')}")
            if args.mode == "dryrun":
                continue
            if sid in done:
                n_skip += 1
                continue
            ok, acc = do_submit(sid, xml_path, "ADD", args.service, pw, log, "register")
            if ok:
                n_ok += 1
                if args.service == "prod" and acc:
                    for r in map_rows:
                        if r["sample_id"] == sid:
                            r["status"], r["biosample_accession"] = "REGISTERED_WGS", acc
            else:
                n_fail += 1
        if args.mode == "submit" and args.service == "prod" and n_ok:
            with open(BIOSAMPLE_MAP, "w") as fh:
                w = csv.DictWriter(fh, fieldnames=["sample_id", "biosample_accession", "status"], delimiter="\t")
                w.writeheader()
                w.writerows(map_rows)
            print(f"{BIOSAMPLE_MAP.name} updated; re-run 00_build_runs_table.py to propagate accessions.")

    else:  # cleanup
        diff_path = WORK / "logs" / "cleanup_diff.tsv"
        diff_rows: list[tuple[str, str, str, str, str]] = []
        summary: Counter = Counter()
        for sid in targets:
            rec = records[sid]
            field = field_row(table, sid)
            if field is None:
                print(f"{sid:<9} NEEDS FIELD DATA - skipped")
                n_fail += 1
                continue
            m = ID_RE.match(sid)
            new, changes = patch_existing(rec, sid, template, field, ctx.get((m.group(1), m.group(2))))
            for tag, old, val in changes:
                diff_rows.append((sid, rec.get("accession"), tag, old, val))
                summary[tag] += 1
            if not changes:
                summary["(no change)"] += 1
                continue
            xml_path = xml_dir / f"{sid}.sample.xml"
            write_xml(new, xml_path)
            if args.mode == "dryrun" or sid in done:
                n_skip += sid in done
                continue
            ok, _ = do_submit(sid, xml_path, "MODIFY", args.service, pw, log, "cleanup")
            n_ok += ok
            n_fail += not ok
        with open(diff_path, "w") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["alias", "accession", "field", "old", "new"])
            w.writerows(diff_rows)
        print(f"\nfields changed (samples affected): {dict(summary.most_common())}")
        print(f"full diff: {diff_path} ({len(diff_rows)} rows)")

    print(f"\n{args.task}: {len(targets)} targets" + ("" if args.mode == "dryrun" else
          f"; {n_ok} OK, {n_fail} failed, {n_skip} skipped on {args.service}"))
    if args.mode == "dryrun":
        print(f"DRY RUN - XML in {xml_dir}, nothing sent.")


if __name__ == "__main__":
    main()
