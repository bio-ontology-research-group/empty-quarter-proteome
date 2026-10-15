#!/usr/bin/env python3
"""Attach child projects to an umbrella project (ENA MODIFY of the umbrella; children keep their own release).

ENA REFUSED converting PRJEB104209 into an umbrella (test receipt 2026-09-03): a NEW umbrella project is
needed, registered last for the PI's approval; the "new umbrella" mode is not written yet.

  python 06_link_umbrella.py --umbrella PRJEB... --children PRJEB106069 PRJEB111735 PRJEB125153 --service test|prod --mode dryrun|submit
"""
from __future__ import annotations

import argparse
import re
import sys
import xml.etree.ElementTree as ET

from ena_common import UMBRELLA, WORK, fetch_record, password, print_receipt, submission_xml, submit


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--umbrella", default=UMBRELLA)
    ap.add_argument("--children", nargs="+", required=True, help="child project accessions (PRJEB...)")
    ap.add_argument("--service", required=True, choices=["test", "prod"])
    ap.add_argument("--mode", required=True, choices=["dryrun", "submit"])
    args = ap.parse_args()
    for c in args.children:
        if not re.fullmatch(r"PRJEB\d+", c):
            sys.exit(f"ERROR: not a project accession: {c}")

    pw = password()
    root = fetch_record("projects", args.umbrella, pw, service=args.service)
    if root is None:
        sys.exit(f"ERROR: could not fetch {args.umbrella} from the {args.service} drop-box")
    cur = root.find("PROJECT") if root.tag != "PROJECT" else root
    if cur is None:
        sys.exit("ERROR: no PROJECT element in the fetched record")

    proj = ET.Element("PROJECT", {"alias": cur.get("alias"), "accession": cur.get("accession", args.umbrella),
                                  "center_name": cur.get("center_name", "")})
    for tag in ("NAME", "TITLE", "DESCRIPTION"):
        el = cur.find(tag)
        if el is not None:
            ET.SubElement(proj, tag).text = el.text
    ET.SubElement(proj, "UMBRELLA_PROJECT")
    rel = ET.SubElement(proj, "RELATED_PROJECTS")
    for c in args.children:
        ET.SubElement(ET.SubElement(rel, "RELATED_PROJECT"), "CHILD_PROJECT", {"accession": c})
    pset = ET.Element("PROJECT_SET")
    pset.append(proj)
    ET.indent(pset)

    xml_dir = WORK / "xml"
    xml_dir.mkdir(parents=True, exist_ok=True)
    pxml = xml_dir / f"umbrella.{args.umbrella}.project.xml"
    sxml = xml_dir / "umbrella.submission.xml"
    ET.ElementTree(pset).write(pxml, encoding="UTF-8", xml_declaration=True)
    sxml.write_text(submission_xml("MODIFY"))
    print(pxml.read_text())
    print(f"current type on {args.service}: "
          f"{'UMBRELLA' if cur.find('UMBRELLA_PROJECT') is not None else 'SEQUENCING'}; "
          f"existing children: {[e.get('accession') for e in cur.iter('CHILD_PROJECT')]}")
    if args.mode == "dryrun":
        print("DRY RUN - nothing sent.")
        return

    receipt = WORK / "receipts" / f"umbrella.{args.umbrella}.{args.service}.receipt.xml"
    ok, info = submit({"SUBMISSION": sxml, "PROJECT": pxml}, pw, args.service, receipt)
    print_receipt(f"umbrella MODIFY ({args.service})", ok, info)
    print(f"receipt: {receipt}")
    if ok:
        after = fetch_record("projects", args.umbrella, pw, service=args.service)
        kids = [e.get("accession") for e in after.iter("CHILD_PROJECT")] if after is not None else "?"
        print(f"re-fetched from {args.service}: type="
              f"{'UMBRELLA' if after is not None and after.find('.//UMBRELLA_PROJECT') is not None else 'SEQUENCING'} "
              f"children={kids}")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
