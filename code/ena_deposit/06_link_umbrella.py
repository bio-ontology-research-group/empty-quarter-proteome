#!/usr/bin/env python3
"""Manage the RAK umbrella project (ENA Webin-65194). Every argument is explicit; no defaults.

  link   - MODIFY an EXISTING umbrella: add child projects (children keep their own data and release)
           and set the umbrella's release date (HOLD). Refuses a target that is not already an umbrella:
           ENA refused converting the regular study PRJEB104209 (test receipt 2026-09-03).
  cancel - CANCEL a private, never-public, childless regular study. NOT for umbrellas: ENA rejects CANCEL on an
           umbrella ("Invalid study accession", test 2026-11-24, it has no ERP); ask the ENA helpdesk instead.

The umbrella is PRJEB107792 (chosen 2026-11-24); PRJEB106047 is the unused duplicate.
Umbrellas have no Webin Portal route and no ERP study accession: manage them only here.

  python 06_link_umbrella.py link --umbrella PRJEB107792 --children PRJEB106069 PRJEB111735 PRJEB125153 \
         --hold 2026-12-03 --service test|prod --mode dryrun|submit
  python 06_link_umbrella.py cancel --project PRJEB106047 --service test|prod --mode dryrun|submit
"""
from __future__ import annotations

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from datetime import date

from ena_common import WORK, fetch_record, password, print_receipt, project_status, submission_xml, submit

ACC = re.compile(r"PRJEB\d+")


def get_project(acc: str, pw: str, service: str) -> ET.Element:
    root = fetch_record("projects", acc, pw, service=service)
    if root is None:
        sys.exit(f"ERROR: could not fetch {acc} from the {service} drop-box")
    cur = root.find("PROJECT") if root.tag != "PROJECT" else root
    if cur is None:
        sys.exit(f"ERROR: no PROJECT element in the fetched record of {acc}")
    return cur


def describe(cur: ET.Element) -> str:
    kind = "UMBRELLA" if cur.find("UMBRELLA_PROJECT") is not None else "REGULAR"
    kids = [e.get("accession") for e in cur.iter("CHILD_PROJECT")]
    return f"type={kind} children={kids}"


def send(sxml_text: str, parts: dict, label: str, acc: str, args, pw: str) -> bool:
    xml_dir = WORK / "xml"
    xml_dir.mkdir(parents=True, exist_ok=True)
    sxml = xml_dir / f"umbrella.{args.cmd}.{acc}.submission.xml"
    sxml.write_text(sxml_text)
    print(sxml_text)
    if args.mode == "dryrun":
        print("DRY RUN - nothing sent.")
        return True
    receipt = WORK / "receipts" / f"umbrella.{args.cmd}.{acc}.{args.service}.receipt.xml"
    ok, info = submit({"SUBMISSION": sxml, **parts}, pw, args.service, receipt)
    print_receipt(f"{label} ({args.service})", ok, info)
    print(f"receipt: {receipt}")
    return ok


def cmd_link(args, pw: str) -> None:
    for c in args.children:
        if not ACC.fullmatch(c):
            sys.exit(f"ERROR: not a project accession: {c}")
    if args.umbrella in args.children:
        sys.exit("ERROR: the umbrella cannot be its own child")
    try:
        if date.fromisoformat(args.hold) <= date.today():
            sys.exit(f"ERROR: hold date {args.hold} is not in the future")
    except ValueError:
        sys.exit(f"ERROR: hold must be YYYY-MM-DD (got {args.hold!r})")

    cur = get_project(args.umbrella, pw, args.service)
    print(f"{args.umbrella} on {args.service} before: {describe(cur)}")
    if cur.find("UMBRELLA_PROJECT") is None:
        sys.exit(f"ERROR: {args.umbrella} is not an umbrella; ENA does not convert regular studies. Stop.")
    for c in args.children:
        child = get_project(c, pw, args.service)
        if child.find("UMBRELLA_PROJECT") is not None:
            sys.exit(f"ERROR: child {c} is itself an umbrella")

    # Umbrella MODIFY: keep alias/title/description exactly; list only the children to add (ENA docs).
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
    pxml = WORK / "xml" / f"umbrella.link.{args.umbrella}.project.xml"
    pxml.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(pset).write(pxml, encoding="UTF-8", xml_declaration=True)
    print(pxml.read_text())

    ok = send(submission_xml("MODIFY", hold_until=args.hold, target=args.umbrella),
              {"PROJECT": pxml}, "umbrella link + hold", args.umbrella, args, pw)
    if args.mode == "submit":
        after = get_project(args.umbrella, pw, args.service)
        st = project_status(args.umbrella, pw, args.service) or {}
        print(f"{args.umbrella} on {args.service} after: {describe(after)} "
              f"status={st.get('releaseStatus')} hold={st.get('holdDate')}")
        print("NOTE: ENA validates CHILD_PROJECT on MODIFY (an unknown child is rejected; test 2026-11-24) but does not "
              "serialise links of a PRIVATE umbrella; confirm children on the public record after release.")
        if not ok:
            sys.exit(1)


def cmd_cancel(args, pw: str) -> None:
    if not ACC.fullmatch(args.project):
        sys.exit(f"ERROR: not a project accession: {args.project}")
    cur = get_project(args.project, pw, args.service)
    st = project_status(args.project, pw, args.service)
    print(f"{args.project} on {args.service}: {describe(cur)} status={st}")
    if st is None:
        sys.exit("ERROR: no status report; refusing to cancel blind")
    if st.get("releaseStatus") != "PRIVATE" or st.get("firstPublic"):
        sys.exit("ERROR: only private, never-public projects can be cancelled. Stop.")
    if any(True for _ in cur.iter("CHILD_PROJECT")):
        sys.exit("ERROR: project has children; CANCEL would propagate. Stop.")
    ok = send(submission_xml("CANCEL", target=args.project), {}, "project cancel", args.project, args, pw)
    if args.mode == "submit":
        print(f"{args.project} on {args.service} after: status={project_status(args.project, pw, args.service)}")
        if not ok:
            sys.exit(1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    lk = sub.add_parser("link", help="add children to an existing umbrella and set its release date")
    lk.add_argument("--umbrella", required=True)
    lk.add_argument("--children", nargs="+", required=True, help="child project accessions (PRJEB...)")
    lk.add_argument("--hold", required=True, help="umbrella release date YYYY-MM-DD")
    cn = sub.add_parser("cancel", help="cancel a private, childless project")
    cn.add_argument("--project", required=True)
    for p in (lk, cn):
        p.add_argument("--service", required=True, choices=["test", "prod"])
        p.add_argument("--mode", required=True, choices=["dryrun", "submit"])
    args = ap.parse_args()
    pw = password()
    {"link": cmd_link, "cancel": cmd_cancel}[args.cmd](args, pw)


if __name__ == "__main__":
    main()
