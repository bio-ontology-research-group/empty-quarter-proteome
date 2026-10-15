#!/usr/bin/env python3
"""Register a sequencing study (used 2026-09-03 for RAK-WGS = PRJEB125153). No defaults: title, name,
hold date, service and mode are all explicit.

  python 02_register_study.py --title "..." --name "..." --hold YYYY-MM-DD --service test|prod --mode dryrun|submit

Outputs: WORK/xml/study.*.xml, WORK/receipts/study.<service>.receipt.xml, WORK/logs/study_accession.<service>.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from xml.sax.saxutils import escape

from ena_common import WORK, password, print_receipt, submission_xml, submit, UMBRELLA

ALIAS = "eq_soil_shotgun_reads"   # submitter alias; must be unique within Webin-65194
CENTER = "KAUST"

# Study description (Alejandra, 2026-09-03): no counts, dates, tools or accessions; sibling-study register.
DESCRIPTION_DRAFT = (
    "The aims of the study are to characterize the structure and function of the soil microbiome of the "
    "Rub' al-Khali (Empty Quarter), Saudi Arabia, the largest contiguous sand desert on Earth. Soil samples "
    "were collected at multiple time points along a transect of more than 1000 km, from surface, deep and "
    "plant root-adjacent soil, together with extensive metadata characterizing the environmental conditions. "
    "Shotgun metagenomic sequencing is used to assemble a catalogue of metagenome-assembled genomes and their "
    "predicted proteins, to quantify the fraction of genes and proteins without detectable homologs in "
    "reference databases, and to relate the taxonomic and functional composition of the communities, "
    "including traits associated with survival under extreme aridity, radiation and salinity, to the measured "
    "environmental gradients. The data complement the amplicon and cultured-genome studies of the same "
    "sampling campaign."
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--title", required=True, help="study title (per Rob)")
    ap.add_argument("--name", required=True, help="short study name (NAME element)")
    ap.add_argument("--description", help="description text; default = DESCRIPTION_DRAFT in this file")
    ap.add_argument("--hold", required=True, help="release date YYYY-MM-DD (HOLD action); explicit, no default")
    ap.add_argument("--service", required=True, choices=["test", "prod"])
    ap.add_argument("--mode", required=True, choices=["dryrun", "submit"])
    args = ap.parse_args()
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", args.hold):
        sys.exit("ERROR: --hold must be YYYY-MM-DD")

    desc = args.description or DESCRIPTION_DRAFT
    project_xml = (
        "<PROJECT_SET>\n"
        f'  <PROJECT alias="{ALIAS}" center_name="{CENTER}">\n'
        f"    <NAME>{escape(args.name)}</NAME>\n"
        f"    <TITLE>{escape(args.title)}</TITLE>\n"
        f"    <DESCRIPTION>{escape(desc)}</DESCRIPTION>\n"
        "    <SUBMISSION_PROJECT>\n      <SEQUENCING_PROJECT/>\n    </SUBMISSION_PROJECT>\n"
        "  </PROJECT>\n</PROJECT_SET>\n"
    )
    xml_dir = WORK / "xml"
    xml_dir.mkdir(parents=True, exist_ok=True)
    pxml = xml_dir / "study.project.xml"
    sxml = xml_dir / "study.submission.xml"
    pxml.write_text(project_xml)
    sxml.write_text(submission_xml("ADD", hold_until=args.hold))
    print(project_xml)
    print(f"written: {pxml}\n         {sxml}")
    if args.mode == "dryrun":
        print("DRY RUN - nothing sent.")
        return

    receipt = WORK / "receipts" / f"study.{args.service}.receipt.xml"
    ok, info = submit({"SUBMISSION": sxml, "PROJECT": pxml}, password(), args.service, receipt)
    print_receipt(f"study ADD ({args.service})", ok, info)
    print(f"receipt: {receipt}")
    if ok:
        log = WORK / "logs" / f"study_accession.{args.service}.tsv"
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "a") as fh:
            w = csv.writer(fh, delimiter="\t")
            for tag, acc, alias in info["accessions"]:
                if acc:
                    w.writerow([tag, acc, alias, args.service, args.hold])
        print(f"accessions appended to {log}")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
