#!/usr/bin/env python3
"""Fold the ENA experiment/run accessions into wgs_runs.tsv (columns ena_experiment, ena_run).

Source of truth = ENA's live lists for the study (report API), cross-checked against the two local logs
(Webin-CLI: logs/webin_submit_status.tsv; XML: logs/runs_xml.prod.tsv). Aliases are webin-reads-<run_alias>.

  python 09_fold_accessions.py
"""
from __future__ import annotations

import csv
import json
import sys

from ena_common import HERE, WORK, curl, password

RUNS = HERE / "wgs_runs.tsv"
REPORT = "https://www.ebi.ac.uk/ena/submit/report/{kind}?format=json&max-results=10000"
STUDY_ERP = "ERP204894"


def live(kind: str, pw: str) -> dict[str, str]:
    rc, body = curl([REPORT.format(kind=kind)], pw)
    if rc != 0 or not body.strip().startswith("["):
        sys.exit(f"ERROR: could not list {kind}")
    return {r["report"]["alias"]: r["report"]["id"] for r in json.loads(body) if r["report"].get("studyId") == STUDY_ERP}


def main() -> None:
    pw = password()
    exps, runs = live("experiments", pw), live("runs", pw)
    local: dict[str, tuple[str, str]] = {}
    for name, cols in (("webin_submit_status.tsv", ("erx", "err")), ("runs_xml.prod.tsv", ("erx", "err"))):
        p = WORK / "logs" / name
        if p.exists():
            for r in csv.DictReader(open(p), delimiter="\t"):
                if r["status"] == "OK":
                    local[r["alias"]] = (r[cols[0]], r[cols[1]])
    with open(RUNS) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows, cols = list(reader), list(reader.fieldnames)
    for c in ("ena_experiment", "ena_run"):
        if c not in cols:
            cols.insert(cols.index("r1_bytes"), c)
    missing, mismatch = [], []
    for r in rows:
        alias = f"webin-reads-{r['run_alias']}"
        r["ena_experiment"], r["ena_run"] = exps.get(alias, ""), runs.get(alias, "")
        if not r["ena_run"]:
            missing.append(r["run_alias"])
        elif r["run_alias"] in local and local[r["run_alias"]] != (r["ena_experiment"], r["ena_run"]):
            mismatch.append((r["run_alias"], local[r["run_alias"]], (r["ena_experiment"], r["ena_run"])))
    with open(RUNS, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} runs; with ENA run accession: {len(rows) - len(missing)}; missing: {missing[:5]}; "
          f"local-vs-ENA mismatches: {mismatch[:3]}")


if __name__ == "__main__":
    main()
