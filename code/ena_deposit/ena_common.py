#!/usr/bin/env python3
"""Shared helpers for the WGS reads deposit (ENA Webin-65194).

Adapted from eqmgc/data/metadata/ena_submission/{07,08}_modify_*.py (2026-09).
Every script in this folder imports from here so that:
  * the Webin password is read ONLY from ~/.webin_pass (never argv, never a file in git);
  * the credential reaches curl through a config on stdin (-K -), never the process list;
  * "test" (wwwdev, rebuilt from production daily at 03:00, discarded next day) and
    "prod" endpoints are selected explicitly - there is no default service.
"""
from __future__ import annotations

import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

WEBIN_USER = "Webin-65194"
PASS_FILE = Path.home() / ".webin_pass"

# Work directory outside git (symlink farm, manifests, receipts, logs). Agreed 2026-09-03.
WORK = Path("/ibex/scratch/projects/c2014/alelopezv/ena_wgs")
HERE = Path(__file__).resolve().parent

SERVICES = {
    "test": "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/",
    "prod": "https://www.ebi.ac.uk/ena/submit/drop-box/",
}

UMBRELLA = "PRJEB104209"          # public; owned by Webin-65194 (account study report, 2026-09-03)
AMPLICON_STUDY = "PRJEB106069"    # private, hold 2026-11-25; owns the 445 reusable BioSamples
RAKPR_STUDY = "PRJEB111735"       # private, hold 2026-11-01


def password() -> str:
    if not PASS_FILE.is_file():
        sys.exit(f"ERROR: {PASS_FILE} not found (Webin password file, mode 600).")
    pw = PASS_FILE.read_text().strip()
    if not pw:
        sys.exit(f"ERROR: {PASS_FILE} is empty.")
    return pw


def service_url(service: str, path: str = "submit/") -> str:
    if service not in SERVICES:
        sys.exit(f"ERROR: service must be one of {list(SERVICES)} (got {service!r}).")
    return SERVICES[service] + path


def curl(args: list[str], pw: str, timeout: int = 300) -> tuple[int, str]:
    """Run curl; the credential goes in as a config file on stdin, not argv."""
    conf = f'user = "{WEBIN_USER}:{pw}"\n'
    p = subprocess.run(["curl", "-s", "-m", str(timeout), "-K", "-", *args],
                       input=conf, capture_output=True, text=True)
    return p.returncode, p.stdout


def fetch_record(kind: str, accession: str, pw: str, service: str = "prod") -> ET.Element | None:
    """GET drop-box/{samples|projects|studies}/{acc}; works for private records with auth."""
    rc, body = curl([service_url(service, f"{kind}/{accession}")], pw)
    if rc != 0 or not body.strip() or "<ErrorDetails" in body:
        return None
    try:
        return ET.fromstring(body)
    except ET.ParseError:
        return None


def submission_xml(action: str, hold_until: str | None = None) -> str:
    """SUBMISSION XML with one action (ADD | MODIFY) and an optional HOLD date (YYYY-MM-DD)."""
    if action not in ("ADD", "MODIFY"):
        sys.exit(f"ERROR: action must be ADD or MODIFY (got {action!r}).")
    hold = f'\n      <ACTION>\n        <HOLD HoldUntilDate="{hold_until}"/>\n      </ACTION>' if hold_until else ""
    return (
        "<SUBMISSION_SET>\n  <SUBMISSION>\n    <ACTIONS>\n      <ACTION>\n"
        f"        <{action}/>\n      </ACTION>{hold}\n    </ACTIONS>\n  </SUBMISSION>\n</SUBMISSION_SET>\n"
    )


def submit(parts: dict[str, Path], pw: str, service: str, receipt_path: Path) -> tuple[bool, dict]:
    """POST XML parts (e.g. {"SUBMISSION": p, "SAMPLE": q}) to the chosen service; save the receipt."""
    args: list[str] = []
    for name, path in parts.items():
        args += ["-F", f"{name}=@{path}"]
    args.append(service_url(service, "submit/"))
    rc, body = curl(args, pw)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(body)
    info = parse_receipt(body)
    return rc == 0 and info["success"], info


def parse_receipt(body: str) -> dict:
    """Receipt -> {success, accessions:[(tag, accession, alias)], messages:[(level, text)]}."""
    out = {"success": False, "accessions": [], "messages": [], "raw": body[:2000]}
    if not body.strip().startswith("<"):
        out["messages"].append(("ERROR", body[:500]))
        return out
    try:
        root = ET.fromstring(body)
    except ET.ParseError as e:
        out["messages"].append(("ERROR", f"unparseable receipt: {e}"))
        return out
    out["success"] = root.get("success") == "true"
    for child in root:
        if child.tag in ("PROJECT", "SAMPLE", "EXPERIMENT", "RUN", "SUBMISSION", "STUDY"):
            out["accessions"].append((child.tag, child.get("accession"), child.get("alias")))
    for msg in root.iter("MESSAGES"):
        for m in msg:
            out["messages"].append((m.tag, (m.text or "").strip()))
    return out


def print_receipt(label: str, ok: bool, info: dict) -> None:
    flag = "OK " if ok else "FAIL"
    accs = " ".join(f"{t}:{a}" for t, a, _ in info["accessions"] if a)
    print(f"[{flag}] {label} {accs}")
    for level, text in info["messages"]:
        print(f"       {level}: {text[:300]}")
