#!/usr/bin/env python3
"""Bounded persistent ws-ts watcher for cross-cluster IBEX SignalP -> Shaheen Stage 2b handoff.

Implements `external-join-final-launch.txt` Item 2:
- Runs on `ws-ts` with project-local state/logs in `science-20261108/watcher/` (no `/tmp`).
- Polls every 90 seconds (within 60-120s window) up to 48 hours (1,920 iterations max).
- Checks IBEX jobs (`55372501`, `55372577`, `55372578`, `55372579`) and Shaheen jobs (`15867464`, `15867465`).
- Aborts immediately if any upstream job enters `FAILED`, `CANCELLED`, `TIMEOUT`, `NODE_FAIL`, `OUT_OF_MEMORY`,
  `PREEMPTED`, or `DependencyNeverSatisfied`.
- Once BOTH IBEX collector `55372579` and Shaheen collector `15867465` have `COMPLETED` (`0:0`):
  1. Downloads only the verified compact SignalP output + manifests from IBEX to `ws-ts`.
  2. Verifies non-empty files and computes SHA-256 checksums.
  3. Stages them to `shaheen:/scratch/project/k10018/eq-science-20261108/results/`.
  4. Submits the full Stage 2b job (`sbatch_stage2b_join_canonical_and_signalp.sh`, without `--allow-pending-signalp`),
     recording an idempotent receipt (`watcher/full_stage2b_submission_receipt.json`) to prevent duplicate submissions.
  5. Monitors the full Stage 2b job to `COMPLETED` and syncs final Stage 2 + Stage 2b results back to `ws-ts`.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

BASE_DIR = Path("/home/leechuck/Documents/papers/empty-quarter/science-20261108")
WATCHER_DIR = BASE_DIR / "watcher"
RESULTS_DIR = BASE_DIR / "results"
STATE_FILE = WATCHER_DIR / "watcher_state.json"
RECEIPT_FILE = WATCHER_DIR / "full_stage2b_submission_receipt.json"

IBEX_DIR = "/ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp"
SHAHEEN_DIR = "/scratch/project/k10018/eq-science-20261108"

IBEX_JOBS = ["55372501", "55372577", "55372578", "55372579"]
IBEX_COLLECTOR_JOB = "55372579"
SHAHEEN_JOBS = ["15867464", "15867465"]
SHAHEEN_COLLECTOR_JOB = "15867465"

POLL_SECONDS = 90
MAX_ITERATIONS = 1920  # 48 hours at 90s

FAILURE_STATES = {
    "FAILED",
    "CANCELLED",
    "TIMEOUT",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "BOOT_FAIL",
    "DEADLINE",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def log(msg: str) -> None:
    print(f"[{now_iso()}] {msg}", flush=True)


def write_state(payload: dict[str, object]) -> None:
    WATCHER_DIR.mkdir(parents=True, exist_ok=True)
    payload["updated_at"] = now_iso()
    payload["pid"] = os.getpid()
    tmp_p = STATE_FILE.with_suffix(".json.tmp")
    tmp_p.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp_p.replace(STATE_FILE)


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def query_cluster_jobs(host: str, job_ids: list[str]) -> tuple[dict[str, list[tuple[str, str]]], list[str]]:
    j_str = ",".join(job_ids)
    cmd = (
        f"sacct -n -X -j {j_str} --format=JobID%25,State%20,ExitCode%12 ; "
        f"echo '---SQUEUE---' ; "
        f"squeue -h -j {j_str} -o '%i|%T|%r' 2>/dev/null || true"
    )
    res = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", host, cmd],
        capture_output=True,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"SSH query to {host} failed (rc={res.returncode}): {res.stderr.strip()}")

    sacct_map: dict[str, list[tuple[str, str]]] = {jid: [] for jid in job_ids}
    squeue_lines: list[str] = []
    in_squeue = False
    for raw in res.stdout.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line == "---SQUEUE---":
            in_squeue = True
            continue
        if not in_squeue:
            parts = line.split()
            if len(parts) >= 3:
                raw_jid, st, ec = parts[0], parts[1].rstrip("+"), parts[2]
                base_jid = raw_jid.split("_")[0]
                if base_jid in sacct_map:
                    sacct_map[base_jid].append((st, ec))
        else:
            squeue_lines.append(line)
    return sacct_map, squeue_lines


def detect_failures(host: str, sacct_map: dict[str, list[tuple[str, str]]], squeue_lines: list[str]) -> str | None:
    for jid, entries in sacct_map.items():
        for st, ec in entries:
            st_base = st.split()[0]
            if st_base in FAILURE_STATES:
                return f"{host} job {jid} entered failure state {st_base} (ExitCode={ec})"
            if st_base == "COMPLETED" and ec != "0:0":
                return f"{host} job {jid} COMPLETED with non-zero ExitCode={ec}"
    for sq in squeue_lines:
        if "DependencyNeverSatisfied" in sq:
            return f"{host} squeue reports DependencyNeverSatisfied: {sq}"
    return None


def is_collector_completed(sacct_map: dict[str, list[tuple[str, str]]], collector_jid: str) -> bool:
    entries = sacct_map.get(collector_jid, [])
    if not entries:
        return False
    return all(st == "COMPLETED" and ec == "0:0" for st, ec in entries)


def transfer_and_stage_signalp_outputs() -> dict[str, str]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    files_to_fetch = [
        (f"{IBEX_DIR}/out/external247_signalp6_predictions.tsv.gz", RESULTS_DIR / "external247_signalp6_predictions.tsv.gz"),
        (f"{IBEX_DIR}/out/external247_signalp6_per_genome_summary.tsv", RESULTS_DIR / "external247_signalp6_per_genome_summary.tsv"),
        (f"{IBEX_DIR}/external_signalp_manifest.tsv", RESULTS_DIR / "external_signalp_manifest.tsv"),
        (f"{IBEX_DIR}/external_signalp_chunks_summary.tsv", RESULTS_DIR / "external_signalp_chunks_summary.tsv"),
    ]
    hashes: dict[str, str] = {}
    for remote_src, local_dst in files_to_fetch:
        log(f"Downloading {remote_src} from IBEX -> {local_dst}")
        subprocess.run(
            ["rsync", "-avz", f"ibex:{remote_src}", str(local_dst)],
            check=True,
        )
        if not local_dst.is_file() or local_dst.stat().st_size == 0:
            raise RuntimeError(f"Downloaded file {local_dst} is missing or empty")
        hashes[local_dst.name] = sha256_file(local_dst)

    log("Staging verified IBEX SignalP outputs and manifest to Shaheen...")
    subprocess.run(
        [
            "rsync",
            "-avz",
            str(RESULTS_DIR / "external247_signalp6_predictions.tsv.gz"),
            str(RESULTS_DIR / "external247_signalp6_per_genome_summary.tsv"),
            str(RESULTS_DIR / "external_signalp_manifest.tsv"),
            str(RESULTS_DIR / "external_signalp_chunks_summary.tsv"),
            f"shaheen:{SHAHEEN_DIR}/results/",
        ],
        check=True,
    )
    return hashes


def submit_full_stage2b_idempotent(hashes: dict[str, str]) -> str:
    if RECEIPT_FILE.is_file():
        rec = json.loads(RECEIPT_FILE.read_text(encoding="utf-8"))
        jid = str(rec.get("shaheen_full_stage2b_job_id", "")).strip()
        if jid:
            log(f"Idempotent receipt already exists for full Stage 2b job {jid}; skipping duplicate sbatch.")
            return jid

    log("Submitting full Stage 2b job on Shaheen (sbatch_stage2b_join_canonical_and_signalp.sh)...")
    res = subprocess.run(
        [
            "ssh",
            "-o",
            "BatchMode=yes",
            "shaheen",
            f"sbatch {SHAHEEN_DIR}/code/sbatch_stage2b_join_canonical_and_signalp.sh",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    out = res.stdout.strip()
    log(f"Shaheen sbatch response: {out}")
    tokens = out.split()
    if len(tokens) < 4 or not tokens[-1].isdigit():
        raise RuntimeError(f"Unexpected sbatch output: {out}")
    job_id = tokens[-1]
    receipt = {
        "submitted_at": now_iso(),
        "shaheen_full_stage2b_job_id": job_id,
        "ibex_collector_job_id": IBEX_COLLECTOR_JOB,
        "shaheen_stage2_collector_job_id": SHAHEEN_COLLECTOR_JOB,
        "staged_file_sha256": hashes,
    }
    RECEIPT_FILE.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return job_id


def main() -> None:
    WATCHER_DIR.mkdir(parents=True, exist_ok=True)
    log(f"Starting ws-ts Stage 2b handoff watcher (PID={os.getpid()}, poll={POLL_SECONDS}s, max_iter={MAX_ITERATIONS})")

    for iteration in range(1, MAX_ITERATIONS + 1):
        try:
            ibex_map, ibex_sq = query_cluster_jobs("ibex", IBEX_JOBS)
            sha_map, sha_sq = query_cluster_jobs("shaheen", SHAHEEN_JOBS)
        except Exception as exc:
            log(f"Transient SSH poll warning on iteration {iteration}: {exc}")
            time.sleep(POLL_SECONDS)
            continue

        ibex_err = detect_failures("ibex", ibex_map, ibex_sq)
        sha_err = detect_failures("shaheen", sha_map, sha_sq)
        if ibex_err or sha_err:
            err_msg = ibex_err or sha_err or "unknown_failure"
            log(f"ABORTING watcher due to upstream failure: {err_msg}")
            write_state(
                {
                    "status": "aborted_on_upstream_job_failure",
                    "iteration": iteration,
                    "error": err_msg,
                    "ibex_sacct": ibex_map,
                    "shaheen_sacct": sha_map,
                }
            )
            sys.exit(1)

        ibex_done = is_collector_completed(ibex_map, IBEX_COLLECTOR_JOB)
        sha_done = is_collector_completed(sha_map, SHAHEEN_COLLECTOR_JOB)

        write_state(
            {
                "status": "polling_upstream_collectors",
                "iteration": iteration,
                "ibex_collector_55372579_completed": ibex_done,
                "shaheen_collector_15867465_completed": sha_done,
                "ibex_squeue": ibex_sq,
                "shaheen_squeue": sha_sq,
            }
        )
        log(
            f"Poll {iteration}/{MAX_ITERATIONS}: "
            f"IBEX_55372579_done={ibex_done} (squeue={len(ibex_sq)}), "
            f"Shaheen_15867465_done={sha_done} (squeue={len(sha_sq)})"
        )

        if ibex_done and sha_done:
            log("Both IBEX collector 55372579 and Shaheen collector 15867465 are COMPLETED! Proceeding to transfer & full Stage 2b submission.")
            hashes = transfer_and_stage_signalp_outputs()
            stage2b_jid = submit_full_stage2b_idempotent(hashes)
            write_state(
                {
                    "status": "submitted_full_stage2b_monitoring",
                    "iteration": iteration,
                    "shaheen_full_stage2b_job_id": stage2b_jid,
                    "staged_file_sha256": hashes,
                }
            )
            # Monitor full Stage 2b job until completion
            for sub_iter in range(1, 240):
                time.sleep(60)
                s2b_map, s2b_sq = query_cluster_jobs("shaheen", [stage2b_jid])
                s2b_err = detect_failures("shaheen", s2b_map, s2b_sq)
                if s2b_err:
                    log(f"Full Stage 2b job {stage2b_jid} failed: {s2b_err}")
                    write_state(
                        {
                            "status": "aborted_on_full_stage2b_failure",
                            "shaheen_full_stage2b_job_id": stage2b_jid,
                            "error": s2b_err,
                        }
                    )
                    sys.exit(1)
                if is_collector_completed(s2b_map, stage2b_jid):
                    log(f"Full Stage 2b job {stage2b_jid} COMPLETED! Syncing final results from Shaheen to ws-ts...")
                    subprocess.run(
                        [
                            "rsync",
                            "-avz",
                            f"shaheen:{SHAHEEN_DIR}/results/task_b_*",
                            f"shaheen:{SHAHEEN_DIR}/results/stage2_summary.json",
                            f"{RESULTS_DIR}/",
                        ],
                        check=True,
                    )
                    write_state(
                        {
                            "status": "completed_full_stage2b",
                            "shaheen_full_stage2b_job_id": stage2b_jid,
                            "staged_file_sha256": hashes,
                        }
                    )
                    log("Handoff watcher completed all tasks successfully.")
                    return
            raise RuntimeError(f"Timed out waiting for full Stage 2b job {stage2b_jid} to complete")

        time.sleep(POLL_SECONDS)

    write_state({"status": "aborted_on_max_iterations_timeout", "iteration": MAX_ITERATIONS})
    sys.exit(1)


if __name__ == "__main__":
    main()
