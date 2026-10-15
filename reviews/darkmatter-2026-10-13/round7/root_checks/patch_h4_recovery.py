from pathlib import Path
import shutil,json
b=Path('reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs')
snap=b/'flash_candidates_v2';snap.mkdir(exist_ok=True)
for p in b.glob('v7_*'):shutil.copy2(p,snap/p.name)
p=b/'v7_crop_and_search_regions.py';s=p.read_text()
s=s.replace('results[name] = {"status": "NOT_ACCESSIBLE_LOCALLY", "path": p}', 'raise FileNotFoundError(f"Required baseline missing: {p}")')
a=s.index('def run_foldseek_search_atomic(');z=s.index('\ndef parse_and_accumulate_hits',a)
new='''def verified_database_identity(db_path):
    """Rehash all 32 frozen PDB database components, excluding its legacy self-entry."""
    pin_path = os.path.join(os.path.dirname(__file__), "pdb_reference_verification.json")
    with open(pin_path) as fh:
        pin = json.load(fh)
    assert pin["all_match"] and len(pin["files"]) == 32
    assert os.path.realpath(db_path) == os.path.realpath(PDB_DB)
    verified = {}
    for row in pin["files"]:
        path = row["path"]
        assert os.path.getsize(path) == row["bytes"], path
        digest = sha256_file(path)
        assert digest == row["expected"], f"Database component changed: {path}"
        verified[path] = digest
    return verified


def run_foldseek_search_atomic(query_dir, db_path, out_tsv, tmp_dir, threads, max_seqs, reuse=False):
    """Fail closed on stale/unsigned evidence; publish output only after subprocess success."""
    sig_path = out_tsv + ".signature.json"
    query_files = sorted(f for f in os.listdir(query_dir) if f.endswith(".pdb"))
    assert query_files, query_dir
    query_hashes = {q[:-4]: sha256_file(os.path.join(query_dir, q)) for q in query_files}
    db_hashes = verified_database_identity(db_path)
    version = subprocess.check_output(["foldseek", "version"], text=True).strip()
    fmt = "query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,fident,pident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend,cigar"
    settings = {"alignment_type": 1, "evalue": 10.0, "sensitivity": 9.5,
                "max_seqs": max_seqs, "threads": threads, "format_output": fmt}
    contract = {"query_id_hashes": query_hashes, "database_components": db_hashes,
                "foldseek_version": version, "settings": settings, "db_path": db_path}
    if os.path.exists(out_tsv) or os.path.exists(sig_path):
        assert reuse, f"Refusing to overwrite existing evidence: {out_tsv}"
        assert os.path.isfile(sig_path) and os.path.isfile(out_tsv), "Unsigned/incomplete cached search"
        with open(sig_path) as fh:
            sig = json.load(fh)
        assert sig.get("status") == "COMPLETED" and sig.get("returncode") == 0
        assert sig.get("contract") == contract, "Search input/version/settings mismatch"
        assert os.path.getsize(out_tsv) == sig["raw_size_bytes"]
        assert sha256_file(out_tsv) == sig["raw_sha256"], "Search output hash mismatch"
        return sig["elapsed_sec"]
    os.makedirs(tmp_dir, exist_ok=True)
    token = f"{os.getpid()}_{time.time_ns()}"
    tmp_out = out_tsv + ".tmp_" + token
    work = os.path.join(tmp_dir, token)
    cmd = ["foldseek", "easy-search", query_dir, db_path, tmp_out, work,
           "--alignment-type", "1", "-e", "10", "-s", "9.5", "--max-seqs", str(max_seqs),
           "--threads", str(threads), "--format-output", fmt]
    print(f"Executing: {' '.join(cmd)}", flush=True)
    started = time.time()
    subprocess.run(cmd, check=True)
    elapsed = round(time.time() - started, 2)
    # Zero-hit outputs are legitimate if the checked subprocess succeeded.
    assert os.path.isfile(tmp_out)
    os.replace(tmp_out, out_tsv)
    sig = {"status": "COMPLETED", "returncode": 0, "contract": contract, "command": cmd,
           "raw_sha256": sha256_file(out_tsv), "raw_size_bytes": os.path.getsize(out_tsv),
           "elapsed_sec": elapsed, "completed_at_utc": datetime.now(timezone.utc).isoformat()}
    with open(sig_path + ".tmp_" + token, "w") as fh:
        json.dump(sig, fh, indent=2)
    os.replace(sig_path + ".tmp_" + token, sig_path)
    print(f"Completed {out_tsv} in {elapsed:.2f}s ({elapsed/60.0:.2f}m)", flush=True)
    return elapsed

'''
s=s[:a]+new+s[z:]
needle='    # Execute search with atomic completion signature\n'
s=s.replace(needle,'''    actual_ids = {f[:-4] for f in os.listdir(rerun_dir) if f.endswith(".pdb")}
    assert actual_ids == set(prob_qids), "Adaptive query directory differs from requested census"
    for qid in prob_qids:
        assert sha256_file(os.path.join(rerun_dir, qid + ".pdb")) == sha256_file(os.path.join(query_dir, qid + ".pdb"))

'''+needle)
p.write_text(s)
shutil.copy2('reviews/darkmatter-2026-10-13/round6/root_checks/pdb_reference_verification.json',b/'pdb_reference_verification.json')
p=b/'v7_recovery_postprocess.py';s=p.read_text().replace('import argparse','import argparse\nimport subprocess')
s=s.replace('generate_sha256_manifest\n','generate_sha256_manifest, verified_database_identity\n')
s=s.replace('def validate_slurm_log_state(log_path, wc_hits_path, crop_hits_path):','def validate_slurm_log_state(log_path, wc_hits_path, crop_hits_path, stderr_path=None):')
s=s.replace('    # Check terminal state\n', '''    if stderr_path is not None:
        assert os.path.isfile(stderr_path), f"Missing Slurm stderr: {stderr_path}"
        with open(stderr_path) as fh:
            log_text += "\\n" + fh.read()

    # Check terminal state
''')
s=s.replace('log_meta = validate_slurm_log_state(log_path, wc_hits_path, crop_hits_path)', '''assert log_path.endswith("h4_v7_53738540.out"), "Recovery contract is for original job53738540"
    state_text = subprocess.check_output(["sacct", "-j", "53738540", "-n", "-X", "-P", "--format=JobID,State,ExitCode"], text=True)
    state_rows = [r.split("|") for r in state_text.splitlines() if r.strip()]
    assert len(state_rows) == 1 and state_rows[0][0] == "53738540"
    assert state_rows[0][1] in ("FAILED", "COMPLETED"), f"Original job not terminal: {state_text}"
    err_path = log_path[:-4] + ".err"
    log_meta = validate_slurm_log_state(log_path, wc_hits_path, crop_hits_path, err_path)
    log_meta.update({"slurm_accounting": state_rows[0], "stderr_path": err_path, "stderr_sha256": sha256_file(err_path)})
    reference_hashes = verified_database_identity(PDB_DB)
    query_hashes = {"whole_chain": {}, "crop": {}}
    with open(manifest_path) as fh:
        for row in csv.DictReader(fh, delimiter="\\t"):
            qid = row["model_id"]
            if qid not in query_hashes["whole_chain"]:
                digest = sha256_file(os.path.join(H4_DIR, "v7_wc_queries", qid + ".pdb"))
                assert digest == row["source_pdb_sha256"], qid
                query_hashes["whole_chain"][qid] = digest
            if row["cropped_pdb"] != "NONE":
                name = row["cropped_pdb"]
                digest = sha256_file(os.path.join(H4_DIR, "v7_crops", name))
                assert digest == row["crop_pdb_sha256"], name
                assert name not in query_hashes["crop"], name
                query_hashes["crop"][name] = digest
    assert len(query_hashes["whole_chain"]) == 5851 and len(query_hashes["crop"]) == 7773
    log_meta["database_components"] = reference_hashes
    log_meta["verified_query_counts"] = {k: len(v) for k,v in query_hashes.items()}''')
p.write_text(s)
