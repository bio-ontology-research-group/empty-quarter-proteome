# RQ/ — research-question records

Each `RQ/<id>/` is a research record, not automatically a manuscript claim. It may contain a current verified
result, a bounded historical/exploratory output, or a plan for work that has not been completed. Read the explicit
`Status` and caveats in `result.md`; folder presence and a historical job number do not imply promotion.

Typical contents are:

- `methods.md` — scope, inputs, tools and inferential boundary;
- `result.md` — the current answer or historical disposition, with evidence and references;
- `run.py` / `run.sh` — a lightweight reproduction entry point where one exists;
- compact result tables or figures, with cluster-scale artifacts represented by checksum-pinned provenance.

Only claims registered in `../reproducibility/claims.tsv`, accepted by `../reproducibility/verify.py`, and present
in `../manuscript/main.tex` are active manuscript claims. Run the complete public harness with:

```bash
bash reproducibility/run_all.sh
```

`RQ/_gen_answered.py` regenerates the curated snapshots and deliberately preserves richer manually maintained
folders. Its statuses use plain descriptions because the older green/yellow/red shorthand conflated “computed”
with “fit for inference.”
