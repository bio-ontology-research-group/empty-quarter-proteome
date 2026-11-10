# Stage 4 v4 Viral Statistical Audit Subsection Pointer (`exp3_viruses/stats_v4`)

The Stage 4 v4 viral statistical audit, compact 591-MAG and candidate-viral 441-sample abundance subsets, 47-site raw independent audit observations, tests, verifier, and Slurm provenance receipts (`15894451`, `15894468`, `15894666`) are packaged in the self-contained [`stats_v4/`](./stats_v4/README.md) subsection directory so as not to conflict with the sibling competitive read-mapping and sequence-followup workflow.

See [`stats_v4/README.md`](./stats_v4/README.md) for full documentation and verification instructions:

```bash
cd /home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/exp3_viruses/stats_v4
uv run python verify_stats_v4.py
uv run pytest -q tests/test_stage4_statistical_audit_v4.py tests/test_stats_v4_package_and_independent_audit.py
```
