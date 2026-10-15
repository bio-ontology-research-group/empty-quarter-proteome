# IBEX login memory incident — 2026-10-14

User reported per-user memory exceeding 20 GB on login509-02-l; notice names systemd PID 1040824. The notice has no timestamp or offending-memory-consumer breakdown.

## Findings

- agy Flash conversation 639411ba-5c55-49dd-a0ee-4642f7d41649 acknowledges direct-SSH Foldseek createsubdb/convert2pdb, standalone g++ compilation, and TM-align execution during the H4 anomaly investigation. These helpers should have used Slurm. Its audit quotes commands and reports a 40-thread Foldseek setting.
- Remote binary/extracted-PDB timestamps independently agree with the reported approximate helper timing (10:05–10:06 local). They do not prove which host executed them or historical memory consumption.
- Root inspected the explicitly named login509-02-l via an internal SSH hop. No heavy process remained. Recorded total RSS was 63,124 KiB including inspection sessions; this is a current process snapshot, not historical peak or cgroup memory. No unrelated process was killed.
- sacct independently confirms the four sensitivity-audit attempts ran on worker gpu510-32, not a login node. The successful job 53746868 recorded MaxRSS 677448K. H3's four active transfer jobs are also on compute nodes.
- Direct login helper execution is established by the agent's self-audit and is a plausible contributor. Exact causation of the 20 GB breach is unproven without administrator accounting. The agent report's asserted cgroup/kill mechanism and memory estimates are not accepted as facts.

## Fix and verification

Added require_slurm_worker.sh locally under reproducibility/darkmatter/round7_analyses and mirrored to the remote p7_analyses_clean root. It rejects login hostnames, requires a numeric Slurm job ID, queries actual RUNNING state, and checks that hostname belongs to the allocation before executing a command.

Actual login invocation exited 64 with `Refusing compute on login host: login509-02-l`. Slurm debug test 53746999 completed exit 0 on gpu510-32, with `Compute allocation verified: job=53746999 host=gpu510-32`.

Updated CLAUDE.md and the prepared Flash collection/Argon H4 review prompts: compilation, database extraction, bulk parsing/checksumming, and helper compute must use the guard inside Slurm; login access is for submission/status/bounded small-file inspection. Explicit tool thread limits must match allocated CPUs. This is a project execution guard, not a cluster-wide restriction: commands bypassing the wrapper are not intercepted.

Existing scientific outputs and queued jobs were not modified. H4 anomaly claims remain provisional pending factual corrections and Argon review. H3 reference retrieval continues (216/373 parts at 07:33 UTC). No public communication or system-wide configuration change was made.
