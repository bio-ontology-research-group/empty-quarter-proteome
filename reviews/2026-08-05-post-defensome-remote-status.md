# Post-defensome dark-matter/MAG execution status

Snapshot: 2026-08-05 00:40 +03
Mode: read-only IBEX inspection; no local heavy computation and no job changes

## Interpretation

The remaining MAG workflows are active or dependency-gated, not silently
stalled. No new phylum, class, or order is currently supported. The existing
5,438 ANI95 representative catalogue has no domain-, phylum-, or class-level
GTDB-Tk R232 assignment gap. A higher-rank novelty claim remains conditional on
the coassembly/rebinning output passing CheckM2, GUNC/contamination checks,
GTDB-Tk R232, marker phylogeny, and placement review.

## Coassembly/rebinning novelty gate

Replacement COMEBin array 49931327 has eight active tasks: 1, 2, 3, 5, 6, 7,
8, and 10. Each had run for approximately 35.5 hours at the snapshot. Their
logs were being updated at inspection time, the jobs were advancing through
representation batches, and their batch steps showed substantial CPU time and
approximately 64--67 GB resident memory. They are therefore computationally
slow but not stalled.

The earlier coassembly branch still has task 49740175_9 active with its GTDB-Tk
jobs dependency-held. Neither branch has yet produced the complete, validated
novel-lineage result needed for manuscript promotion. No downstream novelty
number should be inferred from partial bins.

## CoverM/PICRUSt2 reconstruction

The three declared smoke samples in job 49959297 were recomputed against the
frozen 5,438-representative prebuilt minimap2 index. All three failed the
predeclared 1e-6-percentage-point equivalence criterion relative to historical
profiles:

| Smoke sample task | Maximum genome difference | Values above tolerance | Unmapped difference |
|---:|---:|---:|---:|
| 1 | 0.0532089650 | 1,724 | 0.2254240 |
| 2 | 0.0075929426 | 1,436 | 0.0547950 |
| 3 | 0.0221290000 | 1,476 | 0.1019360 |

This is a successful safety decision, not a pipeline crash: historical and
new profiles must not be mixed. The checksum-gated
`index_equivalence_rejection/EQUIVALENCE_REJECTED` marker is present. Scripts
192--194 implement the prespecified fallback by recomputing every sample with
one frozen index and command, then validating the full cohort.

At the snapshot:

- 280/573 atomic `PROFILE_COMPLETE` bundles existed;
- parent array 49985116 had four active tasks and the remainder limited by its
  `%4` concurrency cap;
- the completed 501--545 and 546--573 array groups were recorded under
  49985130 and 49985772;
- full-cohort validator 49985773 was pending on the remaining remaps;
- obsolete mixed-profile jobs 49959374, 49959377, and 49961508 were cancelled,
  as required after smoke rejection.

The current manuscript's historical C1 values remain explicitly provisional.
Only the four predeclared comparison arms may be promoted after validator
49985773 closes all 573 profiles.

## Other reference and host branches

- MGnify90: all 16 shard records in 49919827 and merge 49919828 completed; the
  bounded current result is already covered by the public reproducibility
  snapshot.
- Virus--host refresh: array 49944250 and downstream jobs 49944251/49944252
  completed; this does not enter the MAG novelty denominator.
- Broad AFDB-v6 CPU search: two shards in 49919824 remain active and merge
  49919825 is dependency-held. Residual large-scale folding on Orix should wait
  for this gate, so only true reference survivors consume GPUs.
- BORG host-map repair: task 50011654_244 remains active; aggregate 50011655 is
  dependency-held. This optional screen is not a blocker for either the frozen
  defensome result or the higher-rank MAG novelty gate.

## Next promotion sequence

1. Let the coherent 573-sample CoverM recomputation and validator finish;
   freeze its hashes and rebuild the four predeclared PICRUSt2 arms.
2. Let the eight COMEBin tasks finish; run the already gated CheckM2,
   GTDB-Tk/GUNC, dereplication, and marker-placement sequence before reading a
   novelty result.
3. Complete the AFDB merge, define the residual structure manifest, and only
   then use Orix GPUs for remaining folds.
4. Integrate compact checksum-closed outputs into the MAG manuscript and rerun
   the complete public harness plus an independent Opus audit.
