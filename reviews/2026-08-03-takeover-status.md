# Dark-matter/MAG manuscript takeover status

Date: 2026-08-03

> **Execution update, 2026-08-05.** The strict three-sample CoverM smoke did
> not stall: it correctly rejected numerical equivalence between historical
> profiles and remapping with the frozen prebuilt index. The checksum-gated
> rejection control is complete, and the predeclared coherent fallback is
> remapping all 573 samples de novo. At the 00:40 +03 snapshot, 280 profiles
> were complete, four tasks from array 49985116 were running, and validator
> 49985773 remained dependency-held. See
> `2026-08-05-post-defensome-remote-status.md`. The eight replacement COMEBin
> tasks in 49931327 also remain active and demonstrably progressing.

This supersedes the execution-status portions of the 2026-07-29 reviews. It
does not replace their scientific critiques.

## Closed blockers

| Blocker from the revision review | Current disposition | Evidence |
|---|---|---|
| No canonical sample/read/MAG universe | Closed | Checksum-verified `canonical_v1` on IBEX: 573 paired metagenomes, 551 metadata samples, 14,426 refined bins (4,225 HQ / 8,209 MQ / 1,992 LQ), and 5,438 HQ/MQ dRep representatives. |
| Misinterpreted GTDB marker-workflow counts and domain-unplaced bins | Closed for the current catalogue | GTDB-Tk R232 census over the 5,438 representatives has no domain-, phylum-, or class-unassigned representative. Eleven lack an order assignment, 142 a family, 1,343 a genus, and 3,755 a species assignment. These are operational rank gaps, not automatic new taxa. |
| Invalid 0.39-versus-0.29 dN/dS comparator | Closed by withdrawal and independent audit | `2026-08-03-dnds-comparator-withdrawal.md` and `2026-08-03-claude-opus-b2-b3-final-audit.md`; active text and figure retain only a bounded dark-family sensitivity, legacy scripts fail by default, the executable verifier passes 54/54 checks, and the final rendered-artifact Opus audit reports PASS. |
| Missing known-fold false-novelty calibration | Closed for the historical 68,527-model structural cohort | Blinded structure-defined controls searched against PDB current through 2026-07-30 and AFDB v6: 101/3,013 (3.35%, Wilson 95% CI 2.77--4.06%) are false dual no-hits at length >=100 aa and pLDDT >=0.70. Candidate-specific distribution matching remains a reporting caveat. |
| Unreproducible PICRUSt2 source/method provenance | Closed | PICRUSt2 2.4.1 source, interpreter, CLI, inputs, and serialization tolerances are checksum-bound; canonical `outputs.sha256` hash `45a5d23d5e7de0d13b380ed8453a93acc531b92830d921aac2b096cfc855ef67`. |
| Stalled COMEBin fork/MKL execution | Closed operationally | Replacement job 49931327 uses single-process data loading while retaining data, model, seed, batch size, and CPU training threads; tasks are actively progressing. |
| Annotation tasks 92--135 rejected after `5_PR` completion | Closed operationally; validation remains gated | `5_PR` converted 44 `48PRr1` catalogue entries from regular files to symlinks. Their resolved byte sizes and sampled SHA-256 values match the frozen inventory, but script 183 compared the symlink object's 96--98-byte size. The frozen script was preserved. Symlink-aware recovery 49963925 follows links for the size check and still requires the original per-FASTA SHA-256; validator 49963934 accepts only the two pinned script hashes and waits for both arrays. |

## Gated analyses still running

| Analysis | Current state | Manuscript consequence |
|---|---|---|
| Current MAG KO reconstruction / PICRUSt2 comparison | Annotation recovery 49963925 and validator 49963934 completed. The strict CoverM smoke 49959297 rejected mixing historical and newly mapped profiles; obsolete mixed-profile chunks 49959374, 49959377, and 49961508 were therefore cancelled. The frozen rejection control now gates a coherent 573-sample full remap through arrays 49985116, 49985130, and 49985772, with validator 49985773 held downstream. | Do not retain any historical C1 number or figure. Promote only the four predeclared arms after full 573-sample validation. |
| Current MGnify90 family search | All 16 shard records and merge 49919828 completed. The current reproducibility snapshot uses the named, versioned, thresholded MGnify90-relative result. | Retain reference-relative language; no global saturation claim. |
| Current AFDB-v6 recurrent-family 3Di search | Thirteen CPU shards 49919824; merge 49919825 held `afterok`. | This, PDB, and MGnify define the residual Orix folding queue. Orix must remain empty until the sequence/structure-reference gate is complete. |
| Novel-lineage recovery by co-assembly/rebinning | COMEBin replacement 49931327 is active; CheckM2/GTDB-Tk downstream is dependency-held. | Existing R232 representatives do not support a new phylum/class claim. Any higher-rank novelty headline depends on recovered HQ/MQ bins and subsequent marker phylogeny, contamination checks, and placement. |
| Virus-host network refresh | Array 49944250 and link/summarize jobs 49944251/49944252 completed. | Relevant chiefly to the mobilome paper; keep it out of MAG novelty denominators. |

## Stable results that can anchor the revision

- The current refined catalogue contains 14,426 bins and 5,438 HQ/MQ species
  representatives. The R232 novelty signal is presently strongest below class:
  11 order gaps, 142 family gaps, 1,343 genus gaps, and 3,755 species gaps.
- The 11 order-gap representatives are all routed through the focus-v3
  adjudication. Ten enter direct cellular groups; one presently fails the
  minimum-alignment gate. These are candidates for phylogenetic evaluation,
  not eleven claimed new orders.
- The current sequence-reference refresh is explicit about universes. In the
  corrected refined historical-no-hit cohort, 43,409/1,902,233 now hit UniRef50
  2026_02; 1,858,824 remain no-hit. Clustering gives 83,219 multi-member,
  >=100-aa families, of which 79,352 span multiple assembly samples. The
  separate MAG-Bakta historical-no-hit cohort must not be combined with this
  denominator.
- The shared MGnify gate contains 298,871 recurrent core family queries and
  1,056 separately labelled defense-neighborhood queries. The primary
  cross-unit core denominator is 287,841; defense queries are excluded from
  core dark-family rates.
- Current PDB/AFDB rescanning of the historical explicit-structure cohort finds
  608 dual current-reference no-hits among 1,875 models with length >=100 aa
  and pLDDT >=0.70. These are calibrated candidate folds, not validated new
  folds; experimental structure is required.
- The bounded B2 result is 22,767 recurrent, alignable dark families with
  median pairwise NG86 dN/dS 0.391, rising only to 0.499 at a >=30-member floor.
  It must not be generalized to all dark proteins or presented as a formal
  neutral/spurious-ORF test.

## Required promotion order

1. Finish and validate current CoverM plus annotation backfill; rebuild the MAG
   KO profiles and run the predeclared PICRUSt2 comparison arms.
2. Finish MGnify90 and AFDB-v6 merges; freeze the exact residual target
   manifest.
3. Fold only that residual manifest on Orix with active GPU-utilization
   monitoring, then rescan and cluster successful models against the blinded
   known-fold controls.
4. Finish co-assembly/rebinning, CheckM2, GTDB-Tk R232, GUNC/contamination
   checks, and marker phylogenies before making any recovered high-rank claim.
5. Promote compact checksum-closed outputs, revise all affected manuscript
   sections together, regenerate figures/claims, and run a final independent
   Opus audit of numbers, semantics, and provenance.
