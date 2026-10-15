# H3 Final Raw-Evidence Independent Audit Report

- **Audit Date**: 2026-10-14T18:56:06Z
- **Slurm Job ID**: 53796328
- **Host**: gpu510-32
- **Verdict**: **`PASS`** (All 7 matched pairs, all 16 .m8 files, all signatures, and full census match)

## 1. Executive Summary

An independent full raw-evidence audit was conducted on the complete Hypothesis 3 (H3) reference-sensitivity analysis artifacts. The audit streamed and verified all raw FAA files, DIAMOND databases, search hit files, and query census sets without invoking production scoring routines. All findings were verified against checksum completion metadata and matched production summary tables with **zero discrepancies**.

## 2. Query Census & Query FASTA Audit

- **Dark Queries (`FAMILIES_FINAL.txt`)**: 106,042 IDs (SHA256: `d3f082f9da305fe11c7ac99f5c7570a45d5066fcca191326c4a5eb019b4b5331`)
- **Control Queries (`CONTROL_FAMILIES.txt`)**: 106,042 IDs (SHA256: `59c572a8a942336bd3117a95994ba66102d370267f94c952f5a671f7f46e9dc1`)
- **Disjointness**: Verified 100% disjoint arms (0 overlapping IDs; combined unique census = 212,084)
- **Query FASTA (`query_all.faa`)**: 212,084 sequences (SHA256: `1b0cda1498dfb08cb56ca224d3c3eebad058e57aeaf1113564ca65e0a087095c`)
- **Query Sequence Integrity**: All 212,084 sequences have non-zero lengths. All FASTA IDs match census IDs exactly with zero extra or missing IDs.

## 3. Matched Reference FAA Pairs Audit (7 Pairs)

| Pair Label | Matched Proteins | Total Residues | 1aa Hist Match | Within-Arm Unique | Shared Seqs | Zero Overlap | Closure Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `harmonized_seed20261014` | 4,210,592 | 1,295,970,055 | PASS | PASS | 22,310 | N/A | PASS |
| `harmonized_seed20261014_no_overlap` | 4,185,525 | 1,289,129,312 | PASS | PASS | 0 | PASS | PASS |
| `harmonized_seed20261015` | 4,143,929 | 1,280,798,499 | PASS | PASS | 17,874 | N/A | PASS |
| `harmonized_seed20261016` | 4,233,572 | 1,302,039,930 | PASS | PASS | 16,373 | N/A | PASS |
| `published_seed20261014` | 4,216,518 | 1,295,529,211 | PASS | PASS | 21,356 | N/A | PASS |
| `published_seed20261015` | 4,156,075 | 1,279,293,753 | PASS | PASS | 17,299 | N/A | PASS |
| `published_seed20261016` | 4,237,831 | 1,301,236,416 | PASS | PASS | 15,560 | N/A | PASS |

## 4. DIAMOND Search Hits Audit (16 .m8 Files)

| Search Run | Cutoff | SMAG Hits M8 SHA256 | SMAG Rows | GEM Hits M8 SHA256 | GEM Rows | Constraints & Signatures |
| :--- | :---: | :--- | :---: | :--- | :---: | :---: |
| `harmonized_seed20261014` | 0 | `30642892286dc0b4...` | 129,737 | `b13914ea9fd3fa22...` | 106,784 | PASS |
| `harmonized_seed20261014_id30` | 30 | `8c3376ad42f1d8a3...` | 128,524 | `518bd56254010eed...` | 105,945 | PASS |
| `harmonized_seed20261014_no_overlap` | 0 | `afe994ac0ff64b97...` | 129,519 | `65bae193dafa2448...` | 106,285 | PASS |
| `harmonized_seed20261015` | 0 | `27ca2dcd4fa1839e...` | 129,455 | `a3f402c9a36c232f...` | 107,957 | PASS |
| `harmonized_seed20261016` | 0 | `33972c42cb94b000...` | 129,742 | `a9bb44bbd23619ab...` | 108,299 | PASS |
| `published_seed20261014` | 0 | `c9255243fc516d7a...` | 129,796 | `d926cdd21473a90e...` | 106,959 | PASS |
| `published_seed20261015` | 0 | `6b5e1350ab84c9f4...` | 129,594 | `a4e2abd0bc42358a...` | 108,073 | PASS |
| `published_seed20261016` | 0 | `e1b95e9e125ac635...` | 129,787 | `4a73f59b4ac833e0...` | 108,422 | PASS |

## 5. Census Contingency & McNemar Replay

All 16 primary census cells (Dark and Control cohorts across 8 search configurations) and all 80 length-stratified strata cells were independently recomputed and compared against production outputs:

| Run Label | Cohort | Total N | Both (a) | SMAG Only (b) | GEM Only (c) | Neither (d) | SMAG % | GEM % | Delta % | Disc Ratio | p-value |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `harmonized_seed20261014` | Dark | 106,042 | 8,015 | 22,049 | 2,384 | 73,594 | 28.35% | 9.81% | +18.54% | 9.25 | `0.000e+00` |
| `harmonized_seed20261014` | Control | 106,042 | 94,277 | 5,396 | 2,108 | 4,261 | 93.99% | 90.89% | +3.10% | 2.56 | `0.000e+00` |
| `harmonized_seed20261014_id30` | Dark | 106,042 | 7,409 | 21,491 | 2,251 | 74,891 | 27.25% | 9.11% | +18.14% | 9.55 | `0.000e+00` |
| `harmonized_seed20261014_id30` | Control | 106,042 | 94,159 | 5,465 | 2,126 | 4,292 | 93.95% | 90.8% | +3.15% | 2.57 | `0.000e+00` |
| `harmonized_seed20261014_no_overlap` | Dark | 106,042 | 7,786 | 22,126 | 2,382 | 73,748 | 28.21% | 9.59% | +18.62% | 9.29 | `0.000e+00` |
| `harmonized_seed20261014_no_overlap` | Control | 106,042 | 94,007 | 5,600 | 2,110 | 4,325 | 93.93% | 90.64% | +3.29% | 2.65 | `0.000e+00` |
| `harmonized_seed20261015` | Dark | 106,042 | 9,000 | 20,823 | 2,972 | 73,247 | 28.12% | 11.29% | +16.83% | 7.01 | `0.000e+00` |
| `harmonized_seed20261015` | Control | 106,042 | 94,014 | 5,618 | 1,971 | 4,439 | 93.96% | 90.52% | +3.44% | 2.85 | `0.000e+00` |
| `harmonized_seed20261016` | Dark | 106,042 | 9,022 | 21,027 | 3,090 | 72,903 | 28.34% | 11.42% | +16.91% | 6.8 | `0.000e+00` |
| `harmonized_seed20261016` | Control | 106,042 | 94,098 | 5,595 | 2,089 | 4,260 | 94.01% | 90.71% | +3.31% | 2.68 | `0.000e+00` |
| `published_seed20261014` | Dark | 106,042 | 8,111 | 21,983 | 2,412 | 73,536 | 28.38% | 9.92% | +18.46% | 9.11 | `0.000e+00` |
| `published_seed20261014` | Control | 106,042 | 94,315 | 5,387 | 2,121 | 4,219 | 94.02% | 90.94% | +3.08% | 2.54 | `2.964e-321` |
| `published_seed20261015` | Dark | 106,042 | 9,046 | 20,870 | 2,982 | 73,144 | 28.21% | 11.34% | +16.87% | 7.0 | `0.000e+00` |
| `published_seed20261015` | Control | 106,042 | 94,088 | 5,590 | 1,957 | 4,407 | 94.0% | 90.57% | +3.43% | 2.86 | `0.000e+00` |
| `published_seed20261016` | Dark | 106,042 | 9,047 | 21,029 | 3,094 | 72,872 | 28.36% | 11.45% | +16.91% | 6.8 | `0.000e+00` |
| `published_seed20261016` | Control | 106,042 | 94,170 | 5,541 | 2,111 | 4,220 | 94.03% | 90.8% | +3.23% | 2.62 | `0.000e+00` |

## 6. Statistical Artifact vs Biological Effect Guidance

- **Binomial Underflow Artifact**: In the primary dark family comparison, discordant pairs are heavily asymmetric (b = 22,049 SMAG-only hits vs c = 2,384 GEM-only hits). The exact binomial test p-value evaluates to ~10^{-324} or smaller, underflowing IEEE-754 double precision float64 to `0.000e+00`. This underflow should be reported descriptively as an artifact of massive discordant sample sizes rather than an anomalous probability.
- **Independent Family Caveat**: The 212,084 query proteins represent dark and control orthogroup representatives. Because queries are biological sequences that may share homology or functional constraints across families, the query-level p-value evaluates sequence-level search sensitivity rather than unlinked independent evolutionary inferences. The descriptive effect sizes — namely the **+18.54% paired gain in dark family sensitivity** and the **9.25x discordant ratio** — provide the appropriate robust summary metrics for manuscript interpretation.

## 7. Scope & Invariants Verified

- [x] Exactly 106,042 unique dark query IDs and 106,042 unique control query IDs (212,084 total, strictly disjoint)
- [x] Query FASTA contains exactly 212,084 positive-length sequences whose IDs match the census union with zero discrepancies
- [x] All seven matched reference FAA pairs have 100% identical exact 1aa length histograms between SMAG and GEM arms
- [x] All seven matched reference FAA pairs have zero duplicate sequences within arms (strict sequence uniqueness)
- [x] Cross-arm sequence overlap matches matching_closure_audit.tsv exactly; zero-overlap is empirically verified for harmonized_seed20261014_no_overlap (0 shared sequences)
- [x] All 14 matched reference FAAs have verified SHA256 hashes matching completion signatures in .dmnd.sig.json and .m8.sig.json
- [x] All 14 DIAMOND .dmnd databases have verified file SHA256 hashes matching .dmnd.sig.json signatures with exit code 0
- [x] All 16 .m8 search files have verified file SHA256 hashes matching .m8.sig.json signatures with exit code 0 and valid parameters
- [x] All rows across all 16 .m8 search hit files conform strictly to 4 columns, finite numeric values, evalue <= 1e-3, and identity cutoff (>= 30 for id30)
- [x] Every query ID in every .m8 file is a strict member of the declared 212,084 census
- [x] Every target ID in every .m8 file is a strict member of the corresponding matched reference FAA
- [x] Every cell in matched_full_census_summary.tsv and matched_length_stratified_results.tsv was independently recomputed from raw .m8 hits and validated with 0 discrepancies
- [x] Preserved binomial test p-value underflow (p=0.000e+00) as an artifact of large discordant sample sizes (e.g., 22,049 vs 2,384 discordant queries); recommend focusing on descriptive effect sizes (discordant ratio ~9.25x, paired delta +18.54%) rather than treating query census p-values as independent family inferences.

