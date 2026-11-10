# Stage 4 Raw (`MAPQ >= 0`) vs Filtered (`MAPQ >= 20`) Focal Viral Coverage QC (`exp3_viruses/coverage_qc`)

## 1. Purpose and Scope

This standalone reproducibility subpackage (`exp3_viruses/coverage_qc`) documents and verifies the bounded Stage 4 raw-versus-filtered (`MAPQ >= 0` vs `MAPQ >= 20`) focal viral coverage quality-control (QC) analysis across all **14 metagenomic samples × 7 unique focal viral targets (98 sample×target rows)**.

- **Single isolated variable**: `samtools coverage` minimum mapping quality threshold (`-q 0` on `<sample>_raw.bam` vs `-q 20` on `<sample>_mapq20.bam`) aligned against the unified competitive reference (`unified_competitive_reference.fna`, `3,173,693` unique contigs / `3,234,560` crosswalk records across `5,438` canonical MAGs [`2,485,945` MAG contigs], `all_plasmids.fna` [`609,160`], `coasm_viruses.fna` [`93,060`], `all_viruses.fna` [`24,871`], `votu_reps.fna` [`19,300`], `all_matched_primary_targets.fna` [`2,187`], and `curated_cellular_controls_verified.fna` [`37`]).
- **Held strictly constant**:
  - Primary alignment filter flags: `--ff 3844` (`UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY` excluded; `0xf04` = `3844`).
  - Reference header identifiers (`VIRUS_TARGET__<entity_id>`) and sequence SHA256s validated in `stage4_reference_validation_receipt.json` (Slurm job `15905724`).
  - Existing competitive BAMs (`<sample>_raw.bam` and `<sample>_mapq20.bam`) and primary `MAPQ >= 20` focal target evaluations (`<sample>_focal_targets_evaluation.tsv`) produced by `stage4_map_sample_competitive.py` (Slurm array `15905725` + `V60PRr2` `inStrain` recovery job `15910797`).
  - Primary Stage 4 `inStrain` estimability gates (`MAPQ >= 20`, `mean_depth >= 10.0x`, `horizontal_breadth_pct >= 50.0%`, `analysis_role == "primary"`, and `source_category_uniquely_viral == True`), which remain **unchanged**.

---

## 2. Slurm Execution Provenance (Shaheen III)

| Slurm Job ID | Job Name | Stage | Partition | Node(s) | State | Exit Code | Elapsed | Window (AST, `+03:00`) | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `15905724` | `exp3_val_ref` | `stage4_validate_competitive_reference` | `workq` | `nid01239` | `COMPLETED` | `0:0` | `00:17:36` | `2026-11-09T20:38:01` – `20:55:37` | Validated `3,173,693` unique reference contigs and all 7 focal target SHA256s |
| `15905725` | `exp3_st4_map` | `stage4_map_sample_competitive_array` | `workq` | `nid[00021,00026-00027,00034,00051-00052,00060,00066,00074,00076,00078,00080,00088-00089]` | `13 completed / 1 failed` (task `1` recovered in `15910797`) | `0:0` (13 tasks), `1:0` (task `1`) | `01:11:23` | `2026-11-09T20:55:41` – `22:07:04` | Produced all 14 `<sample>_raw.bam` & `<sample>_mapq20.bam` files and coverage TSVs; task `1` (`V60PRr2`) failed only at `inStrain` step 4 (`::` header split) |
| `15910797` | `exp3_st4_rec_v60` | `stage4_v60prr2_instrain_recovery` | `workq` | `nid00034` | `COMPLETED` | `0:0` | `00:11:15` | `2026-11-09T23:17:34` – `23:28:49` | Completed `V60PRr2` `inStrain` profile & `V60PRr2_mapping_receipt.json` from existing BAM |
| **`15910231`** | **`exp3_raw_cov_qc`** | **`stage4_raw_vs_mapq20_focal_coverage`** | **`workq`** | **`nid00021`** | **`COMPLETED`** | **`0:0`** | **`00:01:15`** | **`2026-11-09T23:28:54` – `23:30:09`** | **Extracted 98-row raw (`-q 0`) vs filtered (`-q 20`) focal coverage & verified parity against all 14 sample receipts** |

---

## 3. Summary of Actual 98-Row Results (Programmatically Generated from Output TSVs)

### 3.1 Classification Breakdown Across All 98 Sample×Target Rows

| `recruitment_qc_class` | All 98 Rows | Source-Sample Rows (12) | Non-Source Rows (86) | Definition |
| :--- | ---: | ---: | ---: | :--- |
| `primary_mapq20_10x_50pct_eligible` | **1** | **1** | 0 | Meets primary `MAPQ >= 20`, `>=10x` depth, `>=50%` breadth, `analysis_role == "primary"`, and `source_category_uniquely_viral == True` |
| `mapq20_10x_50pct_reached_excluded_by_role_or_ownership_gate` | **1** | **1** | 0 | Meets `MAPQ >= 20`, `>=10x` depth, `>=50%` breadth (`41Sr3` on `41Sr3::41Sr3_56779`), excluded from primary `inStrain` by `analysis_role == "sensitivity_unresolved_virality"` |
| `high_raw_support_10x_50pct_lost_at_mapq20` | **10** | **4** | 6 | Meets `>=10x` depth and `>=50%` breadth at `MAPQ >= 0`, but falls below `10x` / `50%` after `MAPQ >= 20` filtering (`non_estimability_reason = high_raw_support_lost_at_mapq20_exact_reference_attribution_unresolved`) |
| `low_total_read_support_with_mapq20_attrition` | **36** | **5** | 31 | Non-zero raw primary alignments (`raw_numreads > 0`), below `10x` / `50%` even before `MAPQ >= 20` filtering, with additional read loss at `MAPQ >= 20` |
| `measured_zero_assigned_reads` | **50** | **1** | 49 | Exactly `0` raw (`MAPQ >= 0`) and `0` filtered (`MAPQ >= 20`) primary alignments (`NA`/empty retention fractions; never coerced to `0`) |
| **Total** | **98** | **12** | **86** | **14 samples × 7 unique focal viral targets (`48` rows `raw_numreads > 0`, `50` rows `raw_numreads == 0`; `26` rows `filtered_numreads > 0`, `72` rows `filtered_numreads == 0`)** |

### 3.2 Per-Target Summary Across All 14 Samples (`stage4_raw_vs_mapq20_target_summary.tsv`)

| `entity_id` | `analysis_role` | `length_bp` | `host_samples` (`host_mags`) | `virus_samples` | `stage2_cluster_votu_rep_id` (`is_rep`) | `stage2_alias_votu_metadata` | Samples `raw==0` / `raw>0` / `filt>0` | `raw>=10x&50%` / `high_raw_lost` / `primary_elig` | Total Reads Across 14 Samples (`raw` / `filt` / `lost`) |
| :--- | :--- | ---: | :--- | :--- | :--- | :--- | :--- | :--- | ---: |
| `55Dr2::55Dr2_75333` | `primary` | 56,848 | `41PRr2` (`41PRr2__single_SemiBin_97`) | `55Dr2` | `55Dr2::55Dr2_75333` (`1`) | `55Dr2::55Dr2_75333:rep=55Dr2::55Dr2_75333:is_rep=1` | 5 / 9 / 5 | 4 / 4 / 0 | 578,291 / 2,249 / 576,042 |
| `coasm_48_PR::k127_4172870` | `primary` | 56,855 | `V48Dr1` (`V48Dr1__self_SemiBin_3215`) | `coasm_48_PR` | `coasm_48_PR::k127_4172870` (`1`) | `coasm_48_PR::k127_4172870:rep=coasm_48_PR::k127_4172870:is_rep=1` | 6 / 8 / 2 | 0 / 0 / 0 | 273 / 76 / 197 |
| `44Dr2::44Dr2_4002037` | `boundary_sensitive_provisional` | 67,924 | `43PRr3` (`43PRr3__self_SemiBin_469`) | `44Dr2` | `44Dr2::44Dr2_4002037` (`1`) | `44Dr2::44Dr2_4002037:rep=44Dr2::44Dr2_4002037:is_rep=1` | 8 / 6 / 3 | 1 / 1 / 0 | 11,992 / 3,155 / 8,837 |
| `V60PRr2::V60PRr2_715976` | `primary` | 38,593 | `41PRr2` (`41PRr2__single_SemiBin_236`) | `V60PRr2` | `V60PRr2::V60PRr2_715976` (`1`) | `V60PRr2::V60PRr2_715976:rep=V60PRr2::V60PRr2_715976:is_rep=1` | 13 / 1 / 1 | 1 / 0 / 1 | 39,333 / 38,727 / 606 |
| `9Sr2::9Sr2_284237` | `primary` | 56,873 | `41PRr2` (`41PRr2__single_SemiBin_97`) | `9Sr2` | `55Sr3::55Sr3_129891` (`0`) | `9Sr2::9Sr2_284237:rep=55Sr3::55Sr3_129891:is_rep=0` | 6 / 8 / 5 | 4 / 4 / 0 | 582,090 / 7,209 / 574,881 |
| `48Sr2::48Sr2_190170` | `primary` | 38,647 | `S49Dr1` (`S49Dr1__self_SemiBin_1106`) | `48Sr2` | `48PRr2::48PRr2_1185466` (`0`) | `48Sr2::48Sr2_190170:rep=48PRr2::48PRr2_1185466:is_rep=0;coasm_48_S::k127_1210458:rep=coasm_48_S::k127_1210458:is_rep=1` | 6 / 8 / 5 | 1 / 1 / 0 | 8,522 / 858 / 7,664 |
| `41Sr3::41Sr3_56779` | `sensitivity_unresolved_virality` | 196,115 | `41Sr3` (`41Sr3__self_SemiBin_108`) | `41Sr3` | `41Sr3::41Sr3_56779` (`1`) | `41Sr3::41Sr3_56779:rep=41Sr3::41Sr3_56779:is_rep=1` | 6 / 8 / 5 | 1 / 0 / 0 | 134,122 / 130,357 / 3,765 |

### 3.3 Exact 12 Source-Sample Rows (`is_host_source_sample == True` or `is_viral_or_alias_source_sample == True`)

*(Note: `coasm_48_PR::k127_4172870` has `virus_samples = "coasm_48_PR"`, so its single-sample source row in the 14-sample set is `V48Dr1` [`host_source`]; its full 14-sample profile including coassembly-associated samples is shown in Section 3.5.)*

| Row | `sample` (`sample_roles`) | `entity_id` | `analysis_role` | `length_bp` | `stage2_cluster_votu_rep_id` (`is_rep`) | Raw Reads (`-q 0`) | Raw Covbases | Raw Breadth (%) | Raw Samtools (%) | Raw Depth (x) | Raw BaseQ / MAPQ | `MAPQ>=20` Reads (`uniq`/`ambig`) | `MAPQ>=20` Covbases | `MAPQ>=20` Breadth (%) | `MAPQ>=20` Samtools (%) | `MAPQ>=20` Depth (x) | `MAPQ>=20` BaseQ / MAPQ | Lost Reads | Retained Read Fraction | `recruitment_qc_class` |
| ---: | :--- | :--- | :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: | :--- | ---: | :--- | :--- |
| 0 | `41PRr2` (`host_source`) | `55Dr2::55Dr2_75333` | `primary` | 56,848 | `55Dr2::55Dr2_75333` (`1`) | 15 | 1,385 | 2.436321% | 2.4363% | 0.039843 | 35.6 / 0.933 | 0 (0/0) | 0 | 0.000000% | 0.0000% | 0.000000 | 0.0 / 0.0 | 15 | 0.0000% (0.000000) | `low_total_read_support_with_mapq20_attrition` |
| 3 | `41PRr2` (`host_source`) | `V60PRr2::V60PRr2_715976` | `primary` | 38,593 | `V60PRr2::V60PRr2_715976` (`1`) | 0 | 0 | 0.000000% | 0.0000% | 0.000000 | 0.0 / 0.000 | 0 (0/0) | 0 | 0.000000% | 0.0000% | 0.000000 | 0.0 / 0.0 | 0 | NA | `measured_zero_assigned_reads` |
| 4 | `41PRr2` (`host_source`) | `9Sr2::9Sr2_284237` | `primary` | 56,873 | `55Sr3::55Sr3_129891` (`0`) | 9 | 1,011 | 1.777645% | 1.7776% | 0.023737 | 36.7 / 0.444 | 0 (0/0) | 0 | 0.000000% | 0.0000% | 0.000000 | 0.0 / 0.0 | 9 | 0.0000% (0.000000) | `low_total_read_support_with_mapq20_attrition` |
| 10 | `V60PRr2` (`viral_or_alias_source`) | `V60PRr2::V60PRr2_715976` | `primary` | 38,593 | `V60PRr2::V60PRr2_715976` (`1`) | 39,333 | 38,593 | 100.000000% | 100.0000% | 140.550000 | 39.4 / 41.100 | 38,727 (38,727/0) | 38,593 | 100.000000% | 100.0000% | 138.633000 | 39.4 / 41.7 | 606 | 98.4593% (0.984593) | `primary_mapq20_10x_50pct_eligible` |
| 14 | `55Dr2` (`viral_or_alias_source`) | `55Dr2::55Dr2_75333` | `primary` | 56,848 | `55Dr2::55Dr2_75333` (`1`) | 300,087 | 56,848 | 100.000000% | 100.0000% | 787.476000 | 35.9 / 1.150 | 1,266 (1,266/0) | 2,970 | 5.224458% | 5.2245% | 3.335240 | 35.9 / 32.8 | 298,821 | 0.4219% (0.004219) | `high_raw_support_10x_50pct_lost_at_mapq20` |
| 25 | `9Sr2` (`viral_or_alias_source`) | `9Sr2::9Sr2_284237` | `primary` | 56,873 | `55Sr3::55Sr3_129891` (`0`) | 12,065 | 56,860 | 99.977142% | 99.9771% | 30.919000 | 35.6 / 1.510 | 181 (181/0) | 2,124 | 3.734637% | 3.7346% | 0.469854 | 35.7 / 33.6 | 11,884 | 1.5002% (0.015002) | `high_raw_support_10x_50pct_lost_at_mapq20` |
| 40 | `S49Dr1` (`host_source`) | `48Sr2::48Sr2_190170` | `primary` | 38,647 | `48PRr2::48PRr2_1185466` (`0`) | 5,813 | 32,834 | 84.958729% | 84.9587% | 14.462100 | 34.8 / 2.620 | 206 (206/0) | 4,124 | 10.670945% | 10.6709% | 0.518591 | 34.8 / 34.6 | 5,607 | 3.5438% (0.035438) | `high_raw_support_10x_50pct_lost_at_mapq20` |
| 47 | `48Sr2` (`viral_or_alias_source`) | `48Sr2::48Sr2_190170` | `primary` | 38,647 | `48PRr2::48PRr2_1185466` (`0`) | 1,779 | 37,621 | 97.345201% | 97.3452% | 6.924420 | 36.0 / 13.200 | 509 (509/0) | 11,640 | 30.118767% | 30.1188% | 1.981370 | 36.1 / 36.9 | 1,270 | 28.6116% (0.286116) | `low_total_read_support_with_mapq20_attrition` |
| 57 | `V48Dr1` (`host_source`) | `coasm_48_PR::k127_4172870` | `primary` | 56,855 | `coasm_48_PR::k127_4172870` (`1`) | 42 | 1,592 | 2.800106% | 2.8001% | 0.100466 | 38.7 / 2.620 | 2 (2/0) | 101 | 0.177645% | 0.1776% | 0.003553 | 39.8 / 40.0 | 40 | 4.7619% (0.047619) | `low_total_read_support_with_mapq20_attrition` |
| 76 | `41Sr3` (`host_source;viral_or_alias_source`) | `41Sr3::41Sr3_56779` | `sensitivity_unresolved_virality` | 196,115 | `41Sr3::41Sr3_56779` (`1`) | 132,930 | 196,115 | 100.000000% | 100.0000% | 88.592100 | 39.5 / 40.600 | 130,330 (130,330/0) | 196,108 | 99.996431% | 99.9964% | 87.048100 | 39.5 / 41.3 | 2,600 | 98.0441% (0.980441) | `mapq20_10x_50pct_reached_excluded_by_role_or_ownership_gate` |
| 86 | `43PRr3` (`host_source`) | `44Dr2::44Dr2_4002037` | `boundary_sensitive_provisional` | 67,924 | `44Dr2::44Dr2_4002037` (`1`) | 2,923 | 52,887 | 77.862022% | 77.8620% | 6.234790 | 35.7 / 12.400 | 817 (0/817) | 11,343 | 16.699547% | 16.6995% | 1.749880 | 35.5 / 39.9 | 2,106 | 27.9507% (0.279507) | `low_total_read_support_with_mapq20_attrition` |
| 93 | `44Dr2` (`viral_or_alias_source`) | `44Dr2::44Dr2_4002037` | `boundary_sensitive_provisional` | 67,924 | `44Dr2::44Dr2_4002037` (`1`) | 9,048 | 67,879 | 99.933749% | 99.9337% | 19.689800 | 35.8 / 12.100 | 2,337 (0/2,337) | 20,510 | 30.195513% | 30.1955% | 5.067340 | 35.7 / 39.5 | 6,711 | 25.8289% (0.258289) | `high_raw_support_10x_50pct_lost_at_mapq20` |

### 3.4 Exact 6 Non-Source High-Raw-Support Rows Lost at `MAPQ >= 20`

| Row | `sample` (`sample_roles`) | `entity_id` | `length_bp` | Raw Reads (`-q 0`) | Raw Covbases | Raw Breadth (%) | Raw Depth (x) | Raw BaseQ / MAPQ | `MAPQ>=20` Reads (`-q 20`) | `MAPQ>=20` Covbases | `MAPQ>=20` Breadth (%) | `MAPQ>=20` Depth (x) | `MAPQ>=20` BaseQ / MAPQ | Lost Reads | Retained Read Fraction |
| ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | :--- | ---: | :--- |
| 18 | `55Dr2` (`predetermined_coassembly_alias_or_copositive_sample`) | `9Sr2::9Sr2_284237` | 56,873 | 300,866 | 56,873 | 100.000000% | 789.383000 | 35.9 / 1.440 | 3,865 | 3,282 | 5.770752% | 10.158400 | 36.0 / 33.0 | 297,001 | 1.2846% (0.012846) |
| 21 | `9Sr2` (`predetermined_coassembly_alias_or_copositive_sample`) | `55Dr2::55Dr2_75333` | 56,848 | 11,763 | 56,818 | 99.947228% | 30.160100 | 35.7 / 1.120 | 39 | 1,250 | 2.198846% | 0.096943 | 35.8 / 32.3 | 11,724 | 0.3315% (0.003315) |
| 28 | `55Sr3` (`predetermined_coassembly_alias_or_copositive_sample`) | `55Dr2::55Dr2_75333` | 56,848 | 251,534 | 56,848 | 100.000000% | 658.141000 | 36.0 / 1.120 | 938 | 2,499 | 4.395933% | 2.470310 | 35.7 / 33.9 | 250,596 | 0.3729% (0.003729) |
| 32 | `55Sr3` (`predetermined_coassembly_alias_or_copositive_sample`) | `9Sr2::9Sr2_284237` | 56,873 | 254,298 | 56,873 | 100.000000% | 665.158000 | 36.0 / 1.420 | 3,147 | 2,318 | 4.075748% | 8.296010 | 36.0 / 32.9 | 251,151 | 1.2375% (0.012375) |
| 35 | `S49Dr1` (`predetermined_coassembly_alias_or_copositive_sample`) | `55Dr2::55Dr2_75333` | 56,848 | 13,423 | 55,123 | 96.965592% | 22.463500 | 34.7 / 0.820 | 5 | 505 | 0.888334% | 0.008883 | 32.2 / 23.6 | 13,418 | 0.0372% (0.000372) |
| 39 | `S49Dr1` (`predetermined_coassembly_alias_or_copositive_sample`) | `9Sr2::9Sr2_284237` | 56,873 | 13,310 | 54,702 | 96.182723% | 22.288200 | 34.8 / 0.861 | 14 | 809 | 1.422468% | 0.023315 | 35.0 / 27.8 | 13,296 | 0.1052% (0.001052) |

### 3.5 Full 14-Sample Profile for Coassembly Target `coasm_48_PR::k127_4172870` (`length_bp = 56,855`)

| Row | `sample` | `sample_roles` | Raw Reads (`-q 0`) | Raw Covbases | Raw Breadth (%) | Raw Depth (x) | Raw BaseQ / MAPQ | `MAPQ>=20` Reads (`-q 20`) | `MAPQ>=20` Covbases | `MAPQ>=20` Breadth (%) | `MAPQ>=20` Depth (x) | Lost Reads | Retained Read Fraction | `recruitment_qc_class` |
| ---: | :--- | :--- | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: | :--- | :--- |
| 1 | `41PRr2` | `predetermined_coassembly_alias_or_copositive_sample` | 9 | 695 | 1.222408% | 0.023903 | 36.9 / 0.667 | 0 | 0 | 0.000000% | 0.000000 | 9 | 0.0000% (0.000000) | `low_total_read_support_with_mapq20_attrition` |
| 8 | `V60PRr2` | `predetermined_coassembly_alias_or_copositive_sample` | 0 | 0 | 0.000000% | 0.000000 | 0.0 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 0 | NA | `measured_zero_assigned_reads` |
| 15 | `55Dr2` | `predetermined_coassembly_alias_or_copositive_sample` | 2 | 194 | 0.341219% | 0.005312 | 33.6 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 2 | 0.0000% (0.000000) | `low_total_read_support_with_mapq20_attrition` |
| 22 | `9Sr2` | `predetermined_coassembly_alias_or_copositive_sample` | 4 | 206 | 0.362325% | 0.010623 | 36.6 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 4 | 0.0000% (0.000000) | `low_total_read_support_with_mapq20_attrition` |
| 29 | `55Sr3` | `predetermined_coassembly_alias_or_copositive_sample` | 3 | 187 | 0.328907% | 0.007968 | 36.7 / 1.000 | 0 | 0 | 0.000000% | 0.000000 | 3 | 0.0000% (0.000000) | `low_total_read_support_with_mapq20_attrition` |
| 36 | `S49Dr1` | `predetermined_coassembly_alias_or_copositive_sample` | 191 | 5,270 | 9.269194% | 0.314713 | 35.0 / 15.900 | 74 | 2,470 | 4.344385% | 0.122663 | 117 | 38.7435% (0.387435) | `low_total_read_support_with_mapq20_attrition` |
| 43 | `48Sr2` | `predetermined_coassembly_alias_or_copositive_sample` | 14 | 567 | 0.997274% | 0.036971 | 36.3 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 14 | 0.0000% (0.000000) | `low_total_read_support_with_mapq20_attrition` |
| 50 | `F48Sr2` | `predetermined_coassembly_alias_or_copositive_sample` | 0 | 0 | 0.000000% | 0.000000 | 0.0 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 0 | NA | `measured_zero_assigned_reads` |
| 57 | `V48Dr1` | `host_source` | 42 | 1,592 | 2.800106% | 0.100466 | 38.7 / 2.620 | 2 | 101 | 0.177645% | 0.003553 | 40 | 4.7619% (0.047619) | `low_total_read_support_with_mapq20_attrition` |
| 64 | `S48PRr1` | `predetermined_coassembly_alias_or_copositive_sample` | 0 | 0 | 0.000000% | 0.000000 | 0.0 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 0 | NA | `measured_zero_assigned_reads` |
| 71 | `41Sr3` | `predetermined_coassembly_alias_or_copositive_sample` | 8 | 57 | 0.100255% | 0.008020 | 39.8 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 8 | 0.0000% (0.000000) | `low_total_read_support_with_mapq20_attrition` |
| 78 | `10Dr2` | `predetermined_unrelated_control` | 0 | 0 | 0.000000% | 0.000000 | 0.0 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 0 | NA | `measured_zero_assigned_reads` |
| 85 | `43PRr3` | `predetermined_coassembly_alias_or_copositive_sample` | 0 | 0 | 0.000000% | 0.000000 | 0.0 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 0 | NA | `measured_zero_assigned_reads` |
| 92 | `44Dr2` | `predetermined_coassembly_alias_or_copositive_sample` | 0 | 0 | 0.000000% | 0.000000 | 0.0 / 0.000 | 0 | 0 | 0.000000% | 0.000000 | 0 | NA | `measured_zero_assigned_reads` |

### 3.6 Per-Sample Summary Across the 7 Focal Targets (`stage4_raw_vs_mapq20_sample_summary.tsv`)

| `sample` | `slurm_job_id` (`task`) | `raw_total_records` | `raw_primary_mapped_records` | `filtered_high_quality_mapq20_reads` | `filtered_fraction_of_primary` | `focal7_raw_reads_sum` | `focal7_filtered_reads_sum` | `focal7_lost_to_mapq20_reads_sum` | Targets `raw>0` / `raw==0` | Targets `filt>0` | Targets `raw>=10x&50%` | Targets `primary_elig` |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: |
| `41PRr2` | `15908258` (`0`) | 150,592,898 | 82,245,448 | 66,923,270 | 81.3702% (0.813702) | 234 | 40 | 194 | 5 / 2 | 1 | 0 | 0 |
| `V60PRr2` | `15910797` (`not_array`) | 261,298,184 | 165,723,183 | 151,380,724 | 91.3455% (0.913455) | 39,336 | 38,727 | 609 | 2 / 5 | 1 | 1 | 1 |
| `55Dr2` | `15908260` (`2`) | 179,072,038 | 155,033,252 | 118,704,840 | 76.5673% (0.765673) | 601,031 | 5,134 | 595,897 | 5 / 2 | 3 | 2 | 0 |
| `9Sr2` | `15908261` (`3`) | 49,035,820 | 21,914,889 | 13,460,757 | 61.4229% (0.614229) | 23,833 | 220 | 23,613 | 4 / 3 | 2 | 2 | 0 |
| `55Sr3` | `15908262` (`4`) | 161,027,454 | 143,286,072 | 108,918,645 | 76.0148% (0.760148) | 505,904 | 4,086 | 501,818 | 5 / 2 | 3 | 2 | 0 |
| `S49Dr1` | `15908263` (`5`) | 145,943,449 | 64,904,990 | 32,228,172 | 49.6544% (0.496544) | 33,608 | 316 | 33,292 | 5 / 2 | 5 | 3 | 0 |
| `48Sr2` | `15908264` (`6`) | 153,040,178 | 87,781,015 | 62,335,890 | 71.0130% (0.710130) | 3,969 | 512 | 3,457 | 5 / 2 | 3 | 0 | 0 |
| `F48Sr2` | `15908265` (`7`) | 50,542,731 | 32,249,108 | 19,791,999 | 61.3722% (0.613722) | 3 | 0 | 3 | 2 / 5 | 0 | 0 | 0 |
| `V48Dr1` | `15908266` (`8`) | 172,235,818 | 83,510,348 | 63,715,238 | 76.2962% (0.762962) | 1,733 | 109 | 1,624 | 5 / 2 | 3 | 0 | 0 |
| `S48PRr1` | `15908267` (`9`) | 127,841,770 | 49,687,631 | 26,386,402 | 53.1046% (0.531046) | 4 | 0 | 4 | 1 / 6 | 0 | 0 | 0 |
| `41Sr3` | `15908268` (`10`) | 176,762,607 | 122,451,512 | 96,606,924 | 78.8940% (0.788940) | 132,982 | 130,332 | 2,650 | 5 / 2 | 2 | 1 | 0 |
| `10Dr2` | `15908269` (`11`) | 165,457,150 | 52,103,711 | 30,116,021 | 57.8001% (0.578001) | 13 | 1 | 12 | 1 / 6 | 1 | 0 | 0 |
| `43PRr3` | `15908270` (`12`) | 112,957,531 | 55,875,880 | 42,046,005 | 75.2489% (0.752489) | 2,925 | 817 | 2,108 | 2 / 5 | 1 | 0 | 0 |
| `44Dr2` | `15905725` (`13`) | 257,849,373 | 138,133,511 | 96,926,926 | 70.1690% (0.701690) | 9,048 | 2,337 | 6,711 | 1 / 6 | 1 | 1 | 0 |

---

## 4. Scientific Interpretation and Methodological Limits

1. **Authoritative Primary Gates Preserved**:
   - The primary `inStrain` microdiversity gates (`MAPQ >= 20`, `mean_depth >= 10x`, `horizontal_breadth_pct >= 50%`, `analysis_role == "primary"`, `source_category_uniquely_viral == True`) remain the authoritative standard for SNP/microdiversity calling.
2. **Non-Causal Attribution Wording**:
   - For rows in `high_raw_support_10x_50pct_lost_at_mapq20`, `non_estimability_reason` is recorded as `high_raw_support_lost_at_mapq20_exact_reference_attribution_unresolved` with qualification: *"High raw primary-alignment support (>=10x, >=50% breadth at MAPQ>=0) collapses below 10x/50% at MAPQ>=20 (exact-reference attribution unresolved; compatible with competitive ambiguity or weak alignment discrimination)."*
3. **Stage 2 vOTU Context Metadata (Descriptive Only)**:
   - Stage 2 vOTU clustering metadata (`stage2_target_viral_identity_table.tsv`) is joined onto each target (`stage2_cluster_votu_rep_id`, `stage2_is_votu_cluster_rep`, `stage2_alias_votu_metadata`) with `votu_metadata_role = "context_only_not_demonstrated_cause_of_ambiguity; is_votu_cluster_rep=0_is_normal_nonrepresentative_member"`.

---

## 5. Verification and Replay Instructions

### 5.1 Run Standalone 7-Gate Verifier (`verify_coverage_qc.py`)

```bash
uv run python verify_coverage_qc.py
```

### 5.2 Run Unit & Package Test Suite (`pytest`)

```bash
mkdir -p results/.pytest_tmp && uv run pytest -q --basetemp=results/.pytest_tmp tests/test_stage4_raw_vs_mapq20.py tests/test_coverage_qc_package.py && rm -rf results/.pytest_tmp
```

### 5.3 Full End-to-End Bulk-BAM Replay on Shaheen III

Coordinate-sorted competitive BAMs (`<sample>_raw.bam` and `<sample>_mapq20.bam`) reside on Shaheen III under `/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses/results/stage4_competitive_mapping/<sample>/`. To replay raw and filtered `samtools coverage` extraction directly from the BAMs on Shaheen III:

```bash
cd /scratch/project/k10018/eq-mechanism-20261109/exp3_viruses
sbatch code/sbatch_stage4_raw_vs_mapq20.sh
```
