# Genome-Scale Metabolic Network Modelling & Gapfill Ensemble Analysis

## 1. Overview & Experimental Design

Experiment 1 reconstructs, gapfills, and evaluates draft genome-scale metabolic networks across **484 deduplicated genomes** from the Rub' al Khali (Empty Quarter) desert (`EQ_desert`, $n = 237$) and external reference cohorts (`external_comparator`, $n = 247$: 218 GTDB r232 species representatives + 29 MGnify reference genomes). The analysis investigates simulated nutrient dependencies, biomass-template sensitivity, and candidate functional deficiencies in hyper-arid microbial communities under explicit numerical and epistemic boundaries.

The computational workflow consists of five stages:
1. **Model Inventory, Pool Audit & Biomass Template Audit (Stage 1)**: Curation of 484 deduplicated draft SBML models (`exp1_model_inventory_and_qc_484.tsv`, filtered from 538 total candidate archives in `exp1_model_inventory_and_qc_538.tsv`). Audit of the 2,271-reaction cohort-derived candidate pool (`exp1_candidate_reaction_pool_audit.tsv`) confirming 0 reaction ID conflicts (`exp1_candidate_reaction_pool_conflicts.tsv`) and 0 species conflicts (`exp1_candidate_species_conflicts.tsv`) while preserving explicit `unknown_formula_or_charge` ($n = 228$) and `elemental_imbalance` ($n = 1$, excluded) statuses. Audit of the gapseq archaeal (*Methanosarcina barkeri* iAF692) draft biomass template across all 484 models (`exp1_biomass_template_audit_484.tsv`).
2. **Historical Stage 2 Pilot24 & Witness Revalidation (Stage 2)**: Initial validation on a 31-genome pilot (`exp1_pilot24_qc_gate.json`, 806 evaluations) and independent Slurm revalidation of four historical witness models (`focused_release_validation/exp1_four_changed_witnesses_verification.json`). *Note*: The Stage 2 Pilot24 aggregate threshold ($\ge 95\%$ solver completion, 0 timeouts) is preserved strictly as an immutable historical Stage 2 record and is **not** the final scientific acceptance rule for the 484-genome production cohort.
3. **Reproducible Artificial Gene-Masking Controls (Stage 2 Controls, v4)**: Gene-incompleteness masking controls across the 7 external reference genomes in Stage 1 inventory meeting numeric `completeness >= 90%` and `contamination <= 5%` at 10%, 20%, and 30% gene masking $\times$ 3 replicates (**63 evaluations**; Slurm job `15894236`), solved under the exact same `S4_oxic_desert_oa_aa` medium and `B2` biomass as unmasked baselines and independently validated post-hoc (`exp1_masking_controls_7hq_v4_validation_receipt.json`).
4. **Full 484-Genome Production Shard Array & Per-Solution Ledger Verification (Stage 3)**: Parallel execution across **8 Slurm array shards** (`--array=0-7`, Slurm array job `15893066_0..7`) on Shaheen III, evaluating **12,580 planned `(genome, scenario_id, biomass_recipe_id)` cases** (`484 \times 24` core combinations plus `482 \times 2` optional `S0_gapseq_native_predicted` combinations) with binary MILP alternative-gapfill enumeration (`n_Target_solutions = 20` on `B0`/`B2` and `10` on `B1`/`B3`/`B4`, within `25%` weighted cost slack, `6.0` s per MILP solve) and closed-uptake energy-generating cycle (EGC) verification. Across all 8 shards, the full production ledger (`exp1_full484_per_solution_ledger.tsv.gz`, audited in `provenance/exp1_full484_independent_table_audit.json`) records **93,808 solution rows** and **91,806 `accepted_feasible` solutions** across **10,578 feasible ensembles**, alongside **1,936 fixed-medium infeasible cases**, **56 FBA verification-failed cases**, and **10 solver timeouts** (`56 + 10 = 66` explicit unknown/failed cases retained with `NaN` summary metrics). Raw flux vectors are not persisted (`raw_flux_replay_not_claimed`).
5. **Fail-Closed Collection & Candidate Deficiency Homology Search (Stages 3 & 4, v5)**: Versioned collection (`stage3_collect_and_summarize_v4.py`, Slurm job `15894446`) joining per-case solver status and valid solution counts while retaining infeasible and unknown/failed cases strictly as `NaN`. Stage 4 v5 DIAMOND `blastx` homology search (`stage4_candidate_deficiency_search_v5.py`, Slurm job `15895063`) and bounded post-hoc qualification (`posthoc_qualify_stage4_v5_homology_dispositions.py`) evaluating 63 candidate reactions across **3,713 unique `(genome, clean_rxn_id)` pairs** (collapsed from 13,750 scenario rows across 81 MAGs and 47 sample assemblies).

---

## 2. Cohort Composition, Quality Semantics & Candidate Pool Integrity

### Cohort Composition & Quality Boundaries (`exp1_model_inventory_and_qc_484.tsv`)
The evaluation universe consists of **484 unique deduplicated genomes** with **heterogeneous assembly completeness and contamination** (not all 484 genomes or all 247 external comparators are high-quality):
- **Rub' al Khali Desert MAGs (`cohort = EQ_desert`, `source = EQ_MAG`, $n = 237$)**: Uncultivated metagenome-assembled genomes from Rub' al Khali dune, interdune, and sabkha samples. Confirmed **CheckM2** completeness spans `51.11%`–`100.00%` (median `85.42%`, mean `81.63%`) and contamination spans `0.00%`–`9.86%` (median `4.20%`). While `59/237` desert MAGs meet numeric `completeness >= 90%` and `contamination <= 5%`, all 237 carry `quality_tier = UNKNOWN` in Stage 1 inventory because numeric completeness/contamination thresholds alone do not establish formal MIMAG high-quality certification (which additionally requires 5S, 16S, and 23S rRNA genes and $\ge 18$ tRNAs).
- **External Reference Comparators (`cohort = external_comparator`, $n = 247$)**: Comprises **218 GTDB r232 species representatives (`source = GTDB_r232_rep`)** and **29 MGnify reference genomes (`source = MGnify`)** forming a broader **comparative *Thermoplasmatota* reference cohort** (the 247 external genomes do **not** all share taxonomic orders with the `EQ_desert` MAGs). In the Stage 1 inventory table (`exp1_model_inventory_and_qc_484.tsv`), the 29 MGnify genomes carry populated numeric completeness/contamination values (`51.98%`–`99.19%` completeness; `7/29` meet numeric `completeness >= 90%` and `contamination <= 5%` and are labeled `quality_tier = HQ` as the 7 eligible masking-control references, while `22/29` are `MQ_or_unknown`), but the upstream quality-estimation method for those 29 MGnify genomes is **unknown** (despite legacy upstream column naming). The 218 GTDB r232 representatives carry `quality_tier = MQ_or_unknown` in the Stage 1 inventory table, with their known-method quality covariates and shared-order subset reconciled in the separately owned `covariates/` subpackage (defining the **455-genome known-method quality universe** of `237` `EQ_desert` + `218` `GTDB_r232_rep`, within which the **shared-order subset comprises 293 genomes**: `234` `EQ_desert` + `59` `GTDB_r232_rep`). Even for the 7 masking-control reference genomes (`MGYG000519815`, `MGYG000520046`, `MGYG000519748`, `MGYG000510876`, `MGYG000495096`, `MGYG000526805`, `MGYG000523080`), `completeness >= 90%` and `contamination <= 5%` is a numeric selection filter and both formal MIMAG rRNA/tRNA status and upstream quality-estimation method are unknown.

### Candidate Reaction Pool Integrity (`exp1_candidate_reaction_pool_audit.tsv`)
- **Pool Scope (`observed_cohort_538_draft_and_gapseq_gapfilled_reactions_not_full_gapseq_db`)**: 2,271 unique reactions collected across the 538 cohort draft and gapseq-filled models (`1,564` observed in cohort draft models; `707` observed only in cohort gapseq-filled models).
- **Conflict Files**: **0 reaction stoichiometry/compartment conflicts** (`exp1_candidate_reaction_pool_conflicts.tsv`) and **0 species conflicts** (`exp1_candidate_species_conflicts.tsv`). Both conflict files are required by `verify_frozen_outputs.py` (fail-closed if missing).
- **Mass/Charge Balance & Eligibility**: Across the 2,271 reactions, `1,836` are `strict_atom_and_charge_balanced`, `228` have `unknown_formula_or_charge`, `206` are `boundary_or_biomass_exempt`, and `1` has a known `elemental_imbalance` (excluded from gapfilling, leaving **2,270 eligible candidate reactions**).

---

## 3. Biomass Perturbation Recipes (`B0`–`B4`) & Simulated Media Scenarios (`S0`–`S6`)

### Declared Component-Omission Biomass Perturbations (`exp1_sbml_highs_engine.py`, `exp1_biomass_template_audit_484.tsv`)
All 484 draft models use the gapseq archaeal biomass objective (`R_bio1`, `gapseq 2.1.0 dat/biomass/biomass_archaea.json`, derived from *Methanosarcina barkeri* iAF692, Feist et al. 2006). Because this template is a methanogen reference baseline rather than an experimentally validated *Thermoplasmatota* biomass composition, five cumulative component-omission recipes are evaluated:
- **`B0_methanosarcina_template_baseline`**: Unmodified gapseq archaeal draft biomass demanding all **8 *M. barkeri* template trace cofactors** (`n_all8_template_trace_cofactors_demanded = 8` in all 484 models; methanofuran `cpd02431`/`cpd15883` is absent from the gapseq template).
- **`B1_perturbation_omit_8_template_cofactors`**: Component-omission perturbation removing all **8 *M. barkeri* template trace cofactors**:
  - `cpd00649`: Coenzyme $F_{420}$-2
  - `cpd02246`: Coenzyme M
  - `cpd02817`: HTP (7-mercaptoheptanoylthreonine phosphate / coenzyme B)
  - `cpd03425`: Factor $F_{430}$
  - `cpd15862`: Coenzyme $F_{390}$ (adenosine)
  - `cpd15863`: Coenzyme $F_{390}$ (guanosine)
  - `cpd15884`: Tetrahydrosarcinapterin
  - `cpd33626`: Adenosyl factor III
- **`B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420`**: Component-omission perturbation removing the **7 non-$F_{420}$ *M. barkeri* template trace cofactors** (`cpd02246`, `cpd02817`, `cpd03425`, `cpd15862`, `cpd15863`, `cpd15884`, `cpd33626`) while retaining `cpd00649` (Coenzyme $F_{420}$-2, which is used across diverse non-methanogenic archaea and bacteria).
- **`B3_perturbation_omit_template_cofactors_and_surface_glycolipids`**: **Cumulative** omission removing all 8 template trace cofactors in `B1` **plus** `cpd03802` (*sym*-homospermidine), `cpd15880` (glucosaminyl archaetidyl-*myo*-inositol), and the 4 hydroxyphytanyl ether lipids `cpd15817` (`2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycero-1-phosphoethanolamine`), `cpd15818` (`2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycerol-1-phospho-3'-sn-glycerol`), `cpd15819` (`2-O-(3'-hydroxyl)phytanyl-3-O-phytanyl-sn-glycero-1-phospho-myo-inositol`), and `cpd15820` (`2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycero-1-phosphoserine`), while retaining the 4 core diether archaetidyl lipids (`cpd15858`–`cpd15861`) and putrescine (`cpd00118`).
- **`B4_perturbation_core_amino_nucleotide_energy_skeleton`**: **Cumulative** omission removing all 8 template trace cofactors in `B1`, all 11 archaeal ether lipids and polyamines (`cpd15817`–`cpd15820`, `cpd15858`–`cpd15861`, `cpd15880`, `cpd00118` putrescine, `cpd03802` *sym*-homospermidine), and `cpd00155` (glycogen), retaining DNA, RNA, protein (20 amino acids), inorganic ions, and universal redox/energy cofactors (ATP/ADP/AMP, NAD/NADH, NADP/NADPH, CoA, acetyl-CoA, succinyl-CoA, THF).

### Simulated Media Scenarios (`S0`–`S6`)
All media scenarios are **simulated chemical environments or gapseq transporter-predicted media** (not direct desert porewater measurements):
- **`S0_gapseq_native_predicted` (evaluated with `B0`, `B2`)**: Per-genome `<genome>-native_medium.csv` predicted by gapseq. Evaluated for **482/484 genomes**; omitted for the **2 genomes with empty native medium CSVs** (`GCA_016180755.1` and `S22PRr1__self_SemiBin_29`, verified in `provenance/exp1_full484_independent_table_audit.json`).
- **`S1_anoxic_minimal_h2_co2_nh4` & `S2_oxic_minimal_h2_co2_nh4` (evaluated with `B0`, `B2`)**: Simulated minimal autotrophic medium (`BASE_INORGANIC_SALTS`: $\text{H}_2\text{O}$ `cpd00001`, phosphate `cpd00009`, $\text{CO}_2$ `cpd00011`, $\text{NH}_4^+$ `cpd00013`, $\text{SO}_4^{2-}$ `cpd00048`, $\text{H}_2\text{S}$ `cpd00239`, $\text{Na}^+$ `cpd00971`, $\text{K}^+$ `cpd00205`, $\text{Cl}^-$ `cpd00099`, $\text{Mg}^{2+}$ `cpd00254`, $\text{Ca}^{2+}$ `cpd00063`, $\text{Fe}^{2+}$ `cpd10515`, $\text{Fe}^{3+}$ `cpd10516`, $\text{Co}^{2+}$ `cpd00149`, $\text{Ni}^{2+}$ `cpd00244`, $\text{Zn}^{2+}$ `cpd00034`, $\text{Cu}^{2+}$ `cpd00058`, $\text{Mn}^{2+}$ `cpd00030`, $\text{H}^+$ `cpd00067`) plus $\text{H}_2$ (`cpd11640`, $10.0$), without $\text{O}_2$ (`S1`, anoxic) or with $\text{O}_2$ (`cpd00007`, $10.0$; `S2`, oxic).
- **`S3_anoxic_desert_oa_aa` & `S4_oxic_desert_oa_aa` (evaluated with `B0`–`B4`)**: Constructed from each genome's own `<genome>-desert_medium.csv` where non-empty (482/484 genomes, with `cpd00007` $\text{O}_2$ removed in anoxic `S3` and set to $10.0$ in oxic `S4`), or from the fallback simulated desert organic-acid/amino-acid recipe (`BASE_INORGANIC_SALTS` + `DESERT_OA_AA_ADDITIONS`: $\text{H}_2$ `cpd11640`, acetate `cpd00029`, L-lactate `cpd00159`, pyruvate `cpd00020`, formate `cpd00047`, glycerol `cpd00100`, nitrate `cpd00209`, 20 standard amino acids, chorismate `cpd00216`, thiamin `cpd00305`, riboflavin `cpd00220`, $\pm \text{O}_2$) for the 2 genomes with empty `<genome>-desert_medium.csv` (`GCA_016180755.1` and `S22PRr1__self_SemiBin_29`).
- **`S5_anoxic_defined_aa_precursors` & `S6_oxic_defined_aa_precursors` (evaluated with `B0`–`B4`)**: Simulated defined medium (`BASE_INORGANIC_SALTS` + `DEFINED_20AA_AND_PRECURSORS`: all compounds in `DESERT_OA_AA_ADDITIONS` plus putrescine `cpd00118`, niacin `cpd00218`, pantothenate `cpd00644`, folate `cpd00393`, and pyridoxine `cpd00263`), without $\text{O}_2$ (`S5`, anoxic) or with $\text{O}_2$ (`cpd00007`, $10.0$; `S6`, oxic).

Together, the 24 core `(scenario_id, biomass_recipe_id)` combinations across all 484 genomes ($484 \times 24 = 11,616$) plus the 2 optional `S0` combinations across 482 genomes ($482 \times 2 = 964$) define the exact **12,580 planned production cases**.

---

## 4. Stage 3 Ensemble Enumeration Parameters & Full 484-Genome Accounting

### Executed Stage 3 Parameters (`sbatch_exp1_stage3_array.sh`, `stage3_full484_worker.py`, `stage2_pilot24_ensemble.py`)
- **Parallel Sharding**: 8 Slurm array tasks (`--array=0-7`, `--num-tasks 8`, `--workers 32`, `OMP_NUM_THREADS=1`) partitioning the 484 sorted genomes (`61` genomes in shards `0..3`, `60` genomes in shards `4..7`).
- **Target Distinct Solutions**: `--n-target-solutions 20`, setting `n_Target_solutions = 20` for primary recipes `B0` and `B2` and `max(5, 20 // 2) = 10` for perturbation recipes `B1`, `B3`, and `B4`.
- **Weighted Cost Window & Solver Bounds**: Candidate reaction addition weights are assigned from sequence/transporter evidence in `[0.005, 100.0]` (`25.0` for unannotated transport, `100.0` for unannotated internal reactions). Alternative gapfill solutions are enumerated via sequential binary no-good cuts within **`25%` weighted cost slack (`cost_slack_frac = 0.25`, $\le 1.25\times$ anchor cost)** at `min_growth = 0.01` $\text{h}^{-1}$, `mip_rel_gap_tol = 1e-4`, and **`time_limit_per_milp_s = 6.0` s per MILP solve** in Stage 3 (`8.0` s in Stage 2 Pilot24; `15.0` s per MILP solve and `75.0` s ensemble wall cap in Stage 2 v4 masking controls).
- **Per-Solution Numerical & Cycle Verification**: Every candidate MILP solution is independently rechecked by FBA and closed-uptake ATP cycle analysis (`detect_closed_uptake_egc`); acceptance as `accepted_feasible` requires `verified_growth >= 0.01 * 0.999` (`0.00999`), `stoichiometric_max_residual <= 1e-5`, `max_bound_violation <= 1e-5`, `residual_bounds_growth_verified == True`, `cycle_qc_status == "passed_closed_uptake_egc_check"`, and `spurious_atp_cycle == False`. Raw flux vectors are not persisted in the summary or ledger tables (`raw_flux_replay_not_claimed`).
- **Full 484-Genome Case & Solution Ledger Accounting (`exp1_full484_ensemble_summary.tsv.gz`, `exp1_full484_per_solution_ledger.tsv.gz`, `provenance/exp1_full484_independent_table_audit.json`)**:
  - **12,580 evaluated cases**:
    - `10,578` feasible ensembles (`solver_status = heuristic_incumbent_bounded`, `n_distinct_valid_solutions >= 1`)
    - `1,936` fixed-medium infeasible cases (`solver_status = infeasible_under_fixed_medium`, `n_distinct_valid_solutions = 0`, summary metrics strictly `NaN`)
    - `56` FBA verification-failed cases (`solver_status = fba_verification_failed`, `n_distinct_valid_solutions = 0`, summary metrics strictly `NaN`)
    - `10` solver timeouts without incumbent (`solver_status = timeout_no_incumbent`, `n_distinct_valid_solutions = 0`, summary metrics strictly `NaN`)
  - **Explicit `fewer_than_target_reason` breakdown across the 12,580 cases**: `stopped_early_timeout_no_incumbent` (`5,540` cases where a subsequent alternative-cut solve hit the `6.0` s per-MILP limit after $\ge 1$ valid solution was already accepted), `target_cap_reached_more_alternatives_may_exist` (`4,819`), `infeasible_under_fixed_medium` (`1,936`), `stopped_early_fba_verification_failed` (`143`), `stopped_early_time_limit_heuristic_incumbent` (`67`), `fba_verification_failed` (`56`), `timeout_no_incumbent` (`10`), `stopped_early_post_pruning_cost_changed_not_proven_optimal` (`8`), and `stopped_early_optimal` (`1`).
  - **93,808 per-solution ledger rows**: `91,806` `accepted_feasible` solutions + `1,936` `unresolved_infeasible_under_fixed_medium` + `56` `unresolved_fba_verification_failed` + `10` `unresolved_timeout_no_incumbent`.

---

## 5. Reproducible Artificial Masking Controls (v4)

To benchmark gapfill behavior under artificial gene incompleteness:
- **Genomes**: The 7 external reference genomes in Stage 1 inventory meeting numeric `completeness >= 90%` and `contamination <= 5%` (`MGYG000519815`, `MGYG000520046`, `MGYG000519748`, `MGYG000510876`, `MGYG000495096`, `MGYG000526805`, `MGYG000523080`; formal MIMAG rRNA/tRNA status unknown).
- **Design**: 3 masking fractions (`10%`, `20%`, `30%`) $\times$ 3 replicates = **63 evaluations** (`exp1_masking_controls_7hq_v4_validated.tsv`).
- **Matched Medium & Biomass**: Baseline (`0%` mask) and masked models are evaluated under the **exact same `S4_oxic_desert_oa_aa` medium** (verified via canonical SHA-256 JSON hashing in `exp1_masking_controls_media_provenance.json`) and **exact same `B2` biomass** (`B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420`).
- **Deterministic Seeds & Post-Hoc Receipt**: Random seeds are derived from SHA-256 hashes of `exp1_masking_{genome}_{frac}_{rep}`. All 63/63 cases pass post-hoc receipt verification (`exp1_masking_controls_7hq_v4_validation_receipt.json`: `computational_receipt_integrity = PASS`, `control_calibration_status = RESOLVED_ALL_CASES_USABLE`, mean recovery rate `0.168216`, mean alternative fraction `0.767300`, `raw_flux_vector_verification_status = UNKNOWN_NOT_PERSISTED_IN_V4_TABLES`).

---

## 6. Statistical Semantics & Claim Boundaries

1. **Finite Ensemble vs Biochemical Absence**: Gapfill enumeration samples up to `20` (`B0`/`B2`) or `10` (`B1`/`B3`/`B4`) distinct solutions within `25%` weighted cost slack (`6.0` s per MILP solve). A reaction frequency of `0.0` reflects lack of observation within the finite heuristic ensemble (`not_observed_in_finite_ensemble`), NOT proven biological absence or non-requirement.
2. **Biomass Sensitivity Categorization**: Reactions observed across all 5 tested biomass recipes (`B0`–`B4`) are classified as `observed_across_all_tested_biomass_perturbations` (never `directly_verified_robust_candidate`). Reactions required under `B2` supernetwork knockout are labeled as `conditional_b2_essential_candidate`.
3. **No Unwarranted Accuracy Claims**: Masking recovery benchmarks gapfill behavior under the declared 2,270-reaction pool and `S4` medium; candidate weights of `100.0` for masked reactions reflect withheld sequence evidence. Non-recovered reactions reflect lack of essentiality under the simulated medium, not draft annotation error.
4. **Stage 4 Candidate Deficiency Verification Semantics (v5) & Bounded Post-Hoc Qualification**:
   - **Unique `(genome, clean_rxn_id)` Pair Denominator**: The 13,750 scenario-level rows repeat identical genome-level and sample-level sequence alignments across nutrient scenarios (`S0`–`S6`) and collapse to **3,713 unique `(genome, clean_rxn_id)` pairs** across 81 MAGs (47 samples) and 63 reactions (`0` within-pair disposition conflicts). Biological sequence evidence is reported over the **3,713 unique pairs** (`exp1_stage4_v5_posthoc_unique_pair_dispositions.tsv`), never as 13,750 independent observations:
     - `854` broad MAG homolog candidate (`homolog_candidate_detected_in_MAG`)
     - `704` domain/partial MAG hit (`domain_level_or_partial_homolog_candidate_in_MAG`)
     - `137` broad near-contig-edge alignment (`broad_reference_alignment_near_contig_edge`)
     - `68` partial/truncated near-contig-edge alignment (`candidate_truncated_or_split_gene_at_contig_edge`)
     - `396` community broad homolog with unassigned host (`unresolved_sample_community_homolog_host_unassigned`)
     - `37` community domain/partial hit with unassigned host (`unresolved_sample_community_domain_level_hit_host_unassigned`)
     - `1,458` reference unavailable (`not_evaluable_missing_reference_sequences`, 18/63 reactions)
     - `59` with no qualifying MAG or sample community alignment (`tested_no_homolog_detected_conditional_negative_evidence` in raw v5; bounded post-hoc qualified in `exp1_stage4_v5_posthoc_scientific_disposition_summary.tsv`).
   - **Reference Namespace Provenance Boundary**: Generic gapseq crosswalk mappings across all 63 candidate reactions represent **retrieval namespace provenance** (linking ModelSEED `rxn*` IDs to gapseq sequence files in `/project/k10018/sw/db/gapseq_seq/Archaea/{rev,unrev,rxn}/*.fasta`), NOT an experimentally validated reaction-specific enzyme gold standard. `selected_candidate_reference_unique` (`resolved_selected_candidate_reference_unique`) denotes non-sharing among the 63 selected candidate reactions only. No strict function assignments or 59-missing-enzyme claims are derived.
   - **Validated Single-HSP Reference Coverage**: In `diamond blastx` (DNA contigs vs protein references), `qlen` is DNA contig length (bp) and `slen` is reference protein length (aa). Coverage is computed strictly as `(abs(send - sstart) + 1) / slen` in `(0, 1]`, never dividing gapped `alignment_length` by `slen`.
   - **All 59 No-Hit Pairs Involve Only Two Template-Sensitive Reactions (`rxn10457` and `rxn02322`)**:
     1. **`rxn10457` (47 unique `(genome, rxn10457)` pairs / 179 scenario instances)**: All rows are `conditional_b2_essential_candidate` and all vanish when specialized lipids/homospermidine are omitted (`biomass_assumption_sensitive__vanishes_when_specialized_lipids_or_homospermidine_omitted`). The gapseq crosswalk retrieves 5 proteins from `unrev/2.7.7.41.fasta` via `CDPDIGLYSYN-RXN`. Official EC `2.7.7.41` catalyses phosphatidate $\rightarrow$ CDP-diacylglycerol (acyl-lipid chemistry), whereas the ModelSEED target `rxn10457` is `CDP-2-O-(3'-hydroxy)geranyl-3-O-geranyl-sn-glycerol synthase` (archaeal ether-lipid cytidylyltransferase chemistry; official EC `2.7.7.67` is CDP-unsaturated archaeol synthase / `CarS`, itself not automatic proof of exact 3'-hydroxy-substrate activity). Thus the executed search tests gapseq alias-derived reference proteins only and does **NOT** establish absence of archaeal `CarS` or the exact target enzyme (`no_qualifying_alignment_to_alias_derived_acyl_lipid_reference_ec_2_7_7_41_biochemical_scope_unresolved`). Sources: https://iubmb.qmul.ac.uk/enzyme/EC2/7/7/41.html ; https://iubmb.qmul.ac.uk/enzyme/EC2/7/7/67.html ; https://pubmed.ncbi.nlm.nih.gov/25219966/ (Jain et al. 2014).
     2. **`rxn02322` (12 unique `(genome, rxn02322)` pairs / 46 scenario instances)**: All rows are `conditional_b2_essential_candidate`; 44 scenario instances vanish when specialized lipids/homospermidine are omitted and 2 vanish when the core skeleton (`B4`) is used. The target is classical phosphomevalonate kinase (`EC 2.7.4.2`, `ATP:(R)-5-phosphomevalonate phosphotransferase`, 9 reference proteins). Archaea possess known alternative mevalonate (MVA) pathways that bypass classical `EC 2.7.4.2` PMK; non-detection of classical PMK does **NOT** imply IPP/isoprenoid/lipid auxotrophy or a novel dark enzyme, and the actual MVA route in these `JACQPN01` and other MAGs is not yet established (noting `Thermoplasmatota` phylum $\neq$ `Thermoplasmatales` order; `no_qualifying_alignment_to_classical_pmk_ec_2_7_4_2_alternative_archaeal_mva_routes_unresolved`). Sources: https://pmc.ncbi.nlm.nih.gov/articles/PMC11052028/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC4047369/ .
   - **Sample Assembly Input Provenance Addendum**: Job `15895063` recorded the 47 raw DIAMOND output TSV hashes in `stage4_v5_execution_receipts.json` at execution time. The 47 raw input assembly FASTA paths (`/project/k10018/empty-quarter/data/assemblies/<sample>/<sample>_contigs.fa`), byte sizes (`47,555,575,438` bytes total), mtimes, and SHA-256 hashes are provided in `provenance/stage4_v5_posthoc_47_sample_assembly_input_manifest.tsv` (linked post-hoc from the pre-existing Shaheen manifest `stage2_target_assemblies_manifest.tsv` and verified against live filesystem stats with 0 mismatches, while preserving `stage4_v5_execution_receipts.json` unmodified).

---

## 7. Directory Layout & Verification

```
exp1_metabolic/
├── README.md                                         # This documentation
├── verify_frozen_outputs.py                          # 8-gate standalone reproducibility verifier
├── code/                                             # Production pipeline, HiGHS solver engine, and Slurm scripts
│   ├── exp1_sbml_highs_engine.py                     # Custom HiGHS MILP gapfill & FBA solver engine
│   ├── exp1_sbml_highs_engine_v4.py                  # Explicit thread-safe HiGHS solver engine (v4/v5)
│   ├── stage1_inventory_and_qc.py                    # Stage 1 inventory and QC runner
│   ├── stage2_pilot24_ensemble.py                    # Stage 2 Pilot24 & full protocol ensemble runner
│   ├── stage2_focused_release_validation.py          # Stage 2 witness revalidation runner
│   ├── run_masking_controls_7hq.py                   # Stage 2 masking controls runner
│   ├── run_masking_controls_7hq_v4_resumable.py      # Resumable, parallelized masking controls runner (v4)
│   ├── run_masking_controls_7hq_v5_resumable.py      # Future-ready self-contained masking controls runner (v5)
│   ├── validate_masking_controls_v4_outputs.py       # Independent post-hoc v4 masking controls validator
│   ├── stage3_full484_worker.py                      # Stage 3 parallel shard array worker
│   ├── stage3_collect_and_summarize_v3.py            # Stage 3 collector (v3)
│   ├── stage3_collect_and_summarize_v4.py            # Stage 3 collector with validated v4 controls receipt (v4)
│   ├── stage4_candidate_deficiency_search.py         # Stage 4 candidate deficiency search (v5 crosswalk & semantics)
│   ├── stage4_candidate_deficiency_search_v5.py      # Explicit versioned Stage 4 v5 runner
│   ├── posthoc_qualify_stage4_v5_homology_dispositions.py # Bounded post-hoc 3,713-pair scientific-disposition qualifier
│   └── sbatch_*.sh                                   # Production Slurm submission scripts
├── tests/                                            # Comprehensive unit & adversarial test suite
│   ├── test_exp1_metabolic.py                        # Core engine, Stage 1-4, full484 ledger, and post-hoc tests
│   ├── test_masking_controls_v4_validation.py        # Masking controls v4/v5 validation and fail-closed tests
│   └── test_highs_spawn_regression.py                # Multiprocessing spawn + HiGHS OpenMP regression test
├── compact_evidence/                                 # Frozen accepted evidence tables and provenance receipts
│   ├── exp1_biomass_template_audit_484.tsv
│   ├── exp1_candidate_reaction_pool_audit.tsv
│   ├── exp1_candidate_reaction_pool_conflicts.tsv
│   ├── exp1_candidate_species_conflicts.tsv
│   ├── exp1_model_inventory_and_qc_484.tsv
│   ├── exp1_model_inventory_and_qc_538.tsv
│   ├── exp1_pilot24_manifest.tsv
│   ├── exp1_pilot24_qc_gate.json                     # Historical Stage 2 pilot record (not final acceptance gate)
│   ├── exp1_pilot24_supplement_rescue.tsv
│   ├── exp1_masking_controls_7hq_baselines.tsv
│   ├── exp1_masking_controls_media_provenance.json
│   ├── exp1_masking_controls_7hq_v4_validated.tsv
│   ├── exp1_masking_controls_7hq_v4_validation_receipt.json
│   ├── exp1_full484_ensemble_summary.tsv.gz          # Full 484-genome 12,580-case production summary
│   ├── exp1_full484_per_solution_ledger.tsv.gz       # Full 484-genome 93,808-row per-solution ledger (91,806 accepted)
│   ├── exp1_full484_cohort_and_quality_comparisons.tsv
│   ├── exp1_full484_lineage_matched_comparisons_order_family_genus.tsv
│   ├── exp1_stage3_shard_collection_receipts.json
│   ├── exp1_stage4_supported_candidate_deficiencies_summary.tsv
│   ├── exp1_stage4_v5_deficiency_verification_summary.tsv
│   ├── exp1_stage4_v5_posthoc_unique_pair_dispositions.tsv
│   ├── exp1_stage4_v5_posthoc_scientific_disposition_summary.tsv
│   ├── provenance/
│   │   ├── exp1_full484_independent_table_audit.json
│   │   ├── stage4_v5_execution_receipts.json
│   │   ├── stage4_v5_reference_namespace_crosswalk.tsv
│   │   ├── stage4_v5_reference_sequence_ownership.tsv
│   │   ├── stage4_v5_posthoc_47_sample_assembly_input_manifest.tsv
│   │   ├── stage4_v5_posthoc_scientific_disposition_receipt.json
│   │   ├── baseline_15894447_unperformed_namespace_failure/
│   │   └── failed_15894840_inventory_path_typo/
│   └── focused_release_validation/
└── manifests/                                        # Scoped checksums, Slurm receipts, and fetch instructions
    ├── institutional_fetch_instructions.txt          # (Note: sha256_manifest.tsv is scoped to core owned files;
    ├── sha256_manifest.tsv                           #  covariates/, interpretation/, counterfactual_review/ have
    └── slurm_execution_receipts.tsv                  #  independent subpackage manifests)
```

### Verification Commands

```bash
# Run the 8-gate frozen output & scoped manifest verifier
uv run python verify_frozen_outputs.py

# Run the full unit and adversarial test suite
PYTHONPATH=code uv run pytest tests/
```
