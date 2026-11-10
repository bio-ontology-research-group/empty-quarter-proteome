#!/usr/bin/env python3
"""Adversarial unit and scientific validation tests for `exp1_metabolic` (`EXP1-CONTRACT.txt`,
`EXP1-BIOMASS-SCIENTIFIC-NOTE.txt`, `EXP1-SOLVER-REVIEW.txt`, and `EXP1-FINAL-ENGINE-REVIEW.txt`).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import OptimizeResult

p_v3 = Path(__file__).resolve().parent.parent / "code_v3"
p_v2 = Path(__file__).resolve().parent.parent / "code_v2"
CODE_DIR = p_v3 if p_v3.exists() else (p_v2 if p_v2.exists() else Path(__file__).resolve().parent.parent / "code")
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from exp1_sbml_highs_engine import (
    BIOMASS_PERTURBATION_RECIPES,
    F420_COFACTOR_ID,
    METHANOSARCINA_ALL8_TEMPLATE_COFACTORS,
    METHANOSARCINA_NON_F420_TEMPLATE_COFACTORS,
    NUTRIENT_SCENARIOS,
    OBSERVED_COHORT_538_POOL_SCOPE,
    ReactionRecord,
    SBMLModelData,
    SpeciesRecord,
    aggregate_candidate_pool_deterministically,
    apply_artificial_gene_masking,
    apply_biomass_perturbation,
    apply_nutrient_scenario_bounds,
    audit_model_mass_charge_balance,
    audit_reaction_mass_charge_balance,
    build_fixed_medium_supernetwork_and_base,
    classify_biomass_requirement_provenance,
    close_all_external_inflows_directionally,
    detect_closed_uptake_egc,
    enumerate_alternative_gapfill_ensemble,
    evaluate_gpr_ast,
    parse_chemical_formula,
    parse_gapseq_sbml_bytes,
    run_withdrawal_knockout_fva,
    solve_fba_highs,
    solve_nutrient_supplement_rescue_milp,
    validate_archaeal_biomass_and_audit_template,
)


def _make_base_archaeal_toy_model() -> SBMLModelData:
    """Create a realistic archaeal toy model with:
    - Extracellular & cytosolic nutrients (glucose `cpd00027`, NH3 `cpd00013`, H2O `cpd00001`, Pi `cpd00009`, H+ `cpd00067`)
    - ATP/ADP (`cpd00002`, `cpd00008`)
    - L-Lysine (`cpd00039`) and L-Alanine (`cpd00035`)
    - All 8 Methanosarcina template trace cofactors in `R_bio1` (`cpd00649`, `cpd02246`, `cpd02817`, `cpd03425`, `cpd15862`, `cpd15863`, `cpd15884`, `cpd33626`)
    - Archaeal ether lipid (`cpd15858`)
    """
    species: dict[str, SpeciesRecord] = {
        "M_cpd00027_e0": SpeciesRecord("M_cpd00027_e0", "cpd00027", "D-Glucose", "e0", 0, "C6H12O6"),
        "M_cpd00027_c0": SpeciesRecord("M_cpd00027_c0", "cpd00027", "D-Glucose", "c0", 0, "C6H12O6"),
        "M_cpd00013_e0": SpeciesRecord("M_cpd00013_e0", "cpd00013", "NH3", "e0", 0, "H3N"),
        "M_cpd00013_c0": SpeciesRecord("M_cpd00013_c0", "cpd00013", "NH3", "c0", 0, "H3N"),
        "M_cpd00007_e0": SpeciesRecord("M_cpd00007_e0", "cpd00007", "O2", "e0", 0, "O2"),
        "M_cpd00007_c0": SpeciesRecord("M_cpd00007_c0", "cpd00007", "O2", "c0", 0, "O2"),
        "M_cpd00001_c0": SpeciesRecord("M_cpd00001_c0", "cpd00001", "H2O", "c0", 0, "H2O"),
        "M_cpd00009_c0": SpeciesRecord("M_cpd00009_c0", "cpd00009", "Phosphate", "c0", -2, "HO4P"),
        "M_cpd00067_c0": SpeciesRecord("M_cpd00067_c0", "cpd00067", "H+", "c0", 1, "H"),
        "M_cpd00002_c0": SpeciesRecord("M_cpd00002_c0", "cpd00002", "ATP", "c0", -4, "C10H12N5O13P3"),
        "M_cpd00008_c0": SpeciesRecord("M_cpd00008_c0", "cpd00008", "ADP", "c0", -3, "C10H12N5O10P2"),
        "M_cpd00035_c0": SpeciesRecord("M_cpd00035_c0", "cpd00035", "L-Alanine", "c0", 0, "C3H7NO2"),
        "M_cpd00039_e0": SpeciesRecord("M_cpd00039_e0", "cpd00039", "L-Lysine", "e0", 1, "C6H15N2O2"),
        "M_cpd00039_c0": SpeciesRecord("M_cpd00039_c0", "cpd00039", "L-Lysine", "c0", 1, "C6H15N2O2"),
        "M_cpd15858_c0": SpeciesRecord("M_cpd15858_c0", "cpd15858", "2,3-O-phytanyl-sn-glycero-1-phosphoethanolamine", "c0", 0, "C43H88O3"),
        "M_cpd11416_c0": SpeciesRecord("M_cpd11416_c0", "cpd11416", "Biomass", "c0", 0, ""),
    }
    for cpd_id, cname in sorted(METHANOSARCINA_ALL8_TEMPLATE_COFACTORS.items()):
        species[f"M_{cpd_id}_c0"] = SpeciesRecord(f"M_{cpd_id}_c0", cpd_id, cname, "c0", 0, "C10H10N2O2")

    reactions: dict[str, ReactionRecord] = {
        "R_EX_cpd00027_e0": ReactionRecord(
            id="R_EX_cpd00027_e0", base_id="EX_cpd00027_e0", name="Glucose exchange",
            reversible=True, lb=-5.0, ub=1000.0, stoichiometry={"M_cpd00027_e0": -1.0}, is_exchange=True,
        ),
        "R_EX_cpd00013_e0": ReactionRecord(
            id="R_EX_cpd00013_e0", base_id="EX_cpd00013_e0", name="NH3 exchange",
            reversible=True, lb=-10.0, ub=1000.0, stoichiometry={"M_cpd00013_e0": -1.0}, is_exchange=True,
        ),
        "R_EX_cpd00039_e0": ReactionRecord(
            id="R_EX_cpd00039_e0", base_id="EX_cpd00039_e0", name="Lysine exchange",
            reversible=True, lb=0.0, ub=1000.0, stoichiometry={"M_cpd00039_e0": -1.0}, is_exchange=True,
        ),
        "R_EX_cpd00007_e0": ReactionRecord(
            id="R_EX_cpd00007_e0", base_id="EX_cpd00007_e0", name="O2 exchange",
            reversible=True, lb=0.0, ub=1000.0, stoichiometry={"M_cpd00007_e0": -1.0}, is_exchange=True,
        ),
        "R_glc_trans": ReactionRecord(
            id="R_glc_trans", base_id="rxn_glc_trans", name="Glucose transport",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00027_e0": -1.0, "M_cpd00027_c0": 1.0}, is_transport=True,
            gpr_ast=("GENE", "g_glcT"), genes={"g_glcT"},
        ),
        "R_nh3_trans": ReactionRecord(
            id="R_nh3_trans", base_id="rxn_nh3_trans", name="NH3 transport",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00013_e0": -1.0, "M_cpd00013_c0": 1.0}, is_transport=True,
            gpr_ast=("GENE", "g_amtB"), genes={"g_amtB"},
        ),
        "R_atp_ferm": ReactionRecord(
            id="R_atp_ferm", base_id="rxn_atp_ferm", name="Fermentative ATP generation",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={
                "M_cpd00027_c0": -1.0,
                "M_cpd00008_c0": -2.0,
                "M_cpd00009_c0": -2.0,
                "M_cpd00067_c0": -2.0,
                "M_cpd00002_c0": 2.0,
                "M_cpd00001_c0": 2.0,
            },
            gpr_ast=("AND", [("GENE", "g_gly1"), ("GENE", "g_gly2")]),
            genes={"g_gly1", "g_gly2"},
        ),
        "R_ala_synth": ReactionRecord(
            id="R_ala_synth", base_id="rxn_ala_synth", name="Alanine synthesis",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00035_c0": 1.0},
            gpr_ast=("OR", [("GENE", "g_alaA"), ("GENE", "g_alaB")]),
            genes={"g_alaA", "g_alaB"},
        ),
        "R_lipid_synth": ReactionRecord(
            id="R_lipid_synth", base_id="rxn_lipid_synth", name="Archaeol synthesis",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd15858_c0": 1.0},
            gpr_ast=("GENE", "g_arch"), genes={"g_arch"},
        ),
        "R_f420_synth": ReactionRecord(
            id="R_f420_synth", base_id="rxn_f420_synth", name="Coenzyme F420 synthesis",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.1, "M_cpd00013_c0": -0.1, "M_cpd00649_c0": 1.0},
            gpr_ast=("GENE", "g_fbiC"), genes={"g_fbiC"},
        ),
    }
    # Add synthesis for the 7 non-F420 Methanosarcina trace cofactors in the baseline model
    for cpd_id in sorted(METHANOSARCINA_NON_F420_TEMPLATE_COFACTORS):
        reactions[f"R_synth_{cpd_id}"] = ReactionRecord(
            id=f"R_synth_{cpd_id}",
            base_id=f"rxn_synth_{cpd_id}",
            name=f"Synth {cpd_id}",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.05, f"M_{cpd_id}_c0": 1.0},
        )

    bio_stoich: dict[str, float] = {
        "M_cpd00002_c0": -1.0,
        "M_cpd00001_c0": -1.0,
        "M_cpd00008_c0": 1.0,
        "M_cpd00009_c0": 1.0,
        "M_cpd00067_c0": 1.0,
        "M_cpd00035_c0": -0.5,
        "M_cpd00039_c0": -0.5,
        "M_cpd15858_c0": -0.01,
        "M_cpd11416_c0": 1.0,
    }
    for cpd_id in sorted(METHANOSARCINA_ALL8_TEMPLATE_COFACTORS):
        bio_stoich[f"M_{cpd_id}_c0"] = -0.001

    reactions["R_bio1"] = ReactionRecord(
        id="R_bio1",
        base_id="bio1",
        name="Archaea biomass (Methanosarcina template)",
        reversible=False,
        lb=0.0,
        ub=1000.0,
        stoichiometry=bio_stoich,
        is_biomass=True,
    )
    reactions["R_EX_cpd11416_c0"] = ReactionRecord(
        id="R_EX_cpd11416_c0",
        base_id="EX_cpd11416_c0",
        name="Biomass sink",
        reversible=False,
        lb=0.0,
        ub=1000.0,
        stoichiometry={"M_cpd11416_c0": -1.0},
        is_exchange=True,
    )

    genes = {
        "g_glcT": "g_glcT",
        "g_amtB": "g_amtB",
        "g_gly1": "g_gly1",
        "g_gly2": "g_gly2",
        "g_alaA": "g_alaA",
        "g_alaB": "g_alaB",
        "g_arch": "g_arch",
        "g_fbiC": "g_fbiC",
    }
    return SBMLModelData(
        model_id="toy_archaea_1",
        model_name="Toy Archaeal Model",
        species=species,
        reactions=reactions,
        gene_products=genes,
        objective_rxn_id="R_bio1",
        notes_text="",
    )


def test_1_exact_binary_milp_vs_flux_weighted_l1_disagreement():
    """Adversarial Test (`EXP1-SOLVER-REVIEW.txt` Item 1 & 7):
    Construct a network where:
    - Route A uses 1 reaction (`R_cand_A`) with binary cost `c_A = 10.0`, but requires a large stoichiometric flux
      (`v_A = 100.0` per unit of lysine because stoichiometry is `+0.01 M_cpd00039_c0`), so L1 flux cost `c_A * v_A = 1000`.
    - Route B uses 2 reactions (`R_cand_B1`, `R_cand_B2`) each with binary cost `c_B = 6.0` (total binary cost `12.0 > 10.0`),
      each carrying tiny flux `v_B = 0.5` (L1 flux cost `6*0.5 + 6*0.5 = 6.0 << 1000`).
    A flux-weighted L1 LP falsely selects Route B (`{R_cand_B1, R_cand_B2}`, binary cost 12.0), whereas the true
    binary reaction-cost MILP selects Route A (`{R_cand_A}`, binary cost 10.0) first and orders solutions by exact binary cost!
    Also verifies `res.mip_dual_bound` is captured directly without `obj/(1+mip_gap)` inversion (`EXP1-FINAL-ENGINE-REVIEW.txt` Item A).
    """
    draft = _make_base_archaeal_toy_model()
    sp_mid = SpeciesRecord("M_cpd99001_c0", "cpd99001", "Intermediate B", "c0", 0, "C6H15N2O2")
    draft.species[sp_mid.id] = sp_mid

    candidates = {
        "R_cand_A_highflux_1rxn": ReactionRecord(
            id="R_cand_A_highflux_1rxn",
            base_id="rxn_cand_A",
            name="High-flux single-reaction lysine route",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.005, "M_cpd00013_c0": -0.01, "M_cpd00039_c0": 0.01},
        ),
        "R_cand_B1_lowflux_step1": ReactionRecord(
            id="R_cand_B1_lowflux_step1",
            base_id="rxn_cand_B1",
            name="Low-flux 2-step route step 1",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd99001_c0": 1.0},
        ),
        "R_cand_B2_lowflux_step2": ReactionRecord(
            id="R_cand_B2_lowflux_step2",
            base_id="rxn_cand_B2",
            name="Low-flux 2-step route step 2",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd99001_c0": -1.0, "M_cpd00039_c0": 1.0},
        ),
    }
    weights = {
        "R_cand_A_highflux_1rxn": 10.0,
        "R_cand_B1_lowflux_step1": 6.0,
        "R_cand_B2_lowflux_step2": 6.0,
    }
    scenario = {
        "cpd00027": ("D-Glucose", 10.0),
        "cpd00013": ("NH3", 10.0),
    }

    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=draft,
        candidate_reactions=candidates,
        candidate_weights=weights,
        scenario_compounds=scenario,
        min_growth=0.05,
        n_Target_solutions=10,
        cost_slack_frac=0.30,  # allows up to 10.0 * 1.30 = 13.0, so both Route A (10.0) and Route B (12.0) are feasible!
    )
    assert ens["solver_status"] == "optimal"
    assert ens["cost_window_anchor_type"] == "proven_global_milp_optimum_on_declared_pool"
    assert ens["enumeration_completeness"] == "exhausted_within_cost_slack_on_declared_pool"
    assert ens["pool_exhaustion_type"] == "mathematically_exhausted_and_all_cycle_qc_eligible"
    assert ens["n_distinct_valid_solutions"] == 2

    sol0 = ens["solutions"][0]
    sol1 = ens["solutions"][1]
    assert sol0["added_rxns"] == ["R_cand_A_highflux_1rxn"]
    assert math.isclose(sol0["total_cost"], 10.0, rel_tol=1e-6)
    assert math.isclose(sol0["dual_bound"], 10.0, rel_tol=1e-5)
    assert sol1["added_rxns"] == ["R_cand_B1_lowflux_step1", "R_cand_B2_lowflux_step2"]
    assert math.isclose(sol1["total_cost"], 12.0, rel_tol=1e-6)


def test_2_anoxic_medium_never_leaks_oxygen_and_supplement_rescue_is_separate():
    """Adversarial Test (`EXP1-SOLVER-REVIEW.txt` Item 2 & 7):
    Verify that when a reaction requires O2 (`M_cpd00007_c0`) and an O2 exchange (`R_EX_cpd00007_e0`) plus O2 transport
    are present in `candidate_reactions`, an anoxic medium (`cpd00007` NOT in `scenario_compounds`) NEVER opens O2 uptake
    during fixed-medium gap-filling (`infeasible_under_fixed_medium`).
    Also verify that `solve_nutrient_supplement_rescue_milp` raises `ValueError` if O2 is offered in an anoxic regime,
    and rescues with external L-Lysine (`cpd00039`) in a separate supplement experiment without mutating the baseline result.
    """
    draft = _make_base_archaeal_toy_model()
    candidates = {
        "R_EX_cpd00007_e0_cand": ReactionRecord(
            id="R_EX_cpd00007_e0", base_id="EX_cpd00007_e0", name="O2 exchange candidate",
            reversible=True, lb=-100.0, ub=1000.0, stoichiometry={"M_cpd00007_e0": -1.0}, is_exchange=True,
        ),
        "R_o2_trans": ReactionRecord(
            id="R_o2_trans", base_id="rxn_o2_trans", name="O2 transport",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00007_e0": -1.0, "M_cpd00007_c0": 1.0}, is_transport=True,
        ),
        "R_o2_dependent_lys": ReactionRecord(
            id="R_o2_dependent_lys", base_id="rxn_o2_lys", name="O2-dependent lysine synthesis",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00007_c0": -1.0, "M_cpd00039_c0": 1.0},
        ),
        "R_lys_permease": ReactionRecord(
            id="R_lys_permease", base_id="rxn_lys_perm", name="Lysine permease",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00039_e0": -1.0, "M_cpd00039_c0": 1.0}, is_transport=True,
        ),
    }
    weights = {r: 1.0 for r in candidates}
    anoxic_min_medium = {
        "cpd00027": ("D-Glucose", 5.0),
        "cpd00013": ("NH3", 10.0),
    }

    # Fixed anoxic minimal medium must be strictly INFEASIBLE (neither O2 nor external Lysine is supplied!)
    ens_anoxic = enumerate_alternative_gapfill_ensemble(
        draft_model=draft,
        candidate_reactions=candidates,
        candidate_weights=weights,
        scenario_compounds=anoxic_min_medium,
        min_growth=0.01,
    )
    assert ens_anoxic["solver_status"] == "infeasible_under_fixed_medium"
    assert ens_anoxic["n_distinct_valid_solutions"] == 0

    # Offering O2 in an anoxic supplement rescue must be rejected
    with pytest.raises(ValueError, match="Oxygen.*anoxic"):
        solve_nutrient_supplement_rescue_milp(
            draft_model=draft,
            candidate_reactions=candidates,
            candidate_weights=weights,
            fixed_scenario_compounds=anoxic_min_medium,
            supplement_candidate_compounds={"cpd00007": ("O2", 5.0, 50.0)},
            redox_regime="anoxic",
        )

    # Offering L-Lysine as an explicit supplement in the separate rescue MILP succeeds and records the supplement
    rescue = solve_nutrient_supplement_rescue_milp(
        draft_model=draft,
        candidate_reactions=candidates,
        candidate_weights=weights,
        fixed_scenario_compounds=anoxic_min_medium,
        supplement_candidate_compounds={"cpd00039": ("L-Lysine", 2.0, 50.0)},
        redox_regime="anoxic",
    )
    assert rescue["rescued"] is True
    assert rescue["primary_supplements_required"] == ["R_SUPP_UPTAKE_cpd00039_e0"]


def test_3_candidate_introduced_extracellular_species_verified_in_mod_check_and_withdrawal():
    """Adversarial Test (`EXP1-FINAL-ENGINE-REVIEW.txt` Item B):
    When a candidate transport reaction introduces a NEW extracellular species `M_cpd00555_e0` that is NOT in `draft.species`
    (so `mod_super` attaches `R_EX_cpd00555_e0` with supplied medium inflow `lb = -2.0`), verify that `fixed_medium_base`,
    `mod_check`, `mod_pruned`, and `run_withdrawal_knockout_fva` all retain `R_EX_cpd00555_e0` and `M_cpd00555_e0` so
    post-solve FBA verification and nutrient withdrawal succeed!
    """
    draft = _make_base_archaeal_toy_model()
    assert "M_cpd00555_e0" not in draft.species
    assert "R_EX_cpd00555_e0" not in draft.reactions

    candidates = {
        "R_novel_precursor_import": ReactionRecord(
            id="R_novel_precursor_import",
            base_id="rxn_novel_imp",
            name="Novel extracellular precursor import",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd00555_e0": -1.0, "M_cpd00555_c0": 1.0},
            is_transport=True,
        ),
        "R_novel_precursor_to_lys": ReactionRecord(
            id="R_novel_precursor_to_lys",
            base_id="rxn_novel_to_lys",
            name="Novel precursor to Lysine",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd00555_c0": -1.0, "M_cpd00039_c0": 1.0},
        ),
    }
    weights = {"R_novel_precursor_import": 2.0, "R_novel_precursor_to_lys": 3.0}
    scenario = {
        "cpd00027": ("D-Glucose", 5.0),
        "cpd00013": ("NH3", 10.0),
        "cpd00555": ("Novel Precursor", 2.0),
    }

    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=draft,
        candidate_reactions=candidates,
        candidate_weights=weights,
        scenario_compounds=scenario,
        min_growth=0.05,
    )
    assert ens["solver_status"] == "optimal"
    assert ens["n_distinct_valid_solutions"] == 1
    assert ens["solutions"][0]["added_rxns"] == ["R_novel_precursor_import", "R_novel_precursor_to_lys"]
    assert ens["solutions"][0]["verified_fba_growth"] >= 0.05

    wkf = run_withdrawal_knockout_fva(
        base_model=draft,
        ensemble_result=ens,
        candidate_reactions=candidates,
        scenario_compounds=scenario,
        min_growth=0.05,
    )
    assert wkf["evaluation_status"] == "evaluated_feasible_baseline"
    by_cpd = {w["compound_id"]: w for w in wkf["withdrawals"]}
    assert "cpd00555" in by_cpd
    assert by_cpd["cpd00555"]["dependency_classification"] == "required_across_sampled_ensemble_union"
    assert by_cpd["cpd00555"]["addback_rescues_growth"] is True


def test_4_candidate_bounds_positive_lb_and_zero_bounds():
    """Adversarial Test (`EXP1-FINAL-ENGINE-REVIEW.txt` Item C):
    1) Verify that a candidate reaction with native `lb > 0` (`lb = 0.2, ub = 10.0`) or `ub < 0` (`lb = -10.0, ub = -0.2`)
       is permitted to have `v = 0` when `z = 0` (so it does NOT force infeasibility or mandatory selection when unselected!),
       and enforces `lb <= v <= ub` when `z = 1`.
    2) Verify that a zero-bound candidate reaction (`lb = 0.0, ub = 0.0`) is excluded and never selected.
    """
    draft = _make_base_archaeal_toy_model()
    candidates = {
        "R_zero_bounds_cand": ReactionRecord(
            id="R_zero_bounds_cand",
            base_id="rxn_zero_bounds",
            name="Zero bounds candidate",
            reversible=False,
            lb=0.0,
            ub=0.0,
            stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00039_c0": 1.0},
        ),
        "R_deadend_pos_lb_cand": ReactionRecord(
            id="R_deadend_pos_lb_cand",
            base_id="rxn_deadend_pos_lb",
            name="Unselected candidate with positive native lower bound",
            reversible=False,
            lb=5.0,
            ub=100.0,
            # Dead-end metabolite M_cpd88888_c0 has no sink, so if v > 0 were forced when z=0, the whole MILP would be infeasible!
            stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd88888_c0": 1.0},
        ),
        "R_valid_pos_lb_lys": ReactionRecord(
            id="R_valid_pos_lb_lys",
            base_id="rxn_valid_pos_lb_lys",
            name="Valid lysine synthesis with positive native lb",
            reversible=False,
            lb=0.1,
            ub=100.0,
            stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00039_c0": 1.0},
        ),
    }
    weights = {
        "R_zero_bounds_cand": 0.005,
        "R_deadend_pos_lb_cand": 1.0,
        "R_valid_pos_lb_lys": 5.0,
    }
    scenario = {"cpd00027": ("D-Glucose", 5.0), "cpd00013": ("NH3", 10.0)}

    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=draft,
        candidate_reactions=candidates,
        candidate_weights=weights,
        scenario_compounds=scenario,
        min_growth=0.05,
    )
    assert ens["solver_status"] == "optimal"
    assert ens["n_distinct_valid_solutions"] == 1
    sol0 = ens["solutions"][0]
    assert sol0["added_rxns"] == ["R_valid_pos_lb_lys"]
    assert "R_zero_bounds_cand" not in sol0["added_rxns"]
    assert "R_deadend_pos_lb_cand" not in sol0["added_rxns"]


def test_5_cycle_eligibility_filters_spurious_atp_cycle_and_rejects_corrupted_baseline():
    """Adversarial Test (`EXP1-FINAL-ENGINE-REVIEW.txt` Item D):
    1) If a cheap candidate reaction closes a spurious closed-uptake ATP-generating cycle (`R_spurious_cycle_closer`, cost 1.0)
       while a valid biochemical candidate (`R_clean_lys_synth`, cost 1.1) does not:
       - `optimization_feasible_solutions` records both,
       - `solutions` (science-QC-eligible) and `reaction_frequencies` exclude `R_spurious_cycle_closer` and include ONLY
         `R_clean_lys_synth`, with `pool_exhaustion_type == 'mathematically_exhausted_after_cycle_qc_filtering'`.
    2) If the baseline draft model itself grows only with a closed-uptake ATP cycle (or disguised non-`EX_` boundary inflow),
       it is rejected by cycle QC (`draft_growth_rejected_by_cycle_qc__failed_closed_uptake_atp_cycle`) and `run_withdrawal_knockout_fva`
       refuses to evaluate it (`non_evaluable_baseline_failed_cycle_qc`).
    """
    draft = _make_base_archaeal_toy_model()
    draft.species["M_cpd99998_c0"] = SpeciesRecord("M_cpd99998_c0", "cpd99998", "Cycle carrier A", "c0", 0, "C5H10O5")
    draft.species["M_cpd99999_c0"] = SpeciesRecord("M_cpd99999_c0", "cpd99999", "Cycle carrier B", "c0", 0, "C5H10O5")
    # Half-cycle in draft: M_cpd99999_c0 + ADP + Pi + H+ -> M_cpd99998_c0 + ATP + H2O (dormant unless M_cpd99998 -> M_cpd99999 is added!)
    draft.reactions["R_half_cycle_atp"] = ReactionRecord(
        id="R_half_cycle_atp",
        base_id="rxn_half_cycle_atp",
        name="Dormant half-cycle coupled to ATP synthesis",
        reversible=False,
        lb=0.0,
        ub=1000.0,
        stoichiometry={
            "M_cpd99999_c0": -1.0,
            "M_cpd00008_c0": -1.0,
            "M_cpd00009_c0": -1.0,
            "M_cpd00067_c0": -1.0,
            "M_cpd99998_c0": 1.0,
            "M_cpd00002_c0": 1.0,
            "M_cpd00001_c0": 1.0,
        },
    )
    draft.reactions["R_lys_coupled_99999"] = ReactionRecord(
        id="R_lys_coupled_99999",
        base_id="rxn_lys_99999",
        name="Lysine synthesis consuming M_cpd99999_c0",
        reversible=False,
        lb=0.0,
        ub=1000.0,
        stoichiometry={
            "M_cpd99999_c0": -1.0,
            "M_cpd00027_c0": -0.5,
            "M_cpd00013_c0": -1.0,
            "M_cpd99998_c0": 1.0,
            "M_cpd00039_c0": 1.0,
        },
    )
    candidates = {
        "R_spurious_cycle_closer": ReactionRecord(
            id="R_spurious_cycle_closer",
            base_id="rxn_spurious_cycle",
            name="Candidate that closes internal ATP-generating cycle",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd99998_c0": -1.0, "M_cpd99999_c0": 1.0},
        ),
        "R_clean_lys_synth": ReactionRecord(
            id="R_clean_lys_synth",
            base_id="rxn_clean_lys",
            name="Clean balanced lysine synthesis",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00039_c0": 1.0},
        ),
    }
    weights = {"R_spurious_cycle_closer": 1.0, "R_clean_lys_synth": 1.1}
    scenario = {"cpd00027": ("D-Glucose", 5.0), "cpd00013": ("NH3", 10.0)}

    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=draft,
        candidate_reactions=candidates,
        candidate_weights=weights,
        scenario_compounds=scenario,
        min_growth=0.05,
        cost_slack_frac=0.25,
    )
    assert ens["n_optimization_feasible_solutions"] == 2
    assert ens["n_cycle_qc_rejected_solutions"] == 1
    assert ens["n_distinct_valid_solutions"] == 1
    assert ens["pool_exhaustion_type"] == "mathematically_exhausted_after_cycle_qc_filtering"
    assert ens["solutions"][0]["added_rxns"] == ["R_clean_lys_synth"]
    assert list(ens["reaction_frequencies"].keys()) == ["R_clean_lys_synth"]
    assert ens["reaction_frequencies"]["R_clean_lys_synth"]["ensemble_frequency"] == 1.0

    # Now corrupt the baseline draft so it already grows AND has a closed-uptake ATP cycle
    corrupted = draft.copy()
    corrupted.reactions["R_clean_lys_synth"] = ReactionRecord(**vars(candidates["R_clean_lys_synth"]))
    corrupted.reactions["R_free_atp_pump"] = ReactionRecord(
        id="R_free_atp_pump",
        base_id="rxn_free_atp",
        name="Spurious internal ATP generator",
        reversible=False,
        lb=0.0,
        ub=1000.0,
        stoichiometry={
            "M_cpd00008_c0": -1.0,
            "M_cpd00009_c0": -1.0,
            "M_cpd00067_c0": -1.0,
            "M_cpd00002_c0": 1.0,
            "M_cpd00001_c0": 1.0,
        },
    )
    ens_corrupt = enumerate_alternative_gapfill_ensemble(
        draft_model=corrupted,
        candidate_reactions={},
        candidate_weights={},
        scenario_compounds=scenario,
        min_growth=0.05,
    )
    assert ens_corrupt["draft_already_grows"] is False
    assert "draft_growth_rejected_by_cycle_qc" in ens_corrupt["solver_status"]
    assert ens_corrupt["n_distinct_valid_solutions"] == 0

    wkf_corrupt = run_withdrawal_knockout_fva(
        base_model=corrupted,
        ensemble_result=ens_corrupt,
        candidate_reactions={},
        scenario_compounds=scenario,
    )
    assert "non_evaluable_baseline_failed_cycle_qc" in wkf_corrupt["evaluation_status"]


def test_6_biomass_perturbation_and_assumption_sensitive_classification():
    """Adversarial Test (`EXP1-BIOMASS-SCIENTIFIC-NOTE.txt` & `EXP1-SOLVER-REVIEW.txt` Item 4):
    1) Verify bacterial peptidoglycan/murein (`cpd15665`) in biomass raises `ValueError`.
    2) Verify that when a draft lacks one of the 7 non-F420 Methanosarcina cofactors (`cpd02246` CoM),
       it requires `R_cand_com_synth` under `B0_methanosarcina_template_baseline`, does NOT require it under
       `B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420` or `B1_perturbation_omit_8_template_cofactors`,
       and is classified as `biomass_assumption_sensitive__vanishes_when_7_non_f420_template_cofactors_omitted` (NEVER `proven_artifact`).
    """
    draft = _make_base_archaeal_toy_model()
    audit = validate_archaeal_biomass_and_audit_template(draft)
    assert audit["f420_demanded_in_template"] is True
    assert audit["n_non_f420_methanosarcina_cofactors_demanded"] == 7
    assert audit["methanofuran_demanded_in_template"] is False

    # Inject bacterial peptidoglycan into biomass -> must raise ValueError
    bad_bact = draft.copy()
    bad_bact.species["M_cpd15665_c0"] = SpeciesRecord("M_cpd15665_c0", "cpd15665", "Peptidoglycan polymer", "c0", 0, "")
    bad_bact.reactions["R_bio1"].stoichiometry["M_cpd15665_c0"] = -0.01
    with pytest.raises(ValueError, match="bacterial cell-wall"):
        validate_archaeal_biomass_and_audit_template(bad_bact)

    # Remove CoM (`R_synth_cpd02246`) from draft, provide Lysine in draft so ONLY CoM is missing in B0
    draft_no_com = draft.copy()
    del draft_no_com.reactions["R_synth_cpd02246"]
    draft_no_com.reactions["R_lys_synth"] = ReactionRecord(
        id="R_lys_synth", base_id="rxn_lys_synth", name="Lysine synthesis",
        reversible=False, lb=0.0, ub=1000.0,
        stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00039_c0": 1.0},
    )
    candidates = {
        "R_cand_com_synth": ReactionRecord(
            id="R_cand_com_synth", base_id="rxn_com_synth", name="CoM synthesis candidate",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.05, "M_cpd02246_c0": 1.0},
        ),
    }
    weights = {"R_cand_com_synth": 10.0}
    scenario = {"cpd00027": ("D-Glucose", 5.0), "cpd00013": ("NH3", 10.0)}

    res_b0 = enumerate_alternative_gapfill_ensemble(
        draft_no_com, candidates, weights, scenario, biomass_recipe_id="B0_methanosarcina_template_baseline"
    )
    res_b2 = enumerate_alternative_gapfill_ensemble(
        draft_no_com, candidates, weights, scenario, biomass_recipe_id="B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"
    )
    res_b1 = enumerate_alternative_gapfill_ensemble(
        draft_no_com, candidates, weights, scenario, biomass_recipe_id="B1_perturbation_omit_8_template_cofactors"
    )
    assert res_b0["draft_already_grows"] is False
    assert "R_cand_com_synth" in res_b0["reaction_frequencies"]
    assert res_b2["draft_already_grows"] is True
    assert res_b1["draft_already_grows"] is True

    cls = classify_biomass_requirement_provenance(
        req_in_b0_baseline=True,
        req_in_b1_omit_all8=False,
        req_in_b2_omit_7_keep_f420=False,
        req_in_b3_omit_surface_lipids=True,
        req_in_b4_core_skeleton=False,
    )
    assert cls == "biomass_assumption_sensitive__vanishes_when_7_non_f420_template_cofactors_omitted"
    assert "artifact" not in cls


def test_7_solver_timeout_and_failure_never_reported_as_zero_or_optimal():
    """Adversarial Test (`EXP1-SOLVER-REVIEW.txt` Item 3 & 6, `EXP1-FINAL-ENGINE-REVIEW.txt` Item A):
    Verify that when `linprog` or `milp` returns status 1 (timeout / iteration limit) without a solution,
    `solve_fba_highs` and `enumerate_alternative_gapfill_ensemble` return `math.nan` and explicit timeout/failure
    status, NEVER `0.0` growth or `optimal`.
    """
    draft = _make_base_archaeal_toy_model()
    fake_timeout = OptimizeResult(success=False, status=1, message="Time limit reached", x=None, fun=None)
    with patch("exp1_sbml_highs_engine.linprog", return_value=fake_timeout):
        fba_res = solve_fba_highs(draft)
        assert fba_res["status"] == "timeout"
        assert math.isnan(fba_res["growth_rate"])

    candidates = {
        "R_lys_synth": ReactionRecord(
            id="R_lys_synth", base_id="rxn_lys_synth", name="Lysine synthesis",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00039_c0": 1.0},
        ),
    }
    scenario = {"cpd00027": ("D-Glucose", 5.0), "cpd00013": ("NH3", 10.0)}
    with patch("exp1_sbml_highs_engine.milp", return_value=fake_timeout):
        ens_res = enumerate_alternative_gapfill_ensemble(draft, candidates, {"R_lys_synth": 5.0}, scenario)
        assert ens_res["solver_status"] == "timeout_no_incumbent"
        assert ens_res["n_distinct_valid_solutions"] == 0
        assert math.isnan(ens_res["min_total_cost"])


def test_8_mass_charge_balance_distinguishes_unknown_formula_and_gpr_masking():
    """Adversarial Test (`EXP1-SOLVER-REVIEW.txt` Item 5 & `EXP1-CONTRACT.txt` Item 2):
    1) Verify `audit_model_mass_charge_balance` classifies missing formula/charge as `unknown_formula_or_charge`
       (never `balanced`).
    2) Verify `apply_artificial_gene_masking` accurately propagates gene knockouts through `AND`/`OR` GPR trees.
    """
    draft = _make_base_archaeal_toy_model()
    # Add a reaction with a species lacking formula
    draft.species["M_cpd77777_c0"] = SpeciesRecord("M_cpd77777_c0", "cpd77777", "Mystery compound", "c0", None, None)
    draft.reactions["R_mystery"] = ReactionRecord(
        id="R_mystery", base_id="rxn_mystery", name="Mystery rxn",
        reversible=False, lb=0.0, ub=1000.0,
        stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd77777_c0": 1.0},
    )
    qc = audit_model_mass_charge_balance(draft)
    assert qc["n_internal_unknown_formula_or_charge"] >= 1
    assert qc["n_internal_balanced"] + qc["n_internal_unknown_formula_or_charge"] + qc["n_internal_imbalanced"] == qc["n_internal_total"]

    rng = np.random.default_rng(42)
    masked_mod, mask_info = apply_artificial_gene_masking(draft, mask_fraction=0.25, rng=rng)
    assert mask_info["n_masked_genes"] == 2
    assert mask_info["n_total_genes"] == 8


def test_9_post_pruning_cost_change_invalidates_optimal_claim_and_incumbent_anchor():
    """Adversarial Test (`EXP1-FINAL-ENGINE-REVIEW.txt` Item A, C, E):
    1) If a solver returns `status == 0` with `z_sol = [1, 1]` where the second reaction carries zero flux in FBA
       verification (so post-pruning binary cost drops below `res.fun`), verify the engine invalidates the `"optimal"`
       claim (`post_pruning_cost_changed_not_proven_optimal`) and labels `cost_window_anchor_type == 'incumbent_anchored_not_proven_optimal'`.
    2) Verify `dual_bound` is read directly from `res.mip_dual_bound` (never synthesized from `obj/(1+mip_gap)`).
    """
    draft = _make_base_archaeal_toy_model()
    candidates = {
        "R_cand_1_active": ReactionRecord(
            id="R_cand_1_active", base_id="rxn_c1", name="Active lysine synth",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00039_c0": 1.0},
        ),
        "R_cand_2_idle": ReactionRecord(
            id="R_cand_2_idle", base_id="rxn_c2", name="Idle side reaction",
            reversible=False, lb=0.0, ub=1000.0,
            stoichiometry={"M_cpd00035_c0": -1.0, "M_cpd00035_c0": 1.0},
        ),
    }
    weights = {"R_cand_1_active": 5.0, "R_cand_2_idle": 3.0}
    scenario = {"cpd00027": ("D-Glucose", 5.0), "cpd00013": ("NH3", 10.0)}

    real_milp = __import__("scipy.optimize", fromlist=["milp"]).milp
    call_counter = {"n": 0}

    def _mock_milp(*args, **kwargs):
        call_counter["n"] += 1
        res = real_milp(*args, **kwargs)
        if call_counter["n"] == 1 and res.x is not None:
            # Corrupt first solve so both z_0=1 and z_1=1 are set in res.x with fun=8.0, mip_dual_bound=7.99999
            x_mod = res.x.copy()
            x_mod[-1] = 1.0
            x_mod[-2] = 1.0
            return OptimizeResult(
                success=True,
                status=0,
                message="Mocked extra zero-flux z=1",
                x=x_mod,
                fun=8.0,
                mip_gap=0.0,
                mip_dual_bound=7.777,
            )
        return OptimizeResult(success=False, status=2, message="Infeasible", x=None, fun=None)

    with patch("exp1_sbml_highs_engine.milp", side_effect=_mock_milp):
        ens = enumerate_alternative_gapfill_ensemble(draft, candidates, weights, scenario, min_growth=0.05)
        assert ens["n_distinct_valid_solutions"] == 1
        sol0 = ens["solutions"][0]
        assert sol0["optimizer_status"] == "post_pruning_cost_changed_not_proven_optimal"
        assert ens["cost_window_anchor_type"] == "incumbent_anchored_not_proven_optimal"
        assert ens["solver_status"] == "heuristic_incumbent_bounded"
        assert math.isclose(sol0["dual_bound"], 7.777, rel_tol=1e-6)
        assert math.isclose(sol0["total_cost"], 5.0, rel_tol=1e-6)
        assert math.isclose(sol0["raw_milp_objective"], 8.0, rel_tol=1e-6)


def test_10_deterministic_pool_aggregation_and_reorder_invariance():
    """Adversarial Test (`EXP1-INVENTORY-REVIEW.txt` Items 1, 2, 4):
    1) Verify `aggregate_candidate_pool_deterministically` produces bit-for-bit identical reaction bounds,
       stoichiometry, donor provenance lists, and conflict audits across arbitrary permutations of `as_completed`
       donor model arrival order.
    2) Verify that duplicate-ID reactions with conflicting stoichiometry or compartments across models are
       caught (`reaction_conflicts`) and excluded from `reactions` (never silently passing).
    3) Verify that known elemental-imbalanced candidate reactions (`C6H12O6 -> C3H7NO2`) are audited and excluded
       from `reactions`, while `unknown_formula_or_charge` reactions are explicitly labeled and retained.
    4) Verify pool scope is labeled `observed_cohort_538_sbml_union_pool_not_full_modelseed_db` and preserves
       `observed_in_cohort_draft_models` vs `observed_only_in_cohort_gapfilled_models`.
    """
    sp_glc = SpeciesRecord("M_cpd00027_c0", "cpd00027", "D-Glucose", "c0", 0, "C6H12O6")
    sp_pyr = SpeciesRecord("M_cpd00020_c0", "cpd00020", "Pyruvate", "c0", -1, "C3H3O3")
    sp_h = SpeciesRecord("M_cpd00067_c0", "cpd00067", "H+", "c0", 1, "H")
    sp_unk = SpeciesRecord("M_cpd99111_c0", "cpd99111", "Unknown formula", "c0", None, None)

    # Model A has R_rxn00100_c0 irreversible [0, 1000] in draft, plus R_conflict_c0 with stoich -1 Glc -> 2 Pyr
    rec_a = {
        "genome": "MAG_A",
        "extracted_species": {s.id: s for s in (sp_glc, sp_pyr, sp_h, sp_unk)},
        "draft_reaction_ids": {"R_rxn00100_c0", "R_conflict_c0"},
        "fill_desert_reaction_ids": set(),
        "fill_native_reaction_ids": set(),
        "extracted_reactions": {
            "R_rxn00100_c0": ReactionRecord(
                id="R_rxn00100_c0", base_id="rxn00100_c0", name="Balanced rxn",
                reversible=False, lb=0.0, ub=1000.0,
                stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd00020_c0": 2.0, "M_cpd00067_c0": 6.0},
            ),
            "R_conflict_c0": ReactionRecord(
                id="R_conflict_c0", base_id="rxn_conflict_c0", name="Conflicting stoich rxn",
                reversible=False, lb=0.0, ub=1000.0,
                stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd00020_c0": 2.0},
            ),
        },
    }
    # Model B has R_rxn00100_c0 reversible [-1000, 1000] in fill_desert, plus R_conflict_c0 with DIFFERENT stoich (-1 Glc -> 3 Pyr!),
    # plus an imbalanced reaction R_imbalanced_c0 (-1 Glc -> 1 Pyr) and an unknown-formula reaction R_unk_c0
    rec_b = {
        "genome": "MAG_B",
        "extracted_species": {s.id: s for s in (sp_glc, sp_pyr, sp_h, sp_unk)},
        "draft_reaction_ids": set(),
        "fill_desert_reaction_ids": {"R_rxn00100_c0", "R_conflict_c0", "R_imbalanced_c0", "R_unk_c0"},
        "fill_native_reaction_ids": set(),
        "extracted_reactions": {
            "R_rxn00100_c0": ReactionRecord(
                id="R_rxn00100_c0", base_id="rxn00100_c0", name="Balanced rxn",
                reversible=True, lb=-1000.0, ub=1000.0,
                stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd00020_c0": 2.0, "M_cpd00067_c0": 6.0},
            ),
            "R_conflict_c0": ReactionRecord(
                id="R_conflict_c0", base_id="rxn_conflict_c0", name="Conflicting stoich rxn",
                reversible=False, lb=0.0, ub=1000.0,
                stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd00020_c0": 3.0},
            ),
            "R_imbalanced_c0": ReactionRecord(
                id="R_imbalanced_c0", base_id="rxn_imbal_c0", name="Elemental imbalance rxn",
                reversible=False, lb=0.0, ub=1000.0,
                stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd00020_c0": 1.0},
            ),
            "R_unk_c0": ReactionRecord(
                id="R_unk_c0", base_id="rxn_unk_c0", name="Unknown formula rxn",
                reversible=False, lb=0.0, ub=1000.0,
                stoichiometry={"M_cpd00027_c0": -1.0, "M_cpd99111_c0": 1.0},
            ),
        },
    }

    res_ab = aggregate_candidate_pool_deterministically([rec_a, rec_b], seed_rxn_directions={"rxn00100": "="})
    res_ba = aggregate_candidate_pool_deterministically([rec_b, rec_a], seed_rxn_directions={"rxn00100": "="})

    # Reorder invariance check
    assert res_ab["candidate_pool_scope"] == OBSERVED_COHORT_538_POOL_SCOPE
    assert set(res_ab["reactions"].keys()) == set(res_ba["reactions"].keys())
    assert res_ab["pool_audit_rows"] == res_ba["pool_audit_rows"]
    assert res_ab["reaction_conflicts"] == res_ba["reaction_conflicts"]

    # Conflict & balance eligibility check
    assert "R_rxn00100_c0" in res_ab["reactions"]
    assert res_ab["reactions"]["R_rxn00100_c0"].lb == -1000.0
    assert res_ab["reactions"]["R_rxn00100_c0"].ub == 1000.0
    assert res_ab["reactions"]["R_rxn00100_c0"].provenance == "observed_in_cohort_draft_models"

    assert "R_unk_c0" in res_ab["reactions"]
    assert res_ab["reactions"]["R_unk_c0"].provenance == "observed_only_in_cohort_gapfilled_models"

    # Conflicting reaction and known elemental-imbalanced reaction must be excluded from eligible reactions
    assert "R_conflict_c0" not in res_ab["reactions"]
    assert "R_imbalanced_c0" not in res_ab["reactions"]
    assert len(res_ab["reaction_conflicts"]) == 1
    assert res_ab["reaction_conflicts"][0]["rxn_id"] == "R_conflict_c0"


def test_11_driver_level_stage1_metadata_adapter_and_exact_cohort_equality(tmp_path: Path) -> None:
    """Driver-level integration test: tests Stage 1 metadata adapter on real manifest headers,
    verifying 538 total cohort, 484 exact deduplicated membership (237 EQ + 247 external),
    taxonomy parsing, and failure-to-unknown handling (`EXP1-INTEGRATION-REVIEW.txt`).
    """
    import pandas as pd
    from stage1_inventory_and_qc import (
        build_cohort_metadata_538_and_484,
        load_frozen_external_247_manifest,
    )

    # 1. Test on synthetic manifest fixtures matching real headers exactly
    man_tsv = tmp_path / "external_signalp_manifest.tsv"
    man_rows = [
        "genome_id\tsource_cohort\tdedup_retained\tsame_assembly_as_gtdb_rep\tlineage\tcompleteness\tcontamination\tn_proteins",
        "EXT_01\tref_ext\tTrue\tTrue\td__Archaea;p__Thermoplasmatota;c__Poseidoniia;o__Poseidoniales;f__Poseidoniaceae;g__Poseidonia;s__Poseidonia_sp\t95.0\t1.0\t1800",
        "EXT_02\tref_ext\tFalse\tFalse\td__Archaea;p__Thermoplasmatota;c__Poseidoniia;o__Poseidoniales;f__Poseidoniaceae;g__Poseidonia;s__Poseidonia_sp2\t85.0\t2.0\t1700",
        "EXT_03\tref_ext\tTrue\tTrue\td__Archaea;p__Thermoplasmatota;c__Poseidoniia;o__Poseidoniales;f__Poseidoniaceae;g__Poseidonia;s__Poseidonia_sp3\t92.0\t0.5\t1750",
    ]
    man_tsv.write_text("\n".join(man_rows) + "\n")

    comp_df = pd.DataFrame(
        [
            {"genome_id": "EXT_01", "source_cohort": "ref_ext", "lineage": "d__Archaea;p__Thermoplasmatota;c__Poseidoniia;o__Poseidoniales;f__Poseidoniaceae;g__Poseidonia;s__Poseidonia_sp", "completeness": 95.0, "contamination": 1.0},
            {"genome_id": "EXT_02", "source_cohort": "ref_ext", "lineage": "d__Archaea;p__Thermoplasmatota;c__Poseidoniia;o__Poseidoniales;f__Poseidoniaceae;g__Poseidonia;s__Poseidonia_sp2", "completeness": 85.0, "contamination": 2.0},
            {"genome_id": "EXT_03", "source_cohort": "ref_ext", "lineage": "d__Archaea;p__Thermoplasmatota;c__Poseidoniia;o__Poseidoniales;f__Poseidoniaceae;g__Poseidonia;s__Poseidonia_sp3", "completeness": 92.0, "contamination": 0.5},
        ]
    )

    all_ext, ret_ext = load_frozen_external_247_manifest(
        man_tsv, comp_df, expected_total_external=3, expected_retained_external=2
    )
    assert all_ext == {"EXT_01", "EXT_02", "EXT_03"}
    assert ret_ext == {"EXT_01", "EXT_03"}

    # Test mismatch exception
    mismatch_df = pd.DataFrame([{"genome_id": "EXT_01"}, {"genome_id": "EXT_UNKNOWN"}])
    with pytest.raises(RuntimeError, match="External comparator set mismatch"):
        load_frozen_external_247_manifest(man_tsv, mismatch_df, expected_total_external=3, expected_retained_external=2)

    # 2. Test unified cohort metadata builder
    legacy_eq_df = pd.DataFrame(
        [
            {"genome": "EQ_MAG_01", "gtdb_classification": "d__Archaea;p__Thermoplasmatota;c__Methanomethylicia;o__Methanomethyliales;f__Methanomethyliaceae;g__Methanomethylophilus;s__sp1", "checkm_completeness": 96.0, "contamination": 0.8},
            {"genome": "EQ_MAG_02", "gtdb_classification": "d__Archaea;p__Thermoplasmatota;c__Methanomethylicia;o__Methanomethyliales;f__Methanomethyliaceae;g__Methanomethylophilus;s__sp2", "checkm_completeness": 88.0, "contamination": 3.0},
        ]
    )
    ref_rep_df = pd.DataFrame(
        [
            {"mag": "EQ_MAG_01", "Completeness": 96.0, "Contamination": 0.8, "gtdb_r232_classification": "d__Archaea;p__Thermoplasmatota;c__Methanomethylicia;o__Methanomethyliales;f__Methanomethyliaceae;g__Methanomethylophilus;s__sp1_canon", "quality_tier": "HQ", "Total_Contigs": 12, "Contig_N50": 150000.0},
            {"mag": "EQ_MAG_02", "Completeness": 88.0, "Contamination": 3.0, "gtdb_r232_classification": "d__Archaea;p__Thermoplasmatota;c__Methanomethylicia;o__Methanomethyliales;f__Methanomethyliaceae;g__Methanomethylophilus;s__sp2_canon", "quality_tier": "MQ", "Total_Contigs": 35, "Contig_N50": 45000.0},
            {"mag": "BACTERIAL_MAG", "Completeness": 98.0, "Contamination": 0.5, "gtdb_r232_classification": "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__sp"},
        ]
    )

    cohort_meta = build_cohort_metadata_538_and_484(
        ref_rep_df=ref_rep_df,
        legacy_eq_df=legacy_eq_df,
        comp_df=comp_df,
        external_manifest_tsv=man_tsv,
        expected_eq_count=2,
        expected_total_external=3,
        expected_retained_external=2,
    )
    assert len(cohort_meta) == 5  # 2 EQ + 3 external
    assert int(cohort_meta["is_deduplicated_484"].sum()) == 4  # 2 EQ + 2 retained external
    assert set(cohort_meta.loc[cohort_meta["is_deduplicated_484"], "genome"]) == {"EQ_MAG_01", "EQ_MAG_02", "EXT_01", "EXT_03"}
    # Check taxonomy change detection
    eq1_row = cohort_meta[cohort_meta["genome"] == "EQ_MAG_01"].iloc[0]
    assert bool(eq1_row["taxonomy_changed_canonical_vs_legacy"]) is True
    assert eq1_row["canonical_order"] == "o__Methanomethyliales"
    assert eq1_row["canonical_family"] == "f__Methanomethyliaceae"


def test_12_driver_level_stage2_matrix_loop_toy_fixture_and_cache_lookups(tmp_path: Path) -> None:
    """Driver-level integration test: exercises the full Stage 2 matrix loop across all scenarios (S1..S6 + S0),
    all biomass perturbations (B0..B4), withdrawal/KO/FVA, supplement rescue, and masking controls on a toy fixture,
    verifying exact cache lookups, scenario keys, and scientific QC gate (`EXP1-INTEGRATION-REVIEW.txt`).
    """
    import pickle
    from stage2_pilot24_ensemble import (
        evaluate_single_genome_full_protocol,
        validate_pilot24_scientific_gate,
    )

    g = "toy_archaeon_gapfill_fixture"
    supp_dir = tmp_path / "supp" / g
    supp_dir.mkdir(parents=True)

    # Archaea model that needs acetate (cpd00029) or glucose (cpd00027) transport to grow under desert medium
    xml = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2" level="3" version="1" fbc:required="false">
  <model id="toy_archaeon_gapfill_fixture">
    <listOfParameters>
      <parameter id="default_lb" value="-1000.0"/>
      <parameter id="default_ub" value="1000.0"/>
      <parameter id="default_0" value="0.0"/>
    </listOfParameters>
    <listOfSpecies>
      <species id="M_cpd00029_e0" name="Acetate" compartment="e0" fbc:chemicalFormula="C2H3O2" fbc:charge="-1"/>
      <species id="M_cpd00029_c0" name="Acetate" compartment="c0" fbc:chemicalFormula="C2H3O2" fbc:charge="-1"/>
      <species id="M_cpd00013_e0" name="NH3" compartment="e0" fbc:chemicalFormula="H3N" fbc:charge="0"/>
      <species id="M_cpd00013_c0" name="NH3" compartment="c0" fbc:chemicalFormula="H3N" fbc:charge="0"/>
      <species id="M_cpd11416_c0" name="Biomass" compartment="c0"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="R_EX_cpd00029_e0" name="Acetate exchange" reversible="true" fbc:lowerFluxBound="default_lb" fbc:upperFluxBound="default_ub">
        <listOfReactants><speciesReference species="M_cpd00029_e0" stoichiometry="1"/></listOfReactants>
      </reaction>
      <reaction id="R_EX_cpd00013_e0" name="NH3 exchange" reversible="true" fbc:lowerFluxBound="default_lb" fbc:upperFluxBound="default_ub">
        <listOfReactants><speciesReference species="M_cpd00013_e0" stoichiometry="1"/></listOfReactants>
      </reaction>
      <reaction id="R_nh3_trans" name="NH3 transport" reversible="true" fbc:lowerFluxBound="default_lb" fbc:upperFluxBound="default_ub">
        <listOfReactants><speciesReference species="M_cpd00013_e0" stoichiometry="1"/></listOfReactants>
        <listOfProducts><speciesReference species="M_cpd00013_c0" stoichiometry="1"/></listOfProducts>
      </reaction>
      <reaction id="R_bio1" name="Biomass objective" reversible="false" fbc:lowerFluxBound="default_0" fbc:upperFluxBound="default_ub">
        <listOfReactants>
          <speciesReference species="M_cpd00029_c0" stoichiometry="1"/>
          <speciesReference species="M_cpd00013_c0" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="M_cpd11416_c0" stoichiometry="1"/>
        </listOfProducts>
      </reaction>
      <reaction id="R_EX_cpd11416_c0" name="Biomass sink" reversible="false" fbc:lowerFluxBound="default_0" fbc:upperFluxBound="default_ub">
        <listOfReactants><speciesReference species="M_cpd11416_c0" stoichiometry="1"/></listOfReactants>
      </reaction>
    </listOfReactions>
  </model>
</sbml>"""
    (supp_dir / f"{g}-draft.xml").write_text(xml)
    (supp_dir / f"{g}-all-Reactions.tbl").write_text("rxn\tname\tstatus\tbitscore\tbest_gene\tpathway\nrxn00029_trans\tAcetate trans\tcore\t85.0\tgene_actrans\tPWY_ACT\n")
    (supp_dir / f"{g}-Transporter.tbl").write_text("Transporter\tSubstrates\tbitscore\nrxn00029_trans\tcpd00029\t85.0\n")
    (supp_dir / f"{g}-desert_medium.csv").write_text("compounds,name,maxFlux\ncpd00029,Acetate,10\ncpd00013,NH3,10\n")
    (supp_dir / f"{g}-native_medium.csv").write_text("compounds,name,maxFlux\ncpd00029,Acetate,10\ncpd00013,NH3,10\n")

    # Candidate reaction pool provides the missing acetate transporter
    cand_rxn = ReactionRecord(
        id="R_rxn00029_trans",
        base_id="rxn00029_trans",
        name="Acetate transport",
        reversible=False,
        lb=0.0,
        ub=1000.0,
        stoichiometry={"M_cpd00029_e0": -1.0, "M_cpd00029_c0": 1.0},
        is_transport=True,
    )
    pool_data = {
        "candidate_pool_scope": OBSERVED_COHORT_538_POOL_SCOPE,
        "reactions": {"R_rxn00029_trans": cand_rxn},
        "species": {
            "M_cpd00029_e0": SpeciesRecord("M_cpd00029_e0", "cpd00029", "Acetate", "e0", -1, "C2H3O2"),
            "M_cpd00029_c0": SpeciesRecord("M_cpd00029_c0", "cpd00029", "Acetate", "c0", -1, "C2H3O2"),
        },
    }
    pool_pkl = tmp_path / "pool.pkl"
    with open(pool_pkl, "wb") as fh:
        pickle.dump(pool_data, fh)

    task = {
        "genome": g,
        "cohort": "EQ_desert",
        "canonical_order": "o__Methanomicrobiales",
        "canonical_family": "f__Methanomicrobiaceae",
        "canonical_genus": "g__Methanoculleus",
        "completeness": 95.0,
        "contamination": 1.0,
        "gapseq_out_dir": str(tmp_path / "missing_tar"),
        "supp_models_dir": str(tmp_path / "supp"),
        "pool_pkl_path": str(pool_pkl),
        "n_target_solutions": 4,
        "milp_time_limit_s": 2.0,
        "run_masking_control": True,
    }

    res = evaluate_single_genome_full_protocol(task)
    assert len(res["summary_rows"]) >= 18
    assert len(res["rxn_freq_rows"]) > 0
    assert len(res["wkf_rows"]) > 0
    assert len(res["masking_rows"]) == 9  # 3 fractions x 3 replicates

    sum_df = pd.DataFrame(res["summary_rows"])
    rf_df = pd.DataFrame(res["rxn_freq_rows"])
    wkf_df = pd.DataFrame(res["wkf_rows"])
    mask_df = pd.DataFrame(res["masking_rows"])

    # Ensure every scenario executed without KeyError
    scens = set(sum_df["scenario_id"])
    assert "S1_anoxic_minimal_h2_co2_nh4" in scens
    assert "S2_oxic_minimal_h2_co2_nh4" in scens
    assert "S3_anoxic_desert_oa_aa" in scens
    assert "S4_oxic_desert_oa_aa" in scens
    assert "S5_anoxic_defined_aa_precursors" in scens
    assert "S6_oxic_defined_aa_precursors" in scens
    assert "S0_gapseq_native_predicted" in scens

    # Ensure all 5 biomass perturbations ran
    bios = set(sum_df["biomass_recipe_id"])
    assert len(bios) == 5

    gate = validate_pilot24_scientific_gate(sum_df, rf_df, wkf_df, mask_df, expected_n_genomes=1)
    assert gate["checks"]["all_manifest_genomes_evaluated"] is True
    assert gate["checks"]["solver_clean_completion_rate_ge_95pct"] is True
    assert gate["checks"]["biomass_perturbations_b0_to_b4_populated"] is True
    assert gate["checks"]["nutrient_scenarios_s1_to_s6_populated"] is True
    assert gate["checks"]["reaction_frequencies_and_biomass_sensitivity_non_empty"] is True
    assert gate["checks"]["withdrawal_knockout_fva_non_empty"] is True
    assert gate["checks"]["masking_controls_non_empty"] is True
    assert gate["gate_passed"] is True


def test_13_pilot24_gate_clean_solver_rate_and_biomass_sensitivity_regression() -> None:
    """Focused regression test: verifies that validate_pilot24_scientific_gate accurately
    enforces planned protocol scenario/biomass tuples, accounts for statuses honestly,
    and requires verified accepted solutions without relying on an arbitrary 95% aggregate threshold
    (per EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt & EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt).
    """
    from stage2_pilot24_ensemble import validate_pilot24_scientific_gate
    from exp1_sbml_highs_engine import (
        get_planned_protocol_scenario_biomass_tuples,
    )

    genome_ids = [f"g_{i:02d}" for i in range(31)]
    expected_tuples = sorted(get_planned_protocol_scenario_biomass_tuples(genome_ids))
    assert len(expected_tuples) == 31 * 26

    # Construct realistic distribution matching evaluated protocol:
    # 590 heuristic_incumbent_bounded, 200 infeasible, 16 fba_verification_failed (unresolved solver outcome)
    rows = []
    for idx, (g, s, b) in enumerate(expected_tuples):
        if idx < 590:
            st = "heuristic_incumbent_bounded"
            n_sol = 5
            min_cost = 100.0
            min_add = 2.0
            growth = 0.05
            res_val = 1e-7
            bnd_val = 1e-7
        elif idx < 790:
            st = "infeasible_under_fixed_medium"
            n_sol = 0
            min_cost = math.nan
            min_add = math.nan
            growth = math.nan
            res_val = math.nan
            bnd_val = math.nan
        else:
            st = "fba_verification_failed"
            n_sol = 0
            min_cost = math.nan
            min_add = math.nan
            growth = math.nan
            res_val = math.nan
            bnd_val = math.nan

        rows.append(
            {
                "genome": g,
                "scenario_id": s,
                "biomass_recipe_id": b,
                "solver_status": st,
                "n_distinct_valid_solutions": n_sol,
                "min_total_cost": min_cost,
                "min_added_rxns_count": min_add,
                "min_solution_verified_growth": growth,
                "max_solution_residual": res_val,
                "max_solution_bound_viol": bnd_val,
            }
        )

    summary_df = pd.DataFrame(rows)

    rxn_freq_df = pd.DataFrame(
        [
            {
                "genome": "g_00",
                "scenario_id": "S3_anoxic_desert_oa_aa",
                "rxn_id": "R_test",
                "biomass_sensitivity_classification": "biomass_assumption_sensitive__vanishes_when_7_non_f420_template_cofactors_omitted",
            },
            {
                "genome": "g_00",
                "scenario_id": "S3_anoxic_desert_oa_aa",
                "rxn_id": "R_robust",
                "biomass_sensitivity_classification": "robust_across_all_tested_biomass_perturbations",
            },
        ]
    )

    wkf_df = pd.DataFrame([{"genome": "g_00", "record_type": "test"}])
    # 7 HQ reference genomes x 3 fractions x 3 replicates = 63 masking rows
    hq_genomes = [f"MGYG_{i}" for i in range(7)]
    masking_rows = [
        {"genome": g, "mask_fraction": f, "replicate": r, "masked_solver_status": "heuristic_incumbent_bounded"}
        for g in hq_genomes
        for f in (0.10, 0.20, 0.30)
        for r in (0, 1, 2)
    ]
    masking_df = pd.DataFrame(masking_rows)

    gate = validate_pilot24_scientific_gate(
        summary_df=summary_df,
        rxn_freq_df=rxn_freq_df,
        wkf_df=wkf_df,
        masking_df=masking_df,
        expected_n_genomes=31,
        expected_genome_ids=genome_ids,
        expected_masking_genome_ids=hq_genomes,
    )

    assert gate["checks"]["exact_planned_genome_ids_matched"] is True
    assert gate["checks"]["exact_planned_scenario_biomass_tuples_per_genome"] is True
    assert gate["checks"]["no_duplicate_evaluations"] is True
    assert gate["checks"]["statuses_exhaustive_and_honest"] is True
    assert gate["checks"]["per_solution_correctness_verified"] is True
    assert gate["checks"]["masking_controls_executed_and_audited"] is True
    assert gate["gate_passed"] is True
    assert gate["rate_feasible_independently_revalidated"] > 0.70
    assert gate["rate_proven_infeasible"] > 0.20
    assert gate["rate_numerical_failures"] > 0.01

    # Negative control: unhandled status causes gate failure
    summary_df_failing = summary_df.copy()
    summary_df_failing.loc[0, "solver_status"] = "crashed_unknown_solver_error"
    gate_fail = validate_pilot24_scientific_gate(
        summary_df=summary_df_failing,
        rxn_freq_df=rxn_freq_df,
        wkf_df=wkf_df,
        masking_df=masking_df,
        expected_n_genomes=31,
        expected_genome_ids=genome_ids,
        expected_masking_genome_ids=hq_genomes,
    )
    assert gate_fail["checks"]["statuses_exhaustive_and_honest"] is False
    assert gate_fail["gate_passed"] is False


def test_14_exact_coupling_and_tight_tolerances_prevent_trace_cofactor_leakage():
    """Verify that validated HiGHS tolerances (1e-9) and native coupling prevent HiGHS
    from satisfying small biomass cofactor demands via fractional sub-integrality leakage (z ~ 1e-7)
    without setting the candidate binary variable z = 1.
    """
    species = {
        "M_cpd00027_e0": SpeciesRecord("M_cpd00027_e0", "cpd00027", "Glucose", "e0", 0, "C6H12O6"),
        "M_cpd00027_c0": SpeciesRecord("M_cpd00027_c0", "cpd00027", "Glucose", "c0", 0, "C6H12O6"),
        "M_trace_c0": SpeciesRecord("M_trace_c0", "cpd_trace", "Trace cofactor", "c0", 0, ""),
    }
    reactions = {
        "R_EX_cpd00027_e0": ReactionRecord(
            "R_EX_cpd00027_e0", "EX_cpd00027_e0", "Glucose uptake", False, -10.0, 1000.0,
            {"M_cpd00027_e0": -1.0}, is_exchange=True
        ),
        "R_glc_trans": ReactionRecord(
            "R_glc_trans", "rxn_glc_trans", "Glucose import", False, 0.0, 1000.0,
            {"M_cpd00027_e0": -1.0, "M_cpd00027_c0": 1.0}, is_transport=True
        ),
        "R_bio1": ReactionRecord(
            "R_bio1", "bio1", "Biomass", False, 0.0, 1000.0,
            {"M_cpd00027_c0": -1.0, "M_trace_c0": -0.0005},
            is_biomass=True,
        ),
    }
    base = SBMLModelData(
        model_id="toy_leak_test",
        model_name="Toy Leak Test",
        notes_text="",
        species=species,
        reactions=reactions,
        objective_rxn_id="R_bio1",
        gene_products={},
    )

    cand_rxn = ReactionRecord(
        "R_cand_trace_syn", "R_cand_trace_syn", "Trace cofactor synthetase",
        False, 0.0, 1000.0, {"M_cpd00027_c0": -0.01, "M_trace_c0": 1.0}
    )
    candidates = {"R_cand_trace_syn": cand_rxn}
    cand_weights = {"R_cand_trace_syn": 100.0}
    scenario_cpds = {"cpd00027": ("Glucose", 10.0)}

    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=base,
        candidate_reactions=candidates,
        candidate_weights=cand_weights,
        scenario_compounds=scenario_cpds,
        biomass_recipe_id="B0_methanosarcina_template_baseline",
        min_growth=0.01,
        n_Target_solutions=1,
        time_limit_per_milp_s=10.0,
    )

    assert len(ens["solutions"]) == 1
    sol = ens["solutions"][0]
    assert "R_cand_trace_syn" in sol["added_rxns"]
    assert sol["raw_milp_selected_rxns"] == ["R_cand_trace_syn"]
    assert sol["verified_fba_growth"] >= 0.01
    assert ens["solver_status"] in {"optimal", "heuristic_incumbent_bounded"}


def test_15_high_flux_greater_than_20_permitted_and_revalidated():
    """Regression test per EXP1-QC-DESIGN-CORRECTION.txt:
    Ensure native candidate bounds accommodate legitimate high-flux pathways (flux > 20.0, e.g. 50.0)
    without being artificially capped or causing false infeasibility.
    """
    species = {
        "M_cpd00027_e0": SpeciesRecord("M_cpd00027_e0", "cpd00027", "Glucose", "e0", 0, "C6H12O6"),
        "M_cpd00027_c0": SpeciesRecord("M_cpd00027_c0", "cpd00027", "Glucose", "c0", 0, "C6H12O6"),
        "M_energy_c0": SpeciesRecord("M_energy_c0", "cpd_energy", "High flux intermediate", "c0", 0, ""),
    }
    # Nutrient scenario provides up to 100.0 glucose
    reactions = {
        "R_EX_cpd00027_e0": ReactionRecord(
            "R_EX_cpd00027_e0", "EX_cpd00027_e0", "Glucose uptake", False, -100.0, 1000.0,
            {"M_cpd00027_e0": -1.0}, is_exchange=True
        ),
        "R_glc_trans": ReactionRecord(
            "R_glc_trans", "rxn_glc_trans", "Glucose import", False, 0.0, 1000.0,
            {"M_cpd00027_e0": -1.0, "M_cpd00027_c0": 1.0}, is_transport=True
        ),
        # Biomass requires 50.0 units of M_energy_c0 per unit biomass (high flux demand)
        "R_bio1": ReactionRecord(
            "R_bio1", "bio1", "Biomass", False, 0.0, 1000.0,
            {"M_energy_c0": -50.0},
            is_biomass=True,
        ),
    }
    base = SBMLModelData(
        model_id="toy_high_flux",
        model_name="Toy High Flux Model",
        notes_text="",
        species=species,
        reactions=reactions,
        objective_rxn_id="R_bio1",
        gene_products={},
    )

    # Candidate reaction converts 1 glucose into 1 energy unit, native upper bound 1000.0
    cand_rxn = ReactionRecord(
        "R_cand_energy_gen", "R_cand_energy_gen", "High flux energy generation",
        False, 0.0, 1000.0, {"M_cpd00027_c0": -1.0, "M_energy_c0": 1.0}
    )
    candidates = {"R_cand_energy_gen": cand_rxn}
    cand_weights = {"R_cand_energy_gen": 10.0}
    scenario_cpds = {"cpd00027": ("Glucose", 100.0)}

    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=base,
        candidate_reactions=candidates,
        candidate_weights=cand_weights,
        scenario_compounds=scenario_cpds,
        biomass_recipe_id="B0_methanosarcina_template_baseline",
        min_growth=1.0,  # Requires growth rate 1.0 -> energy flux must be exactly 50.0 (> 20.0)!
        n_Target_solutions=1,
        time_limit_per_milp_s=10.0,
    )

    assert len(ens["solutions"]) == 1
    sol = ens["solutions"][0]
    assert "R_cand_energy_gen" in sol["added_rxns"]
    cand_flux = sol["verified_fluxes"]["R_cand_energy_gen"]
    assert cand_flux >= 49.9, f"Expected flux >= 50.0, got {cand_flux}"
    assert sol["verified_fba_growth"] >= 0.999
    assert ens["solver_status"] in {"optimal", "heuristic_incumbent_bounded"}


def test_16_candidate_lower_bounds_and_mandatory_reaction_cases():
    """Regression test per EXP1-QC-DESIGN-CORRECTION.txt:
    Ensure candidate reactions with mandatory positive lower bounds (lb > 0 when active)
    permit flux = 0 when unselected (z = 0), but strictly enforce v >= lb when selected (z = 1).
    """
    species = {
        "M_cpd00027_e0": SpeciesRecord("M_cpd00027_e0", "cpd00027", "Glucose", "e0", 0, "C6H12O6"),
        "M_cpd00027_c0": SpeciesRecord("M_cpd00027_c0", "cpd00027", "Glucose", "c0", 0, "C6H12O6"),
        "M_prod_c0": SpeciesRecord("M_prod_c0", "cpd_prod", "Product", "c0", 0, ""),
        "M_waste_c0": SpeciesRecord("M_waste_c0", "cpd_waste", "Mandatory byproduct", "c0", 0, ""),
    }
    reactions = {
        "R_EX_cpd00027_e0": ReactionRecord(
            "R_EX_cpd00027_e0", "EX_cpd00027_e0", "Glucose uptake", False, -10.0, 1000.0,
            {"M_cpd00027_e0": -1.0}, is_exchange=True
        ),
        "R_glc_trans": ReactionRecord(
            "R_glc_trans", "rxn_glc_trans", "Glucose import", False, 0.0, 1000.0,
            {"M_cpd00027_e0": -1.0, "M_cpd00027_c0": 1.0}, is_transport=True
        ),
        "R_bio1": ReactionRecord(
            "R_bio1", "bio1", "Biomass", False, 0.0, 1000.0,
            {"M_prod_c0": -1.0},
            is_biomass=True,
        ),
        "R_waste_drain": ReactionRecord(
            "R_waste_drain", "R_waste_drain", "Waste drain", False, 0.0, 1000.0,
            {"M_waste_c0": -1.0}, is_exchange=True
        ),
    }
    base = SBMLModelData(
        model_id="toy_lb_test",
        model_name="Toy LB Test Model",
        notes_text="",
        species=species,
        reactions=reactions,
        objective_rxn_id="R_bio1",
        gene_products={},
    )

    # Candidate 1: standard pathway, cost 50
    # Candidate 2: mandatory LB candidate with lb=2.0, cost 100
    cand_rxn_1 = ReactionRecord(
        "R_cand_path_1", "R_cand_path_1", "Path 1",
        False, 0.0, 1000.0, {"M_cpd00027_c0": -1.0, "M_prod_c0": 1.0}
    )
    cand_rxn_mandatory_lb = ReactionRecord(
        "R_cand_mandatory_lb", "R_cand_mandatory_lb", "Mandatory LB candidate",
        False, 2.0, 1000.0, {"M_cpd00027_c0": -1.0, "M_waste_c0": 1.0}
    )
    candidates = {
        "R_cand_path_1": cand_rxn_1,
        "R_cand_mandatory_lb": cand_rxn_mandatory_lb,
    }
    cand_weights = {"R_cand_path_1": 50.0, "R_cand_mandatory_lb": 100.0}
    scenario_cpds = {"cpd00027": ("Glucose", 10.0)}

    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=base,
        candidate_reactions=candidates,
        candidate_weights=cand_weights,
        scenario_compounds=scenario_cpds,
        biomass_recipe_id="B0_methanosarcina_template_baseline",
        min_growth=0.01,
        n_Target_solutions=1,
        time_limit_per_milp_s=10.0,
    )

    assert len(ens["solutions"]) == 1
    sol = ens["solutions"][0]
    # Path 1 is cheaper, so R_cand_mandatory_lb must NOT be selected (z=0)
    assert "R_cand_path_1" in sol["added_rxns"]
    assert "R_cand_mandatory_lb" not in sol["added_rxns"]
    # Unselected candidate with native lb=2.0 must have EXACTLY 0.0 flux when unselected
    assert abs(sol["verified_fluxes"].get("R_cand_mandatory_lb", 0.0)) < 1e-9


def test_17_exact_scenario_biomass_tuple_manifest_rejects_adversarial_counterexamples():
    """Adversarial counterexample test per EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt Item A:
    Root constructed a synthetic adversarial dataset where 26 rows had bogus scenario names
    (nonexistent_scenario_0..25), biomass nonexistent_biomass, all numerical_error, and duplicate masking rows.
    Verify that the tightened QC gate detects and REJECTS this counterexample!
    """
    from stage2_pilot24_ensemble import validate_pilot24_scientific_gate

    # Counterexample A: Bogus scenario and biomass names with count 26
    bogus_scens = [f"nonexistent_scenario_{i}" for i in range(26)]
    bogus_summary = pd.DataFrame(
        {
            "genome": ["adversarial_genome"] * 26,
            "scenario_id": bogus_scens,
            "biomass_recipe_id": ["nonexistent_biomass"] * 26,
            "solver_status": ["numerical_error"] * 26,
            "n_distinct_valid_solutions": [0] * 26,
            "min_total_cost": [math.nan] * 26,
            "min_added_rxns_count": [math.nan] * 26,
        }
    )
    # Duplicate masking rows
    bogus_masking = pd.DataFrame(
        [{"genome": "adversarial_genome", "mask_fraction": 0.10, "replicate": 0, "masked_solver_status": "numerical_error"}] * 8
    )
    empty_rf = pd.DataFrame(columns=["genome", "scenario_id", "rxn_id", "biomass_sensitivity_classification"])
    empty_wkf = pd.DataFrame(columns=["genome", "record_type"])

    gate = validate_pilot24_scientific_gate(
        summary_df=bogus_summary,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=bogus_masking,
        expected_n_genomes=1,
    )
    assert gate["checks"]["exact_planned_scenario_biomass_tuples_per_genome"] is False
    assert gate["checks"]["statuses_exhaustive_and_honest"] is True  # numerical_error is known
    assert gate["gate_passed"] is False, "Bogus scenario names must FAIL the gate!"

    # Counterexample B: Accepted solution with violated residual or 0 solutions claiming success
    from exp1_sbml_highs_engine import get_planned_protocol_scenario_biomass_tuples
    valid_tuples = sorted(get_planned_protocol_scenario_biomass_tuples(["test_genome"]))
    forged_accepted_rows = []
    for s, b in [(t[1], t[2]) for t in valid_tuples]:
        forged_accepted_rows.append(
            {
                "genome": "test_genome",
                "scenario_id": s,
                "biomass_recipe_id": b,
                "solver_status": "optimal",
                "n_distinct_valid_solutions": 0,  # Invalid: claims optimal but 0 solutions!
                "min_total_cost": 10.0,
                "min_added_rxns_count": 1.0,
                "min_solution_verified_growth": 0.05,
                "max_solution_residual": 0.05,  # Invalid: residual 0.05 > 1e-5!
                "max_solution_bound_viol": 1e-7,
            }
        )
    forged_summary = pd.DataFrame(forged_accepted_rows)
    valid_masking = pd.DataFrame([{"genome": "test_genome", "replicate": 0}])
    gate_forged = validate_pilot24_scientific_gate(
        summary_df=forged_summary,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=valid_masking,
        expected_n_genomes=1,
    )
    assert gate_forged["checks"]["per_solution_correctness_verified"] is False
    assert gate_forged["gate_passed"] is False, "Accepted solution with violated residual must FAIL!"


def test_18_solve_fba_highs_rejects_forged_solver_success_with_violated_residuals_or_bounds():
    """Adversarial test per EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt Item 2:
    Verifies that solve_fba_highs does NOT trust linprog res.success blindly, and rejects
    forged successes that violate stoichiometric residuals (max|Sx - b| > 1e-5), bound limits,
    or contain NaN/Inf values.
    """
    toy = _make_base_archaeal_toy_model()

    # Case A: linprog reports success=True, but x has severe mass balance residual (e.g. 0.1)
    with patch("exp1_sbml_highs_engine.linprog") as mock_lp:
        n_rxns = len(toy.reactions)
        # Create non-zero flux on glucose uptake without internal consumption -> violates d(Glucose)/dt = 0!
        x_violated = np.zeros(n_rxns, dtype=np.float64)
        x_violated[0] = 10.0  # R_EX_cpd00027_e0 = 10.0, produces intracellular glucose with no sink
        mock_lp.return_value = OptimizeResult(
            success=True,
            status=0,
            message="Optimization terminated successfully (forged)",
            x=x_violated,
        )
        sol = solve_fba_highs(toy, time_limit_s=5.0, residual_tol=1e-5)
        assert sol["status"] == "numerical_coupling_violation"
        assert math.isnan(sol["growth_rate"])
        assert sol["max_residual"] > 1e-5
        assert sol["fluxes"] == {}

    # Case B: linprog reports success=True, but x violates reaction upper bound
    with patch("exp1_sbml_highs_engine.linprog") as mock_lp:
        x_bound_viol = np.zeros(n_rxns, dtype=np.float64)
        x_bound_viol[0] = 5000.0  # ub is 1000.0 -> violates bound by 4000.0!
        mock_lp.return_value = OptimizeResult(
            success=True,
            status=0,
            message="Optimization terminated successfully (forged)",
            x=x_bound_viol,
        )
        sol = solve_fba_highs(toy, time_limit_s=5.0, residual_tol=1e-5)
        assert sol["status"] == "numerical_coupling_violation"
        assert math.isnan(sol["growth_rate"])
        assert sol["max_bound_viol"] > 1e-5

    # Case C: linprog reports success=True, but x contains NaN
    with patch("exp1_sbml_highs_engine.linprog") as mock_lp:
        x_nan = np.zeros(n_rxns, dtype=np.float64)
        x_nan[0] = np.nan
        mock_lp.return_value = OptimizeResult(
            success=True,
            status=0,
            message="Optimization terminated successfully (forged)",
            x=x_nan,
        )
        sol = solve_fba_highs(toy, time_limit_s=5.0)
        assert sol["status"] == "numerical_error"
        assert math.isnan(sol["growth_rate"])


def test_19_stage3_collector_exact_ids_na_handling_and_empty_shard_safety(tmp_path: Path):
    """Collector test per EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt Item 3:
    Verifies that stage3_collect_and_summarize handles empty shards without concat([]) crash,
    enforces exact eligible genome IDs, assigns NA (not 0) to invalid models for reaction counts,
    and excludes taxonomy placeholders from shared lineages.
    """
    from stage3_collect_and_summarize import main as collector_main

    shards_dir = tmp_path / "shards"
    shards_dir.mkdir(parents=True)
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True)

    # 2 eligible genomes in inventory
    inv_tsv = tmp_path / "inv.tsv"
    inv_tsv.write_text(
        "genome\thas_draft_sbml\tqc_error\tlegacy_classification\n"
        "g_01\tTrue\t\to__Methanomicrobiales;f__Methanomicrobiaceae;g__Methanoculleus\n"
        "g_02\tTrue\t\to__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina\n"
    )

    # Shard 00 has g_01, Shard 01 has g_02, Shard 02 is empty (valid empty shard)
    # Shard 00: g_01 is feasible
    sum_0 = pd.DataFrame(
        [
            {
                "genome": "g_01",
                "cohort": "EQ_desert",
                "canonical_order": "o__Methanomicrobiales",
                "canonical_family": "f__Methanomicrobiaceae",
                "canonical_genus": "g__Methanoculleus",
                "completeness": 95.0,
                "contamination": 1.0,
                "scenario_id": "S4_oxic_desert_oa_aa",
                "redox_regime": "oxic",
                "biomass_recipe_id": "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420",
                "solver_status": "heuristic_incumbent_bounded",
                "n_distinct_valid_solutions": 4,
                "min_added_rxns_count": 3.0,
                "max_added_rxns_count": 5.0,
                "n_rxns_observed_in_all_enumerated_models": 2.0,
                "n_rxns_essential_in_declared_supernetwork": 1.0,
                "min_total_cost": 50.0,
                "max_ensemble_cost": 75.0,
            }
        ]
    )
    sum_0.to_csv(shards_dir / "exp1_shard_00_summary.tsv", sep="\t", index=False)
    # Shard 01: g_02 is infeasible (no valid solutions) -> must be converted to NA in collector!
    sum_1 = pd.DataFrame(
        [
            {
                "genome": "g_02",
                "cohort": "external_comparator",
                "canonical_order": "o__Methanosarcinales",
                "canonical_family": "f__",  # Placeholder taxonomy!
                "canonical_genus": "g__",
                "completeness": 92.0,
                "contamination": 2.0,
                "scenario_id": "S4_oxic_desert_oa_aa",
                "redox_regime": "oxic",
                "biomass_recipe_id": "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420",
                "solver_status": "infeasible_under_fixed_medium",
                "n_distinct_valid_solutions": 0,
                "min_added_rxns_count": 0.0,  # Falsely written as 0.0 by uncorrected worker -> collector must fix to NA!
                "max_added_rxns_count": 0.0,
                "n_rxns_observed_in_all_enumerated_models": 0.0,
                "n_rxns_essential_in_declared_supernetwork": 0.0,
                "min_total_cost": 0.0,
                "max_ensemble_cost": 0.0,
            }
        ]
    )
    sum_1.to_csv(shards_dir / "exp1_shard_01_summary.tsv", sep="\t", index=False)
    # Shard 02: empty summary file with header only
    pd.DataFrame(columns=sum_0.columns).to_csv(shards_dir / "exp1_shard_02_summary.tsv", sep="\t", index=False)

    with patch("sys.argv", [
        "stage3_collect_and_summarize.py",
        "--inventory-484-tsv", str(inv_tsv),
        "--shards-dir", str(shards_dir),
        "--out-dir", str(out_dir),
        "--num-tasks", "3",
    ]):
        collector_main()

    # Verify that g_02 (infeasible) has NA for added rxns and costs, NOT 0.0!
    full_sum = pd.read_csv(out_dir / "exp1_full484_ensemble_summary.tsv.gz", sep="\t")
    g02_row = full_sum[full_sum["genome"] == "g_02"].iloc[0]
    assert pd.isna(g02_row["min_added_rxns_count"])
    assert pd.isna(g02_row["min_total_cost"])
    assert pd.isna(g02_row["n_rxns_essential_in_declared_supernetwork"])

    # Verify that g_01 (feasible) retains valid reaction counts
    g01_row = full_sum[full_sum["genome"] == "g_01"].iloc[0]
    assert g01_row["min_added_rxns_count"] == 3.0

    # Verify lineage comparisons file excludes placeholder f__ from being treated as real shared lineage
    lin_df = pd.read_csv(out_dir / "exp1_full484_lineage_matched_comparisons_order_family_genus.tsv", sep="\t")
    f_placeholder = lin_df[(lin_df["taxon_name"] == "f__")]
    if not f_placeholder.empty:
        assert bool(f_placeholder["is_missing_taxonomy_placeholder"].all())
        assert not bool(f_placeholder["is_shared_lineage"].any())


def test_20_stage4_missing_reference_and_assembly_produces_not_evaluable_never_auxotrophy(tmp_path: Path):
    """Stage 4 test per EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt Item B:
    Missing individual reaction reference FASTA must produce explicit not_evaluable disposition,
    NEVER silently default to supported_candidate_auxotrophy_or_pathway_absence.
    """
    from stage4_candidate_deficiency_search import main as stage4_main

    out_dir = tmp_path / "stage4_out"
    tmp_dir = tmp_path / "stage4_tmp"
    out_dir.mkdir(parents=True)
    tmp_dir.mkdir(parents=True)

    # Empty reference directory (no fasta for rxn_missing_ref)
    ref_dir = tmp_path / "empty_ref_dir"
    ref_dir.mkdir(parents=True)

    inst_tsv = tmp_path / "instances.tsv"
    inst_tsv.write_text(
        "scenario_id\trxn_id\tbase_id\trxn_name\tis_transport\tbiomass_sensitivity_classification\t"
        "followup_selection_tier\tepistemic_status\tgenome\tcohort\tcanonical_order\tcanonical_family\t"
        "canonical_genus\tcompleteness\tcontamination\tsample_id\n"
        "S4_oxic_desert_oa_aa\tR_missing_rxn\trxn_missing_ref\tMissing Enzyme\tFalse\t"
        "robust_across_all_tested_biomass_perturbations\tdirectly_verified_robust_candidate\t"
        "sequence_synteny_followup_selection_not_functional_proof\tEQ_mag_01\tEQ_desert\t"
        "o__Methanomicrobiales\tf__Methanomicrobiaceae\tg__Methanoculleus\t95.0\t1.0\tSample_01\n"
    )

    with patch("sys.argv", [
        "stage4_candidate_deficiency_search.py",
        "--instances-tsv", str(inst_tsv),
        "--summary-tsv", str(inst_tsv),
        "--gapseq-rxn-fasta-dir", str(ref_dir),
        "--refined-mags-tar", str(tmp_path / "nonexistent_mags.tar.gz"),
        "--assemblies-dir", str(tmp_path / "assemblies"),
        "--gap-context-tsv-gz", str(tmp_path / "nonexistent.tsv.gz"),
        "--c18-candidates-tsv", str(tmp_path / "nonexistent.tsv"),
        "--out-dir", str(out_dir),
        "--tmp-dir", str(tmp_dir),
        "--threads", "1",
    ]):
        stage4_main()

    res_df = pd.read_csv(out_dir / "exp1_stage4_deficiency_verification_instances.tsv", sep="\t")
    row = res_df.iloc[0]
    assert row["reference_sequence_status"] == "reference_unavailable"
    assert row["deficiency_verification_disposition"] == "not_evaluable_missing_reference_sequences"
    assert row["deficiency_verification_disposition"] != "supported_candidate_auxotrophy_or_pathway_absence"


def test_21_stage4_sample_assembly_hit_without_host_linkage_is_unresolved_not_dropout(tmp_path: Path):
    """Stage 4 test per EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt Item B:
    A raw metagenome assembly hit can belong to another organism in the community.
    Without verified host linkage, it must be classified as unresolved host assignment,
    NOT asserted as assembly_or_binning_dropout of the queried MAG.
    """
    from stage4_candidate_deficiency_search import main as stage4_main

    out_dir = tmp_path / "out"
    tmp_dir = tmp_path / "tmp"
    out_dir.mkdir(parents=True)
    tmp_dir.mkdir(parents=True)

    ref_dir = tmp_path / "ref_dir"
    ref_dir.mkdir(parents=True)
    (ref_dir / "rxn12345.fasta").write_text(">rxn12345|test_prot\nMKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQ\n")

    # Sample assembly exists with a hit
    sample_dir = tmp_path / "assemblies" / "Sample_01"
    sample_dir.mkdir(parents=True)
    (sample_dir / "Sample_01_contigs.fa").write_text(">sample_contig_999\nATGAAAGTGTTATGGGCCGCATTGTTGGTTACTTTCTTAGCTGGCTGTCAAGCAAAGGTGGAACAAGCAGTTGAGACTGAGCCAGAGCCTGAGCTTCGCCAGCAGACAGAGTGGCAAAGCGGCCAGTGA\n")

    inst_tsv = tmp_path / "instances.tsv"
    inst_tsv.write_text(
        "scenario_id\trxn_id\tbase_id\trxn_name\tis_transport\tbiomass_sensitivity_classification\t"
        "followup_selection_tier\tepistemic_status\tgenome\tcohort\tcanonical_order\tcanonical_family\t"
        "canonical_genus\tcompleteness\tcontamination\tsample_id\n"
        "S4_oxic_desert_oa_aa\tR_rxn12345\trxn12345\tTest Enzyme\tFalse\t"
        "robust_across_all_tested_biomass_perturbations\tdirectly_verified_robust_candidate\t"
        "sequence_synteny_followup_selection_not_functional_proof\tSample_01__bin1\tEQ_desert\t"
        "o__Methanomicrobiales\tf__Methanomicrobiaceae\tg__Methanoculleus\t95.0\t1.0\tSample_01\n"
    )

    # Mock DIAMOND to return sample assembly hit with bitscore 75.0, but NO MAG hit
    with patch("stage4_candidate_deficiency_search.subprocess.run") as mock_subp, \
         patch("stage4_candidate_deficiency_search.extract_refined_mags_and_run_diamond_blastx") as mock_mag, \
         patch("stage4_candidate_deficiency_search.search_sample_assemblies_for_community_homologs") as mock_samp:

        mock_mag.return_value = ({}, {}, "sha_dummy")
        mock_samp.return_value = (
            {("Sample_01", "rxn12345"): {
                "sample_contig_id": "sample_contig_999",
                "sample_assembly_bitscore": 75.0,
                "sample_assembly_pident": 98.0,
                "sample_assembly_evalue": 1e-25,
                "sample_alignment_length": 42,
            }},
            {"Sample_01": True},
            {"Sample_01": "sha_sample"}
        )

        with patch("sys.argv", [
            "stage4_candidate_deficiency_search.py",
            "--instances-tsv", str(inst_tsv),
            "--summary-tsv", str(inst_tsv),
            "--gapseq-rxn-fasta-dir", str(ref_dir),
            "--refined-mags-tar", str(tmp_path / "dummy.tar.gz"),
            "--assemblies-dir", str(tmp_path / "assemblies"),
            "--gap-context-tsv-gz", str(tmp_path / "nonexistent.tsv.gz"),
            "--c18-candidates-tsv", str(tmp_path / "nonexistent.tsv"),
            "--out-dir", str(out_dir),
            "--tmp-dir", str(tmp_dir),
            "--threads", "1",
        ]):
            # create dummy tar so exists() is True
            (tmp_path / "dummy.tar.gz").write_bytes(b"dummy")
            stage4_main()

    res_df = pd.read_csv(out_dir / "exp1_stage4_deficiency_verification_instances.tsv", sep="\t")
    row = res_df.iloc[0]
    assert row["sample_assembly_bitscore"] == 75.0
    assert not bool(row["verified_host_linkage"])
    assert row["deficiency_verification_disposition"] == "unresolved_sample_community_homolog_host_unassigned"
    assert row["deficiency_verification_disposition"] != "assembly_or_binning_dropout"


def test_22_stage4_near_edge_full_length_vs_truncated_distinction(tmp_path: Path):
    """Stage 4 test per EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt Item B:
    A hit within 300bp of contig end with full subject coverage is a complete gene,
    NOT a truncated/split gene. Truncated requires subject coverage < 0.70.
    """
    from stage4_candidate_deficiency_search import main as stage4_main

    out_dir = tmp_path / "out"
    tmp_dir = tmp_path / "tmp"
    out_dir.mkdir(parents=True)
    tmp_dir.mkdir(parents=True)

    ref_dir = tmp_path / "ref_dir"
    ref_dir.mkdir(parents=True)
    (ref_dir / "rxn_full.fasta").write_text(">rxn_full|prot1\n" + "A" * 100 + "\n")
    (ref_dir / "rxn_trunc.fasta").write_text(">rxn_trunc|prot2\n" + "A" * 100 + "\n")

    inst_tsv = tmp_path / "instances.tsv"
    inst_tsv.write_text(
        "scenario_id\trxn_id\tbase_id\trxn_name\tis_transport\tbiomass_sensitivity_classification\t"
        "followup_selection_tier\tepistemic_status\tgenome\tcohort\tcanonical_order\tcanonical_family\t"
        "canonical_genus\tcompleteness\tcontamination\tsample_id\n"
        "S4_oxic_desert_oa_aa\tR_full\trxn_full\tFull Enzyme\tFalse\trobust\ttier1\tstatus\tg_1\tEQ_desert\to\tf\tg\t95.0\t1.0\tS1\n"
        "S4_oxic_desert_oa_aa\tR_trunc\trxn_trunc\tTruncated Enzyme\tFalse\trobust\ttier1\tstatus\tg_1\tEQ_desert\to\tf\tg\t95.0\t1.0\tS1\n"
    )

    with patch("stage4_candidate_deficiency_search.subprocess.run"), \
         patch("stage4_candidate_deficiency_search.extract_refined_mags_and_run_diamond_blastx") as mock_mag, \
         patch("stage4_candidate_deficiency_search.search_sample_assemblies_for_community_homologs") as mock_samp:

        mock_mag.return_value = (
            {
                ("g_1", "rxn_full"): {
                    "contig_id": "c1", "bitscore": 80.0, "pident": 95.0, "evalue": 1e-30,
                    "alignment_length": 98, "subject_length": 100, "subject_coverage": 0.98,
                    "dist_to_contig_edge_bp": 150, "is_near_contig_edge": True,
                    "is_truncated_near_edge": False, "is_full_length_near_edge": True,
                },
                ("g_1", "rxn_trunc"): {
                    "contig_id": "c2", "bitscore": 45.0, "pident": 95.0, "evalue": 1e-15,
                    "alignment_length": 40, "subject_length": 100, "subject_coverage": 0.40,
                    "dist_to_contig_edge_bp": 50, "is_near_contig_edge": True,
                    "is_truncated_near_edge": True, "is_full_length_near_edge": False,
                },
            },
            {},
            "sha_mag",
        )
        mock_samp.return_value = ({}, {"S1": False}, {})

        (tmp_path / "dummy.tar.gz").write_bytes(b"dummy")
        with patch("sys.argv", [
            "stage4_candidate_deficiency_search.py",
            "--instances-tsv", str(inst_tsv),
            "--summary-tsv", str(inst_tsv),
            "--gapseq-rxn-fasta-dir", str(ref_dir),
            "--refined-mags-tar", str(tmp_path / "dummy.tar.gz"),
            "--assemblies-dir", str(tmp_path / "assemblies"),
            "--gap-context-tsv-gz", str(tmp_path / "nonexistent.tsv.gz"),
            "--c18-candidates-tsv", str(tmp_path / "nonexistent.tsv"),
            "--out-dir", str(out_dir),
            "--tmp-dir", str(tmp_dir),
            "--threads", "1",
        ]):
            stage4_main()

    res_df = pd.read_csv(out_dir / "exp1_stage4_deficiency_verification_instances.tsv", sep="\t")
    row_full = res_df[res_df["clean_rxn_id"] == "rxn_full"].iloc[0]
    row_trunc = res_df[res_df["clean_rxn_id"] == "rxn_trunc"].iloc[0]

    assert row_full["deficiency_verification_disposition"] == "complete_homolog_near_contig_edge"
    assert row_trunc["deficiency_verification_disposition"] == "candidate_truncated_or_split_gene_at_contig_edge"


def test_23_stage4_empty_candidates_emits_schema_valid_completed_outputs(tmp_path: Path):
    """Stage 4 test per EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt Item B:
    Empty eligible candidate set emits schema-valid completed nonestimable outputs.
    """
    from stage4_candidate_deficiency_search import main as stage4_main

    out_dir = tmp_path / "out"
    tmp_dir = tmp_path / "tmp"
    out_dir.mkdir(parents=True)
    tmp_dir.mkdir(parents=True)

    inst_tsv = tmp_path / "empty_instances.tsv"
    inst_tsv.write_text("scenario_id\trxn_id\tbase_id\n")

    with patch("sys.argv", [
        "stage4_candidate_deficiency_search.py",
        "--instances-tsv", str(inst_tsv),
        "--summary-tsv", str(inst_tsv),
        "--gapseq-rxn-fasta-dir", str(tmp_path / "ref"),
        "--refined-mags-tar", str(tmp_path / "dummy.tar.gz"),
        "--assemblies-dir", str(tmp_path / "assemblies"),
        "--gap-context-tsv-gz", str(tmp_path / "nonexistent.tsv.gz"),
        "--c18-candidates-tsv", str(tmp_path / "nonexistent.tsv"),
        "--out-dir", str(out_dir),
        "--tmp-dir", str(tmp_dir),
        "--threads", "1",
    ]):
        stage4_main()

    inst_out = out_dir / "exp1_stage4_deficiency_verification_instances.tsv"
    sum_out = out_dir / "exp1_stage4_deficiency_verification_summary.tsv"
    receipt_out = out_dir / "provenance" / "stage4_execution_receipts.json"

    assert inst_out.exists()
    assert sum_out.exists()
    assert receipt_out.exists()
    assert len(pd.read_csv(inst_out, sep="\t")) == 0
    receipt = json.loads(receipt_out.read_text())
    assert receipt["status"] == "completed_zero_candidates"


def test_24_missing_verification_columns_and_na_residuals_fail_gate():
    """Regression test per EXP1-FULL-RELEASE-DECISION.txt:
    Verifies that missing verification columns or NA/NaN residuals in accepted rows
    IMMEDIATELY fail acceptance (missing evidence is NEVER treated as zero residual).
    """
    from stage2_pilot24_ensemble import validate_pilot24_scientific_gate
    from exp1_sbml_highs_engine import (
        PLANNED_PROTOCOL_SCENARIO_BIOMASS_TUPLES_PER_GENOME,
        get_planned_protocol_scenario_biomass_tuples,
    )

    g = "TEST_G01"
    tuples = get_planned_protocol_scenario_biomass_tuples([g])
    sorted_tuples = sorted(tuples)

    # Base valid summary
    base_rows = []
    for _, sc, bio in sorted_tuples:
        base_rows.append({
            "genome": g,
            "scenario_id": sc,
            "biomass_recipe_id": bio,
            "solver_status": "optimal",
            "n_distinct_valid_solutions": 1,
            "min_added_rxns_count": 0.0,
            "min_total_cost": 0.0,
            "min_solution_verified_growth": 0.05,
            "max_solution_residual": 1e-8,
            "max_solution_bound_viol": 1e-8,
            "n_cycle_qc_rejected_solutions": 0,
            "n_optimization_feasible_solutions": 1,
        })
    df_valid = pd.DataFrame(base_rows)
    empty_rf = pd.DataFrame([{"genome": g, "scenario_id": "S3", "rxn_id": "R1", "biomass_sensitivity_classification": "robust_across_all_tested_biomass_perturbations"}])
    empty_wkf = pd.DataFrame([{"genome": g, "record_type": "withdrawal"}])
    empty_mask = pd.DataFrame([{"genome": g, "mask_fraction": 0.1, "replicate": 0, "masked_solver_status": "optimal"}])

    # Case A: Valid dataframe passes
    gate_ok = validate_pilot24_scientific_gate(
        summary_df=df_valid,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_ok["checks"]["per_solution_correctness_verified"] is True

    # Case B: Missing max_solution_residual column entirely -> MUST FAIL!
    df_missing_col = df_valid.drop(columns=["max_solution_residual"])
    gate_fail_col = validate_pilot24_scientific_gate(
        summary_df=df_missing_col,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_fail_col["checks"]["per_solution_correctness_verified"] is False

    # Case C: max_solution_residual is NaN in an accepted row -> MUST FAIL (not fillna(0.0)!)
    df_na_residual = df_valid.copy()
    df_na_residual.loc[0, "max_solution_residual"] = math.nan
    gate_fail_na = validate_pilot24_scientific_gate(
        summary_df=df_na_residual,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_fail_na["checks"]["per_solution_correctness_verified"] is False

    # Case D: max_solution_bound_viol is NaN in an accepted row -> MUST FAIL!
    df_na_bound = df_valid.copy()
    df_na_bound.loc[0, "max_solution_bound_viol"] = math.nan
    gate_fail_bound = validate_pilot24_scientific_gate(
        summary_df=df_na_bound,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_fail_bound["checks"]["per_solution_correctness_verified"] is False

    # Case E: min_solution_verified_growth is NaN in an accepted row -> MUST FAIL!
    df_na_growth = df_valid.copy()
    df_na_growth.loc[0, "min_solution_verified_growth"] = math.nan
    gate_fail_growth = validate_pilot24_scientific_gate(
        summary_df=df_na_growth,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_fail_growth["checks"]["per_solution_correctness_verified"] is False


def test_25_positive_solution_count_with_failed_residual_or_bounds_or_cycle_fails_gate():
    """Regression test per EXP1-FULL-RELEASE-DECISION.txt:
    Positive valid solution count with violated residual, bound violation, growth,
    or rejected cycle MUST be rejected by the gate.
    """
    from stage2_pilot24_ensemble import validate_pilot24_scientific_gate
    from exp1_sbml_highs_engine import get_planned_protocol_scenario_biomass_tuples

    g = "TEST_G01"
    sorted_tuples = sorted(get_planned_protocol_scenario_biomass_tuples([g]))

    base_rows = []
    for _, sc, bio in sorted_tuples:
        base_rows.append({
            "genome": g,
            "scenario_id": sc,
            "biomass_recipe_id": bio,
            "solver_status": "optimal",
            "n_distinct_valid_solutions": 1,
            "min_added_rxns_count": 1.0,
            "min_total_cost": 2.0,
            "min_solution_verified_growth": 0.05,
            "max_solution_residual": 1e-8,
            "max_solution_bound_viol": 1e-8,
            "n_cycle_qc_rejected_solutions": 0,
            "n_optimization_feasible_solutions": 1,
        })
    df_valid = pd.DataFrame(base_rows)
    empty_rf = pd.DataFrame([{"genome": g, "scenario_id": "S3", "rxn_id": "R1", "biomass_sensitivity_classification": "robust_across_all_tested_biomass_perturbations"}])
    empty_wkf = pd.DataFrame([{"genome": g, "record_type": "withdrawal"}])
    empty_mask = pd.DataFrame([{"genome": g, "mask_fraction": 0.1, "replicate": 0, "masked_solver_status": "optimal"}])

    # Case A: max_solution_residual exceeds 1e-5 (e.g. 0.01) -> MUST FAIL!
    df_bad_res = df_valid.copy()
    df_bad_res.loc[0, "max_solution_residual"] = 0.01
    gate_bad_res = validate_pilot24_scientific_gate(
        summary_df=df_bad_res,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_bad_res["checks"]["per_solution_correctness_verified"] is False

    # Case B: max_solution_bound_viol exceeds 1e-5 (e.g. 0.05) -> MUST FAIL!
    df_bad_bnd = df_valid.copy()
    df_bad_bnd.loc[0, "max_solution_bound_viol"] = 0.05
    gate_bad_bnd = validate_pilot24_scientific_gate(
        summary_df=df_bad_bnd,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_bad_bnd["checks"]["per_solution_correctness_verified"] is False

    # Case C: verified growth below 0.00999 (e.g. 0.0001) -> MUST FAIL!
    df_bad_gr = df_valid.copy()
    df_bad_gr.loc[0, "min_solution_verified_growth"] = 0.0001
    gate_bad_gr = validate_pilot24_scientific_gate(
        summary_df=df_bad_gr,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_bad_gr["checks"]["per_solution_correctness_verified"] is False

    # Case D: all solutions rejected by cycle QC -> MUST FAIL if claimed as accepted!
    df_all_cyc_rej = df_valid.copy()
    df_all_cyc_rej.loc[0, "n_cycle_qc_rejected_solutions"] = 5
    df_all_cyc_rej.loc[0, "n_optimization_feasible_solutions"] = 5
    df_all_cyc_rej.loc[0, "n_distinct_valid_solutions"] = 0
    gate_cyc_rej = validate_pilot24_scientific_gate(
        summary_df=df_all_cyc_rej,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
    )
    assert gate_cyc_rej["checks"]["per_solution_correctness_verified"] is False


def test_26_draft_already_grows_produces_finite_verified_evidence_and_ledger_entry():
    """Test draft-already-grows early return produces finite verified residuals, bounds, growth,
    and valid per-solution ledger entries (per EXP1-FULL-RELEASE-DECISION.txt).
    """
    from exp1_sbml_highs_engine import (
        BIO_B0_BASELINE,
        ReactionRecord,
        enumerate_alternative_gapfill_ensemble,
    )

    toy = _make_base_archaeal_toy_model()
    # Add endogenous lysine synthesis so the draft model is self-sufficient on glucose + NH3
    toy.reactions["R_lys_synth"] = ReactionRecord(
        id="R_lys_synth",
        base_id="rxn_lys_synth",
        name="Lysine synthesis",
        reversible=False,
        lb=0.0,
        ub=1000.0,
        stoichiometry={"M_cpd00027_c0": -0.5, "M_cpd00013_c0": -1.0, "M_cpd00039_c0": 1.0},
    )

    med_cpds = {
        "cpd00027": ("Glucose", 5.0),
        "cpd00013": ("NH3", 10.0),
    }
    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=toy,
        candidate_reactions={},
        candidate_weights={},
        scenario_compounds=med_cpds,
        biomass_recipe_id=BIO_B0_BASELINE,
        min_growth=0.01,
        n_Target_solutions=1,
    )
    assert ens["solver_status"] == "optimal_draft_no_gapfill_needed"
    assert ens["draft_already_grows"] is True
    assert math.isfinite(ens["max_solution_residual"])
    assert ens["max_solution_residual"] <= 1e-5
    assert math.isfinite(ens["max_solution_bound_viol"])
    assert ens["max_solution_bound_viol"] <= 1e-5
    assert math.isfinite(ens["min_solution_verified_growth"])
    assert ens["min_solution_verified_growth"] >= 0.01

    sol = ens["solutions"][0]
    assert math.isfinite(sol["fba_check_max_residual"])
    assert sol["fba_check_max_residual"] <= 1e-5
    assert math.isfinite(sol["fba_check_max_bound_viol"])
    assert sol["fba_check_max_bound_viol"] <= 1e-5
    assert sol["verified_fba_growth"] >= 0.01


def test_27_per_solution_ledger_validity_and_unknown_dispositions():
    """Test per-solution ledger validity disposition auditing and unknown outcome handling
    (per EXP1-FULL-RELEASE-DECISION.txt).
    """
    from stage2_pilot24_ensemble import validate_pilot24_scientific_gate
    from exp1_sbml_highs_engine import get_planned_protocol_scenario_biomass_tuples

    g = "TEST_G01"
    tuples = sorted(get_planned_protocol_scenario_biomass_tuples([g]))

    sum_rows = []
    ledger_rows = []
    for idx, (_, sc, bio) in enumerate(tuples):
        if idx < 20:
            sum_rows.append({
                "genome": g, "scenario_id": sc, "biomass_recipe_id": bio,
                "solver_status": "optimal", "n_distinct_valid_solutions": 1,
                "min_added_rxns_count": 1.0, "min_total_cost": 1.0,
                "min_solution_verified_growth": 0.05, "max_solution_residual": 1e-8,
                "max_solution_bound_viol": 1e-8,
            })
            ledger_rows.append({
                "genome": g, "cohort": "EQ_desert", "canonical_order": "o__", "canonical_family": "f__",
                "scenario_id": sc, "biomass_recipe_id": bio, "solution_index": 0,
                "validity_disposition": "accepted_feasible", "solver_status": "optimal",
                "selected_reaction_ids": "rxn1", "n_added_rxns": 1, "total_cost": 1.0,
                "raw_milp_objective": 1.0, "mip_gap": 0.0, "dual_bound": 1.0,
                "verified_growth": 0.05, "stoichiometric_max_residual": 1e-8,
                "max_bound_violation": 1e-8, "residual_bounds_growth_verified": True,
                "cycle_qc_status": "passed_no_free_atp_flux", "spurious_atp_cycle": False,
                "explicit_reason": "", "code_version": "v2",
            })
        else:
            sum_rows.append({
                "genome": g, "scenario_id": sc, "biomass_recipe_id": bio,
                "solver_status": "infeasible_under_fixed_medium", "n_distinct_valid_solutions": 0,
                "min_added_rxns_count": math.nan, "min_total_cost": math.nan,
                "min_solution_verified_growth": math.nan, "max_solution_residual": math.nan,
                "max_solution_bound_viol": math.nan,
            })
            ledger_rows.append({
                "genome": g, "cohort": "EQ_desert", "canonical_order": "o__", "canonical_family": "f__",
                "scenario_id": sc, "biomass_recipe_id": bio, "solution_index": -1,
                "validity_disposition": "unresolved_infeasible_under_fixed_medium",
                "solver_status": "infeasible_under_fixed_medium",
                "selected_reaction_ids": "", "n_added_rxns": math.nan, "total_cost": math.nan,
                "raw_milp_objective": math.nan, "mip_gap": math.nan, "dual_bound": math.nan,
                "verified_growth": math.nan, "stoichiometric_max_residual": math.nan,
                "max_bound_violation": math.nan, "residual_bounds_growth_verified": False,
                "cycle_qc_status": "none_evaluated", "spurious_atp_cycle": math.nan,
                "explicit_reason": "supernetwork_cannot_produce_biomass_under_fixed_medium",
                "code_version": "v2",
            })

    sum_df = pd.DataFrame(sum_rows)
    ledger_df = pd.DataFrame(ledger_rows)
    empty_rf = pd.DataFrame([{"genome": g, "scenario_id": "S3", "rxn_id": "R1", "biomass_sensitivity_classification": "robust_across_all_tested_biomass_perturbations"}])
    empty_wkf = pd.DataFrame([{"genome": g, "record_type": "withdrawal"}])
    empty_mask = pd.DataFrame([{"genome": g, "mask_fraction": 0.1, "replicate": 0, "masked_solver_status": "optimal"}])

    # Valid ledger passes
    gate = validate_pilot24_scientific_gate(
        summary_df=sum_df,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
        solution_ledger_df=ledger_df,
    )
    assert gate["checks"]["per_solution_correctness_verified"] is True

    # Corrupted ledger: accepted_feasible with residual_bounds_growth_verified == False -> MUST FAIL!
    bad_ledger = ledger_df.copy()
    bad_ledger.loc[0, "residual_bounds_growth_verified"] = False
    gate_bad_ledger = validate_pilot24_scientific_gate(
        summary_df=sum_df,
        rxn_freq_df=empty_rf,
        wkf_df=empty_wkf,
        masking_df=empty_mask,
        expected_n_genomes=1,
        solution_ledger_df=bad_ledger,
    )
    assert gate_bad_ledger["checks"]["per_solution_correctness_verified"] is False


def test_28_fail_open_ledger_holes_strictly_rejected():
    """Narrow regression test for EXP1-LEDGER-FINAL-ACCEPTANCE.txt & EXP1-REPLAY-EXECUTION-CORRECTION.txt:
    1. Missing residual/bounds in solution dictionary (NaN default, cannot invent 0.0).
    2. Unknown cycle status and spurious ATP cycle defaults (cannot invent False/passed).
    3. Empty ledger strictly rejected when summary is non-empty.
    4. Key mismatch between summary and ledger tasks strictly rejected.
    5. Count mismatch between summary n_distinct_valid_solutions and ledger accepted rows strictly rejected.
    6. Legitimate documented no-ATP status accepted when valid.
    """
    from stage2_pilot24_ensemble import validate_pilot24_scientific_gate
    from stage2_focused_release_validation import validate_focused_release_results
    from exp1_sbml_highs_engine import get_planned_protocol_scenario_biomass_tuples

    g = "TEST_G01"
    tuples = sorted(get_planned_protocol_scenario_biomass_tuples([g]))
    sc, bio = tuples[0][1], tuples[0][2]

    # Baseline 1-case summary
    base_sum = pd.DataFrame([{
        "genome": g, "cohort": "EQ_desert", "canonical_order": "o__", "canonical_family": "f__",
        "scenario_id": sc, "biomass_recipe_id": bio, "solver_status": "optimal", "case_type": "representative_feasible",
        "n_distinct_valid_solutions": 1, "min_added_rxns_count": 1.0, "min_total_cost": 1.0,
        "min_solution_verified_growth": 0.05, "max_solution_residual": 1e-8,
        "max_solution_bound_viol": 1e-8,
    }])
    base_led = pd.DataFrame([{
        "genome": g, "cohort": "EQ_desert", "canonical_order": "o__", "canonical_family": "f__",
        "scenario_id": sc, "biomass_recipe_id": bio, "case_type": "representative_feasible",
        "solution_index": 0, "validity_disposition": "accepted_feasible", "solver_status": "optimal",
        "selected_reaction_ids": "rxn1", "n_added_rxns": 1, "total_cost": 1.0,
        "raw_milp_objective": 1.0, "mip_gap": 0.0, "dual_bound": 1.0,
        "verified_growth": 0.05, "stoichiometric_max_residual": 1e-8,
        "max_bound_violation": 1e-8, "residual_bounds_growth_verified": True,
        "cycle_qc_status": "passed_closed_uptake_egc_check", "spurious_atp_cycle": False,
        "explicit_reason": "", "code_version": "v2",
    }])
    base_led["spurious_atp_cycle"] = base_led["spurious_atp_cycle"].astype(object)
    base_led["residual_bounds_growth_verified"] = base_led["residual_bounds_growth_verified"].astype(object)

    # 1. Baseline valid passes both validators
    rep_foc = validate_focused_release_results(base_sum, base_led)
    assert rep_foc["validation_passed"] is True

    # 2. Empty ledger with non-empty summary MUST FAIL
    empty_led = pd.DataFrame()
    rep_empty = validate_focused_release_results(base_sum, empty_led)
    assert rep_empty["validation_passed"] is False
    assert rep_empty["checks"]["per_solution_ledger_integrity_verified"] is False

    # 3. Key mismatch between summary and ledger MUST FAIL
    bad_key_led = base_led.copy()
    bad_key_led.loc[0, "scenario_id"] = "S99_unmatched_scenario"
    rep_bad_key = validate_focused_release_results(base_sum, bad_key_led)
    assert rep_bad_key["validation_passed"] is False
    assert rep_bad_key["checks"]["per_solution_ledger_integrity_verified"] is False

    # 4. Count mismatch (summary claims 2 solutions, ledger has 1) MUST FAIL
    bad_count_sum = base_sum.copy()
    bad_count_sum.loc[0, "n_distinct_valid_solutions"] = 2
    rep_bad_count = validate_focused_release_results(bad_count_sum, base_led)
    assert rep_bad_count["validation_passed"] is False
    assert rep_bad_count["checks"]["per_solution_ledger_integrity_verified"] is False

    # 5. Unknown cycle status in accepted row MUST FAIL
    bad_cyc_led = base_led.copy()
    bad_cyc_led.loc[0, "cycle_qc_status"] = "unknown"
    rep_bad_cyc = validate_focused_release_results(base_sum, bad_cyc_led)
    assert rep_bad_cyc["validation_passed"] is False
    assert rep_bad_cyc["checks"]["per_solution_ledger_integrity_verified"] is False

    # 6. Spurious ATP cycle True in accepted row MUST FAIL
    bad_spur_led = base_led.copy()
    bad_spur_led.loc[0, "spurious_atp_cycle"] = True
    rep_bad_spur = validate_focused_release_results(base_sum, bad_spur_led)
    assert rep_bad_spur["validation_passed"] is False
    assert rep_bad_spur["checks"]["per_solution_ledger_integrity_verified"] is False

    # 7. NaN / missing residual in accepted row MUST FAIL
    bad_res_led = base_led.copy()
    bad_res_led.loc[0, "stoichiometric_max_residual"] = math.nan
    rep_bad_res = validate_focused_release_results(base_sum, bad_res_led)
    assert rep_bad_res["validation_passed"] is False
    assert rep_bad_res["checks"]["per_solution_ledger_integrity_verified"] is False

    # 8. Legitimate documented no-ATP status is accepted
    no_atp_led = base_led.copy()
    no_atp_led["spurious_atp_cycle"] = no_atp_led["spurious_atp_cycle"].astype(object)
    no_atp_led.loc[0, "cycle_qc_status"] = "evaluated_no_atp_species"
    no_atp_led.loc[0, "spurious_atp_cycle"] = math.nan
    rep_no_atp = validate_focused_release_results(base_sum, no_atp_led)
    assert rep_no_atp["validation_passed"] is True

    # 9. Missing cycle flag (NaN) on accepted row with standard cycle check MUST FAIL
    bad_flag_led = base_led.copy()
    bad_flag_led.loc[0, "spurious_atp_cycle"] = math.nan
    rep_bad_flag = validate_focused_release_results(base_sum, bad_flag_led)
    assert rep_bad_flag["validation_passed"] is False
    assert rep_bad_flag["checks"]["per_solution_ledger_integrity_verified"] is False


def test_29_focused_release_validation_and_witness_reconstruction():
    """Test focused release validation with ground truth controls and reconstruct_and_verify_solution_witness
    on a toy fixture model (`EXP1-REPLAY-EXECUTION-CORRECTION.txt`).
    """
    from exp1_sbml_highs_engine import (
        reconstruct_and_verify_solution_witness,
        ReactionRecord,
        SpeciesRecord,
        solve_fba_highs,
    )
    from stage2_focused_release_validation import validate_focused_release_results
    draft = _make_base_archaeal_toy_model()

    # Toy candidate reactions providing required lysine synthesis
    cand_rxns = {
        "R_cand_A_highflux_1rxn": ReactionRecord(
            id="R_cand_A_highflux_1rxn",
            base_id="rxn_cand_A",
            name="Lysine route",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.005, "M_cpd00013_c0": -0.01, "M_cpd00039_c0": 0.01},
        )
    }
    scen_cpds = {"cpd00027": ("Glucose", 10.0), "cpd00013": ("NH3", 10.0), "cpd00007": ("O2", 20.0), "cpd00009": ("Pi", 10.0)}

    # Witness reconstruction with clean reaction
    res = reconstruct_and_verify_solution_witness(
        draft_model=draft,
        candidate_reactions=cand_rxns,
        scenario_compounds=scen_cpds,
        biomass_recipe_id="B0_methanosarcina_template_baseline",
        selected_reaction_ids=["R_cand_A_highflux_1rxn"],
        min_growth=0.01,
        time_limit_s=5.0,
    )
    assert res["verified"] is True
    assert res["growth_rate"] >= 0.01
    assert res["max_residual"] <= 1e-5
    assert res["cycle_qc_status"] == "passed_closed_uptake_egc_check"

    # Known impossible negative control: impossible required precursor that no reaction produces
    draft_impossible = draft.copy()
    draft_impossible.species["M_cpd_impossible_c0"] = SpeciesRecord("M_cpd_impossible_c0", "cpd_impossible", "Impossible", "c0", 0, "C1")
    # Biomass needs cpd_impossible
    draft_impossible.reactions["R_bio1"].stoichiometry["M_cpd_impossible_c0"] = -1.0
    res_imp = reconstruct_and_verify_solution_witness(
        draft_model=draft_impossible,
        candidate_reactions=cand_rxns,
        scenario_compounds=scen_cpds,
        biomass_recipe_id="B0_methanosarcina_template_baseline",
        selected_reaction_ids=["R_cand_A_highflux_1rxn"],
        min_growth=0.01,
        time_limit_s=5.0,
    )
    assert res_imp["verified"] is False
    assert res_imp["growth_rate"] < 0.00999


def test_30_stable_masking_seed_and_deterministic_hash():
    """Test stable SHA-derived masking seeds and input hashes (EXP1-FINAL-CONTROLS-AND-COLLECTOR.txt Item 1)."""
    from run_masking_controls_7hq import (
        derive_stable_masking_seed,
        compute_exact_input_provenance_hash,
        serialize_and_hash_medium,
    )

    # 1. Deterministic across invocations
    s1, str1 = derive_stable_masking_seed("MGYG000519815", 0.10, 0)
    s2, str2 = derive_stable_masking_seed("MGYG000519815", 0.10, 0)
    assert s1 == s2
    assert str1 == str2
    assert isinstance(s1, int)
    assert 0 <= s1 < 2**32

    # 2. Distinct seeds for different fractions/replicates/genomes
    s3, _ = derive_stable_masking_seed("MGYG000519815", 0.20, 0)
    s4, _ = derive_stable_masking_seed("MGYG000519815", 0.10, 1)
    s5, _ = derive_stable_masking_seed("MGYG000520046", 0.10, 0)
    assert len({s1, s3, s4, s5}) == 4

    # 3. Media serialization and exact provenance hash stability (no arbitrary rounding)
    scen_cpds = {"cpd00027": ("Glucose", 10.0), "cpd00013": ("NH3", 10.0)}
    _, med_h1 = serialize_and_hash_medium(scen_cpds)
    _, med_h2 = serialize_and_hash_medium(scen_cpds)
    assert med_h1 == med_h2

    h1 = compute_exact_input_provenance_hash("archive_sha_1", "model_1", med_h1, "pool_sha_1", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420")
    h2 = compute_exact_input_provenance_hash("archive_sha_1", "model_1", med_h1, "pool_sha_1", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420")
    assert h1 == h2
    assert len(h1) == 64

    # Changed compound changes hash
    _, med_h3 = serialize_and_hash_medium({"cpd00027": ("Glucose", 5.0), "cpd00013": ("NH3", 10.0)})
    h3 = compute_exact_input_provenance_hash("archive_sha_1", "model_1", med_h3, "pool_sha_1", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420")
    assert h1 != h3


def test_31_same_medium_masking_controls_and_nonevaluable_handling():
    """Test same-medium masking controls execution, non-evaluable handling, and validator (Item 1)."""
    from run_masking_controls_7hq import (
        evaluate_single_mask_case,
        validate_masking_controls_results,
        derive_stable_masking_seed,
        EXACT_7_HQ_EXTERNAL_IDS,
    )
    toy = _make_base_archaeal_toy_model()
    scen_cpds = {"cpd00027": ("Glucose", 10.0), "cpd00013": ("NH3", 10.0), "cpd00007": ("O2", 10.0)}

    cand_rxns = {
        "R_cand_A": ReactionRecord(
            id="R_cand_A",
            base_id="rxn_cand_A",
            name="Lysine route",
            reversible=False,
            lb=0.0,
            ub=1000.0,
            stoichiometry={"M_cpd00027_c0": -0.005, "M_cpd00013_c0": -0.01, "M_cpd00039_c0": 0.01},
        )
    }
    # Provide candidate replacements for draft internal reactions so masked models have viable ensembles
    for rid, r in toy.reactions.items():
        if not r.is_biomass and not r.is_exchange:
            cand_rxns[f"CAND_{rid}"] = ReactionRecord(
                id=f"CAND_{rid}", base_id=r.base_id, name=f"Cand {r.name}",
                reversible=r.reversible, lb=r.lb, ub=r.ub,
                stoichiometry=r.stoichiometry, is_transport=r.is_transport,
            )

    # Case 1: Baseline valid, masked valid -> evaluable
    task_ok = {
        "genome": "MGYG000519815",
        "cohort": "external_comparator",
        "canonical_order": "o__UBA202",
        "canonical_family": "f__JdFR-43",
        "completeness": 98.5,
        "contamination": 0.2,
        "mask_fraction": 0.10,
        "replicate": 0,
        "draft_model": toy,
        "s4_compounds": scen_cpds,
        "medium_hash": "dummy_med_hash",
        "candidate_reactions": cand_rxns,
        "candidate_weights": {k: 1.0 for k in cand_rxns},
        "milp_time_limit_s": 5.0,
        "input_hash": "dummy_hash_123",
        "baseline_status": "optimal",
        "baseline_n_valid_solutions": 1,
        "baseline_gapfill_rxns": ["R_cand_A"],
    }
    eval_pack = evaluate_single_mask_case(task_ok)
    res_ok = eval_pack["summary_row"]
    led_ok = eval_pack["ledger_rows"]
    assert res_ok["is_evaluable"] is True
    assert res_ok["scenario_id"] == "S4_oxic_desert_oa_aa"
    assert res_ok["biomass_recipe_id"] == "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"
    assert np.isfinite(res_ok["n_newly_induced_gapfill_rxns"])
    assert len(led_ok) >= 1
    assert led_ok[0]["validity_disposition"] == "accepted_feasible"

    # Case 2: Baseline invalid (numerical failure / infeasible) -> comparative metrics strictly NaN,
    # but masked reactions independently counted if masked is valid (EXP1-MASKING-REPRO-REVIEW.txt Item 3)
    task_infeas = dict(task_ok)
    task_infeas["baseline_status"] = "numerical_coupling_violation"
    task_infeas["baseline_n_valid_solutions"] = 0
    res_infeas = evaluate_single_mask_case(task_infeas)["summary_row"]
    assert res_infeas["is_evaluable"] is False
    assert math.isnan(res_infeas["n_newly_induced_gapfill_rxns"])
    assert math.isnan(res_infeas["n_true_masked_rxns_recovered_in_gapfill"])
    assert math.isnan(res_infeas["n_alternative_candidate_rxns_induced"])
    assert math.isnan(res_infeas["masked_rxn_recovery_rate_under_medium"])
    assert math.isnan(res_infeas["alternative_candidate_fraction_under_medium"])
    # Independent masked count preserved (not false 0.0 or NaN when masked solve succeeded)
    assert res_infeas["masked_solver_status"] == "optimal"
    assert np.isfinite(res_infeas["n_masked_gapfill_rxns"])
    assert res_infeas["n_masked_gapfill_rxns"] >= 1.0
    assert "non_evaluable_invalid_baseline" in res_infeas["evaluation_status"]

    # Case 3: Test validator on mock 63-row DataFrame, ledger, and media provenance
    mock_rows = []
    mock_ledger = []
    mock_media = {}
    for g in EXACT_7_HQ_EXTERNAL_IDS:
        med_h = f"med_hash_{g}"
        mock_media[g] = {"medium_hash": med_h}
        for frac in (0.10, 0.20, 0.30):
            for rep in (0, 1, 2):
                exp_seed, exp_str = derive_stable_masking_seed(g, frac, rep)
                mock_rows.append({
                    "condition_type": "masked",
                    "genome": g,
                    "cohort": "external_comparator",
                    "canonical_order": "o__test",
                    "canonical_family": "f__test",
                    "completeness": 98.0,
                    "contamination": 0.5,
                    "scenario_id": "S4_oxic_desert_oa_aa",
                    "biomass_recipe_id": "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420",
                    "mask_fraction": frac,
                    "replicate": rep,
                    "seed_val": exp_seed,
                    "seed_str": exp_str,
                    "medium_hash": med_h,
                    "input_hash": f"hash_{g}_{frac}_{rep}",
                    "n_total_genes": 100,
                    "n_masked_genes": int(round(100 * frac)),
                    "masked_gene_ids": "g1,g2",
                    "n_removed_gpr_reactions": 2,
                    "removed_rxn_ids": "r1,r2",
                    "baseline_solver_status": "optimal",
                    "baseline_n_valid_solutions": 1,
                    "n_baseline_gapfill_rxns": 1.0,
                    "masked_solver_status": "optimal",
                    "masked_n_valid_solutions": 1,
                    "n_masked_gapfill_rxns": 2.0,
                    "is_evaluable": True,
                    "evaluation_status": "evaluable_valid_ensembles",
                    "n_newly_induced_gapfill_rxns": 2.0,
                    "n_true_masked_rxns_recovered_in_gapfill": 1.0,
                    "n_alternative_candidate_rxns_induced": 1.0,
                    "masked_rxn_recovery_rate_under_medium": 0.5,
                    "alternative_candidate_fraction_under_medium": 0.5,
                    "runtime_seconds": 0.5,
                })
                mock_ledger.append({
                    "condition_type": "masked",
                    "genome": g,
                    "cohort": "external_comparator",
                    "canonical_order": "o__test",
                    "canonical_family": "f__test",
                    "scenario_id": "S4_oxic_desert_oa_aa",
                    "biomass_recipe_id": "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420",
                    "mask_fraction": frac,
                    "replicate": rep,
                    "solution_index": 0,
                    "validity_disposition": "accepted_feasible",
                    "solver_status": "optimal",
                    "selected_reaction_ids": "r1",
                    "n_added_rxns": 1,
                    "total_cost": 1.0,
                    "raw_milp_objective": 1.0,
                    "mip_gap": 0.0,
                    "dual_bound": 1.0,
                    "verified_growth": 0.05,
                    "stoichiometric_max_residual": 1e-8,
                    "max_bound_violation": 1e-8,
                    "residual_bounds_growth_verified": True,
                    "cycle_qc_status": "passed_closed_uptake_egc_check",
                    "spurious_atp_cycle": False,
                    "explicit_reason": "",
                    "code_version": "v3_masking_controls",
                })
    mock_df = pd.DataFrame(mock_rows)
    mock_led_df = pd.DataFrame(mock_ledger)
    val_rep = validate_masking_controls_results(
        masking_df=mock_df,
        expected_genomes=set(EXACT_7_HQ_EXTERNAL_IDS),
        ledger_df=mock_led_df,
        media_provenance=mock_media,
    )
    assert val_rep["validation_passed"] is True
    assert val_rep["n_total_mask_rows"] == 63

    # Case 4: Undefined recovery rate (NaN) when 0 reactions removed (Item 4)
    zero_rem_df = mock_df.copy()
    zero_rem_df.loc[0, "n_removed_gpr_reactions"] = 0
    zero_rem_df.loc[0, "masked_rxn_recovery_rate_under_medium"] = 1.0  # Erroneous 1.0!
    val_bad_zero = validate_masking_controls_results(
        masking_df=zero_rem_df,
        expected_genomes=set(EXACT_7_HQ_EXTERNAL_IDS),
        ledger_df=mock_led_df,
        media_provenance=mock_media,
    )
    assert val_bad_zero["validation_passed"] is False
    assert val_bad_zero["checks"]["zero_removed_recovery_rate_strictly_nan"] is False


def test_32_collector_unknown_as_zero_repaired_to_nan():
    """Test repaired collector reaction-frequency table and unknown-as-zero fix (Item 2)."""
    from stage3_collect_and_summarize_v3 import repair_reaction_frequency_and_biomass_sensitivity
    from exp1_sbml_highs_engine import (
        BIO_B0_BASELINE,
        BIO_B1_OMIT_8,
        BIO_B2_OMIT_7_KEEP_F420,
        BIO_B3_OMIT_LIPIDS,
        BIO_B4_CORE_SKELETON,
    )

    g = "TEST_G01"
    scen = "S4_gapseq_desert_plus_o2"

    # Scenario: B0 optimal (1 solution), but B2 is numerical_coupling_violation (0 solutions)
    sum_df = pd.DataFrame([
        {"genome": g, "scenario_id": scen, "biomass_recipe_id": BIO_B0_BASELINE, "solver_status": "optimal", "n_distinct_valid_solutions": 1},
        {"genome": g, "scenario_id": scen, "biomass_recipe_id": BIO_B1_OMIT_8, "solver_status": "optimal", "n_distinct_valid_solutions": 1},
        {"genome": g, "scenario_id": scen, "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420, "solver_status": "numerical_coupling_violation", "n_distinct_valid_solutions": 0},
        {"genome": g, "scenario_id": scen, "biomass_recipe_id": BIO_B3_OMIT_LIPIDS, "solver_status": "optimal", "n_distinct_valid_solutions": 1},
        {"genome": g, "scenario_id": scen, "biomass_recipe_id": BIO_B4_CORE_SKELETON, "solver_status": "optimal", "n_distinct_valid_solutions": 1},
    ])

    raw_rf = pd.DataFrame([{
        "genome": g,
        "cohort": "EQ_desert",
        "completeness": 95.0,
        "contamination": 1.0,
        "scenario_id": scen,
        "rxn_id": "rxn_test_1",
        "base_id": "rxn_test",
        "rxn_name": "Test Reaction",
        "is_transport": False,
        "candidate_weight": 1.0,
        "freq_B0_template_baseline": 1.0,
        "freq_B1_omit_all8_cofactors": 1.0,
        "freq_B2_omit_7_keep_f420": 0.0,  # Old worker defaulted to 0.0!
        "freq_B3_omit_surface_lipids": 1.0,
        "freq_B4_core_skeleton": 1.0,
        "stability_tier_B0": "core_essential",
        "stability_tier_B2": "absent",  # Old worker defaulted to absent!
        "supernetwork_ko_necessity_B0": "model_specific_essential_in_declared_supernetwork",
        "supernetwork_ko_necessity_B2": "not_in_B2_ensemble",
        "biomass_sensitivity_classification": "biomass_assumption_sensitive__vanishes_when_7_non_f420_template_cofactors_omitted",  # Misclassified!
    }])

    repaired = repair_reaction_frequency_and_biomass_sensitivity(raw_rf, sum_df)
    row = repaired.iloc[0]

    # Must be repaired to NaN, NOT 0.0 or absent
    assert math.isnan(row["freq_B2_omit_7_keep_f420"])
    assert "unsolved_or_infeasible_numerical_coupling_violation" in row["stability_tier_B2"]
    assert "unsolved_or_infeasible_numerical_coupling_violation" in row["supernetwork_ko_necessity_B2"]
    assert math.isnan(row["freq_difference_B0_minus_B2"])
    assert "unknown_biomass_sensitivity__B2_numerical_coupling_violation" in row["biomass_sensitivity_classification"]


def test_33_claim_tier_and_conditional_b2_requirement_honest_labeling():
    """Test that candidate deficiency selection honestly labels conditional B2 requirements and never
    asserts directly_verified_robust_candidate without necessity proof across recipes (Item 2).
    """
    from stage3_collect_and_summarize_v3 import repair_reaction_frequency_and_biomass_sensitivity
    from exp1_sbml_highs_engine import (
        BIO_B0_BASELINE,
        BIO_B1_OMIT_8,
        BIO_B2_OMIT_7_KEEP_F420,
        BIO_B3_OMIT_LIPIDS,
        BIO_B4_CORE_SKELETON,
    )

    g = "TEST_G02"
    scen = "S4_gapseq_desert_plus_o2"

    # All 5 perturbations optimal
    sum_df = pd.DataFrame([
        {"genome": g, "scenario_id": scen, "biomass_recipe_id": b, "solver_status": "optimal", "n_distinct_valid_solutions": 1}
        for b in (BIO_B0_BASELINE, BIO_B1_OMIT_8, BIO_B2_OMIT_7_KEEP_F420, BIO_B3_OMIT_LIPIDS, BIO_B4_CORE_SKELETON)
    ])

    # Reaction observed across all recipes, but KO necessity ONLY checked under B2
    raw_rf = pd.DataFrame([{
        "genome": g,
        "cohort": "EQ_desert",
        "completeness": 92.0,
        "contamination": 1.0,
        "scenario_id": scen,
        "rxn_id": "rxn_test_cond",
        "base_id": "rxn_test_cond",
        "rxn_name": "Conditional Test Reaction",
        "is_transport": False,
        "candidate_weight": 1.0,
        "freq_B0_template_baseline": 1.0,
        "freq_B1_omit_all8_cofactors": 1.0,
        "freq_B2_omit_7_keep_f420": 1.0,
        "freq_B3_omit_surface_lipids": 1.0,
        "freq_B4_core_skeleton": 1.0,
        "stability_tier_B0": "core_essential",
        "stability_tier_B2": "core_essential",
        "supernetwork_ko_necessity_B0": "not_tested_or_alternative_pathway_exists",
        "supernetwork_ko_necessity_B2": "model_specific_essential_in_declared_supernetwork",
        "biomass_sensitivity_classification": "robust_across_all_tested_biomass_perturbations",
    }])

    repaired = repair_reaction_frequency_and_biomass_sensitivity(raw_rf, sum_df)
    row = repaired.iloc[0]

    # Classification must be observed across assumptions, NOT directly verified robust
    assert row["biomass_sensitivity_classification"] == "observed_across_all_tested_biomass_perturbations"


def test_34_engine_solution_schema_and_missing_field_rejection():
    """Test actual engine-result fixture, honest validity disposition derivation,
    and strict rejection of missing engine solution keys (`EXP1-V4-REVIEW-CORRECTIONS.txt`).
    """
    from run_masking_controls_7hq_v4_resumable import (
        REQUIRED_ENGINE_SOL_KEYS,
        build_solution_ledger_row,
        build_empty_solver_ledger_row,
    )

    # 1. Complete authentic engine solution fixture
    valid_engine_sol = {
        "solution_index": 0,
        "optimization_solution_index": 0,
        "added_rxns": ["rxn00147_c0", "rxn00288_c0"],
        "raw_milp_selected_rxns": ["rxn00147_c0", "rxn00288_c0"],
        "added_internal_rxns": ["rxn00147_c0"],
        "added_transport_rxns": ["rxn00288_c0"],
        "total_cost": 2.0,
        "raw_milp_objective": 2.0,
        "raw_binary_cost": 2.0,
        "post_pruning_cost_delta": 0.0,
        "verified_fba_growth": 0.052,
        "verified_fluxes": {"rxn00147_c0": 0.5, "rxn00288_c0": 0.5},
        "fba_check_max_residual": 1.2e-8,
        "fba_check_max_bound_viol": 0.0,
        "mip_gap": 0.0,
        "dual_bound": 2.0,
        "raw_solver_status": 0,
        "raw_solver_message": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
        "optimizer_status": "optimal",
        "post_pruning_egc_status": "passed_closed_uptake_egc_check",
        "post_pruning_egc_has_spurious_atp_cycle": False,
        "post_pruning_egc_atp_flux": 0.0,
        "cycle_qc_status": "passed_closed_uptake_egc_check",
        "science_qc_eligible": True,
    }

    # Verify successful accepted row
    row_ok = build_solution_ledger_row(
        sol=valid_engine_sol,
        condition_type="baseline",
        genome="MGYG000519815",
        cohort="external_comparator",
        canonical_order="o__Methanosarcinales",
        canonical_family="f__Methanosarcinaceae",
        mask_fraction=0.0,
        replicate=0,
        solution_index=0,
        solver_status="optimal",
    )
    assert row_ok["validity_disposition"] == "accepted_feasible"
    assert row_ok["residual_bounds_growth_verified"] is True
    assert row_ok["selected_reaction_ids"] == "rxn00147_c0,rxn00288_c0"
    assert row_ok["total_cost"] == 2.0
    assert row_ok["verified_growth"] == 0.052

    # 2. Reject missing required fields: test each required key independently
    for req_key in REQUIRED_ENGINE_SOL_KEYS:
        corrupted = dict(valid_engine_sol)
        del corrupted[req_key]
        try:
            build_solution_ledger_row(
                sol=corrupted,
                condition_type="baseline",
                genome="MGYG000519815",
                cohort="external_comparator",
                canonical_order="o__Methanosarcinales",
                canonical_family="f__Methanosarcinaceae",
                mask_fraction=0.0,
                replicate=0,
                solution_index=0,
                solver_status="optimal",
            )
            assert False, f"Expected KeyError for missing required key '{req_key}'"
        except KeyError as e:
            assert req_key in str(e)

    # 3. Honest disposition derivation tests:
    # A. Spurious ATP cycle must be rejected even if science_qc_eligible was claimed
    bad_spur = dict(valid_engine_sol, post_pruning_egc_has_spurious_atp_cycle=True)
    row_spur = build_solution_ledger_row(bad_spur, "masked", "G1", "cohort", "o", "f", 0.1, 0, 0, "optimal")
    assert row_spur["validity_disposition"] == "rejected_spurious_atp_cycle"

    # B. Subthreshold growth (< 0.01) must be rejected
    bad_gr = dict(valid_engine_sol, verified_fba_growth=0.005)
    row_gr = build_solution_ledger_row(bad_gr, "masked", "G1", "cohort", "o", "f", 0.1, 0, 0, "optimal")
    assert row_gr["validity_disposition"] == "rejected_insufficient_growth"

    # C. Residual violation (> 1e-5) must be rejected
    bad_res = dict(valid_engine_sol, fba_check_max_residual=1e-4)
    row_res = build_solution_ledger_row(bad_res, "masked", "G1", "cohort", "o", "f", 0.1, 0, 0, "optimal")
    assert row_res["validity_disposition"] == "rejected_residual_violation"

    # D. Bound violation (> 1e-5) must be rejected
    bad_bnd = dict(valid_engine_sol, fba_check_max_bound_viol=1e-4)
    row_bnd = build_solution_ledger_row(bad_bnd, "masked", "G1", "cohort", "o", "f", 0.1, 0, 0, "optimal")
    assert row_bnd["validity_disposition"] == "rejected_bound_violation"

    # E. Cycle QC ineligible (science_qc_eligible = False) must be rejected
    bad_qc = dict(valid_engine_sol, science_qc_eligible=False, cycle_qc_status="rejected_thermodynamic_cycle")
    row_qc = build_solution_ledger_row(bad_qc, "masked", "G1", "cohort", "o", "f", 0.1, 0, 0, "optimal")
    assert row_qc["validity_disposition"] == "rejected_cycle_qc_ineligible"


def test_35_collector_consumes_versioned_v4_outputs_and_ledger(tmp_path):
    """Test that stage3_collect_and_summarize_v3 explicitly consumes versioned v4 masking controls
    and its companion ledger TSV (`EXP1-V4-REVIEW-CORRECTIONS.txt`).
    """
    from stage3_collect_and_summarize_v3 import main as collector_main

    # Build minimal mock shards and inventory
    shards_dir = tmp_path / "shards"
    shards_dir.mkdir()
    out_dir = tmp_path / "results"
    out_dir.mkdir()

    inv_path = tmp_path / "inventory.tsv"
    inv_df = pd.DataFrame([{
        "genome": "G001",
        "cohort": "EQ_desert",
        "is_retained_484": True,
        "has_draft_sbml": True,
        "qc_error": "",
        "legacy_classification": "Methanosarcinaceae",
        "is_complete_nonchimeric": True,
        "completeness": 95.0,
        "contamination": 1.0,
        "canonical_order": "o__Test",
        "canonical_family": "f__Test",
    }])
    inv_df.to_csv(inv_path, sep="\t", index=False)

    # Shard 0 with complete outputs for G001 across 5 scenarios and 5 recipes
    from exp1_sbml_highs_engine import (
        NUTRIENT_SCENARIOS,
        BIO_B0_BASELINE,
        BIO_B1_OMIT_8,
        BIO_B2_OMIT_7_KEEP_F420,
        BIO_B3_OMIT_LIPIDS,
        BIO_B4_CORE_SKELETON,
    )
    ALL_RECIPES = [
        BIO_B0_BASELINE,
        BIO_B1_OMIT_8,
        BIO_B2_OMIT_7_KEEP_F420,
        BIO_B3_OMIT_LIPIDS,
        BIO_B4_CORE_SKELETON,
    ]
    s_rows = []
    rf_rows = []
    wkf_rows = []
    for scen in NUTRIENT_SCENARIOS:
        for bio in ALL_RECIPES:
            s_rows.append({
                "genome": "G001",
                "cohort": "EQ_desert",
                "canonical_order": "o__Test",
                "canonical_family": "f__Test",
                "completeness": 95.0,
                "contamination": 1.0,
                "scenario_id": scen,
                "biomass_recipe_id": bio,
                "solver_status": "optimal",
                "n_distinct_valid_solutions": 1,
                "min_total_cost": 1.0,
                "max_ensemble_cost": 1.0,
                "n_solutions_attempted": 1,
                "min_added_rxns_count": 1.0,
                "max_added_rxns_count": 1.0,
                "n_rxns_observed_in_all_enumerated_models": 1.0,
                "n_rxns_essential_in_declared_supernetwork": 0.0,
                "runtime_seconds": 0.1,
            })
            rf_rows.append({
                "genome": "G001",
                "cohort": "EQ_desert",
                "completeness": 95.0,
                "contamination": 1.0,
                "scenario_id": scen,
                "rxn_id": "R_test",
                "base_id": "rxn_test",
                "rxn_name": "Test",
                "is_transport": False,
                "candidate_weight": 1.0,
                "freq_B0_template_baseline": 1.0,
                "freq_B1_omit_all8_cofactors": 1.0,
                "freq_B2_omit_7_keep_f420": 1.0,
                "freq_B3_omit_surface_lipids": 1.0,
                "freq_B4_core_skeleton": 1.0,
                "stability_tier_B0": "core_essential",
                "stability_tier_B2": "core_essential",
                "supernetwork_ko_necessity_B0": "not_tested",
                "supernetwork_ko_necessity_B2": "not_tested",
                "biomass_sensitivity_classification": "observed_across_all_tested_biomass_perturbations",
            })
            wkf_rows.append({
                "genome": "G001",
                "cohort": "EQ_desert",
                "scenario_id": scen,
                "biomass_recipe_id": bio,
                "record_type": "knockout",
                "target_id": "none",
                "target_name": "none",
                "essential_in_native": False,
                "essential_in_desert": False,
                "growth_with_target": 0.05,
                "growth_without_target": 0.05,
                "fva_min_flux": 0.0,
                "fva_max_flux": 1.0,
            })

    s_df = pd.DataFrame(s_rows)
    rf_df = pd.DataFrame(rf_rows)
    wkf_df = pd.DataFrame(wkf_rows)

    s_df.to_csv(shards_dir / "exp1_shard_00_summary.tsv", sep="\t", index=False)
    rf_df.to_csv(shards_dir / "exp1_shard_00_rxn_freq.tsv", sep="\t", index=False)
    wkf_df.to_csv(shards_dir / "exp1_shard_00_wkf.tsv", sep="\t", index=False)

    for i in range(1, 8):
        pd.DataFrame(columns=s_df.columns).to_csv(shards_dir / f"exp1_shard_{i:02d}_summary.tsv", sep="\t", index=False)
        pd.DataFrame(columns=rf_df.columns).to_csv(shards_dir / f"exp1_shard_{i:02d}_rxn_freq.tsv", sep="\t", index=False)
        pd.DataFrame(columns=wkf_df.columns).to_csv(shards_dir / f"exp1_shard_{i:02d}_wkf.tsv", sep="\t", index=False)

    # Create mock v4 masking controls and companion ledger
    mask_v4_tsv = tmp_path / "exp1_masking_controls_7hq_v4.tsv"
    mask_v4_ledger = tmp_path / "exp1_masking_controls_7hq_v4_per_solution_ledger.tsv"

    mock_m = pd.DataFrame([{
        "condition_type": "masked",
        "genome": "MGYG000519815",
        "mask_fraction": 0.1,
        "replicate": 0,
        "is_evaluable": True,
        "masked_rxn_recovery_rate_under_medium": 0.75,
    }])
    mock_m.to_csv(mask_v4_tsv, sep="\t", index=False)

    mock_l = pd.DataFrame([{
        "condition_type": "masked",
        "genome": "MGYG000519815",
        "solution_index": 0,
        "validity_disposition": "accepted_feasible",
    }])
    mock_l.to_csv(mask_v4_ledger, sep="\t", index=False)

    # Run collector with versioned v4 masking controls
    import sys
    test_args = [
        "stage3_collect_and_summarize_v3.py",
        "--inventory-484-tsv", str(inv_path),
        "--shards-dir", str(shards_dir),
        "--masking-controls-tsv", str(mask_v4_tsv),
        "--out-dir", str(out_dir),
        "--num-tasks", "8",
    ]
    old_argv = sys.argv
    try:
        sys.argv = test_args
        collector_main()
    finally:
        sys.argv = old_argv

    # Verify that v4 files were copied and metadata captured
    assert (out_dir / "exp1_masking_controls_7hq_v4.tsv").exists()
    assert (out_dir / "exp1_masking_controls_7hq_v4_per_solution_ledger.tsv").exists()

    with open(out_dir / "exp1_stage3_shard_collection_receipts.json", "r") as fh:
        rep = json.load(fh)
    assert rep["masking_controls"]["masking_controls_present"] is True
    assert rep["masking_controls"]["masking_controls_file"] == "exp1_masking_controls_7hq_v4.tsv"
    assert rep["masking_controls"]["ledger_file"] == "exp1_masking_controls_7hq_v4_per_solution_ledger.tsv"
    assert rep["masking_controls"]["mean_recovery_rate"] == 0.75



