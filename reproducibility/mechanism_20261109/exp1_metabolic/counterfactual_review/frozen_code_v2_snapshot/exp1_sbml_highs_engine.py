#!/usr/bin/env python3
"""Core SBML Level 3 FBC v2 parser, QC auditor, Methanosarcina biomass perturbation engine,
and HiGHS binary MILP (`scipy.optimize.milp`) alternative-gap ensemble solver for Experiment 1 (`exp1_metabolic`).

Implements all scientific & methodological requirements of:
- `mechanism-20261109/EXP1-CONTRACT.txt`
- `mechanism-20261109/EXP1-BIOMASS-SCIENTIFIC-NOTE.txt`
- `mechanism-20261109/EXP1-SOLVER-REVIEW.txt`
- `coordinated-revision-20261108/mechanism-experiments-20261109.txt`
"""
from __future__ import annotations

import csv
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy import sparse
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

# =============================================================================
# 1. Explicit Constants, Provenance, Biomass Perturbations, and Scenarios
# =============================================================================

METHANOSARCINA_BIOMASS_PROVENANCE = (
    "gapseq 2.1.0 dat/biomass/biomass_archaea.json derived from Methanosarcina barkeri iAF692 "
    "(Feist et al. 2006, doi:10.1038/msb4100046); template baseline for comparison, NOT validated "
    "Thermoplasmatota biomass."
)

# 7 Methanosarcina template cofactors associated with methanogenesis / Methanosarcina-specific C1 metabolism:
METHANOSARCINA_NON_F420_TEMPLATE_COFACTORS: dict[str, str] = {
    "cpd02246": "coenzyme m",
    "cpd02817": "HTP (7-mercaptoheptanoylthreonine phosphate / coenzyme B)",
    "cpd03425": "Factor F430",
    "cpd15862": "coenzyme F390 (adenosine)",
    "cpd15863": "coenzyme F390 (guanosine)",
    "cpd15884": "tetrahydrosarcinapterin",
    "cpd33626": "adenosyl factor III",
}

# Coenzyme F420-2 is in the M. barkeri template but is NOT exclusive to methanogens (occurs across diverse archaea/bacteria):
F420_COFACTOR_ID = "cpd00649"
F420_COFACTOR_NAME = "Coenzyme F420-2"

# All 8 special trace cofactors in the M. barkeri template Cofactors group beyond universal NAD/NADP/CoA/THF/ATP:
METHANOSARCINA_ALL8_TEMPLATE_COFACTORS: dict[str, str] = {
    F420_COFACTOR_ID: F420_COFACTOR_NAME,
    **METHANOSARCINA_NON_F420_TEMPLATE_COFACTORS,
}

# Checked empirically in R_bio1 (absent from gapseq biomass_archaea.json):
METHANOFURAN_COMPOUND_IDS: dict[str, str] = {
    "cpd02431": "Methanofuran",
    "cpd15883": "Methanofuran-b",
}

ARCHAEAL_ETHER_LIPIDS_AND_POLYAMINES: dict[str, str] = {
    "cpd15817": "2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycero-l-phosphoethanolamine",
    "cpd15818": "2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycerol-1-phospho-3'-sn-glycerol",
    "cpd15819": "2-O-(3'-hydroxyl)phytanyl-3-O-phytanyl-sn-glycero-1-phospho-myo-inositol",
    "cpd15820": "2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycero-l-phosphoserine",
    "cpd15858": "2,3-O-phytanyl-sn-glycero-1-phosphoethanolamine",
    "cpd15859": "2,3-di-O-phytanyl-sn-glycerol-1-phospho-3'-sn-glycerol",
    "cpd15860": "2,3-O-phytanyl-sn-glycero-1-phospho-myo-inositol",
    "cpd15861": "2,3-O-phytanyl-sn-glycero-1-phosphoserine",
    "cpd15880": "glucosaminyl archaetidyl-myo-inositol",
    "cpd00118": "Putrescine",
    "cpd03802": "sym-Homospermidine",
}

# Bacterial cell-wall / peptidoglycan / LPS markers that MUST NEVER appear in an archaeal biomass objective:
BACTERIAL_BIOMASS_FORBIDDEN_CPDS: dict[str, str] = {
    "cpd15665": "Peptidoglycan polymer (n subunits)",
    "cpd15666": "Peptidoglycan polymer (n-1 subunits)",
    "cpd15529": "Lipopolysaccharide",
    "cpd15667": "Undecaprenyl-diphospho-N-acetylmuramoyl-pentapeptide",
    "cpd02229": "Bactoprenyl diphosphate",
    "cpd15386": "Teichoic acid",
}

# Declared component-omission sensitivity perturbations (`EXP1-BIOMASS-SCIENTIFIC-NOTE.txt`):
BIOMASS_PERTURBATION_RECIPES: dict[str, dict[str, Any]] = {
    "B0_methanosarcina_template_baseline": {
        "description": "Unmodified gapseq archaeal biomass (Feist 2006 Methanosarcina barkeri iAF692 template).",
        "epistemic_role": "template_baseline_for_comparison_not_ground_truth",
        "omit_cpd_ids": set(),
    },
    "B1_perturbation_omit_8_template_cofactors": {
        "description": (
            "Declared component-omission sensitivity perturbation removing the 8 M. barkeri template trace cofactors "
            "(cpd00649 F420-2, cpd02246 CoM, cpd02817 HTP/CoB, cpd03425 F430, cpd15862 F390-A, cpd15863 F390-G, "
            "cpd15884 tetrahydrosarcinapterin, cpd33626 adenosyl factor III). Not an experimentally validated recipe."
        ),
        "epistemic_role": "component_omission_perturbation_not_experimentally_validated_biomass",
        "omit_cpd_ids": set(METHANOSARCINA_ALL8_TEMPLATE_COFACTORS.keys()),
    },
    "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420": {
        "description": (
            "Declared component-omission sensitivity perturbation removing the 7 non-F420 M. barkeri template cofactors "
            "(CoM, HTP/CoB, F430, F390-A, F390-G, tetrahydrosarcinapterin, adenosyl factor III) while retaining "
            "Coenzyme F420-2 (cpd00649, which is not exclusive to methanogens). Not an experimentally validated recipe."
        ),
        "epistemic_role": "component_omission_perturbation_not_experimentally_validated_biomass",
        "omit_cpd_ids": set(METHANOSARCINA_NON_F420_TEMPLATE_COFACTORS.keys()),
    },
    "B3_perturbation_omit_template_cofactors_and_surface_glycolipids": {
        "description": (
            "Declared component-omission sensitivity perturbation removing the 8 template trace cofactors plus "
            "sym-Homospermidine (cpd03802), glucosaminyl archaetidyl-myo-inositol (cpd15880), and hydroxyphytanyl "
            "lipids (cpd15817..cpd15820), retaining core diether archaetidyl lipids (cpd15858..cpd15861). "
            "Not an experimentally validated recipe."
        ),
        "epistemic_role": "component_omission_perturbation_not_experimentally_validated_biomass",
        "omit_cpd_ids": set(METHANOSARCINA_ALL8_TEMPLATE_COFACTORS.keys())
        | {"cpd03802", "cpd15880", "cpd15817", "cpd15818", "cpd15819", "cpd15820"},
    },
    "B4_perturbation_core_amino_nucleotide_energy_skeleton": {
        "description": (
            "Declared component-omission sensitivity perturbation retaining DNA, RNA, Protein (20 amino acids), "
            "and universal redox/energy cofactors (ATP/ADP/AMP, NAD/NADH, NADP/NADPH, CoA, Acetyl-CoA, Succinyl-CoA, THF). "
            "Not an experimentally validated recipe."
        ),
        "epistemic_role": "component_omission_perturbation_not_experimentally_validated_biomass",
        "omit_cpd_ids": set(METHANOSARCINA_ALL8_TEMPLATE_COFACTORS.keys())
        | set(ARCHAEAL_ETHER_LIPIDS_AND_POLYAMINES.keys())
        | {"cpd00155"},
    },
}

BASE_INORGANIC_SALTS: dict[str, tuple[str, float]] = {
    "cpd00001": ("H2O", 100.0),
    "cpd00009": ("Phosphate", 10.0),
    "cpd00011": ("CO2", 20.0),
    "cpd00013": ("Ammonium", 10.0),
    "cpd00048": ("Sulfate", 10.0),
    "cpd00239": ("H2S", 2.0),
    "cpd00971": ("Sodium", 10.0),
    "cpd00205": ("K+", 10.0),
    "cpd00099": ("Cl-", 10.0),
    "cpd00254": ("Mg2+", 10.0),
    "cpd00063": ("Ca2+", 10.0),
    "cpd10515": ("Fe2+", 10.0),
    "cpd10516": ("Fe3+", 1.0),
    "cpd00149": ("Cobalt", 2.0),
    "cpd00244": ("Nickel", 2.0),
    "cpd00034": ("Zn2+", 2.0),
    "cpd00058": ("Cu2+", 2.0),
    "cpd00030": ("Mn2+", 2.0),
    "cpd00067": ("H+", 100.0),
}

DESERT_OA_AA_ADDITIONS: dict[str, tuple[str, float]] = {
    "cpd11640": ("H2", 5.0),
    "cpd00029": ("Acetate", 2.0),
    "cpd00159": ("L-Lactate", 2.0),
    "cpd00020": ("Pyruvate", 1.0),
    "cpd00047": ("Formate", 2.0),
    "cpd00100": ("Glycerol", 1.0),
    "cpd00209": ("Nitrate", 10.0),
    # 20 standard amino acids and cofactors present in gapseq desert_medium.csv
    "cpd00033": ("Glycine", 0.5),
    "cpd00035": ("L-Alanine", 0.5),
    "cpd00023": ("L-Glutamate", 0.5),
    "cpd00051": ("L-Arginine", 0.1),
    "cpd00132": ("L-Asparagine", 0.1),
    "cpd00041": ("L-Aspartate", 0.1),
    "cpd00084": ("L-Cysteine", 0.1),
    "cpd00053": ("L-Glutamine", 0.1),
    "cpd00119": ("L-Histidine", 0.1),
    "cpd00322": ("L-Isoleucine", 0.1),
    "cpd00107": ("L-Leucine", 0.1),
    "cpd00039": ("L-Lysine", 0.1),
    "cpd00060": ("L-Methionine", 0.1),
    "cpd00066": ("L-Phenylalanine", 0.1),
    "cpd00129": ("L-Proline", 0.1),
    "cpd00054": ("L-Serine", 0.1),
    "cpd00161": ("L-Threonine", 0.1),
    "cpd00065": ("L-Tryptophan", 0.1),
    "cpd00069": ("L-Tyrosine", 0.1),
    "cpd00156": ("L-Valine", 0.1),
    "cpd00216": ("Chorismate", 0.1),
    "cpd00305": ("Thiamin", 0.025),
    "cpd00220": ("Riboflavin", 0.025),
}

DEFINED_20AA_AND_PRECURSORS: dict[str, tuple[str, float]] = {
    **DESERT_OA_AA_ADDITIONS,
    "cpd00118": ("Putrescine", 0.1),
    "cpd00218": ("Niacin", 0.05),
    "cpd00644": ("Pantothenate", 0.05),
    "cpd00393": ("Folate", 0.05),
    "cpd00263": ("Pyridoxine", 0.05),
}

SCEN_S0_NATIVE = "S0_gapseq_native_predicted"
SCEN_S1_ANOXIC_MINIMAL = "S1_anoxic_minimal_h2_co2_nh4"
SCEN_S2_OXIC_MINIMAL = "S2_oxic_minimal_h2_co2_nh4"
SCEN_S3_ANOXIC_DESERT = "S3_anoxic_desert_oa_aa"
SCEN_S4_OXIC_DESERT = "S4_oxic_desert_oa_aa"
SCEN_S5_ANOXIC_DEFINED = "S5_anoxic_defined_aa_precursors"
SCEN_S6_OXIC_DEFINED = "S6_oxic_defined_aa_precursors"

BIO_B0_BASELINE = "B0_methanosarcina_template_baseline"
BIO_B1_OMIT_8 = "B1_perturbation_omit_8_template_cofactors"
BIO_B2_OMIT_7_KEEP_F420 = "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"
BIO_B3_OMIT_LIPIDS = "B3_perturbation_omit_template_cofactors_and_surface_glycolipids"
BIO_B4_CORE_SKELETON = "B4_perturbation_core_amino_nucleotide_energy_skeleton"

NUTRIENT_SCENARIOS: dict[str, dict[str, Any]] = {
    SCEN_S1_ANOXIC_MINIMAL: {
        "redox_regime": "anoxic",
        "epistemic_status": "exploratory_chemical_scenario_not_measured_porewater",
        "compounds": {**BASE_INORGANIC_SALTS, "cpd11640": ("H2", 10.0)},
    },
    SCEN_S2_OXIC_MINIMAL: {
        "redox_regime": "oxic",
        "epistemic_status": "exploratory_chemical_scenario_not_measured_porewater",
        "compounds": {**BASE_INORGANIC_SALTS, "cpd11640": ("H2", 10.0), "cpd00007": ("O2", 10.0)},
    },
    SCEN_S3_ANOXIC_DESERT: {
        "redox_regime": "anoxic",
        "epistemic_status": "exploratory_chemical_scenario_not_measured_porewater",
        "compounds": {**BASE_INORGANIC_SALTS, **DESERT_OA_AA_ADDITIONS},
    },
    SCEN_S4_OXIC_DESERT: {
        "redox_regime": "oxic",
        "epistemic_status": "exploratory_chemical_scenario_not_measured_porewater",
        "compounds": {**BASE_INORGANIC_SALTS, **DESERT_OA_AA_ADDITIONS, "cpd00007": ("O2", 10.0)},
    },
    SCEN_S5_ANOXIC_DEFINED: {
        "redox_regime": "anoxic",
        "epistemic_status": "exploratory_chemical_scenario_not_measured_porewater",
        "compounds": {**BASE_INORGANIC_SALTS, **DEFINED_20AA_AND_PRECURSORS},
    },
    SCEN_S6_OXIC_DEFINED: {
        "redox_regime": "oxic",
        "epistemic_status": "exploratory_chemical_scenario_not_measured_porewater",
        "compounds": {**BASE_INORGANIC_SALTS, **DEFINED_20AA_AND_PRECURSORS, "cpd00007": ("O2", 10.0)},
    },
}

PLANNED_PROTOCOL_SCENARIO_BIOMASS_TUPLES_PER_GENOME: tuple[tuple[str, str], ...] = (
    (SCEN_S0_NATIVE, BIO_B0_BASELINE),
    (SCEN_S0_NATIVE, BIO_B2_OMIT_7_KEEP_F420),
    (SCEN_S1_ANOXIC_MINIMAL, BIO_B0_BASELINE),
    (SCEN_S1_ANOXIC_MINIMAL, BIO_B2_OMIT_7_KEEP_F420),
    (SCEN_S2_OXIC_MINIMAL, BIO_B0_BASELINE),
    (SCEN_S2_OXIC_MINIMAL, BIO_B2_OMIT_7_KEEP_F420),
    (SCEN_S3_ANOXIC_DESERT, BIO_B0_BASELINE),
    (SCEN_S3_ANOXIC_DESERT, BIO_B1_OMIT_8),
    (SCEN_S3_ANOXIC_DESERT, BIO_B2_OMIT_7_KEEP_F420),
    (SCEN_S3_ANOXIC_DESERT, BIO_B3_OMIT_LIPIDS),
    (SCEN_S3_ANOXIC_DESERT, BIO_B4_CORE_SKELETON),
    (SCEN_S4_OXIC_DESERT, BIO_B0_BASELINE),
    (SCEN_S4_OXIC_DESERT, BIO_B1_OMIT_8),
    (SCEN_S4_OXIC_DESERT, BIO_B2_OMIT_7_KEEP_F420),
    (SCEN_S4_OXIC_DESERT, BIO_B3_OMIT_LIPIDS),
    (SCEN_S4_OXIC_DESERT, BIO_B4_CORE_SKELETON),
    (SCEN_S5_ANOXIC_DEFINED, BIO_B0_BASELINE),
    (SCEN_S5_ANOXIC_DEFINED, BIO_B1_OMIT_8),
    (SCEN_S5_ANOXIC_DEFINED, BIO_B2_OMIT_7_KEEP_F420),
    (SCEN_S5_ANOXIC_DEFINED, BIO_B3_OMIT_LIPIDS),
    (SCEN_S5_ANOXIC_DEFINED, BIO_B4_CORE_SKELETON),
    (SCEN_S6_OXIC_DEFINED, BIO_B0_BASELINE),
    (SCEN_S6_OXIC_DEFINED, BIO_B1_OMIT_8),
    (SCEN_S6_OXIC_DEFINED, BIO_B2_OMIT_7_KEEP_F420),
    (SCEN_S6_OXIC_DEFINED, BIO_B3_OMIT_LIPIDS),
    (SCEN_S6_OXIC_DEFINED, BIO_B4_CORE_SKELETON),
)


def get_planned_protocol_scenario_biomass_tuples(
    genome_ids: Iterable[str],
) -> set[tuple[str, str, str]]:
    """Return the exact set of (genome, scenario_id, biomass_recipe_id) triples mandated by the
    planned protocol (`EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt` Item A).
    Each genome must evaluate exactly these 26 declared triples, preventing arbitrary or corrupted
    scenario names from passing gate checks by count alone.
    """
    tuples: set[tuple[str, str, str]] = set()
    for g in genome_ids:
        for s_id, b_id in PLANNED_PROTOCOL_SCENARIO_BIOMASS_TUPLES_PER_GENOME:
            tuples.add((str(g), s_id, b_id))
    return tuples


# =============================================================================
# 2. Data Structures for SBML Models and GPR Trees
# =============================================================================

@dataclass
class SpeciesRecord:
    id: str
    cpd_id: str
    name: str
    compartment: str
    charge: int | None
    formula: str | None


@dataclass
class ReactionRecord:
    id: str
    base_id: str
    name: str
    reversible: bool
    lb: float
    ub: float
    stoichiometry: dict[str, float]
    gpr_ast: Any | None = None
    genes: set[str] = field(default_factory=set)
    ec_numbers: list[str] = field(default_factory=list)
    is_exchange: bool = False
    is_demand: bool = False
    is_biomass: bool = False
    is_transport: bool = False
    provenance: str = "draft_sequence_supported"


@dataclass
class SBMLModelData:
    model_id: str
    model_name: str
    notes_text: str
    species: dict[str, SpeciesRecord]
    reactions: dict[str, ReactionRecord]
    objective_rxn_id: str
    gene_products: dict[str, str]

    def copy(self) -> SBMLModelData:
        return SBMLModelData(
            model_id=self.model_id,
            model_name=self.model_name,
            notes_text=self.notes_text,
            species={k: SpeciesRecord(**vars(v)) for k, v in self.species.items()},
            reactions={
                k: ReactionRecord(
                    id=v.id,
                    base_id=v.base_id,
                    name=v.name,
                    reversible=v.reversible,
                    lb=v.lb,
                    ub=v.ub,
                    stoichiometry=dict(v.stoichiometry),
                    gpr_ast=v.gpr_ast,
                    genes=set(v.genes),
                    ec_numbers=list(v.ec_numbers),
                    is_exchange=v.is_exchange,
                    is_demand=v.is_demand,
                    is_biomass=v.is_biomass,
                    is_transport=v.is_transport,
                    provenance=v.provenance,
                )
                for k, v in self.reactions.items()
            },
            objective_rxn_id=self.objective_rxn_id,
            gene_products=dict(self.gene_products),
        )


# =============================================================================
# 3. SBML Level 3 FBC v2 Parser & GPR Evaluator
# =============================================================================

_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
_VALID_FORMULA_RE = re.compile(r"^(?:[A-Z][a-z]?\d*)+$")


def parse_chemical_formula(formula_str: str | None) -> dict[str, int] | None:
    """Parse chemical formula into element -> atom count.
    Returns None if formula is missing, 'null', 'NA', empty, or contains generic R/X/polymer tokens.
    """
    if formula_str is None:
        return None
    s = formula_str.strip()
    if not s or s.lower() in {"null", "na", "none", "unknown", "noformula"}:
        return None
    if not _VALID_FORMULA_RE.match(s):
        return None
    counts: dict[str, int] = {}
    for elem, num_s in _FORMULA_TOKEN_RE.findall(s):
        if elem in {"R", "X", "Z"}:
            return None
        n = int(num_s) if num_s else 1
        counts[elem] = counts.get(elem, 0) + n
    return counts if counts else None


def _strip_tag(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


def _parse_gpr_element(elem: ET.Element) -> tuple[Any, set[str]]:
    tag = _strip_tag(elem.tag)
    if tag == "geneProductRef":
        gp = (
            elem.attrib.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct")
            or elem.attrib.get("fbc:geneProduct")
            or elem.attrib.get("geneProduct")
            or ""
        )
        return ("GENE", gp), ({gp} if gp else set())
    if tag in {"and", "or"}:
        children = []
        genes: set[str] = set()
        for ch in elem:
            ast_c, g_c = _parse_gpr_element(ch)
            if ast_c is not None:
                children.append(ast_c)
                genes |= g_c
        return (tag.upper(), children), genes
    for ch in elem:
        return _parse_gpr_element(ch)
    return None, set()


def evaluate_gpr_ast(ast: Any, active_genes: set[str]) -> bool:
    """Evaluate GPR boolean AST given the set of unmasked (active) gene product IDs."""
    if ast is None:
        return True
    kind, payload = ast
    if kind == "GENE":
        return payload in active_genes
    if kind == "AND":
        return all(evaluate_gpr_ast(c, active_genes) for c in payload)
    if kind == "OR":
        return any(evaluate_gpr_ast(c, active_genes) for c in payload)
    raise ValueError(f"Unknown GPR AST node kind: {kind}")


def is_single_sided_boundary_reaction(rxn: ReactionRecord) -> bool:
    """Return True if reaction has only reactants (all stoich < 0) or only products (all stoich > 0)
    or is explicitly flagged as exchange/demand.
    """
    if rxn.is_biomass:
        return False
    if rxn.is_exchange or rxn.is_demand:
        return True
    non_zero = [v for v in rxn.stoichiometry.values() if abs(v) > 1e-12]
    if not non_zero:
        return True
    has_neg = any(v < 0 for v in non_zero)
    has_pos = any(v > 0 for v in non_zero)
    return not (has_neg and has_pos)


def parse_gapseq_medium_csv(csv_text: str) -> dict[str, tuple[str, float]]:
    """Parse a gapseq medium CSV or TSV definition into a dictionary of
    `{cpd_id: (compound_name, max_uptake)}`.
    """
    media: dict[str, tuple[str, float]] = {}
    if not csv_text or not csv_text.strip():
        return media

    lines = [ln.strip() for ln in csv_text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        return media

    first = lines[0]
    delim = "\t" if "\t" in first else ("," if "," in first else ";")
    reader = csv.reader(lines, delimiter=delim)
    header = [h.strip().lower() for h in next(reader, [])]

    cpd_idx = -1
    name_idx = -1
    flux_idx = -1

    for idx, col in enumerate(header):
        if col in {"compounds", "compound", "id", "cpd", "cpd_id", "metabolite"}:
            cpd_idx = idx
        elif col in {"name", "compound_name", "metabolite_name"}:
            name_idx = idx
        elif col in {"maxflux", "max_flux", "flux", "uptake", "rate", "max_uptake"}:
            flux_idx = idx

    rows_to_process: list[list[str]] = []
    if cpd_idx == -1 and any(re.match(r"^cpd\d+", h, re.IGNORECASE) for h in header):
        rows_to_process.append(header)
        cpd_idx = 0
        if len(header) > 2 and re.match(r"^-?\d+(?:\.\d+)?$", header[2]):
            name_idx = 1
            flux_idx = 2
        elif len(header) > 1 and re.match(r"^-?\d+(?:\.\d+)?$", header[1]):
            flux_idx = 1
            name_idx = -1

    for row in reader:
        rows_to_process.append(row)

    for row in rows_to_process:
        if not row:
            continue
        cid = ""
        cname = ""
        flux = 10.0
        if cpd_idx != -1 and cpd_idx < len(row):
            cid = row[cpd_idx].strip()
        else:
            for c in row:
                m = re.search(r"\b(cpd\d+)\b", c)
                if m:
                    cid = m.group(1)
                    break
        if not cid.startswith("cpd"):
            continue

        if name_idx != -1 and name_idx < len(row):
            cname = row[name_idx].strip()
        if not cname:
            cname = cid

        if flux_idx != -1 and flux_idx < len(row):
            try:
                flux = abs(float(row[flux_idx].strip()))
            except ValueError:
                flux = 10.0
        media[cid] = (cname, flux)

    return media


def parse_gapseq_sbml_bytes(xml_bytes: bytes, provenance_default: str = "draft_sequence_supported") -> SBMLModelData:
    """Parse gapseq SBML Level 3 Version 2 + FBC v2 XML bytes into SBMLModelData."""
    root = ET.fromstring(xml_bytes)
    model_el = None
    for ch in root:
        if _strip_tag(ch.tag) == "model":
            model_el = ch
            break
    if model_el is None:
        raise ValueError("SBML root contains no <model> element")

    model_id = model_el.attrib.get("id", "unknown_model")
    model_name = model_el.attrib.get("name", model_id)
    notes_text = ""

    params: dict[str, float] = {"default_0": 0.0, "default_lb": -1000.0, "default_ub": 1000.0}
    species: dict[str, SpeciesRecord] = {}
    reactions: dict[str, ReactionRecord] = {}
    gene_products: dict[str, str] = {}
    objective_rxn_id = "R_bio1"

    fbc_ns = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}"

    for sec in model_el:
        sec_tag = _strip_tag(sec.tag)
        if sec_tag == "notes":
            notes_text = "".join(sec.itertext()).strip()
        elif sec_tag == "listOfParameters":
            for p in sec:
                if _strip_tag(p.tag) == "parameter":
                    pid = p.attrib["id"]
                    params[pid] = float(p.attrib["value"])
        elif sec_tag == "listOfSpecies":
            for sp in sec:
                if _strip_tag(sp.tag) != "species":
                    continue
                sid = sp.attrib["id"]
                sname = sp.attrib.get("name", sid)
                comp = sp.attrib.get("compartment", "c0")
                ch_raw = sp.attrib.get(f"{fbc_ns}charge") or sp.attrib.get("fbc:charge")
                charge: int | None = None
                if ch_raw is not None and ch_raw.strip() not in {"", "NA", "null"}:
                    try:
                        charge = int(float(ch_raw))
                    except ValueError:
                        charge = None
                form_raw = sp.attrib.get(f"{fbc_ns}chemicalFormula") or sp.attrib.get("fbc:chemicalFormula")
                m_cpd = re.match(r"^M_(cpd\d+)_([cep]0)$", sid)
                cpd_id = m_cpd.group(1) if m_cpd else sid
                species[sid] = SpeciesRecord(
                    id=sid,
                    cpd_id=cpd_id,
                    name=sname,
                    compartment=comp,
                    charge=charge,
                    formula=form_raw.strip() if form_raw else None,
                )
        elif sec_tag == "listOfObjectives":
            for obj in sec:
                for lfo in obj:
                    for fo in lfo:
                        if _strip_tag(fo.tag) == "fluxObjective":
                            rxn_ref = (
                                fo.attrib.get(f"{fbc_ns}reaction")
                                or fo.attrib.get("fbc:reaction")
                                or fo.attrib.get("reaction")
                            )
                            if rxn_ref:
                                objective_rxn_id = rxn_ref
        elif sec_tag == "listOfGeneProducts":
            for gp in sec:
                if _strip_tag(gp.tag) == "geneProduct":
                    gid = gp.attrib.get(f"{fbc_ns}id") or gp.attrib.get("fbc:id") or gp.attrib.get("id") or ""
                    glabel = gp.attrib.get(f"{fbc_ns}label") or gp.attrib.get("fbc:label") or gid
                    if gid:
                        gene_products[gid] = glabel
        elif sec_tag == "listOfReactions":
            for rxn in sec:
                if _strip_tag(rxn.tag) != "reaction":
                    continue
                rid = rxn.attrib["id"]
                rname = rxn.attrib.get("name", rid)
                rev = rxn.attrib.get("reversible", "false").lower() == "true"
                lb_ref = (
                    rxn.attrib.get(f"{fbc_ns}lowerFluxBound")
                    or rxn.attrib.get("fbc:lowerFluxBound")
                    or "default_lb"
                )
                ub_ref = (
                    rxn.attrib.get(f"{fbc_ns}upperFluxBound")
                    or rxn.attrib.get("fbc:upperFluxBound")
                    or "default_ub"
                )
                lb = float(params.get(lb_ref, -1000.0))
                ub = float(params.get(ub_ref, 1000.0))
                stoich: dict[str, float] = {}
                gpr_ast = None
                genes: set[str] = set()
                ecs: list[str] = []

                for sub in rxn:
                    stag = _strip_tag(sub.tag)
                    if stag == "listOfReactants":
                        for sr in sub:
                            if _strip_tag(sr.tag) == "speciesReference":
                                sp_id = sr.attrib["species"]
                                st_val = float(sr.attrib.get("stoichiometry", "1"))
                                stoich[sp_id] = stoich.get(sp_id, 0.0) - st_val
                    elif stag == "listOfProducts":
                        for sr in sub:
                            if _strip_tag(sr.tag) == "speciesReference":
                                sp_id = sr.attrib["species"]
                                st_val = float(sr.attrib.get("stoichiometry", "1"))
                                stoich[sp_id] = stoich.get(sp_id, 0.0) + st_val
                    elif stag == "geneProductAssociation":
                        gpr_ast, genes = _parse_gpr_element(sub)
                    elif stag == "annotation":
                        for elem in sub.iter():
                            res_attr = elem.attrib.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource", "")
                            if "ec-code:" in res_attr:
                                ecs.append(res_attr.split("ec-code:", 1)[-1])

                base_id = rid[2:] if rid.startswith("R_") else rid
                base_id = re.sub(r"_[cep]0$", "", base_id)
                is_ex = rid.startswith("R_EX_") or base_id.startswith("EX_")
                is_dm = rid.startswith("R_DM_") or base_id.startswith("DM_") or rid.startswith("R_SINK_")
                is_bio = rid == "R_bio1" or base_id == "bio1" or rid == objective_rxn_id
                comps_touched = {species[s].compartment for s in stoich if s in species}
                is_trans = (len(comps_touched) > 1) and not is_ex and not is_bio

                reactions[rid] = ReactionRecord(
                    id=rid,
                    base_id=base_id,
                    name=rname,
                    reversible=rev,
                    lb=lb,
                    ub=ub,
                    stoichiometry=stoich,
                    gpr_ast=gpr_ast,
                    genes=genes,
                    ec_numbers=ecs,
                    is_exchange=is_ex,
                    is_demand=is_dm,
                    is_biomass=is_bio,
                    is_transport=is_trans,
                    provenance=provenance_default,
                )

    return SBMLModelData(
        model_id=model_id,
        model_name=model_name,
        notes_text=notes_text,
        species=species,
        reactions=reactions,
        objective_rxn_id=objective_rxn_id,
        gene_products=gene_products,
    )


def parse_gapseq_sbml_file(path: Path, provenance_default: str = "draft_sequence_supported") -> SBMLModelData:
    return parse_gapseq_sbml_bytes(path.read_bytes(), provenance_default=provenance_default)


# =============================================================================
# 4. Reaction Mass/Charge Balance, Biomass Audit, and Closed-Boundary EGC Check
# =============================================================================

def audit_reaction_mass_charge_balance(rxn: ReactionRecord, species_map: dict[str, SpeciesRecord]) -> dict[str, Any]:
    """Audit elemental mass and charge balance of a single reaction.
    CRITICAL RULE (`EXP1-CONTRACT.txt` Item 2): Unknown formula or charge is NEVER treated as balanced!
    """
    if rxn.is_biomass or is_single_sided_boundary_reaction(rxn):
        return {
            "rxn_id": rxn.id,
            "balance_status": "boundary_or_biomass_exempt",
            "mass_balanced": None,
            "charge_balanced": None,
            "net_charge": None,
        }

    atom_net: dict[str, float] = {}
    net_charge = 0.0
    formula_known_all = True
    charge_known_all = True

    for sp_id, coeff in rxn.stoichiometry.items():
        if abs(coeff) < 1e-12:
            continue
        sp = species_map.get(sp_id)
        if sp is None:
            formula_known_all = False
            charge_known_all = False
            break
        parsed_f = parse_chemical_formula(sp.formula)
        if parsed_f is None:
            formula_known_all = False
        else:
            for elem, cnt in parsed_f.items():
                atom_net[elem] = atom_net.get(elem, 0.0) + coeff * cnt
        if sp.charge is None:
            charge_known_all = False
        else:
            net_charge += coeff * float(sp.charge)

    if not formula_known_all or not charge_known_all:
        return {
            "rxn_id": rxn.id,
            "balance_status": "unknown_formula_or_charge",
            "mass_balanced": None,
            "charge_balanced": None,
            "net_charge": None,
        }

    heavy_imbalance = any(abs(v) > 1e-5 for k, v in atom_net.items() if k != "H")
    h_imbalance = abs(atom_net.get("H", 0.0)) > 1e-5
    q_imbalance = abs(net_charge) > 1e-5

    if not heavy_imbalance and not h_imbalance and not q_imbalance:
        status = "strict_atom_and_charge_balanced"
    elif not heavy_imbalance and (h_imbalance or q_imbalance):
        status = "heavy_atom_balanced_proton_or_charge_offset"
    else:
        status = "elemental_imbalance"

    return {
        "rxn_id": rxn.id,
        "balance_status": status,
        "mass_balanced": not heavy_imbalance and not h_imbalance,
        "charge_balanced": not q_imbalance,
        "net_charge": float(net_charge),
    }


def audit_model_mass_charge_balance(model: SBMLModelData) -> dict[str, int]:
    """Summarize reaction mass and charge balance across all internal/transport reactions of `model`,
    strictly separating `unknown_formula_or_charge` from balanced and imbalanced counts.
    """
    n_total = 0
    n_strict_bal = 0
    n_heavy_bal_h_offset = 0
    n_unknown = 0
    n_imbalanced = 0
    for rxn in model.reactions.values():
        res = audit_reaction_mass_charge_balance(rxn, model.species)
        st = res["balance_status"]
        if st == "boundary_or_biomass_exempt":
            continue
        n_total += 1
        if st == "strict_atom_and_charge_balanced":
            n_strict_bal += 1
        elif st == "heavy_atom_balanced_proton_or_charge_offset":
            n_heavy_bal_h_offset += 1
        elif st == "unknown_formula_or_charge":
            n_unknown += 1
        else:
            n_imbalanced += 1
    return {
        "n_internal_total": n_total,
        "n_internal_strict_balanced": n_strict_bal,
        "n_internal_heavy_atom_balanced_h_offset": n_heavy_bal_h_offset,
        "n_internal_balanced": n_strict_bal + n_heavy_bal_h_offset,
        "n_internal_unknown_formula_or_charge": n_unknown,
        "n_internal_imbalanced": n_imbalanced,
    }


def validate_archaeal_biomass_and_audit_template(
    model: SBMLModelData,
    rxn_evidence_map: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Validate that the biomass objective is archaeal (never bacterial default) and audit
    Methanosarcina barkeri template demands vs sequence-supported capacities (`EXP1-BIOMASS-SCIENTIFIC-NOTE.txt`).
    """
    if model.objective_rxn_id not in model.reactions:
        raise ValueError(f"Model {model.model_id} missing objective reaction {model.objective_rxn_id}")
    bio_rxn = model.reactions[model.objective_rxn_id]
    rname_low = bio_rxn.name.lower()
    if "bacteria" in rname_low or "gram" in rname_low or "murein" in rname_low:
        raise ValueError(f"Invalid bacterial default biomass objective in {model.model_id}: {bio_rxn.name!r}")

    reactant_cpds: dict[str, float] = {}
    product_cpds: dict[str, float] = {}
    for sp_id, coeff in bio_rxn.stoichiometry.items():
        sp = model.species.get(sp_id)
        cpd = sp.cpd_id if sp else re.sub(r"^M_|_[cep]0$", "", sp_id)
        if coeff < 0:
            reactant_cpds[cpd] = abs(float(coeff))
        elif coeff > 0:
            product_cpds[cpd] = float(coeff)

    if not reactant_cpds:
        raise ValueError(f"Biomass reaction {bio_rxn.id} in {model.model_id} has zero reactants")

    forbidden_found = sorted(set(reactant_cpds.keys()) & set(BACTERIAL_BIOMASS_FORBIDDEN_CPDS.keys()))
    if forbidden_found:
        raise ValueError(
            f"Invalid bacterial cell-wall/peptidoglycan components {forbidden_found} found in biomass of {model.model_id}"
        )

    all8_present = sorted(set(reactant_cpds.keys()) & set(METHANOSARCINA_ALL8_TEMPLATE_COFACTORS.keys()))
    non_f420_7_present = sorted(set(reactant_cpds.keys()) & set(METHANOSARCINA_NON_F420_TEMPLATE_COFACTORS.keys()))
    f420_present = F420_COFACTOR_ID in reactant_cpds
    methanofuran_present = sorted(set(reactant_cpds.keys()) & set(METHANOFURAN_COMPOUND_IDS.keys()))
    ether_lipids_present = sorted(set(reactant_cpds.keys()) & set(ARCHAEAL_ETHER_LIPIDS_AND_POLYAMINES.keys()))

    # Check sequence-supported production or utilization in draft model
    draft_produced_cpds: set[str] = set()
    draft_utilized_cpds: set[str] = set()
    for r in model.reactions.values():
        if r.is_biomass or is_single_sided_boundary_reaction(r):
            continue
        for sp_id, coeff in r.stoichiometry.items():
            sp = model.species.get(sp_id)
            cpd = sp.cpd_id if sp else re.sub(r"^M_|_[cep]0$", "", sp_id)
            if coeff > 0 or (r.reversible and coeff != 0):
                draft_produced_cpds.add(cpd)
            if coeff < 0 or (r.reversible and coeff != 0):
                draft_utilized_cpds.add(cpd)

    supported_in_draft = sorted(c for c in all8_present if c in draft_produced_cpds)
    unsupported_in_draft = sorted(c for c in all8_present if c not in draft_produced_cpds)

    return {
        "model_id": model.model_id,
        "biomass_rxn_id": bio_rxn.id,
        "biomass_rxn_name": bio_rxn.name,
        "biomass_provenance": METHANOSARCINA_BIOMASS_PROVENANCE,
        "n_biomass_reactants": len(reactant_cpds),
        "n_biomass_products": len(product_cpds),
        "gam_atp_stoichiometry": float(reactant_cpds.get("cpd00002", 0.0)),
        "f420_demanded_in_template": f420_present,
        "f420_sequence_supported_synthesis_in_draft": F420_COFACTOR_ID in draft_produced_cpds,
        "f420_sequence_supported_use_in_draft": F420_COFACTOR_ID in draft_utilized_cpds,
        "n_non_f420_methanosarcina_cofactors_demanded": len(non_f420_7_present),
        "n_all8_template_trace_cofactors_demanded": len(all8_present),
        "template_trace_cofactors_demanded_csv": ",".join(all8_present),
        "methanofuran_demanded_in_template": len(methanofuran_present) > 0,
        "n_template_trace_cofactors_supported_in_draft": len(supported_in_draft),
        "n_template_trace_cofactors_unsupported_in_draft": len(unsupported_in_draft),
        "template_trace_cofactors_unsupported_in_draft_csv": ",".join(unsupported_in_draft),
        "n_archaeal_ether_lipids_and_polyamines_demanded": len(ether_lipids_present),
    }


def apply_biomass_perturbation(model: SBMLModelData, recipe_id: str) -> SBMLModelData:
    """Return a copy of `model` with the specified declared component-omission biomass perturbation applied."""
    if recipe_id not in BIOMASS_PERTURBATION_RECIPES:
        raise ValueError(f"Unknown biomass recipe_id: {recipe_id}")
    validate_archaeal_biomass_and_audit_template(model)
    mod_copy = model.copy()
    omit_cpds: set[str] = BIOMASS_PERTURBATION_RECIPES[recipe_id]["omit_cpd_ids"]
    if not omit_cpds:
        return mod_copy

    bio_rxn = mod_copy.reactions[mod_copy.objective_rxn_id]
    new_stoich: dict[str, float] = {}
    for sp_id, coeff in bio_rxn.stoichiometry.items():
        sp = mod_copy.species.get(sp_id)
        cpd = sp.cpd_id if sp else re.sub(r"^M_|_[cep]0$", "", sp_id)
        if coeff < 0 and cpd in omit_cpds:
            continue
        new_stoich[sp_id] = coeff

    if not any(v < 0 for v in new_stoich.values()):
        raise ValueError(f"Biomass perturbation {recipe_id} removed all reactants from {model.model_id}")
    bio_rxn.stoichiometry = new_stoich
    return mod_copy


def classify_biomass_requirement_provenance(
    req_in_b0_baseline: bool,
    req_in_b1_omit_all8: bool,
    req_in_b2_omit_7_keep_f420: bool,
    req_in_b3_omit_surface_lipids: bool,
    req_in_b4_core_skeleton: bool,
) -> str:
    """Classify dependence of a gap-fill or nutrient requirement on biomass formulation (`EXP1-BIOMASS-SCIENTIFIC-NOTE.txt` & `EXP1-SOLVER-REVIEW.txt` Item 4).
    CRITICAL RULE: If a requirement disappears when a biomass component is omitted, label it
    `biomass_assumption_sensitive` (with sub-detail), NEVER `proven_artifact` or `proven_non_auxotrophy`!
    """
    if req_in_b0_baseline and req_in_b1_omit_all8 and req_in_b2_omit_7_keep_f420 and req_in_b3_omit_surface_lipids and req_in_b4_core_skeleton:
        return "robust_across_all_tested_biomass_perturbations"
    if req_in_b0_baseline and not (req_in_b1_omit_all8 and req_in_b2_omit_7_keep_f420 and req_in_b3_omit_surface_lipids and req_in_b4_core_skeleton):
        if not req_in_b2_omit_7_keep_f420:
            return "biomass_assumption_sensitive__vanishes_when_7_non_f420_template_cofactors_omitted"
        if req_in_b2_omit_7_keep_f420 and not req_in_b1_omit_all8:
            return "biomass_assumption_sensitive__vanishes_when_f420_omitted"
        if req_in_b1_omit_all8 and not req_in_b3_omit_surface_lipids:
            return "biomass_assumption_sensitive__vanishes_when_specialized_lipids_or_homospermidine_omitted"
        return "biomass_assumption_sensitive__vanishes_in_core_skeleton_perturbation"
    if not req_in_b0_baseline and any([req_in_b1_omit_all8, req_in_b2_omit_7_keep_f420, req_in_b3_omit_surface_lipids, req_in_b4_core_skeleton]):
        return "biomass_assumption_sensitive__appears_only_in_perturbed_biomass"
    return "not_required"


# =============================================================================
# 5. HiGHS Linear Programming & Closed-Boundary EGC Check
# =============================================================================

def build_sparse_stoichiometric_matrix(
    model: SBMLModelData,
) -> tuple[list[str], list[str], sparse.csc_matrix, np.ndarray, np.ndarray, int]:
    """Assemble sparse stoichiometric matrix S (n_species x n_reactions) and bound arrays."""
    sp_ids = sorted(model.species.keys())
    rxn_ids = sorted(model.reactions.keys())
    sp_idx = {s: i for i, s in enumerate(sp_ids)}
    rxn_idx = {r: j for j, r in enumerate(rxn_ids)}

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    lb = np.zeros(len(rxn_ids), dtype=np.float64)
    ub = np.zeros(len(rxn_ids), dtype=np.float64)

    for j, rid in enumerate(rxn_ids):
        rxn = model.reactions[rid]
        if rxn.lb > rxn.ub + 1e-12:
            raise ValueError(f"Invalid flux bounds on {rid}: lb={rxn.lb} > ub={rxn.ub}")
        lb[j] = rxn.lb
        ub[j] = rxn.ub
        for sid, coeff in rxn.stoichiometry.items():
            if abs(coeff) < 1e-14:
                continue
            if sid not in sp_idx:
                sp_idx[sid] = len(sp_ids)
                sp_ids.append(sid)
            rows.append(sp_idx[sid])
            cols.append(j)
            vals.append(float(coeff))

    S = sparse.csc_matrix((vals, (rows, cols)), shape=(len(sp_ids), len(rxn_ids)), dtype=np.float64)
    obj_j = rxn_idx[model.objective_rxn_id]
    return sp_ids, rxn_ids, S, lb, ub, obj_j


def close_all_external_inflows_directionally(model: SBMLModelData) -> SBMLModelData:
    r"""Close ALL external inflow across EVERY single-sided boundary reaction (exchange, demand, sink, or custom source),
    regardless of reaction naming (`EXP1-SOLVER-REVIEW.txt` Item 5).
    - If a boundary reaction has negative stoichiometry ($A \to \emptyset$), inflow is $v < 0$, so set `lb = max(0.0, lb)`.
    - If a boundary reaction has positive stoichiometry ($\emptyset \to A$), inflow is $v > 0$, so set `ub = min(0.0, ub)`.
    """
    mod = model.copy()
    for rid, rxn in mod.reactions.items():
        if rxn.is_biomass:
            continue
        if is_single_sided_boundary_reaction(rxn):
            non_zero = [v for v in rxn.stoichiometry.values() if abs(v) > 1e-12]
            if not non_zero:
                rxn.lb = 0.0
                rxn.ub = 0.0
            elif all(v < 0 for v in non_zero):
                # Standard sink/exchange convention: M -> ø. Close uptake (v < 0), allow secretion (v >= 0).
                rxn.lb = max(0.0, rxn.lb)
                rxn.ub = max(rxn.ub, 0.0)
            elif all(v > 0 for v in non_zero):
                # Source convention: ø -> M. Close inflow (v > 0), allow outflow (v <= 0).
                rxn.ub = min(0.0, rxn.ub)
                rxn.lb = min(rxn.lb, 0.0)
    return mod


def apply_nutrient_scenario_bounds(
    model: SBMLModelData,
    compounds_max_uptake: dict[str, tuple[str, float]],
) -> SBMLModelData:
    """Apply chemically explicit uptake bounds to exchange reactions (`R_EX_cpd*_e0`).
    First directionally closes all single-sided boundary inflows, then opens ONLY the exchanges
    specified in `compounds_max_uptake`. Validates that exchange stoichiometry is `-1` (`M_e0 -> ø`).
    """
    mod = model.copy()
    existing_ex_by_cpd: dict[str, str] = {}
    for rid, rxn in mod.reactions.items():
        if rxn.is_exchange:
            if len(rxn.stoichiometry) != 1:
                raise ValueError(f"Exchange reaction {rid} must have 1 species, found {len(rxn.stoichiometry)}")
            sp_id, coeff = next(iter(rxn.stoichiometry.items()))
            if coeff >= 0:
                raise ValueError(
                    f"Exchange reaction {rid} has positive stoichiometry ({coeff}) on {sp_id}; "
                    "expected -1 (M_e0 -> ø) so negative flux = uptake."
                )
            m = re.match(r"^R_EX_(cpd\d+)_e0$", rid) or re.match(r"^EX_(cpd\d+)_e0$", rid)
            if m:
                existing_ex_by_cpd[m.group(1)] = rid
            else:
                sp = mod.species.get(sp_id)
                if sp:
                    existing_ex_by_cpd[sp.cpd_id] = rid

    # Directionally close all boundary inflows first
    mod = close_all_external_inflows_directionally(mod)
    for rid, rxn in mod.reactions.items():
        if rxn.is_exchange:
            rxn.lb = 0.0
            rxn.ub = max(rxn.ub, 1000.0)

    # Open only supplied scenario compounds
    for cpd_id, (_cname, max_inflow) in compounds_max_uptake.items():
        if max_inflow < 0:
            raise ValueError(f"Medium max_inflow must be >= 0 for {cpd_id}, got {max_inflow}")
        ex_rid = existing_ex_by_cpd.get(cpd_id)
        if ex_rid is not None:
            mod.reactions[ex_rid].lb = -float(max_inflow)
            mod.reactions[ex_rid].ub = 1000.0
        else:
            sp_e0 = f"M_{cpd_id}_e0"
            if sp_e0 in mod.species:
                new_rid = f"R_EX_{cpd_id}_e0"
                mod.reactions[new_rid] = ReactionRecord(
                    id=new_rid,
                    base_id=f"EX_{cpd_id}_e0",
                    name=f"{_cname} Exchange",
                    reversible=True,
                    lb=-float(max_inflow),
                    ub=1000.0,
                    stoichiometry={sp_e0: -1.0},
                    is_exchange=True,
                    provenance="scenario_exchange_boundary",
                )
    return mod


def solve_fba_highs(
    model: SBMLModelData,
    time_limit_s: float = 15.0,
    feasibility_tol: float = 1e-9,
    residual_tol: float = 1e-5,
) -> dict[str, Any]:
    """Solve standard FBA (maximize objective reaction) using HiGHS LP (`scipy.optimize.linprog`).
    CRITICAL RULE (`EXP1-SOLVER-REVIEW.txt` & `EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt` Item 2):
    - Passes tightened primal feasibility tolerance (1e-9).
    - Independently checks unthresholded solution vector x for finite values and NaN/Inf growth.
    - Independently verifies stoichiometric mass balance max|Sx - b| <= residual_tol (1e-5).
    - Independently verifies bound compliance max(lb - x, x - ub) <= residual_tol (1e-5).
    - Solver failure or residual violation is NEVER reported as optimal or 0.0 growth; `growth_rate` is `math.nan`.
    - Preserves unthresholded diagnostics before flux sparsification.
    """
    sp_ids, rxn_ids, S, lb, ub, obj_j = build_sparse_stoichiometric_matrix(model)
    c = np.zeros(len(rxn_ids), dtype=np.float64)
    c[obj_j] = -1.0  # linprog minimizes c^T v
    b_eq = np.zeros(len(sp_ids), dtype=np.float64)

    res = linprog(
        c=c,
        A_eq=S,
        b_eq=b_eq,
        bounds=list(zip(lb, ub)),
        method="highs",
        options={
            "time_limit": float(time_limit_s),
            "presolve": True,
            "primal_feasibility_tolerance": float(feasibility_tol),
        },
    )
    raw_st = int(getattr(res, "status", -1))
    if not res.success or res.x is None:
        if raw_st == 2:
            return {
                "status": "infeasible",
                "raw_status": raw_st,
                "message": str(getattr(res, "message", "")),
                "growth_rate": 0.0,
                "fluxes": {},
                "max_residual": 0.0,
                "max_bound_viol": 0.0,
                "feasibility_tol": feasibility_tol,
                "residual_tol": residual_tol,
            }
        status_map = {
            1: "timeout",
            3: "unbounded",
            4: "numerical_error",
        }
        return {
            "status": status_map.get(raw_st, "solver_failed"),
            "raw_status": raw_st,
            "message": str(getattr(res, "message", "")),
            "growth_rate": math.nan,
            "fluxes": {},
            "max_residual": math.nan,
            "max_bound_viol": math.nan,
            "feasibility_tol": feasibility_tol,
            "residual_tol": residual_tol,
        }

    x_arr = np.asarray(res.x, dtype=np.float64)

    # 1. Finite checking on unthresholded x
    if not np.all(np.isfinite(x_arr)):
        return {
            "status": "numerical_error",
            "raw_status": raw_st,
            "message": "Solution vector contains non-finite values (NaN/Inf)",
            "growth_rate": math.nan,
            "fluxes": {},
            "max_residual": math.nan,
            "max_bound_viol": math.nan,
            "feasibility_tol": feasibility_tol,
            "residual_tol": residual_tol,
        }

    growth = float(x_arr[obj_j])
    if math.isnan(growth) or math.isinf(growth):
        return {
            "status": "numerical_error",
            "raw_status": raw_st,
            "message": "Objective growth rate is non-finite (NaN/Inf)",
            "growth_rate": math.nan,
            "fluxes": {},
            "max_residual": math.nan,
            "max_bound_viol": math.nan,
            "feasibility_tol": feasibility_tol,
            "residual_tol": residual_tol,
        }

    # 2. Independently compute full unthresholded x residual max|Sx| and bound violations
    residuals = S.dot(x_arr) - b_eq
    max_res = float(np.max(np.abs(residuals))) if len(residuals) > 0 else 0.0

    lb_viol = np.maximum(0.0, lb - x_arr)
    ub_viol = np.maximum(0.0, x_arr - ub)
    max_bound_viol = float(max(np.max(lb_viol) if len(lb_viol) > 0 else 0.0, np.max(ub_viol) if len(ub_viol) > 0 else 0.0))

    if max_res > residual_tol or max_bound_viol > residual_tol:
        return {
            "status": "numerical_coupling_violation",
            "raw_status": raw_st,
            "message": f"Independent residual validation failed: max|Sx|={max_res:.2e}, max_bound_viol={max_bound_viol:.2e} (tol={residual_tol:.2e})",
            "growth_rate": math.nan,
            "fluxes": {},
            "max_residual": max_res,
            "max_bound_viol": max_bound_viol,
            "feasibility_tol": feasibility_tol,
            "residual_tol": residual_tol,
        }

    flux_map = {rid: float(x_arr[j]) for j, rid in enumerate(rxn_ids) if abs(x_arr[j]) > 1e-9}
    return {
        "status": "optimal",
        "raw_status": raw_st,
        "message": str(res.message),
        "growth_rate": growth,
        "fluxes": flux_map,
        "max_residual": max_res,
        "max_bound_viol": max_bound_viol,
        "feasibility_tol": feasibility_tol,
        "residual_tol": residual_tol,
    }


def detect_closed_uptake_egc(model: SBMLModelData, time_limit_s: float = 10.0) -> dict[str, Any]:
    """Check whether the model can generate ATP or redox potential when ALL external boundary inflows
    (including non-`EX_` sources/sinks) are directionally closed (`EXP1-SOLVER-REVIEW.txt` Item 5).
    """
    mod = close_all_external_inflows_directionally(model)
    for rid, rxn in mod.reactions.items():
        if rxn.is_biomass:
            rxn.lb = 0.0
            rxn.ub = 0.0

    atp_sp = {
        "M_cpd00002_c0": -1.0,
        "M_cpd00001_c0": -1.0,
        "M_cpd00008_c0": 1.0,
        "M_cpd00009_c0": 1.0,
        "M_cpd00067_c0": 1.0,
    }
    if not all(s in mod.species for s in atp_sp):
        return {
            "egc_check_status": "evaluated_no_atp_species",
            "has_spurious_atp_cycle": False,
            "closed_uptake_atp_dissipation_flux": 0.0,
        }

    mod_atp = mod.copy()
    mod_atp.reactions["R_EGC_ATP_TEST"] = ReactionRecord(
        id="R_EGC_ATP_TEST",
        base_id="EGC_ATP_TEST",
        name="Closed-uptake ATP dissipation test",
        reversible=False,
        lb=0.0,
        ub=1000.0,
        stoichiometry=atp_sp,
        provenance="egc_test",
    )
    mod_atp.objective_rxn_id = "R_EGC_ATP_TEST"
    sol_atp = solve_fba_highs(mod_atp, time_limit_s=time_limit_s)
    if sol_atp["status"] != "optimal":
        if sol_atp["status"] == "infeasible":
            return {
                "egc_check_status": "evaluated_infeasible_zero_flux",
                "has_spurious_atp_cycle": False,
                "closed_uptake_atp_dissipation_flux": 0.0,
            }
        return {
            "egc_check_status": f"unknown_{sol_atp['status']}",
            "has_spurious_atp_cycle": None,
            "closed_uptake_atp_dissipation_flux": math.nan,
        }

    atp_flux = float(sol_atp["growth_rate"])
    return {
        "egc_check_status": "evaluated_optimal",
        "has_spurious_atp_cycle": atp_flux > 1e-5,
        "closed_uptake_atp_dissipation_flux": atp_flux,
    }


# =============================================================================
# 6. Exact Binary Reaction-Cost MILP (`scipy.optimize.milp`) & No-Good Cuts
# =============================================================================

def parse_genome_reaction_and_transporter_evidence(
    reactions_tbl_text: str,
    transporter_tbl_text: str,
    high_evi_bs: float = 200.0,
    min_core_bs: float = 50.0,
    dummy_weight: float = 100.0,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    """Parse `<genome>-all-Reactions.tbl` and `<genome>-Transporter.tbl` into per-reaction and per-exchange
    evidence records and gapseq-compatible base costs (`w in [0.005, 100.0]`).
    """
    rxn_evi: dict[str, dict[str, Any]] = {}
    lines_r = [ln for ln in reactions_tbl_text.splitlines() if ln.strip() and not ln.startswith("#")]
    if lines_r:
        reader = csv.DictReader(lines_r, delimiter="\t")
        for row in reader:
            dbhits = (row.get("dbhit") or "").strip().split()
            status = (row.get("status") or "").strip()
            pwy_status = (row.get("pathway.status") or "").strip()
            spont = (row.get("spont") or "").strip().upper() == "TRUE"
            bs_raw = (row.get("bitscore") or "").strip()
            try:
                bs = float(bs_raw) if bs_raw else 0.0
            except ValueError:
                bs = 0.0
            gene_hit = ((row.get("stitle") or "").strip().split() or [""])[0]
            ec = (row.get("ec") or "").strip()
            pwy = (row.get("pathway") or "").strip()

            eff_bs = bs
            if spont or (status in {"bad_blast", "no_blast"} and pwy_status in {"full", "threshold", "keyenzyme"}):
                eff_bs = max(eff_bs, high_evi_bs)
            bs_clamped = min(eff_bs, high_evi_bs)
            if eff_bs > 0:
                w = (bs_clamped - high_evi_bs) * ((0.005 - dummy_weight) / (high_evi_bs - min_core_bs)) + 0.005
                w = float(min(max(w, 0.005), dummy_weight))
            else:
                w = float(dummy_weight)

            for r_id in dbhits:
                r_id = r_id.strip()
                if not r_id.startswith("rxn"):
                    continue
                prev = rxn_evi.get(r_id)
                if prev is None or eff_bs > prev["effective_bitscore"] or (status == "good_blast" and prev["status"] != "good_blast"):
                    rxn_evi[r_id] = {
                        "rxn_id": r_id,
                        "status": status if status else ("spontaneous" if spont else "no_blast"),
                        "bitscore": bs,
                        "effective_bitscore": eff_bs,
                        "weight": w,
                        "best_gene": gene_hit,
                        "ec": ec,
                        "pathway": pwy,
                        "pathway_status": pwy_status,
                        "is_core_candidate": eff_bs >= min_core_bs,
                    }

    trans_evi: dict[str, dict[str, Any]] = {}
    lines_t = [ln for ln in transporter_tbl_text.splitlines() if ln.strip() and not ln.startswith("#")]
    if lines_t:
        reader_t = csv.DictReader(lines_t, delimiter="\t")
        for row in reader_t:
            exid = (row.get("exid") or "").strip()
            reas = [x.strip() for x in (row.get("rea") or "").split(",") if x.strip()]
            bs_raw = (row.get("bitscore") or "").strip()
            try:
                bs = float(bs_raw) if bs_raw else 0.0
            except ValueError:
                bs = 0.0
            gene_hit = ((row.get("stitle") or "").strip().split() or [""])[0]
            tc = (row.get("tc") or "").strip()
            sub = (row.get("sub_gapseq") or row.get("sub") or "").strip()

            if exid:
                prev_ex = trans_evi.get(exid)
                if prev_ex is None or bs > prev_ex["bitscore"]:
                    trans_evi[exid] = {
                        "target_id": exid,
                        "kind": "exchange",
                        "bitscore": bs,
                        "best_gene": gene_hit,
                        "tc": tc,
                        "substrate": sub,
                        "weight": 0.05 if bs >= min_core_bs else 25.0,
                    }
            for r_id in reas:
                if not r_id.startswith("rxn"):
                    continue
                prev_r = trans_evi.get(r_id)
                if prev_r is None or bs > prev_r["bitscore"]:
                    trans_evi[r_id] = {
                        "target_id": r_id,
                        "kind": "transport_rxn",
                        "bitscore": bs,
                        "best_gene": gene_hit,
                        "tc": tc,
                        "substrate": sub,
                        "weight": 0.05 if bs >= min_core_bs else 25.0,
                    }

    return rxn_evi, trans_evi


def _classify_ensemble_frequency_label(freq: float) -> str:
    """Classify ensemble inclusion frequency using exact non-overclaiming terminology (`EXP1-SOLVER-REVIEW.txt` Item 3)."""
    if freq >= 0.999999:
        return "observed_in_all_enumerated_models"
    if freq >= 0.80:
        return "high_ensemble_stability"
    if freq >= 0.20:
        return "intermediate_ensemble_stability"
    return "low_ensemble_stability"


def build_fixed_medium_supernetwork_and_base(
    mod_base: SBMLModelData,
    candidate_reactions: dict[str, ReactionRecord],
    scenario_compounds: dict[str, tuple[str, float]],
) -> tuple[SBMLModelData, SBMLModelData, list[str]]:
    """Construct both:
    1) `fixed_medium_base`: `mod_base` augmented with ALL permitted fixed-medium boundary exchanges and extracellular
       species introduced by `candidate_reactions` (with supplied scenario compounds opened to `-max_inflow` and
       unsupplied compounds strictly closed at `lb = 0.0`), but NO unselected internal/transport candidate reactions.
    2) `mod_super`: `fixed_medium_base` plus all eligible non-zero-bound internal/transport candidate reactions.
    3) `cand_rids`: sorted list of eligible candidate reaction IDs in `mod_super`.

    Fixes `EXP1-FINAL-ENGINE-REVIEW.txt` Item B & C:
    - Verification models (`mod_check`, `mod_pruned`) and downstream withdrawal/KO/FVA models start from
      `fixed_medium_base.copy()` so candidate-introduced `e0` boundary exchanges are never omitted.
    - Candidate reactions with zero-width zero bounds (`abs(lb) <= 1e-12 and abs(ub) <= 1e-12`) cannot carry flux
      and are excluded from `cand_rids` so they can never be selected as phantom additions.
    """
    fixed_medium_base = mod_base.copy()
    mod_super = mod_base.copy()
    draft_rxn_ids = set(mod_base.reactions.keys())
    draft_base_ids = {r.base_id for r in mod_base.reactions.values()}

    cand_rids: list[str] = []
    for cid, crxn in sorted(candidate_reactions.items()):
        if cid in draft_rxn_ids or crxn.base_id in draft_base_ids:
            continue
        if crxn.is_biomass:
            continue

        # Check if this is a single-sided boundary/exchange reaction
        if is_single_sided_boundary_reaction(crxn) and crxn.provenance != "explicit_candidate_nutrient_supplement_not_baseline_medium":
            if crxn.is_exchange:
                m_ex = re.match(r"^R_EX_(cpd\d+)_e0$", crxn.id) or re.match(r"^EX_(cpd\d+)_e0$", crxn.base_id)
                cpd_ex = m_ex.group(1) if m_ex else ""
                for sp_id in crxn.stoichiometry:
                    if sp_id not in mod_super.species:
                        m_sp = re.match(r"^M_(cpd\d+)_([cep]0)$", sp_id)
                        sp_rec = SpeciesRecord(
                            id=sp_id,
                            cpd_id=m_sp.group(1) if m_sp else sp_id,
                            name=sp_id,
                            compartment=m_sp.group(2) if m_sp else "e0",
                            charge=None,
                            formula=None,
                        )
                        mod_super.species[sp_id] = sp_rec
                        fixed_medium_base.species[sp_id] = SpeciesRecord(**vars(sp_rec))
                rcopy = ReactionRecord(**vars(crxn))
                if cpd_ex in scenario_compounds:
                    rcopy.lb = -float(scenario_compounds[cpd_ex][1])
                    rcopy.ub = 1000.0
                else:
                    rcopy.lb = 0.0
                    rcopy.ub = 1000.0
                mod_super.reactions[rcopy.id] = rcopy
                fixed_medium_base.reactions[rcopy.id] = ReactionRecord(**vars(rcopy))
                draft_rxn_ids.add(rcopy.id)
            continue

        # Exclude zero-bound candidate reactions (cannot carry flux, never allow as phantom z=1 addition)
        if abs(crxn.lb) <= 1e-12 and abs(crxn.ub) <= 1e-12:
            continue

        # Internal or multi-compartment transport candidate reaction
        for sp_id in crxn.stoichiometry:
            if sp_id not in mod_super.species:
                m_sp = re.match(r"^M_(cpd\d+)_([cep]0)$", sp_id)
                sp_rec = SpeciesRecord(
                    id=sp_id,
                    cpd_id=m_sp.group(1) if m_sp else sp_id,
                    name=sp_id,
                    compartment=m_sp.group(2) if m_sp else "c0",
                    charge=None,
                    formula=None,
                )
                mod_super.species[sp_id] = sp_rec
                if m_sp and m_sp.group(2) == "e0":
                    fixed_medium_base.species[sp_id] = SpeciesRecord(**vars(sp_rec))
                    cpd_e0 = m_sp.group(1)
                    ex_id = f"R_EX_{cpd_e0}_e0"
                    if ex_id not in mod_super.reactions:
                        inflow_lb = -float(scenario_compounds[cpd_e0][1]) if cpd_e0 in scenario_compounds else 0.0
                        ex_rec = ReactionRecord(
                            id=ex_id,
                            base_id=f"EX_{cpd_e0}_e0",
                            name=f"{cpd_e0} Exchange",
                            reversible=True,
                            lb=inflow_lb,
                            ub=1000.0,
                            stoichiometry={sp_id: -1.0},
                            is_exchange=True,
                            provenance="fixed_medium_boundary",
                        )
                        mod_super.reactions[ex_id] = ex_rec
                        fixed_medium_base.reactions[ex_id] = ReactionRecord(**vars(ex_rec))
                        draft_rxn_ids.add(ex_id)

        rcopy = ReactionRecord(**vars(crxn))
        mod_super.reactions[rcopy.id] = rcopy
        cand_rids.append(rcopy.id)

    return fixed_medium_base, mod_super, cand_rids


def _evaluate_cycle_qc_eligibility(egc_res: dict[str, Any]) -> tuple[bool, str]:
    """Classify cycle QC result (`EXP1-FINAL-ENGINE-REVIEW.txt` Item D).
    A model is science-QC-eligible ONLY when closed-uptake EGC check completed cleanly and found NO spurious ATP cycle.
    """
    st = str(egc_res.get("egc_check_status", "unknown"))
    has_cyc = egc_res.get("has_spurious_atp_cycle", None)
    if has_cyc is True:
        return False, "failed_closed_uptake_atp_cycle"
    if has_cyc is False and st in {"evaluated_optimal", "evaluated_infeasible_zero_flux", "evaluated_no_atp_species"}:
        return True, "passed_closed_uptake_egc_check"
    return False, f"unknown_cycle_check_{st}"


def reconstruct_and_verify_solution_witness(
    draft_model: SBMLModelData,
    candidate_reactions: dict[str, ReactionRecord],
    scenario_compounds: dict[str, tuple[str, float]],
    biomass_recipe_id: str = "B0_methanosarcina_template_baseline",
    selected_reaction_ids: list[str] | set[str] | None = None,
    min_growth: float = 0.01,
    time_limit_s: float = 10.0,
) -> dict[str, Any]:
    """Reconstruct an exact witness model from a base draft, scenario bounds, biomass recipe,
    and a specific set of selected reactions, then rigorously verify FBA feasibility, residual,
    bound violations, and closed-uptake EGC cycle QC (`EXP1-REPLAY-EXECUTION-CORRECTION.txt`).
    """
    selected_rids = set(selected_reaction_ids or [])
    mod_base = apply_biomass_perturbation(draft_model, biomass_recipe_id)
    mod_base = apply_nutrient_scenario_bounds(mod_base, scenario_compounds)
    fixed_medium_base, mod_super, _ = build_fixed_medium_supernetwork_and_base(
        mod_base, candidate_reactions, scenario_compounds
    )
    mod_check = fixed_medium_base.copy()
    for rid in sorted(selected_rids):
        if rid in mod_super.reactions:
            crxn = mod_super.reactions[rid]
            for sp_id in crxn.stoichiometry:
                if sp_id not in mod_check.species:
                    mod_check.species[sp_id] = SpeciesRecord(**vars(mod_super.species[sp_id]))
            mod_check.reactions[rid] = ReactionRecord(**vars(crxn))

    fba_check = solve_fba_highs(mod_check, time_limit_s=time_limit_s)
    egc_res = detect_closed_uptake_egc(mod_check, time_limit_s=5.0)
    cyc_eligible, cyc_status = _evaluate_cycle_qc_eligibility(egc_res)

    res_val = float(fba_check.get("max_residual", math.nan))
    bnd_val = float(fba_check.get("max_bound_viol", math.nan))
    gr_val = float(fba_check.get("growth_rate", math.nan)) if fba_check.get("status") == "optimal" else math.nan

    verified = bool(
        fba_check.get("status") == "optimal"
        and math.isfinite(gr_val) and gr_val >= min_growth * 0.999
        and math.isfinite(res_val) and 0.0 <= res_val <= 1e-5
        and math.isfinite(bnd_val) and 0.0 <= bnd_val <= 1e-5
        and cyc_eligible
    )
    return {
        "verified": verified,
        "fba_status": fba_check.get("status"),
        "growth_rate": gr_val,
        "max_residual": res_val,
        "max_bound_viol": bnd_val,
        "cycle_qc_status": cyc_status,
        "has_spurious_atp_cycle": egc_res.get("has_spurious_atp_cycle"),
        "science_qc_eligible": cyc_eligible,
        "model": mod_check,
    }


def enumerate_alternative_gapfill_ensemble(
    draft_model: SBMLModelData,
    candidate_reactions: dict[str, ReactionRecord],
    candidate_weights: dict[str, float],
    scenario_compounds: dict[str, tuple[str, float]],
    biomass_recipe_id: str = "B0_methanosarcina_template_baseline",
    min_growth: float = 0.01,
    n_Target_solutions: int = 20,
    cost_slack_frac: float = 0.25,
    candidate_pool_scope: str = "declared_candidate_reaction_pool",
    time_limit_per_milp_s: float = 8.0,
    mip_rel_gap_tol: float = 1e-4,
) -> dict[str, Any]:
    """Enumerate distinct minimum-cost and near-minimum-cost gap-filling reaction sets using
    EXACT binary MILP (`scipy.optimize.milp`) and exact binary no-good cuts
    (`EXP1-SOLVER-REVIEW.txt` Items 1-5 & `EXP1-FINAL-ENGINE-REVIEW.txt` Items A-E).
    """
    mod_base = apply_biomass_perturbation(draft_model, biomass_recipe_id)
    mod_base = apply_nutrient_scenario_bounds(mod_base, scenario_compounds)
    base_egc = detect_closed_uptake_egc(mod_base, time_limit_s=5.0)
    base_cycle_eligible, base_cycle_status = _evaluate_cycle_qc_eligibility(base_egc)

    # 1. Check if sequence-supported draft already grows on the fixed medium
    draft_fba = solve_fba_highs(mod_base, time_limit_s=time_limit_per_milp_s)
    if draft_fba["status"] not in {"optimal", "infeasible"}:
        return {
            "solver_status": f"draft_fba_{draft_fba['status']}",
            "optimizer_status": draft_fba["status"],
            "raw_solver_status": draft_fba["raw_status"],
            "raw_solver_message": draft_fba["message"],
            "candidate_pool_scope": candidate_pool_scope,
            "enumeration_completeness": "not_evaluated_solver_failure",
            "pool_exhaustion_type": "not_evaluated_solver_failure",
            "cost_window_anchor_type": "not_evaluated_solver_failure",
            "draft_already_grows": False,
            "draft_growth_rate": math.nan,
            "biomass_recipe_id": biomass_recipe_id,
            "min_total_cost": math.nan,
            "max_ensemble_cost": math.nan,
            "first_solve_mip_gap": math.nan,
            "first_solve_dual_bound": math.nan,
            "n_solutions_attempted": 1,
            "n_optimization_feasible_solutions": 0,
            "n_cycle_qc_rejected_solutions": 0,
            "n_distinct_valid_solutions": 0,
            "fewer_than_target_reason": f"draft_fba_{draft_fba['status']}",
            "solutions": [],
            "optimization_feasible_solutions": [],
            "reaction_frequencies": {},
        }

    draft_growth = float(draft_fba["growth_rate"]) if draft_fba["status"] == "optimal" else 0.0
    if draft_fba["status"] == "optimal" and draft_growth >= min_growth:
        draft_sol_entry = {
            "solution_index": 0,
            "added_rxns": [],
            "raw_milp_selected_rxns": [],
            "added_internal_rxns": [],
            "added_transport_rxns": [],
            "total_cost": 0.0,
            "raw_milp_objective": 0.0,
            "raw_binary_cost": 0.0,
            "post_pruning_cost_delta": 0.0,
            "verified_fba_growth": draft_growth,
            "verified_fluxes": draft_fba.get("fluxes", {}),
            "fba_check_max_residual": float(draft_fba.get("max_residual", 0.0)),
            "fba_check_max_bound_viol": float(draft_fba.get("max_bound_viol", 0.0)),
            "mip_gap": 0.0,
            "dual_bound": 0.0,
            "raw_solver_status": int(draft_fba["raw_status"]),
            "raw_solver_message": str(draft_fba["message"]),
            "optimizer_status": "optimal",
            "post_pruning_egc_status": base_egc["egc_check_status"],
            "post_pruning_egc_has_spurious_atp_cycle": base_egc["has_spurious_atp_cycle"],
            "post_pruning_egc_atp_flux": base_egc["closed_uptake_atp_dissipation_flux"],
            "cycle_qc_status": base_cycle_status,
            "science_qc_eligible": base_cycle_eligible,
        }
        if not base_cycle_eligible:
            # Item D: Baseline grows with cycle (or unknown cycle check) is NOT acceptable growth support!
            return {
                "solver_status": f"draft_growth_rejected_by_cycle_qc__{base_cycle_status}",
                "optimizer_status": "optimal",
                "raw_solver_status": int(draft_fba["raw_status"]),
                "raw_solver_message": str(draft_fba["message"]),
                "candidate_pool_scope": candidate_pool_scope,
                "enumeration_completeness": f"rejected_baseline_cycle_qc__{base_cycle_status}",
                "pool_exhaustion_type": "rejected_baseline_cycle_qc",
                "cost_window_anchor_type": "rejected_baseline_cycle_qc",
                "draft_already_grows": False,
                "draft_growth_rate": draft_growth,
                "biomass_recipe_id": biomass_recipe_id,
                "min_total_cost": math.nan,
                "max_ensemble_cost": math.nan,
                "first_solve_mip_gap": 0.0,
                "first_solve_dual_bound": 0.0,
                "max_solution_residual": float(draft_fba.get("max_residual", 0.0)),
                "max_solution_bound_viol": float(draft_fba.get("max_bound_viol", 0.0)),
                "min_solution_verified_growth": draft_growth,
                "n_solutions_attempted": 1,
                "n_optimization_feasible_solutions": 1,
                "n_cycle_qc_rejected_solutions": 1,
                "n_distinct_valid_solutions": 0,
                "fewer_than_target_reason": f"draft_growth_rejected_by_cycle_qc__{base_cycle_status}",
                "solutions": [],
                "optimization_feasible_solutions": [draft_sol_entry],
                "reaction_frequencies": {},
            }
        return {
            "solver_status": "optimal_draft_no_gapfill_needed",
            "optimizer_status": "optimal",
            "raw_solver_status": int(draft_fba["raw_status"]),
            "raw_solver_message": str(draft_fba["message"]),
            "candidate_pool_scope": candidate_pool_scope,
            "enumeration_completeness": "draft_sufficient_zero_additions",
            "pool_exhaustion_type": "draft_sufficient_zero_additions",
            "cost_window_anchor_type": "zero_cost_draft_sufficient",
            "draft_already_grows": True,
            "draft_growth_rate": draft_growth,
            "biomass_recipe_id": biomass_recipe_id,
            "min_total_cost": 0.0,
            "max_ensemble_cost": 0.0,
            "first_solve_mip_gap": 0.0,
            "first_solve_dual_bound": 0.0,
            "max_solution_residual": float(draft_fba.get("max_residual", 0.0)),
            "max_solution_bound_viol": float(draft_fba.get("max_bound_viol", 0.0)),
            "min_solution_verified_growth": draft_growth,
            "n_solutions_attempted": 1,
            "n_optimization_feasible_solutions": 1,
            "n_cycle_qc_rejected_solutions": 0,
            "n_distinct_valid_solutions": 1,
            "fewer_than_target_reason": "draft_model_grows_without_additions",
            "solutions": [draft_sol_entry],
            "optimization_feasible_solutions": [draft_sol_entry],
            "reaction_frequencies": {},
        }

    # 2. Build frozen fixed_medium_base and mod_super (Fixes Item B & C)
    fixed_medium_base, mod_super, cand_rids = build_fixed_medium_supernetwork_and_base(
        mod_base=mod_base,
        candidate_reactions=candidate_reactions,
        scenario_compounds=scenario_compounds,
    )

    # 3. Check supernetwork feasibility under the fixed medium (relaxing candidate bounds to include 0 when z_k = 0)
    mod_super_lp = mod_super.copy()
    for crid in cand_rids:
        mod_super_lp.reactions[crid].lb = min(0.0, mod_super_lp.reactions[crid].lb)
        mod_super_lp.reactions[crid].ub = max(0.0, mod_super_lp.reactions[crid].ub)
    super_fba = solve_fba_highs(mod_super_lp, time_limit_s=time_limit_per_milp_s)
    if super_fba["status"] == "infeasible" or (super_fba["status"] == "optimal" and float(super_fba["growth_rate"]) < min_growth):
        return {
            "solver_status": "infeasible_under_fixed_medium",
            "optimizer_status": "infeasible",
            "raw_solver_status": int(super_fba["raw_status"]),
            "raw_solver_message": str(super_fba["message"]),
            "candidate_pool_scope": candidate_pool_scope,
            "enumeration_completeness": "proven_infeasible_under_fixed_medium_on_declared_pool",
            "pool_exhaustion_type": "proven_infeasible_supernetwork",
            "cost_window_anchor_type": "infeasible_no_anchor",
            "draft_already_grows": False,
            "draft_growth_rate": draft_growth,
            "biomass_recipe_id": biomass_recipe_id,
            "min_total_cost": math.nan,
            "max_ensemble_cost": math.nan,
            "first_solve_mip_gap": math.nan,
            "first_solve_dual_bound": math.nan,
            "n_solutions_attempted": 1,
            "n_optimization_feasible_solutions": 0,
            "n_cycle_qc_rejected_solutions": 0,
            "n_distinct_valid_solutions": 0,
            "fewer_than_target_reason": "infeasible_under_fixed_medium",
            "solutions": [],
            "optimization_feasible_solutions": [],
            "reaction_frequencies": {},
        }
    if super_fba["status"] != "optimal":
        return {
            "solver_status": f"supernetwork_fba_{super_fba['status']}",
            "optimizer_status": super_fba["status"],
            "raw_solver_status": int(super_fba["raw_status"]),
            "raw_solver_message": str(super_fba["message"]),
            "candidate_pool_scope": candidate_pool_scope,
            "enumeration_completeness": "not_evaluated_solver_failure",
            "pool_exhaustion_type": "not_evaluated_solver_failure",
            "cost_window_anchor_type": "not_evaluated_solver_failure",
            "draft_already_grows": False,
            "draft_growth_rate": draft_growth,
            "biomass_recipe_id": biomass_recipe_id,
            "min_total_cost": math.nan,
            "max_ensemble_cost": math.nan,
            "first_solve_mip_gap": math.nan,
            "first_solve_dual_bound": math.nan,
            "n_solutions_attempted": 1,
            "n_optimization_feasible_solutions": 0,
            "n_cycle_qc_rejected_solutions": 0,
            "n_distinct_valid_solutions": 0,
            "fewer_than_target_reason": f"supernetwork_fba_{super_fba['status']}",
            "solutions": [],
            "optimization_feasible_solutions": [],
            "reaction_frequencies": {},
        }

    # 4. Assemble exact MILP: x = [v (n_rxn), z (n_cand)]
    sp_ids, rxn_ids, S, lb_v, ub_v, obj_j = build_sparse_stoichiometric_matrix(mod_super)
    n_sp = len(sp_ids)
    n_rxn = len(rxn_ids)
    rxn_idx = {r: j for j, r in enumerate(rxn_ids)}

    cand_rids = [r for r in cand_rids if r in rxn_idx]
    n_cand = len(cand_rids)
    if n_cand == 0:
        return {
            "solver_status": "infeasible_empty_candidate_pool",
            "optimizer_status": "infeasible",
            "raw_solver_status": 2,
            "raw_solver_message": "No eligible candidate reactions",
            "candidate_pool_scope": candidate_pool_scope,
            "enumeration_completeness": "proven_infeasible_empty_candidate_pool",
            "pool_exhaustion_type": "proven_infeasible_empty_candidate_pool",
            "cost_window_anchor_type": "infeasible_no_anchor",
            "draft_already_grows": False,
            "draft_growth_rate": draft_growth,
            "biomass_recipe_id": biomass_recipe_id,
            "min_total_cost": math.nan,
            "max_ensemble_cost": math.nan,
            "first_solve_mip_gap": math.nan,
            "first_solve_dual_bound": math.nan,
            "n_solutions_attempted": 1,
            "n_optimization_feasible_solutions": 0,
            "n_cycle_qc_rejected_solutions": 0,
            "n_distinct_valid_solutions": 0,
            "fewer_than_target_reason": "infeasible_empty_candidate_pool",
            "solutions": [],
            "optimization_feasible_solutions": [],
            "reaction_frequencies": {},
        }

    cand_v_cols = np.array([rxn_idx[r] for r in cand_rids], dtype=np.int64)
    cand_idx_map = {r: k for k, r in enumerate(cand_rids)}

    # Enforce minimum biomass growth on v[obj_j]
    lb_v = lb_v.copy()
    lb_v[obj_j] = max(lb_v[obj_j], float(min_growth))

    # Candidate base costs c_z (strictly positive)
    c_z = np.zeros(n_cand, dtype=np.float64)
    for k, rid in enumerate(cand_rids):
        robj = mod_super.reactions[rid]
        w = float(candidate_weights.get(rid, candidate_weights.get(robj.base_id, 100.0)))
        c_z[k] = max(w, 0.005)

    # Item C fix: variable bounds on candidate fluxes v_j must permit 0 when z_k = 0 even if native lb_j > 0 or ub_j < 0,
    # while coupling constraints v_j - ub_j * z_k <= 0 and -v_j + lb_j * z_k <= 0 enforce native [lb_j, ub_j] when z_k = 1!
    ub_cand_native = ub_v[cand_v_cols].copy()
    lb_cand_native = lb_v[cand_v_cols].copy()
    lb_v_var = lb_v.copy()
    ub_v_var = ub_v.copy()
    lb_v_var[cand_v_cols] = np.minimum(0.0, lb_cand_native)
    ub_v_var[cand_v_cols] = np.maximum(0.0, ub_cand_native)

    c_milp = np.concatenate([np.zeros(n_rxn, dtype=np.float64), c_z])
    integrality = np.concatenate([np.zeros(n_rxn, dtype=np.int8), np.ones(n_cand, dtype=np.int8)])
    lb_x = np.concatenate([lb_v_var, np.zeros(n_cand, dtype=np.float64)])
    ub_x = np.concatenate([ub_v_var, np.ones(n_cand, dtype=np.float64)])

    # Steady-state constraint: [S, 0] x == 0
    A_eq = sparse.hstack([S, sparse.csc_matrix((n_sp, n_cand), dtype=np.float64)], format="csc")
    lb_eq = np.zeros(n_sp, dtype=np.float64)
    ub_eq = np.zeros(n_sp, dtype=np.float64)

    # Exact Big-M coupling constraints using native bounds:
    # 1) v_j - ub_cand_native[k] * z_k <= 0  (so v_j <= 0 when z_k=0, and v_j <= ub_j when z_k=1)
    # 2) -v_j + lb_cand_native[k] * z_k <= 0 (so v_j >= 0 when z_k=0, and v_j >= lb_j when z_k=1)
    k_idx = np.arange(n_cand, dtype=np.int64)
    z_cols = n_rxn + k_idx
    row_1 = k_idx
    row_2 = n_cand + k_idx
    coup_rows = np.concatenate([row_1, row_1, row_2, row_2])
    coup_cols = np.concatenate([cand_v_cols, z_cols, cand_v_cols, z_cols])
    coup_vals = np.concatenate([
        np.ones(n_cand, dtype=np.float64),
        -ub_cand_native,
        -np.ones(n_cand, dtype=np.float64),
        lb_cand_native,
    ])
    A_coup = sparse.csc_matrix((coup_vals, (coup_rows, coup_cols)), shape=(2 * n_cand, n_rxn + n_cand), dtype=np.float64)
    lb_coup = np.full(2 * n_cand, -np.inf, dtype=np.float64)
    ub_coup = np.zeros(2 * n_cand, dtype=np.float64)

    A_base = sparse.vstack([A_eq, A_coup], format="csc")
    lb_base = np.concatenate([lb_eq, lb_coup])
    ub_base = np.concatenate([ub_eq, ub_coup])

    extra_rows_list: list[sparse.csc_matrix] = []
    extra_lb_list: list[float] = []
    extra_ub_list: list[float] = []

    def _add_linear_cut(cols: np.ndarray, vals: np.ndarray, lb_val: float, ub_val: float) -> None:
        r_mat = sparse.csc_matrix(
            (vals, (np.zeros(len(cols), dtype=np.int64), cols)),
            shape=(1, n_rxn + n_cand),
            dtype=np.float64,
        )
        extra_rows_list.append(r_mat)
        extra_lb_list.append(float(lb_val))
        extra_ub_list.append(float(ub_val))

    def _add_binary_nogood_cut(selected_k_indices: set[int]) -> None:
        # Exact binary no-good cut: sum_{k in S} z_k - sum_{k not in S} z_k <= |S| - 1
        coeffs = np.full(n_cand, -1.0, dtype=np.float64)
        for k_sel in selected_k_indices:
            coeffs[k_sel] = 1.0
        _add_linear_cut(z_cols, coeffs, -np.inf, float(len(selected_k_indices) - 1))

    seen_sets: set[frozenset[str]] = set()
    optimization_feasible_solutions: list[dict[str, Any]] = []
    science_qc_eligible_solutions: list[dict[str, Any]] = []
    n_attempted = 0
    first_cost_anchor: float | None = None
    first_mip_gap: float = math.nan
    first_dual_bound: float = math.nan
    first_raw_status: int = -1
    first_raw_message: str = ""
    cost_window_anchor_type: str = "not_set"
    all_solves_proven_optimal = True
    terminated_by_infeasibility = False
    last_solver_status = "not_run"

    max_milp_solves = max(int(n_Target_solutions) + 10, 28)
    for solve_idx in range(max_milp_solves):
        if len(science_qc_eligible_solutions) >= n_Target_solutions:
            break
        n_attempted += 1
        if extra_rows_list:
            A_curr = sparse.vstack([A_base, *extra_rows_list], format="csc")
            lb_curr = np.concatenate([lb_base, np.array(extra_lb_list, dtype=np.float64)])
            ub_curr = np.concatenate([ub_base, np.array(extra_ub_list, dtype=np.float64)])
        else:
            A_curr = A_base
            lb_curr = lb_base
            ub_curr = ub_base

        res = milp(
            c=c_milp,
            integrality=integrality,
            bounds=Bounds(lb_x, ub_x),
            constraints=LinearConstraint(A_curr, lb_curr, ub_curr),
            options={
                "time_limit": float(time_limit_per_milp_s),
                "mip_rel_gap": float(mip_rel_gap_tol),
                "mip_feasibility_tolerance": 1e-9,
                "primal_feasibility_tolerance": 1e-9,
            },
        )
        raw_st = int(getattr(res, "status", -1))
        raw_msg = str(getattr(res, "message", ""))
        raw_gap_attr = getattr(res, "mip_gap", None)
        mip_gap_val = float(raw_gap_attr) if raw_gap_attr is not None else math.nan
        raw_db_attr = getattr(res, "mip_dual_bound", None)
        dual_bnd = float(raw_db_attr) if raw_db_attr is not None else math.nan

        if solve_idx == 0:
            first_raw_status = raw_st
            first_raw_message = raw_msg

        if raw_st == 2:
            terminated_by_infeasibility = True
            last_solver_status = "infeasible"
            break

        if res.x is None:
            all_solves_proven_optimal = False
            last_solver_status = "timeout_no_incumbent" if raw_st == 1 else f"solver_failed_status_{raw_st}"
            break

        z_sol = res.x[n_rxn:]
        v_sol = res.x[:n_rxn]
        obj_val = float(res.fun) if getattr(res, "fun", None) is not None else float(np.dot(c_milp, res.x))

        # Item C fix: binary no-good cut uses binary vector z_sol ONLY, with numerical integrality & coupling checked
        max_z_integrality_err = float(np.max(np.abs(z_sol - np.round(z_sol)))) if n_cand > 0 else 0.0
        raw_selected_k = {int(k) for k in range(n_cand) if z_sol[k] >= 0.5}
        unselected_k = [k for k in range(n_cand) if k not in raw_selected_k]
        max_unselected_flux_leak = (
            float(np.max(np.abs(v_sol[cand_v_cols[unselected_k]]))) if unselected_k else 0.0
        )
        numerical_coupling_valid = (max_z_integrality_err <= 1e-4) and (max_unselected_flux_leak <= 1e-5)

        # Always add binary no-good cut on raw_selected_k so the MILP never repeats the same binary vector
        _add_binary_nogood_cut(raw_selected_k)

        if not numerical_coupling_valid:
            all_solves_proven_optimal = False
            last_solver_status = "numerical_coupling_violation"
            continue

        # Item A fix: prove optimal ONLY when raw_st == 0 and actual res.mip_dual_bound & res.mip_gap meet tolerance
        bound_gap_ok = (
            not math.isnan(mip_gap_val)
            and not math.isnan(dual_bnd)
            and mip_gap_val <= mip_rel_gap_tol + 1e-6
            and (obj_val - dual_bnd) <= max(1e-6, (mip_rel_gap_tol + 1e-6) * max(1.0, abs(obj_val)))
        )
        if raw_st == 0 and bound_gap_ok:
            opt_stat = "optimal"
        elif raw_st == 1:
            opt_stat = "time_limit_heuristic_incumbent"
            all_solves_proven_optimal = False
        else:
            opt_stat = f"suboptimal_status_{raw_st}"
            all_solves_proven_optimal = False

        # Item B fix: verify on fixed_medium_base + selected candidate reactions (never omitting candidate e0 boundary exchanges)
        raw_selected_rids = sorted(cand_rids[k] for k in raw_selected_k)
        mod_check = fixed_medium_base.copy()
        for rid in raw_selected_rids:
            crxn = mod_super.reactions[rid]
            for sp_id in crxn.stoichiometry:
                if sp_id not in mod_check.species:
                    mod_check.species[sp_id] = SpeciesRecord(**vars(mod_super.species[sp_id]))
            mod_check.reactions[rid] = ReactionRecord(**vars(crxn))

        fba_check = solve_fba_highs(mod_check, time_limit_s=time_limit_per_milp_s)
        if fba_check["status"] != "optimal" or float(fba_check["growth_rate"]) < min_growth * 0.999:
            all_solves_proven_optimal = False
            last_solver_status = "fba_verification_failed"
            continue

        pruned_rids = sorted(
            rid for rid in raw_selected_rids if abs(fba_check["fluxes"].get(rid, 0.0)) > 1e-8
        )
        pruned_k = {cand_idx_map[r] for r in pruned_rids}
        if pruned_k != raw_selected_k:
            _add_binary_nogood_cut(pruned_k)
            mod_pruned = fixed_medium_base.copy()
            for rid in pruned_rids:
                crxn = mod_super.reactions[rid]
                for sp_id in crxn.stoichiometry:
                    if sp_id not in mod_pruned.species:
                        mod_pruned.species[sp_id] = SpeciesRecord(**vars(mod_super.species[sp_id]))
                mod_pruned.reactions[rid] = ReactionRecord(**vars(crxn))
            fba_pruned = solve_fba_highs(mod_pruned, time_limit_s=time_limit_per_milp_s)
            if fba_pruned["status"] == "optimal" and float(fba_pruned["growth_rate"]) >= min_growth * 0.999:
                fba_check = fba_pruned
                mod_check = mod_pruned
            else:
                pruned_rids = raw_selected_rids
                pruned_k = raw_selected_k

        raw_binary_cost = float(sum(c_z[k] for k in raw_selected_k))
        pruned_binary_cost = float(sum(c_z[k] for k in pruned_k))
        # Item E fix: if post-pruning cost changes unexpectedly beyond tolerance, optimum claim is invalid
        cost_change_after_pruning = abs(pruned_binary_cost - obj_val)
        if cost_change_after_pruning > max(1e-6, (mip_rel_gap_tol + 1e-6) * max(1.0, abs(obj_val))):
            opt_stat = "post_pruning_cost_changed_not_proven_optimal"
            all_solves_proven_optimal = False

        last_solver_status = opt_stat

        if first_cost_anchor is None:
            first_mip_gap = mip_gap_val
            first_dual_bound = dual_bnd
            if opt_stat == "optimal":
                # Item E: anchor near-optimal pool to raw proven MILP objective on declared pool
                first_cost_anchor = obj_val
                cost_window_anchor_type = "proven_global_milp_optimum_on_declared_pool"
            else:
                first_cost_anchor = pruned_binary_cost
                cost_window_anchor_type = "incumbent_anchored_not_proven_optimal"
            max_allowed_cost = first_cost_anchor * (1.0 + float(cost_slack_frac)) + 1e-7
            _add_linear_cut(z_cols, c_z, -np.inf, max_allowed_cost)

        cand_key = frozenset(pruned_rids)
        if cand_key in seen_sets:
            continue
        seen_sets.add(cand_key)

        # Item D fix: evaluate closed-uptake EGC cycle QC and separate optimization-feasible vs science-QC-eligible solutions
        sol_egc = detect_closed_uptake_egc(mod_check, time_limit_s=5.0)
        cyc_eligible, cyc_status = _evaluate_cycle_qc_eligibility(sol_egc)

        internal_added = [r for r in pruned_rids if not mod_super.reactions[r].is_transport]
        transport_added = [r for r in pruned_rids if mod_super.reactions[r].is_transport]
        sol_record = {
            "solution_index": len(science_qc_eligible_solutions) if cyc_eligible else -1,
            "optimization_solution_index": len(optimization_feasible_solutions),
            "added_rxns": pruned_rids,
            "raw_milp_selected_rxns": raw_selected_rids,
            "added_internal_rxns": internal_added,
            "added_transport_rxns": transport_added,
            "total_cost": pruned_binary_cost,
            "raw_milp_objective": obj_val,
            "raw_binary_cost": raw_binary_cost,
            "post_pruning_cost_delta": cost_change_after_pruning,
            "verified_fba_growth": float(fba_check["growth_rate"]),
            "verified_fluxes": fba_check["fluxes"],
            "fba_check_max_residual": float(fba_check.get("max_residual", 0.0)),
            "fba_check_max_bound_viol": float(fba_check.get("max_bound_viol", 0.0)),
            "mip_gap": mip_gap_val,
            "dual_bound": dual_bnd,
            "raw_solver_status": raw_st,
            "raw_solver_message": raw_msg,
            "optimizer_status": opt_stat,
            "post_pruning_egc_status": sol_egc["egc_check_status"],
            "post_pruning_egc_has_spurious_atp_cycle": sol_egc["has_spurious_atp_cycle"],
            "post_pruning_egc_atp_flux": sol_egc["closed_uptake_atp_dissipation_flux"],
            "cycle_qc_status": cyc_status,
            "science_qc_eligible": cyc_eligible,
        }
        optimization_feasible_solutions.append(sol_record)
        if cyc_eligible:
            science_qc_eligible_solutions.append(sol_record)

    n_opt_feasible = len(optimization_feasible_solutions)
    n_distinct = len(science_qc_eligible_solutions)
    n_cyc_rejected = n_opt_feasible - n_distinct

    if n_distinct == 0:
        if n_opt_feasible > 0:
            fail_reason = "all_optimization_feasible_solutions_rejected_by_cycle_qc"
            enum_comp = "rejected_all_solutions_failed_or_unknown_cycle_qc"
            pool_exh = (
                "mathematically_exhausted_all_rejected_by_cycle_qc"
                if (all_solves_proven_optimal and terminated_by_infeasibility)
                else "sampled_all_rejected_by_cycle_qc"
            )
        else:
            fail_reason = last_solver_status if last_solver_status != "not_run" else "no_verified_solution"
            enum_comp = "no_solution_found"
            pool_exh = "no_solution_found"
        return {
            "solver_status": fail_reason,
            "optimizer_status": last_solver_status,
            "raw_solver_status": first_raw_status,
            "raw_solver_message": first_raw_message,
            "candidate_pool_scope": candidate_pool_scope,
            "enumeration_completeness": enum_comp,
            "pool_exhaustion_type": pool_exh,
            "cost_window_anchor_type": cost_window_anchor_type,
            "draft_already_grows": False,
            "draft_growth_rate": draft_growth,
            "biomass_recipe_id": biomass_recipe_id,
            "min_total_cost": math.nan,
            "max_ensemble_cost": math.nan,
            "first_solve_mip_gap": first_mip_gap,
            "first_solve_dual_bound": first_dual_bound,
            "n_solutions_attempted": n_attempted,
            "n_optimization_feasible_solutions": n_opt_feasible,
            "n_cycle_qc_rejected_solutions": n_cyc_rejected,
            "n_distinct_valid_solutions": 0,
            "fewer_than_target_reason": fail_reason,
            "solutions": [],
            "optimization_feasible_solutions": optimization_feasible_solutions,
            "reaction_frequencies": {},
        }

    # Item D & E: explicit mathematical vs cycle-QC-sampled pool exhaustion classification
    if all_solves_proven_optimal and terminated_by_infeasibility:
        if n_cyc_rejected == 0:
            enum_comp = "exhausted_within_cost_slack_on_declared_pool"
            pool_exh = "mathematically_exhausted_and_all_cycle_qc_eligible"
            fewer_reason = "exhausted_all_feasible_sets_within_cost_slack_on_declared_pool"
        else:
            enum_comp = "exhausted_within_cost_slack_after_cycle_qc_filtering"
            pool_exh = "mathematically_exhausted_after_cycle_qc_filtering"
            fewer_reason = f"exhausted_after_filtering_{n_cyc_rejected}_cycle_qc_ineligible_solutions"
    elif n_distinct >= n_Target_solutions:
        enum_comp = (
            "bounded_sample_at_target_cap_not_exhausted"
            if all_solves_proven_optimal
            else "heuristic_incumbent_sample_not_proven_optimal"
        )
        pool_exh = (
            "cycle_qc_eligible_target_cap_sampled_not_exhausted"
            if n_cyc_rejected > 0
            else "target_cap_sampled_not_mathematically_exhausted"
        )
        fewer_reason = "target_cap_reached_more_alternatives_may_exist"
    else:
        enum_comp = "heuristic_or_time_limited_sample_not_exhausted"
        pool_exh = "incomplete_or_heuristic_sample_not_exhausted"
        fewer_reason = f"stopped_early_{last_solver_status}"

    # Item D: reaction frequencies are computed ONLY among cycle-QC-eligible solutions (`science_qc_eligible_solutions`)
    rxn_counts: dict[str, int] = {}
    for sol in science_qc_eligible_solutions:
        for r in sol["added_rxns"]:
            rxn_counts[r] = rxn_counts.get(r, 0) + 1

    rxn_freqs: dict[str, dict[str, Any]] = {}
    for r, cnt in sorted(rxn_counts.items()):
        freq = float(cnt) / float(n_distinct)
        rxn_obj = mod_super.reactions[r]
        w = float(c_z[cand_idx_map[r]])

        # Direct supernetwork KO test
        mod_ko_super = mod_super_lp.copy()
        mod_ko_super.reactions[r].lb = 0.0
        mod_ko_super.reactions[r].ub = 0.0
        ko_sol = solve_fba_highs(mod_ko_super, time_limit_s=time_limit_per_milp_s)
        if ko_sol["status"] == "infeasible" or (ko_sol["status"] == "optimal" and float(ko_sol["growth_rate"]) < min_growth):
            super_necessity = "model_specific_essential_in_declared_supernetwork"
            ko_growth = float(ko_sol["growth_rate"]) if ko_sol["status"] == "optimal" else 0.0
        elif ko_sol["status"] == "optimal":
            super_necessity = "bypassable_in_declared_supernetwork"
            ko_growth = float(ko_sol["growth_rate"])
        else:
            super_necessity = f"unknown_ko_solver_{ko_sol['status']}"
            ko_growth = math.nan

        rxn_freqs[r] = {
            "rxn_id": r,
            "base_id": rxn_obj.base_id,
            "name": rxn_obj.name,
            "is_exchange": rxn_obj.is_exchange,
            "is_transport": rxn_obj.is_transport,
            "candidate_weight": w,
            "n_solutions_present": cnt,
            "n_distinct_solutions_total": n_distinct,
            "n_optimization_feasible_solutions_total": n_opt_feasible,
            "n_cycle_qc_rejected_solutions_total": n_cyc_rejected,
            "ensemble_frequency": freq,
            "stability_tier": _classify_ensemble_frequency_label(freq),
            "supernetwork_ko_necessity": super_necessity,
            "supernetwork_ko_growth": ko_growth,
            "candidate_pool_scope": candidate_pool_scope,
            "enumeration_completeness": enum_comp,
            "pool_exhaustion_type": pool_exh,
            "cost_window_anchor_type": cost_window_anchor_type,
        }

    costs = [s["total_cost"] for s in science_qc_eligible_solutions]
    if all_solves_proven_optimal and cost_window_anchor_type == "proven_global_milp_optimum_on_declared_pool":
        overall_opt_status = "optimal" if n_cyc_rejected == 0 else "optimal_after_cycle_qc_filtering"
    else:
        overall_opt_status = "heuristic_incumbent_bounded"
    return {
        "solver_status": overall_opt_status,
        "optimizer_status": overall_opt_status,
        "raw_solver_status": first_raw_status,
        "raw_solver_message": first_raw_message,
        "candidate_pool_scope": candidate_pool_scope,
        "enumeration_completeness": enum_comp,
        "pool_exhaustion_type": pool_exh,
        "cost_window_anchor_type": cost_window_anchor_type,
        "draft_already_grows": False,
        "draft_growth_rate": draft_growth,
        "biomass_recipe_id": biomass_recipe_id,
        "min_total_cost": float(min(costs)),
        "max_ensemble_cost": float(max(costs)),
        "first_solve_mip_gap": first_mip_gap,
        "first_solve_dual_bound": first_dual_bound,
        "max_solution_residual": float(max((s.get("fba_check_max_residual", 0.0) for s in science_qc_eligible_solutions), default=math.nan)),
        "max_solution_bound_viol": float(max((s.get("fba_check_max_bound_viol", 0.0) for s in science_qc_eligible_solutions), default=math.nan)),
        "min_solution_verified_growth": float(min((s["verified_fba_growth"] for s in science_qc_eligible_solutions), default=math.nan)),
        "n_solutions_attempted": n_attempted,
        "n_optimization_feasible_solutions": n_opt_feasible,
        "n_cycle_qc_rejected_solutions": n_cyc_rejected,
        "n_distinct_valid_solutions": n_distinct,
        "fewer_than_target_reason": fewer_reason,
        "solutions": science_qc_eligible_solutions,
        "optimization_feasible_solutions": optimization_feasible_solutions,
        "reaction_frequencies": rxn_freqs,
    }


def solve_nutrient_supplement_rescue_milp(
    draft_model: SBMLModelData,
    candidate_reactions: dict[str, ReactionRecord],
    candidate_weights: dict[str, float],
    fixed_scenario_compounds: dict[str, tuple[str, float]],
    supplement_candidate_compounds: dict[str, tuple[str, float, float]],
    biomass_recipe_id: str = "B0_methanosarcina_template_baseline",
    redox_regime: str = "anoxic",
    min_growth: float = 0.01,
    time_limit_s: float = 10.0,
) -> dict[str, Any]:
    """Explicit, SEPARATE exploratory nutrient-supplement rescue MILP (`EXP1-SOLVER-REVIEW.txt` Item 2 & 6).
    Never relabels the fixed baseline medium; tests whether adding binary-penalized external nutrient
    supplements (`supplement_candidate_compounds`: `cpd_id -> (name, max_uptake, supplement_cost)`)
    can rescue a model that is infeasible on `fixed_scenario_compounds`.
    Strictly forbids O2 (`cpd00007`) as a supplement when `redox_regime == 'anoxic'`.
    """
    if redox_regime == "anoxic" and "cpd00007" in supplement_candidate_compounds:
        raise ValueError("Oxygen (cpd00007) cannot be supplied as a candidate supplement in an anoxic regime")

    # Build augmented candidate reaction pool including explicit supplement exchanges
    aug_candidates = dict(candidate_reactions)
    aug_weights = dict(candidate_weights)
    aug_scenario = dict(fixed_scenario_compounds)

    mod_test = apply_biomass_perturbation(draft_model, biomass_recipe_id)
    mod_test = apply_nutrient_scenario_bounds(mod_test, fixed_scenario_compounds)

    # Represent each candidate supplement as an explicit candidate import reaction `R_SUPP_IMPORT_<cpd>_e0` (`ø -> M_<cpd>_e0`)
    # Wait: to preserve single-sided boundary closure on `mod_base`, we model supplement import as a candidate reaction
    # bounded by `[0, max_uptake]` with binary indicator `z_supp`!
    for cpd_id, (cname, max_up, supp_cost) in sorted(supplement_candidate_compounds.items()):
        if cpd_id in fixed_scenario_compounds:
            continue
        sp_e0 = f"M_{cpd_id}_e0"
        supp_rid = f"R_SUPP_UPTAKE_{cpd_id}_e0"
        aug_candidates[supp_rid] = ReactionRecord(
            id=supp_rid,
            base_id=f"SUPP_UPTAKE_{cpd_id}_e0",
            name=f"Candidate supplement uptake for {cname} ({cpd_id})",
            reversible=False,
            lb=0.0,
            ub=float(max_up),
            stoichiometry={sp_e0: 1.0, f"M_{cpd_id}_c0": 0.0} if f"M_{cpd_id}_c0" in draft_model.species else {sp_e0: 1.0},
            is_exchange=False,
            is_transport=True,
            provenance="explicit_candidate_nutrient_supplement_not_baseline_medium",
        )
        aug_weights[supp_rid] = float(supp_cost)

    res = enumerate_alternative_gapfill_ensemble(
        draft_model=draft_model,
        candidate_reactions=aug_candidates,
        candidate_weights=aug_weights,
        scenario_compounds=aug_scenario,
        biomass_recipe_id=biomass_recipe_id,
        min_growth=min_growth,
        n_Target_solutions=5,
        candidate_pool_scope="explicit_supplement_rescue_experiment_not_baseline_medium",
        time_limit_per_milp_s=time_limit_s,
    )
    supplements_used: list[str] = []
    if res["solutions"]:
        for r in res["solutions"][0]["added_rxns"]:
            if r.startswith("R_SUPP_UPTAKE_"):
                supplements_used.append(r)
    return {
        "experiment_type": "separate_exploratory_nutrient_supplement_rescue_not_baseline_medium",
        "rescued": len(res["solutions"]) > 0,
        "solver_status": res["solver_status"],
        "primary_supplements_required": supplements_used,
        "ensemble_result": res,
    }


# =============================================================================
# 7. Nutrient Withdrawal / Add-Back Rescue, Knockouts, & FVA Counterfactuals
# =============================================================================

def run_withdrawal_knockout_fva(
    base_model: SBMLModelData,
    ensemble_result: dict[str, Any],
    candidate_reactions: dict[str, ReactionRecord],
    scenario_compounds: dict[str, tuple[str, float]],
    biomass_recipe_id: str = "B0_methanosarcina_template_baseline",
    min_growth: float = 0.01,
    time_limit_s: float = 10.0,
) -> dict[str, Any]:
    """Run nutrient withdrawal/add-back rescue, reaction knockouts, and FVA on a model + ensemble (`EXP1-CONTRACT.txt` Item 3).
    CRITICAL RULE: If baseline model has zero growth (`< 1e-6`), failed solver, or fails closed-uptake cycle QC,
    return a `non_evaluable_*` status and NEVER classify a zero-growth or cycle-corrupted model as an auxotroph!
    """
    mod_scenario = apply_biomass_perturbation(base_model, biomass_recipe_id)
    mod_scenario = apply_nutrient_scenario_bounds(mod_scenario, scenario_compounds)
    fixed_medium_base, _mod_super_unused, _cand_unused = build_fixed_medium_supernetwork_and_base(
        mod_base=mod_scenario,
        candidate_reactions=candidate_reactions,
        scenario_compounds=scenario_compounds,
    )

    solutions = [
        s for s in ensemble_result.get("solutions", [])
        if s.get("science_qc_eligible", True) is True
    ]
    if not solutions:
        if int(ensemble_result.get("n_cycle_qc_rejected_solutions", 0)) > 0:
            return {
                "evaluation_status": "non_evaluable_baseline_failed_cycle_qc",
                "baseline_growth_rate": math.nan,
                "withdrawals": [],
                "knockouts": [],
                "fva": [],
            }
        base_sol = solve_fba_highs(fixed_medium_base, time_limit_s=time_limit_s)
        if base_sol["status"] not in {"optimal", "infeasible"}:
            return {
                "evaluation_status": f"non_evaluable_solver_{base_sol['status']}",
                "baseline_growth_rate": math.nan,
                "withdrawals": [],
                "knockouts": [],
                "fva": [],
            }
        g0 = float(base_sol["growth_rate"]) if base_sol["status"] == "optimal" else 0.0
        if g0 < 1e-6:
            return {
                "evaluation_status": "non_evaluable_zero_baseline_growth",
                "baseline_growth_rate": 0.0,
                "withdrawals": [],
                "knockouts": [],
                "fva": [],
            }
        base_egc = detect_closed_uptake_egc(fixed_medium_base, time_limit_s=5.0)
        cyc_ok, cyc_st = _evaluate_cycle_qc_eligibility(base_egc)
        if not cyc_ok:
            return {
                "evaluation_status": f"non_evaluable_baseline_failed_cycle_qc__{cyc_st}",
                "baseline_growth_rate": math.nan,
                "withdrawals": [],
                "knockouts": [],
                "fva": [],
            }
        union_added: list[str] = []
        primary_added: list[str] = []
    else:
        primary_added = list(solutions[0]["added_rxns"])
        union_set: set[str] = set()
        for s in solutions:
            union_set.update(s["added_rxns"])
        union_added = sorted(union_set)

    def _build_augmented(added_list: list[str]) -> SBMLModelData:
        m = fixed_medium_base.copy()
        for rid in added_list:
            crxn = candidate_reactions[rid]
            for sp_id in crxn.stoichiometry:
                if sp_id not in m.species:
                    m_sp = re.match(r"^M_(cpd\d+)_([cep]0)$", sp_id)
                    m.species[sp_id] = SpeciesRecord(
                        id=sp_id,
                        cpd_id=m_sp.group(1) if m_sp else sp_id,
                        name=sp_id,
                        compartment=m_sp.group(2) if m_sp else "c0",
                        charge=None,
                        formula=None,
                    )
            m.reactions[rid] = ReactionRecord(**vars(crxn))
        return m

    mod_primary = _build_augmented(primary_added)
    mod_union = _build_augmented(union_added)

    prim_egc = detect_closed_uptake_egc(mod_primary, time_limit_s=5.0)
    prim_cyc_ok, prim_cyc_st = _evaluate_cycle_qc_eligibility(prim_egc)
    if not prim_cyc_ok:
        return {
            "evaluation_status": f"non_evaluable_baseline_failed_cycle_qc__{prim_cyc_st}",
            "baseline_growth_rate": math.nan,
            "withdrawals": [],
            "knockouts": [],
            "fva": [],
        }

    sol0 = solve_fba_highs(mod_primary, time_limit_s=time_limit_s)
    if sol0["status"] != "optimal":
        return {
            "evaluation_status": f"non_evaluable_baseline_{sol0['status']}",
            "baseline_growth_rate": math.nan if sol0["status"] != "infeasible" else 0.0,
            "withdrawals": [],
            "knockouts": [],
            "fva": [],
        }
    g0 = float(sol0["growth_rate"])
    if g0 < 1e-6:
        return {
            "evaluation_status": "non_evaluable_zero_baseline_growth",
            "baseline_growth_rate": g0,
            "withdrawals": [],
            "knockouts": [],
            "fva": [],
        }

    # 1. Nutrient withdrawal + single-metabolite add-back rescue on all actively imported exchange reactions
    withdrawals: list[dict[str, Any]] = []
    inorganic_core = {"cpd00001", "cpd00009", "cpd00067", "cpd00971", "cpd00205", "cpd00099", "cpd00254", "cpd00063"}
    for rid, flux_val in sorted(sol0["fluxes"].items()):
        rxn = mod_primary.reactions.get(rid)
        if rxn is None or not rxn.is_exchange or flux_val >= -1e-6:
            continue
        m_ex = re.match(r"^R_EX_(cpd\d+)_e0$", rid) or re.match(r"^EX_(cpd\d+)_e0$", rxn.base_id)
        cpd_id = m_ex.group(1) if m_ex else rxn.base_id
        if cpd_id in inorganic_core:
            continue

        mod_w_prim = mod_primary.copy()
        orig_lb = mod_w_prim.reactions[rid].lb
        mod_w_prim.reactions[rid].lb = 0.0
        sol_w_prim = solve_fba_highs(mod_w_prim, time_limit_s=time_limit_s)
        gw_prim = float(sol_w_prim["growth_rate"]) if sol_w_prim["status"] in {"optimal", "infeasible"} else math.nan

        mod_w_uni = mod_union.copy()
        mod_w_uni.reactions[rid].lb = 0.0
        sol_w_uni = solve_fba_highs(mod_w_uni, time_limit_s=time_limit_s)
        gw_uni = float(sol_w_uni["growth_rate"]) if sol_w_uni["status"] in {"optimal", "infeasible"} else math.nan

        mod_addback = mod_w_prim.copy()
        mod_addback.reactions[rid].lb = orig_lb
        sol_addback = solve_fba_highs(mod_addback, time_limit_s=time_limit_s)
        g_addback = float(sol_addback["growth_rate"]) if sol_addback["status"] in {"optimal", "infeasible"} else math.nan

        if math.isnan(gw_prim) or math.isnan(gw_uni):
            dep_class = "unknown_solver_failure"
        elif gw_prim < 1e-4 and gw_uni < 1e-4:
            dep_class = "required_across_sampled_ensemble_union"
        elif gw_prim < 1e-4 and gw_uni >= min_growth * 0.5:
            dep_class = "model_dependent_bypassed_in_sampled_ensemble_union"
        else:
            dep_class = "non_essential_stimulatory_uptake"

        withdrawals.append(
            {
                "exchange_rxn_id": rid,
                "compound_id": cpd_id,
                "compound_name": rxn.name,
                "baseline_uptake_flux": float(-flux_val),
                "baseline_growth": g0,
                "withdrawn_growth_primary": gw_prim,
                "withdrawn_growth_ensemble_union": gw_uni,
                "addback_rescue_growth": g_addback,
                "addback_rescues_growth": bool(
                    not math.isnan(gw_prim)
                    and not math.isnan(g_addback)
                    and gw_prim < 1e-4
                    and g_addback >= min_growth * 0.99
                ),
                "dependency_classification": dep_class,
            }
        )

    # 2. Targeted reaction knockouts of primary gap-filled additions inside the sampled ensemble-union network
    knockouts: list[dict[str, Any]] = []
    for rid in primary_added:
        mod_ko_uni = mod_union.copy()
        mod_ko_uni.reactions[rid].lb = 0.0
        mod_ko_uni.reactions[rid].ub = 0.0
        sol_ko = solve_fba_highs(mod_ko_uni, time_limit_s=time_limit_s)
        g_ko = float(sol_ko["growth_rate"]) if sol_ko["status"] in {"optimal", "infeasible"} else math.nan
        if math.isnan(g_ko):
            ko_cls = f"unknown_solver_{sol_ko['status']}"
            ess_flag = None
        elif g_ko < 1e-4:
            ko_cls = "essential_in_sampled_ensemble_union"
            ess_flag = True
        else:
            ko_cls = "bypassable_in_sampled_ensemble_union"
            ess_flag = False
        knockouts.append(
            {
                "rxn_id": rid,
                "base_id": mod_union.reactions[rid].base_id,
                "baseline_growth": g0,
                "knockout_growth_in_ensemble_union": g_ko,
                "essential_across_ensemble_union": ess_flag,
                "knockout_class": ko_cls,
            }
        )

    # 3. Flux Variability Analysis (FVA) at >= 95% optimal growth on mod_union for primary_added reactions
    fva_rows: list[dict[str, Any]] = []
    if primary_added:
        sp_ids, rxn_ids, S, lb, ub, obj_j = build_sparse_stoichiometric_matrix(mod_union)
        rxn_idx = {r: j for j, r in enumerate(rxn_ids)}
        lb_fva = lb.copy()
        lb_fva[obj_j] = max(lb_fva[obj_j], 0.95 * g0)
        b_eq = np.zeros(len(sp_ids), dtype=np.float64)
        bounds_fva = list(zip(lb_fva, ub))

        for rid in primary_added[:30]:
            j_target = rxn_idx[rid]
            c_min = np.zeros(len(rxn_ids), dtype=np.float64)
            c_min[j_target] = 1.0
            r_min = linprog(c=c_min, A_eq=S, b_eq=b_eq, bounds=bounds_fva, method="highs", options={"time_limit": 5.0})
            v_min = float(r_min.x[j_target]) if (r_min.success and r_min.x is not None) else math.nan

            c_max = np.zeros(len(rxn_ids), dtype=np.float64)
            c_max[j_target] = -1.0
            r_max = linprog(c=c_max, A_eq=S, b_eq=b_eq, bounds=bounds_fva, method="highs", options={"time_limit": 5.0})
            v_max = float(r_max.x[j_target]) if (r_max.success and r_max.x is not None) else math.nan

            if math.isnan(v_min) or math.isnan(v_max):
                strictly_active = None
            else:
                strictly_active = bool(v_min > 1e-6 or v_max < -1e-6)
            fva_rows.append(
                {
                    "rxn_id": rid,
                    "fva_min_flux_at_95pct_opt": v_min,
                    "fva_max_flux_at_95pct_opt": v_max,
                    "fva_strictly_required": strictly_active,
                }
            )

    return {
        "evaluation_status": "evaluated_feasible_baseline",
        "baseline_growth_rate": g0,
        "withdrawals": withdrawals,
        "knockouts": knockouts,
        "fva": fva_rows,
    }


# =============================================================================
# 8. Artificial Genome-Incompleteness Masking Control (`EXP1-CONTRACT.txt` Item 2)
# =============================================================================

def apply_artificial_gene_masking(
    draft_model: SBMLModelData,
    mask_fraction: float,
    rng: np.random.Generator,
) -> tuple[SBMLModelData, dict[str, Any]]:
    """Randomly mask `mask_fraction` (e.g. 0.10, 0.20, 0.30) of sequence-supported genes in a high-quality
    draft model, propagate through GPR boolean trees (`AND`/`OR`), and return the pruned draft model
    plus the exact list of knocked-out reactions (which move into the unannotated candidate pool).
    """
    if not (0.0 < mask_fraction < 1.0):
        raise ValueError(f"mask_fraction must be in (0, 1), got {mask_fraction}")
    all_genes = sorted(draft_model.gene_products.keys())
    if not all_genes:
        gset: set[str] = set()
        for r in draft_model.reactions.values():
            gset |= r.genes
        all_genes = sorted(gset)

    n_total = len(all_genes)
    n_mask = max(1, int(round(n_total * mask_fraction)))
    masked_genes = set(rng.choice(all_genes, size=min(n_mask, n_total), replace=False).tolist())
    active_genes = set(all_genes) - masked_genes

    masked_model = draft_model.copy()
    removed_rxns: list[str] = []
    for rid, rxn in list(masked_model.reactions.items()):
        if rxn.is_biomass or is_single_sided_boundary_reaction(rxn):
            continue
        if rxn.gpr_ast is not None and rxn.genes:
            if not evaluate_gpr_ast(rxn.gpr_ast, active_genes):
                removed_rxns.append(rid)
                del masked_model.reactions[rid]

    return masked_model, {
        "mask_fraction": float(mask_fraction),
        "n_total_genes": n_total,
        "n_masked_genes": len(masked_genes),
        "n_removed_gpr_reactions": len(removed_rxns),
        "removed_rxn_ids": sorted(removed_rxns),
    }


# =============================================================================
# 9. Deterministic Candidate Reaction Pool Aggregation & Conflict Audit
#    (`EXP1-INVENTORY-REVIEW.txt`)
# =============================================================================

OBSERVED_COHORT_538_POOL_SCOPE = "observed_cohort_538_sbml_union_pool_not_full_modelseed_db"


def _normalize_stoichiometry_signature(stoich: dict[str, float]) -> tuple[tuple[str, float], ...]:
    """Return a canonical sorted tuple `(species_id, rounded_coeff)` for non-zero stoichiometric coefficients."""
    items = [(sp_id, round(float(coeff), 6)) for sp_id, coeff in stoich.items() if abs(float(coeff)) > 1e-12]
    items.sort(key=lambda x: x[0])
    return tuple(items)


def _extract_compartment_signature(stoich: dict[str, float]) -> tuple[str, ...]:
    comps: set[str] = set()
    for sp_id, coeff in stoich.items():
        if abs(float(coeff)) <= 1e-12:
            continue
        m = re.search(r"_([cep]0)$", sp_id)
        comps.add(m.group(1) if m else "c0")
    return tuple(sorted(comps))


def aggregate_candidate_pool_deterministically(
    per_genome_extractions: list[dict[str, Any]],
    seed_rxn_directions: dict[str, str] | None = None,
    seed_metabolite_metadata: dict[str, tuple[str | None, int | None]] | None = None,
    excluded_base_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Deterministically aggregate candidate reactions and species across all cohort models (`EXP1-INVENTORY-REVIEW.txt`):
    - Invariant to worker completion order (`as_completed`) by sorting donor records by `genome` first.
    - Tracks exact donor provenance (`donor_draft_models`, `donor_fill_desert_models`, `donor_fill_native_models`,
      and `pool_source_tier`: `'observed_in_cohort_draft_models'` vs `'observed_only_in_cohort_gapfilled_models'`).
    - Verifies stoichiometry + compartment equality across all donor models for each reaction ID; if any donor model
      carries a conflicting stoichiometry or compartment signature for the same reaction ID, records the conflict in
      `reaction_conflicts` and excludes the conflicting reaction from the eligible candidate pool (`eligible_reactions`)
      so duplicate-ID mismatches NEVER silently pass.
    - Derives canonical reaction direction and flux bounds from frozen gapseq `seed_reactions_corrected.tsv`
      (`seed_rxn_directions`: `'='` / `'<=>'` -> `[-1000, 1000]`, `'>'` / `'=>'` -> `[0, 1000]`, `'<'` / `'<='` -> `[-1000, 0]`)
      when available, or else from the deterministic conservative union of observed donor bounds across all sorted donors.
    - Audits mass/charge balance of every candidate reaction (`audit_reaction_mass_charge_balance`); known-imbalanced
      internal candidate reactions (`balance_status == 'elemental_imbalance'`) are flagged `eligible_for_mechanism_gapfill = False`
      and excluded from `eligible_reactions` (while `unknown_formula_or_charge` is explicitly labeled and retained).
    """
    seed_dirs = seed_rxn_directions or {}
    seed_mets = seed_metabolite_metadata or {}
    excl_ids = excluded_base_ids or set()

    # Sort deterministically by genome ID so arrival order from as_completed never affects aggregation
    sorted_records = sorted(per_genome_extractions, key=lambda r: str(r["genome"]))

    canonical_species: dict[str, SpeciesRecord] = {}
    species_conflicts: list[dict[str, Any]] = []

    # 1. Deterministic species aggregation
    for rec in sorted_records:
        g = str(rec["genome"])
        sp_map: dict[str, SpeciesRecord] = rec.get("extracted_species", {})
        for sid in sorted(sp_map.keys()):
            srec = sp_map[sid]
            cpd_clean = srec.cpd_id
            seed_f, seed_q = seed_mets.get(cpd_clean, (None, None))
            f_val = seed_f if seed_f is not None else srec.formula
            q_val = seed_q if seed_q is not None else srec.charge

            if sid not in canonical_species:
                canonical_species[sid] = SpeciesRecord(
                    id=sid,
                    cpd_id=cpd_clean,
                    name=srec.name,
                    compartment=srec.compartment,
                    charge=q_val,
                    formula=f_val,
                )
            else:
                prev_sp = canonical_species[sid]
                if prev_sp.formula is None and f_val is not None:
                    prev_sp.formula = f_val
                elif prev_sp.formula is not None and f_val is not None and prev_sp.formula != f_val:
                    species_conflicts.append(
                        {
                            "species_id": sid,
                            "conflict_type": "formula_mismatch",
                            "prev_value": prev_sp.formula,
                            "new_value": f_val,
                            "donor_genome": g,
                        }
                    )
                if prev_sp.charge is None and q_val is not None:
                    prev_sp.charge = q_val
                elif prev_sp.charge is not None and q_val is not None and int(prev_sp.charge) != int(q_val):
                    species_conflicts.append(
                        {
                            "species_id": sid,
                            "conflict_type": "charge_mismatch",
                            "prev_value": prev_sp.charge,
                            "new_value": q_val,
                            "donor_genome": g,
                        }
                    )

    # 2. Group reaction observations deterministically by reaction ID
    rxn_obs: dict[str, dict[str, Any]] = {}
    for rec in sorted_records:
        g = str(rec["genome"])
        draft_rids: set[str] = set(rec.get("draft_reaction_ids", set()))
        des_rids: set[str] = set(rec.get("fill_desert_reaction_ids", set()))
        nat_rids: set[str] = set(rec.get("fill_native_reaction_ids", set()))
        rxn_map: dict[str, ReactionRecord] = rec.get("extracted_reactions", {})

        for rid in sorted(rxn_map.keys()):
            rrec = rxn_map[rid]
            if rrec.is_biomass:
                continue
            clean_base = rrec.base_id.split("_")[0] if rrec.base_id.startswith("rxn") else rrec.base_id
            if rrec.base_id in excl_ids or clean_base in excl_ids:
                continue

            stoich_sig = _normalize_stoichiometry_signature(rrec.stoichiometry)
            comp_sig = _extract_compartment_signature(rrec.stoichiometry)

            entry = rxn_obs.get(rid)
            if entry is None:
                entry = {
                    "template_rxn": ReactionRecord(**vars(rrec)),
                    "stoich_sig": stoich_sig,
                    "comp_sig": comp_sig,
                    "has_conflict": False,
                    "conflicts": [],
                    "donor_draft": [],
                    "donor_desert": [],
                    "donor_native": [],
                    "observed_min_lb": float(rrec.lb),
                    "observed_max_ub": float(rrec.ub),
                    "observed_bounds_set": {(round(float(rrec.lb), 4), round(float(rrec.ub), 4))},
                }
                rxn_obs[rid] = entry
            else:
                if stoich_sig != entry["stoich_sig"] or comp_sig != entry["comp_sig"]:
                    entry["has_conflict"] = True
                    entry["conflicts"].append(
                        {
                            "rxn_id": rid,
                            "donor_genome": g,
                            "expected_stoich": str(entry["stoich_sig"]),
                            "observed_stoich": str(stoich_sig),
                            "expected_compartments": ",".join(entry["comp_sig"]),
                            "observed_compartments": ",".join(comp_sig),
                        }
                    )
                entry["observed_min_lb"] = min(entry["observed_min_lb"], float(rrec.lb))
                entry["observed_max_ub"] = max(entry["observed_max_ub"], float(rrec.ub))
                entry["observed_bounds_set"].add((round(float(rrec.lb), 4), round(float(rrec.ub), 4)))

            if rid in draft_rids:
                entry["donor_draft"].append(g)
            if rid in des_rids:
                entry["donor_desert"].append(g)
            if rid in nat_rids:
                entry["donor_native"].append(g)

    eligible_reactions: dict[str, ReactionRecord] = {}
    pool_audit_rows: list[dict[str, Any]] = []
    reaction_conflicts: list[dict[str, Any]] = []

    for rid in sorted(rxn_obs.keys()):
        info = rxn_obs[rid]
        tmpl: ReactionRecord = info["template_rxn"]
        clean_base = tmpl.base_id.split("_")[0] if tmpl.base_id.startswith("rxn") else tmpl.base_id
        reaction_conflicts.extend(info["conflicts"])

        # Canonical direction & bounds resolution (`EXP1-INVENTORY-REVIEW.txt` Item 2):
        # 1) Exchange reactions: standard boundary convention (`M_e0 -> ø`, `lb = 0.0, ub = 1000.0`, opened only by scenario)
        # 2) If clean_base is in frozen gapseq `seed_reactions_corrected.tsv` (`seed_dirs`) AND all donors agree on bounds,
        #    or if we take the conservative union of frozen gapseq direction and observed donor directions:
        if tmpl.is_exchange:
            canon_lb, canon_ub = 0.0, 1000.0
            dir_source = "exchange_boundary_standard_closed_uptake"
        elif clean_base in seed_dirs and not tmpl.is_transport:
            s_dir = seed_dirs[clean_base]
            if s_dir in {"=", "<=>"}:
                canon_lb, canon_ub = -1000.0, 1000.0
            elif s_dir in {"<", "<="}:
                # Note: if donor SBML models also allow forward flux (`observed_max_ub > 0`), preserve conservative union
                canon_lb = -1000.0
                canon_ub = 1000.0 if info["observed_max_ub"] > 1e-9 else 0.0
            else:
                canon_lb = -1000.0 if info["observed_min_lb"] < -1e-9 else 0.0
                canon_ub = 1000.0
            dir_source = f"gapseq_seed_reactions_corrected_tsv({s_dir})_with_donor_union"
        else:
            canon_lb = -1000.0 if info["observed_min_lb"] < -1e-9 else 0.0
            canon_ub = 1000.0 if info["observed_max_ub"] > 1e-9 else 0.0
            dir_source = "deterministic_cohort_donor_union"

        canon_rxn = ReactionRecord(
            id=tmpl.id,
            base_id=tmpl.base_id,
            name=tmpl.name,
            reversible=(canon_lb < 0.0 and canon_ub > 0.0),
            lb=float(canon_lb),
            ub=float(canon_ub),
            stoichiometry=dict(info["stoich_sig"]),
            gpr_ast=None,
            genes=set(),
            is_exchange=tmpl.is_exchange,
            is_biomass=False,
            is_transport=tmpl.is_transport,
            provenance=(
                "observed_in_cohort_draft_models"
                if len(info["donor_draft"]) > 0
                else "observed_only_in_cohort_gapfilled_models"
            ),
        )

        bal = audit_reaction_mass_charge_balance(canon_rxn, canonical_species)
        bal_status = str(bal["balance_status"])

        # Known elemental-imbalanced internal reactions or stoichiometry-conflicting reactions are NOT eligible for mechanism gap-filling!
        has_stoich_conflict = bool(info["has_conflict"])
        is_known_imbalanced = (bal_status == "elemental_imbalance") and not canon_rxn.is_exchange
        is_zero_bounds = (abs(canon_rxn.lb) <= 1e-12 and abs(canon_rxn.ub) <= 1e-12)
        eligible = (not has_stoich_conflict) and (not is_known_imbalanced) and (not is_zero_bounds)

        if eligible:
            eligible_reactions[rid] = canon_rxn

        pool_audit_rows.append(
            {
                "rxn_id": rid,
                "base_id": canon_rxn.base_id,
                "clean_rxn_id": clean_base,
                "name": canon_rxn.name,
                "is_exchange": canon_rxn.is_exchange,
                "is_transport": canon_rxn.is_transport,
                "canonical_lb": canon_rxn.lb,
                "canonical_ub": canon_rxn.ub,
                "direction_policy_source": dir_source,
                "n_distinct_donor_bound_pairs": len(info["observed_bounds_set"]),
                "has_stoichiometry_or_compartment_conflict": has_stoich_conflict,
                "balance_status": bal_status,
                "net_charge": bal["net_charge"],
                "eligible_for_mechanism_gapfill": eligible,
                "exclusion_reason": (
                    "stoichiometry_or_compartment_conflict"
                    if has_stoich_conflict
                    else ("known_elemental_imbalance" if is_known_imbalanced else ("zero_bounds" if is_zero_bounds else ""))
                ),
                "pool_source_tier": canon_rxn.provenance,
                "n_donor_draft_models": len(info["donor_draft"]),
                "n_donor_fill_desert_models": len(info["donor_desert"]),
                "n_donor_fill_native_models": len(info["donor_native"]),
                "example_donor_draft_models_csv": ",".join(info["donor_draft"][:5]),
                "example_donor_fill_desert_models_csv": ",".join(info["donor_desert"][:5]),
            }
        )

    return {
        "candidate_pool_scope": OBSERVED_COHORT_538_POOL_SCOPE,
        "reactions": eligible_reactions,
        "species": canonical_species,
        "pool_audit_rows": pool_audit_rows,
        "reaction_conflicts": reaction_conflicts,
        "species_conflicts": species_conflicts,
    }

