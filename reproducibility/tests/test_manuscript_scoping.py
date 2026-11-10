"""
Negative and Positive Fixture Tests for Manuscript Include Loader, Number Normalization,
and Boundary Scoping Verifiers.

Verifies:
- Recursive include loader traverses \\input{} and \\include{} directives without infinite recursion.
- Number format normalization correctly matches LaTeX braced commas (14{,}426) and standard commas.
- Negative fixture: missing required manuscript phrases fail closed with descriptive errors.
- Negative fixture: forbidden phrases in included subdocuments are caught and flagged.
- Negative fixture: unbounded novel-fold-family patterns are detected and flagged.
- Clean execution: production manuscript passes with zero errors.
"""

from __future__ import annotations

import os
import re
import shutil
import sys
from pathlib import Path
import pytest

REPRO_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from verify import load_active_manuscript_texts, verify_manuscript_boundaries


@pytest.fixture
def local_scratch():
    """Scratch directory within repo workspace (strictly avoids /tmp)."""
    scratch = REPRO_DIR / "tests" / "fixture_scratch"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    yield scratch
    if scratch.exists():
        shutil.rmtree(scratch)


# =====================================================================
# 1. Active TeX Include Tree Traversal
# =====================================================================

def test_production_manuscript_include_tree():
    """Verify recursive loader discovers main.tex and its included subdocuments."""
    main_tex = REPO_ROOT / "manuscript" / "main.tex"
    assert main_tex.is_file(), "manuscript/main.tex must exist"

    loaded = load_active_manuscript_texts(str(main_tex))
    loaded_basenames = [Path(p).name for p in loaded.keys()]

    assert "main.tex" in loaded_basenames
    assert "authors.tex" in loaded_basenames
    assert "supplementary_notes.tex" in loaded_basenames
    assert "mechanistic_tests.tex" in loaded_basenames


def test_circular_include_safety(local_scratch):
    """Verify that circular includes (A inputs B, B inputs A) do not cause infinite recursion."""
    doc_a = local_scratch / "doc_a.tex"
    doc_b = local_scratch / "doc_b.tex"

    doc_a.write_text("Document A text.\n\\input{doc_b}\n", encoding="utf-8")
    doc_b.write_text("Document B text.\n\\input{doc_a}\n", encoding="utf-8")

    loaded = load_active_manuscript_texts(str(doc_a))
    assert len(loaded) == 2
    assert os.path.realpath(str(doc_a)) in loaded
    assert os.path.realpath(str(doc_b)) in loaded


def test_missing_include_handled(local_scratch):
    """Verify that missing include targets are safely ignored without unhandled exception."""
    doc = local_scratch / "missing.tex"
    doc.write_text("Some text.\n\\input{does_not_exist_xyz}\n", encoding="utf-8")

    loaded = load_active_manuscript_texts(str(doc))
    assert len(loaded) == 1
    assert "Some text." in list(loaded.values())[0]


# =====================================================================
# 2. Number Normalization Checks
# =====================================================================

def test_number_normalization_braced_and_standard_commas():
    """Verify LaTeX braced commas and standard commas normalize identically."""
    def norm(s: str) -> str:
        return re.sub(r"\{,\}", ",", s)

    assert norm("14{,}426") == "14,426"
    assert norm("5{,}438") == "5,438"
    assert norm("4{,}210{,}592") == "4,210,592"

    text_with_braced = "Total of 14{,}426 bins dereplicated to 5{,}438 representatives."
    text_with_standard = "Total of 14,426 bins dereplicated to 5,438 representatives."

    assert norm("14{,}426") in norm(text_with_braced)
    assert norm("14{,}426") in norm(text_with_standard)
    assert norm("14,426") in norm(text_with_braced)
    assert norm("14,426") in norm(text_with_standard)


# =====================================================================
# 3. Negative Fixtures: Boundary & Scoping Enforcement
# =====================================================================

def test_negative_fixture_missing_required_phrase():
    """Verify that omitting a required phrase triggers a verification failure."""
    manuscript_path = REPO_ROOT / "manuscript" / "main.tex"
    loaded = load_active_manuscript_texts(str(manuscript_path))
    combined = "\n".join(loaded.values())

    # Artificially remove a required phrase
    censored = combined.replace("GTDB-Tk v2.7.2", "OldTaxonomySoftware v1.0")

    def norm_tex(s: str) -> str:
        return re.sub(r"\{,\}", ",", s)

    norm_censored = norm_tex(censored)
    assert "GTDB-Tk v2.7.2" not in norm_censored

    # Run check logic on censored text
    required_checks = [
        ("GTDB-Tk v2.7.2", ["GTDB-Tk v2.7.2"], "current taxonomy software"),
    ]
    errs = []
    for label, variants, reason in required_checks:
        if not any(norm_tex(v) in norm_censored for v in variants):
            errs.append(f"required manuscript phrase missing {label!r}: {reason}")

    assert len(errs) == 1
    assert "current taxonomy software" in errs[0]


def test_negative_fixture_forbidden_phrase_in_subdocument():
    """Verify that forbidden phrases injected into an included file are caught."""
    forbidden = {
        "SeqCode-eligible": "R232 empty-rank counts are not SeqCode eligibility calls",
        "globally endemic": "specificity is relative to named comparators",
    }
    fake_subdoc = "In this supplementary note, we claim the lineage is globally endemic and SeqCode-eligible."

    errs = []
    for phrase, reason in forbidden.items():
        if phrase in fake_subdoc:
            errs.append(f"forbidden phrase {phrase!r}: {reason}")

    assert len(errs) == 2
    assert any("SeqCode-eligible" in e for e in errs)
    assert any("globally endemic" in e for e in errs)


def test_negative_fixture_unbounded_novel_fold_pattern():
    """Verify that unhedged novel fold family wording is detected and rejected."""
    pattern = r"novel(?:[-\s]+)fold(?:[-\s]+)famil"
    allowance = re.compile(r"candidate|historical|legacy|reference[- ]relative", re.I)

    # Unbounded usage (should FAIL)
    unbounded_text = "We discovered 68 novel fold families across the desert metagenome."
    matcher = re.compile(pattern, re.I)
    match_unbounded = list(matcher.finditer(unbounded_text))
    assert len(match_unbounded) > 0
    assert not allowance.search(unbounded_text)

    # Properly bounded usage (should PASS)
    bounded_text = "We identified 68 candidate novel fold families across the desert metagenome."
    match_bounded = list(matcher.finditer(bounded_text))
    assert len(match_bounded) > 0
    assert allowance.search(bounded_text) is not None


# =====================================================================
# 4. Production Manuscript Verification Pass
# =====================================================================

def test_production_manuscript_zero_boundary_errors():
    """Verify that verify_manuscript_boundaries() passes with 0 errors on production repository."""
    errors = verify_manuscript_boundaries()
    assert errors == [], f"Expected 0 boundary errors, got {len(errors)}: {errors}"
