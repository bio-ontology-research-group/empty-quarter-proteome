# Independent Claude Code Opus follow-up audit — pre-final repair

Date: 2026-08-04

Scope: read-only audit of the MAG/dark-matter manuscript, RQ records, current
structure-calibration artifacts, claims verifier, and reproducibility deposit after
the first repair pass. The agent was launched with `claude -p --model opus
--effort high`; it reported model `claude-opus-5`, 67 turns, 42,617 output tokens,
and no file modifications.

## Verdict at audit time: FAIL (1 blocker, 4 major, 10 minor)

### Blocker

- **BL-2 — the corrected deposit and manuscript submodule were not committed.**
  A clean committed checkout therefore still resolved to the superseded manuscript
  and omitted required compact artifacts and curated tables.

### Major findings

- **MJ-A — candidate no-hit decomposition was incomplete.** The manuscript did not
  split no-row retrieval failures from retrieved-below-threshold cases before and
  after the confidence/length gate.
- **MJ-B — ungated known-fold control rates were omitted.** This hid the strong
  confidence dependence of the retrieval/error floor.
- **MJ-C — legacy RQ prose retained overclaims and stale R226/H3 statements.** These
  contradicted the bounded R232 manuscript even though the headline text had been
  repaired.
- **MJ-D — build instructions named Biber while the manuscript actually uses
  BibTeX.** The documented clean-build path was therefore wrong.

## Disposition in the consolidated revision

- The manuscript submodule is committed first, and the parent commit records that
  exact pointer together with every compact artifact and required curated input.
- Both TM-score normalisations now report the no-row/retrieved-below-threshold split
  for candidates and controls, gated and ungated. Ungated control errors are
  440/3,684 (11.94%) for query-normalised TM and 236/3,684 (6.41%) for
  alignment-normalised TM; the qualified rates are 101/3,013 (3.35%) and
  68/3,013 (2.26%).
- Active RQ prose is bounded to named references and the current R232 catalogue;
  recursive claim-language guards cover `results/**/*.md`, `RQ/**/*.md`, and the
  RQ generator.
- `reproducibility/build_manuscript.sh` and all documentation use
  Latexmk/BibTeX and fail on unresolved citations or references.
- Historical structure counts now close through a checksum-pinned exact
  claim/value/source ledger. Figure F5c uses a deterministic, checksum-pinned
  15,469-MAG taxonomy snapshot instead of a workstation-only workbench symlink.
- The verifier recomputes dN/dS to the reported three decimals, reads MGnify counts
  from the deposited summary table rather than a hard-coded constant, regenerates
  all eleven figures, and builds the manuscript.

This report records the independent audit that drove the second repair pass. A
separate final Opus audit is run against the committed clean-checkout state.
