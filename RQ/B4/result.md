# RQ-B4 — Is the vanadium → alternative-nitrogenase link real?

**Status:** ✅ answered at the family level (negative) — via RQ-B3; a finer marker-abundance test remains optional 🟡

## Answer
**No support.** Two independent lines of evidence say the vanadium→alternative-nitrogenase idea is not borne out:
1. **Vanadium is null in the family-level PERMANOVA** of dark-proteome composition (V block pseudo-F≈1.03, p=0.325; `results/darkmatter_ecology/phase_b_summary.txt`).
2. **The systematic gene↔geochemistry map (RQ-B3) finds that V recruits *zero* dark families** at FDR q<0.05 — as do the other candidate metalloenzyme cofactors Ni, Cu, Co (Mo recruits only 177, all positive). The dark proteome is structured by the **water/salt-chemistry gradient** (Ca 4066, S 3964, Cl 3594, salinity 3179 families), not by trace-metal cofactors. See `RQ/B3/gene_env_summary.tsv`.

So at the resolution of dark-family abundance, there is **no vanadium–metabolism coupling**, consistent with vanadium-nitrogenase being a minor/contextual rather than geochemically-driven trait here.

**Remaining (optional, 🟡):** a targeted test of `vnf`/`anf` *marker* abundance vs measured `V` across the 251 XRF sites would confirm this directly at the nitrogenase-gene level (needs per-sample vnf/anf abundance extraction from the workbench; ~2 d). Expected to remain negative.

**Evidence / reproduce:** `RQ/B3/run.py` (V column → 0 significant families); `results/darkmatter_ecology/phase_b_summary.txt` (PERMANOVA vanadium block).

**References:** the alternative-nitrogenase context; see also `RQ/B3/result.md`.
