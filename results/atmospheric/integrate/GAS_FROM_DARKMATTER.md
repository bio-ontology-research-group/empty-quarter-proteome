# How much functional dark matter is recoverable as atmospheric trace-gas scavenging? — ≈ none

Generated 2026-06-16 (`scripts/ibex/82h_darkmatter_recover_gas.sh`, jobs 47563356 strict / 47563441 broad).
Searches the dark-matter protein cohorts against the atmospheric-enzyme references to ask whether any
*divergent* coxL/hydB that reference annotation (EggNOG/UniRef) missed can be rescued by targeted search:
HydDB DIAMOND (group, incl. 1h/1l) + HydDB NiFe/FeFe HMMs (divergent detector) + Cordero form-I CoxL refs +
AYRGAGR active-site motif. Length-stratified (a real hydrogenase LSU ≈500-600 aa, CoxL ≈650-800 aa).

## Result — across both dark tiers
| cohort | size | hydrogenase (any) | HMM-supported | high-affinity 1h/1l | form-I CoxL | AYRGAGR | **high-conf total** |
|---|---|---|---|---|---|---|---|
| **reallydark** (no UniRef50 **and** no AFDB) | 468,211 | 1 | 0 | 0 | 0 | 0 | **0** |
| **isdark** (no UniRef50; fold may be known) | 563,619 | 6 | 1 | 0 | 0 | 0 | **1** |

**≈ 0% of the functional dark matter is recoverable as a trace-gas enzyme.** Even the non-tautological tests
that detect *divergent* enzymes — the HydDB **HMM** and the **AYRGAGR motif** — return essentially zero
(1 borderline divergent hydrogenase in the entire 0.56 M sequence-novel set; **no** high-affinity 1h/1l, **no**
form-I CODH, **no** active-site motif).

## Interpretation — two independent axes
The atmospheric guild and the structural dark matter are **disjoint**:
1. **"Life on air" runs on conventional, reference-annotatable enzymes.** The 13,675 coxL + 6,602 hydB were found
   *because* they have homology (EggNOG/KO). The trace-gas machinery is **not** novel — it is well-conserved and
   reference-represented even in this never-sampled biome. (Consistent with the operon result: gas operons are
   characterizable — `OPERON_DARKMATTER.md`.)
2. **The dark matter is doing something else.** It is not cryptic/divergent gas metabolism. Whatever functions the
   ~66% EC-unassignable dark families encode, atmospheric trace-gas scavenging is not among them.

## Why this is useful for the paper
- It **cleanly separates the two headline strands** (atmospheric guild vs functional dark matter) so neither is
  conflated with the other — each stands on its own evidence.
- It is a **negative reference-bias result with a twist**: one might expect an extreme, unsampled desert to harbour
  divergent trace-gas enzymes invisible to reference databases. It does **not** — the gas enzymes are conserved.
  Reference bias in this biome bites elsewhere (the ~66% unassignable dark matter), not in the trace-gas guild.

Files: `dark_recover_gas_summary.txt` / `.tsv` (reallydark), `isdark_recover_gas_summary.txt` / `.tsv` (isdark).
