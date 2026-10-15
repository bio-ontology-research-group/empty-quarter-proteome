# Community metabolic interdependence in the Empty Quarter (pilot, 16 communities)

GEMs: gapseq on HQ single-sample MAGs. **No commercial MILP solver is available on this cluster (CPLEX/Gurobi access-restricted), so SMETANA *detailed* mode (MILP) is intractable with GLPK at these community sizes.** We therefore report (i) SMETANA *global* MRO (computed) and (ii) an LP-based FBA cross-feeding-potential analysis for the guild-as-feeder test.

## (1) Resource overlap (SMETANA global MRO)

MRO median 0.69 (range 0.55-0.78) across 16 communities — moderate-to-high metabolic resource overlap (a competition signal typical of oligotrophic resource sharing). MIP returned n/a (non-interacting baseline infeasible for communities containing non-growing gap-filled models).

## (2) Are the atmospheric guild the metabolic donors? (FBA cross-feeding potential)

Donor potential = # (receiver,metabolite) organic cross-feeds a MAG can supply on a complete environment. Guild MAGs (n=27) vs non-guild (n=289):
- donor potential: guild median 232.0 vs non-guild 243.0; one-sided MWU (guild>non-guild) p=0.694.
- donor FRACTION: guild 0.503 vs non-guild 0.495; one-sided MWU p=0.387.

## CONTROL — is it a MAG-incompleteness artefact?

receiver potential vs completeness Spearman rho=0.03 (n=316); donor vs completeness rho=0.06. All MAGs are HQ (>=90%), and |rho|<~0.3 would support a biological (not completeness-driven) signal.

## Caveats
- Complete-environment FBA => an UPPER BOUND on secretion/cross-feeding potential; a defined desert-minimal medium would tighten it.
- Potential (capability) cross-feeding, not realised flux or MILP-scored SMETANA cross-feeding.
- GEMs see only decodable metabolism (dark proteome invisible); pilot = 16/547 samples.
- Commercial-solver access would enable full SMETANA detailed + scale-up (flagged for revision).
