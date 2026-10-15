# Historical registry reconciliation

The historical verifier was actually run. Its numerical/source checks report:
62 COMPUTE pass,72 FROZEN pass,6 TABLE sources present and25 PROVENANCE entries
ledger-pinned. It still exits nonzero on31text-boundary gates. These results do
not validate the current manuscript's numbers automatically: some passed recipes
refer to older cohorts. No gates were disabled to manufacture a pass.

Repository guidance in `CLAUDE.md` reserves changes to `claims.tsv` for coordinated
rebinding and gives `notes/CORRECTIONS.md` precedence in the interim. The review's
49 current bindings, exact inputs and regression guards are implemented separately
in `reproducibility/darkmatter/`. The following is the concrete reconciliation
needed for the historical gate layer.

| Old gate group | Current disposition/evidence |
|---|---|
| RQ/SPINE_DOSSIER novel-floor wording (2 gates) | Superseded planning document, not the current submitted manuscript; preserve historical status rather than importing its claims. |
| dRep12,433-input phrase | Methods reports12,434eligible bins and one internal drop, yielding5,438representatives; wording-level gate needs rebinding. |
| GTDB-Tk2.7.2 phrase | Current Methods states v2.7.2; the literal matcher rejects the inserted v. |
| Structural608/1,875;97/1,875;74no-row;534retrieved;74/23 (5 gates) | Superseded candidate cohort. Current all-qualified counts2,838, query failures1,031, aligned failures95, no-row22, retrieved failures1,009/73; raw-hit audit53679366. |
| Control101/3,013;68/3,013 (2 gates) | Current values agree, but prose uses “of” rather than slash. No-row66separately identified. |
| Unfiltered41,127/16,386 and440/3,684/236/3,684 (4 gates) | Historical unfiltered-cohort claims are not in the narrowed manuscript; retain as historical provenance, not active text requirements. |
| Low-MSA order candidate | Current Methods excludes the eleventh candidate with46%marker coverage from adjudication. Same boundary, different wording. |
|31.41PICRUSt2diagnostic | Not reported in the current paired-profile comparison; current recall/correlation figures are bound to retained per-sample tables. |
|184named species | Separate named-label and assigned-rank rules; current census187assigned-through-species and122without an empty/digit-bearing label. Do not conflate the rules. |
|151,400historical tier snapshot | Historical non-closure, no current claim on that basis. |
| Named MGnify90/GEM/five-desert wording | Current comparison includes SPIRE and names catalogues explicitly; bind the actual search cohort. |
|51named phyla | Current text correctly says51GTDBphylum labels, without treating every label as a named/validated taxon. |
|2,839,030dark count;22.2%length fraction | Historical dark-set basis; current genome-resolved4,292,446,65.8%under100aa and56.4%partial. |
|4,482median KO count;56.7%no-KO | Current text uses the protein-denominator ceiling55.6%; historical median not reported. |
|45.6%occupancy;21.4%reference no-hit | Superseded scope/search. Current combined-catalogue no-hit30,479/106,042=28.74%; occupancy reported by61locations and length stratum. |
|64.5%terpene | Not an active quantitative claim in the narrowed draft. |
|median0.391dN/dS | Current97,200-family median0.4253, with frame/ascertainment/divergence qualifications and flagged-family calibration. |
|23/332pocket and positive-selection warning (2 gates) | Historical pocket/selection block is outside this paper's current scope. |

This document records the proposed old-to-current reconciliation; it is not an
assertion that the canonical registry has been rebound or that every manuscript
claim is independently reproduced. The compact revision checks remain explicitly
scoped, and the original31failures are retained in `legacy_verify.log`.
