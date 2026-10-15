# RQ-B3 — Result: which variables select which dark protein families

## Headline answer
The Empty Quarter dark-protein repertoire is organised by a **water/salt–chemistry gradient,
not by metal cofactors**. Salinity chemistry (Ca, S, Cl, the salinity_proxy, Na) and aridity
recruit thousands of dark families each; the trace-metal cofactor axes that motivated the
"life-on-air" metalloenzyme hypotheses — **V, Ni, Cu, Co — recruit zero** dark families, and Mo
only a modest 177. For every salt/aridity axis the response is **overwhelmingly negative**:
rising stress suppresses most of the dark functional repertoire and positively recruits only a
small specialist set headed by phosphorus-scavenging (alkaline phosphatase) and osmolyte
phosphatases. The associations are mostly genuine chemistry: **57–84 % survive conditioning on
geography**.

## Are the indicator files the full result?
**No.** Each `indicators/<VAR>_top.tsv` holds exactly **40 rows** — the top hits only (the
generating script writes `TOP_N = 40`). The full per-axis q<0.05 counts live in the narrative
`INDICATORS.md`/`PARTIAL.md`. `run.py` therefore recomputes the full per-family Spearman from
the raw abundance matrix; its q<0.05 counts **reproduce `INDICATORS.md` exactly**, and its
space-partial counts and retained % **reproduce `PARTIAL.md` exactly** — so the synthesis below
is verified, not merely transcribed.

## Per-variable map (q<0.05; full table in `gene_env_summary.tsv`)

| variable | n_assoc q05 | n_pos | n_neg | q10 | strongest family (signed r, label) |
|---|---:|---:|---:|---:|---|
| Ca | **4066** | 323 | 3743 | 4973 | −0.46 PWY0-1337 (fatty-acid β-ox) |
| S | **3964** | 416 | 3548 | 4951 | +0.43 alkaline phosphatase |
| Cl | **3594** | 407 | 3187 | 4568 | +0.42 alkaline phosphatase |
| salinity_proxy | **3179** | 432 | 2747 | 4105 | −0.43 PRA-succinocarboxamide synthase |
| Fe | 1644 | 66 | 1578 | 2562 | −0.36 PWY0-1337 (fatty-acid β-ox) |
| K | 955 | 419 | 536 | 1652 | −0.37 PWY-6728 / +0.34 |
| Na | 705 | 128 | 577 | 1228 | −0.34 allophanate hydrolase |
| aridity_index | 491 | 58 | 433 | 832 | −0.36 adenosine nucleosidase |
| Mo | 177 | 177 | 0 | 324 | +0.30 PRGAM formyltransferase |
| Zn | 1 | 1 | 0 | 2 | +0.29 aspartate aminotransferase |
| Mn | 1 | 1 | 0 | 1 | +0.29 PWY-5941 |
| V | 0 | 0 | 0 | 0 | — |
| Ni | 0 | 0 | 0 | 0 | — |
| Cu | 0 | 0 | 0 | 0 | — |
| Co | 0 | 0 | 0 | 0 | — |

## Cross-variable synthesis

**Which axes recruit the most / fewest families.** The desert master variables dominate:
calcium (4,066), sulfur (3,964), chloride (3,594) and the salinity_proxy (3,179) each select a
quarter to a third of all tested families, and these four are strongly collinear (they recruit
overlapping sets — counts are not independent). Iron (1,644) and the macronutrient cations
K (955) / Na (705) form a middle tier; **aridity_index** recruits a more focused 491. At the
bottom the **trace-metal cofactor axes are essentially flat — V, Ni, Cu, Co all zero, Zn and Mn
one each** — and **Mo** a modest 177. The vanadium–nitrogenase hypothesis therefore gets **no
per-family support**: no dark family tracks the vanadium gradient at any FDR threshold.

**Direction — selection is mostly a filter, not a recruiter.** For every salt/Ca/aridity/Fe
axis the great majority of associations are **negative** (Ca 92 %, Fe 96 %, S/Cl ≈ 89 %,
salinity 86 %, aridity 88 % of significant families have r<0). Rising salinity, calcite and
aridity act as an **environmental filter that suppresses most of the dark repertoire**, while a
minority is positively recruited. The exceptions are Mo (100 % positive) and K (balanced,
44 % positive) — the macronutrient/cofactor axes behave as enrichers rather than filters.

**What the recruited specialists are.** The top positively-recruited families at S, Cl and
salinity_proxy are the **same family — an alkaline phosphatase** (EC 3.1.3.1, r ≈ +0.42–0.43),
i.e. phosphorus scavenging rises with salinity; the curated **osmolyte phosphatases** (e.g.
glucosylglycerol-P phosphatase, which also tops the host_Tribulus contrast) recur among the
positives, with 34–44 osmolyte-EC families significant on each salt axis. This is the expected
desiccation/salt-solute signature.

**How interpretable.** ~98–100 % of associated families carry a decoded best EC
(`frac_with_ec` ≈ 0.99 on every axis), so **"EC-less dark" is negligible** — the C·V·T decoder
assigns a plausible EC to almost all of them. The honest limit is decoder *confidence*
(calibrated 0.22–0.43), not missing labels: the recruited families read as plausible
desiccation / salt / phosphorus / β-oxidation enzymes, but each individual EC is a hypothesis.

**Geography vs chemistry.** Conditioning on lat+lon, the master-variable associations are
largely robust — **84 % (Cl), 79 % (S), 72 % (Ca), 71 % (aridity)** and **57 % (salinity_proxy)**
of marginal hits remain significant. Salinity is the most spatially confounded axis (largest
geographic component); chloride and sulfur the least. The gene–environment map is therefore
predominantly real environmental chemistry, with salinity carrying the heaviest spatial load.

## Interpretation in context
This is the genome-level realisation of the classic biogeographic result that **accessory gene
content assorts by niche and that nitrogen and water availability drive core-vs-accessory
functional diversity in soil** (2019_core_accessory_drivers_soil.pdf): the Empty Quarter's
"accessory" stratum — its dark protein families — partitions cleanly along the water/salt
chemistry axis rather than along trace-metal cofactors. The dominance of a salt/aridity filter
with a small osmolyte/phosphate-scavenging specialist set is consistent with the dryland
literature framing **aridity as the master stressor structuring desert microbiomes**
(2024_dryland_microbiomes_desertification_ISME_wrae056.pdf).

## Caveats
- Indicator `*_top.tsv` files are top-40 only; all full counts here are recomputed (and match
  `INDICATORS.md`/`PARTIAL.md`).
- EC labels are low-confidence decoder calls — plausible, not certain.
- Ca/Cl/S/salinity/Na are collinear → overlapping family sets; per-axis counts not independent.
- Relative-abundance correlations are associative, not causal; XRF axes use 251 samples,
  aridity 456; families present in <25 subset samples are untested.

## Files
- `gene_env_summary.tsv` — unified per-variable map (counts, pos/neg, top families, partial).
- `fig_gene_env.pdf` / `.png` — panel a recruitment per axis, panel b strongest family per axis.
- `run.py` — reproducible (`python3 RQ/B3/run.py`); `methods.md` — full method.
