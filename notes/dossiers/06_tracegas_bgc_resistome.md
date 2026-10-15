# Dossier 06 — Trace gases (§2.2) · Biosynthesis, resistome and the biosafety control (§2.4)

Sources read: `paper/main.tex` (§2.2, §2.4, M.4, M.5, M.9, M.11, `fig:lifeonair`/`fig:bgc`/`fig:amr` captions,
Box 1, fenced `tab:atmo`), `paper/supplementary.tex` (Supplementary Note 1 + Tables S1/S2),
`notes/CATALOGUE_NUMBERS.md`, `notes/CORRECTIONS.md`, `notes/PIPELINE_GUIDE.md` §6.3/§§2.2–2.4 (rows carrying
its **V** = verified-against-run-script mark), `notes/DISCUSSION_NOTES.md`, `notes/PAPER_STRUCTURE_PROPOSED.md`,
`resources/literature/{READING_NOTES,citation_index,CITATION_AUDIT}.md`, `RQ/F2/result.md`,
`code/pathogen_census.py`, `reproducibility/claims.tsv`, the workbench scripts themselves
(`M/06_atmospheric/lifeonair_clean/scripts/lifeonair_1..4`, `M/04_functional/repko_clean/scripts/repko_1..3`,
plus that folder's `README.md`, `ENERGETICS.md`, `classify/`, `results/` and job logs), and the repo mirrors
`results/{bgc_census,bgc_product,amr_class,stress_class}_summary_5438.tsv`, `results/vfdb/`,
`results/atmospheric/` (15,469 snapshot). M = `/ibex/scratch/projects/c2014/rob/metagenomics`.

## Questions (both)

- §2.2 (RQ D1–D3): does the annotated fraction encode *bona-fide high-affinity* atmospheric trace-gas machinery rather than its wastebasket homologues; are the carriers abundant; and is atmospheric supply energetically sufficient to sustain them?
- §2.4 (RQ F1 + the A5 paragraph): what chemical-interaction repertoire — biosynthetic and defensive — does the catalogue encode, is the resistome geogenic or anthropogenic, and is the community depleted enough of human-associated organisms to serve as **the contamination control the rest of the paper leans on**?
- Shared payload: both sections are "the known parts" argument — familiar metabolism at unfamiliar scale, with §2.4 doubling as the biosafety control feeding §3.1's specificity claims.

## Basis

| | §2.2 trace gas | §2.4 BGC / resistome |
|---|---|---|
| genome set | 5,438 species reps (dRep 3.6.2, ANI95) | same 5,438 |
| marker input | `rep_ko_v2.tsv` → protein sequences from the Bakta proteomes | antiSMASH: rep genome FASTA · AMRFinderPlus: rep proteomes (Bakta `.faa`) |
| completion | 5 `_sub` reps dropped post-typing (marker on contigs refinement removed) | 5,438/5,438 both arrays, logs clean |
| abundance cohort | 548 profiled soil samples (CoverM v0.7.0) | same 548 |
| taxonomy | — | GTDB-Tk 2.7.2 / **R232** for the pathogen census |
| workbench | `M/06_atmospheric/lifeonair_clean/` | `M/06_metabolism/{bgc,resistome,pathway_check}_clean/` + `M/06_resistome_safety/vfdb/` (PI, reused) |

**`rep_ko_v2` and what it fixed.** §2.2 (and §2.1–§2.5) run on `M/04_functional/repko_clean/results/rep_ko_v2.tsv`,
never the PI's `M/04_functional/refs/rep_ko.tsv`, which is preserved unmodified. Root cause of the defects:
annotation ran on the **raw, pre-refinement** HQ+MQ bins (15,469 rows) while the catalogue is built from
**post-refinement** bins, and refinement both created `_sub` bins and promoted LQ bins. Three fixes:
(i) 44 reps named `48PRr1__SemiBin_N` were unannotated because the resolver splits on a `single_`/`self_`
token — "the annotation existed all along"; (ii) 80 reps whose raw bin was LQ pre-refinement plus 44 `_sub`
reps whose parent was also unannotated (168 unannotated total → **0**), re-annotated with the identical stack
(Bakta 1.12.0, DB 6.0; eggNOG-mapper 2.1.13 executable against 2.1.12 reference data, both quoted, not
harmonised); (iii) **966 `_sub` reps had inherited the PARENT bin's annotation** — a superset; on one measured
rep 422 of 1,521 KOs (28%) were spurious — switched to contig-subset KOs (1,010/1,010 corrected).
Result: 5,438 rows, **0 empty rows** (the PI's table has 168), 11,395 distinct KOs (his 11,357), 4,580 rows
identical, **690 reps shed ≥1 spurious KO (8,134 instances)**, 168 reps gained 195,004. Effect on §2.2:
guild 1,822 → **1,860**; whole-community abundance 18.0% → **18.4%**. M.4 documents fix (iii) obliquely
("Representatives emitted by refinement as strict contig subsets … were annotated by restricting the parent's
gene calls to their own contigs"); the paper never names `rep_ko_v2` or the superseded table.

**The two denominators** (M.5, used identically by both sections): *whole-community* = share of all of a
sample's reads, treating organisms outside the catalogue as non-carriers — a conservative lower bound;
*mapped-community* = share of the reads that map to any representative. Median unmapped ≈ 51%, so the two
differ by roughly a factor of 2.4. Detection is CoverM's never-passed `--min-covered-fraction` default of 10%,
"an inherited default rather than a chosen threshold".

## How it ran

**§2.2 step 1 — crude KO census.** `lifeonair_1_census.sh`, awk over `rep_ko_v2.tsv`: hydB = **K06281**,
coxL = **K03520**. The script labels itself a *"CRUDE KO CENSUS … FIRST PASS only: KO presence over-counts,
because it does not distinguish"* the groups. The reported counts are **not** this census.

**Step 2 — HydDB subgroup assignment.** `lifeonair_2_type_guild.sh`, DIAMOND 2.0.15
(`$CONDA/eggnog/bin/diamond`): `blastp -d $REFS/hyddb --ultra-sensitive --max-target-seqs 1 --outfmt 6
qseqid sseqid pident qcovhsp evalue bitscore`. Reference = the GreeningLab `HydDB_all_hydrogenases.faa`
(`--depth 1` clone of 2026-06-16), whose headers carry the class in field 3 of the pipe split
(`>WP_005035184.1|Holophaga_foetida|[FeFe]_Group_A1`); class = that field of the single best hit.
High-affinity: `if any(t in grp.lower() for t in ("1h","1l")) and "nife" in grp.lower()`.
⚠ **No identity, coverage, bitscore or e-value floor is applied** — only DIAMOND's implicit `--evalue 0.001`
gates whether a hit exists. Measured from `classify/hydB_diamond.tsv`: group 1h spans 53.8–98.8% identity
(min query coverage 42.2%), **group 1l spans 37.8–92.7% identity with minimum query coverage 20.5%**.
Current class table: 1l 1,124 · 1h 445 · 2a 140 · 1f 11 · 1d 11 · 1c 7 · rest ≤3 each.

**Step 2b — form-I coxL: the rule is OR, and M.9 is right.** DIAMOND against **three** Cordero 2019 RefSeq
form-I references only (WP_013913730.1 *Afipia carboxidovorans*, WP_003892166.1 *Mycolicibacterium smegmatis*,
WP_006067999.1 *Natronorubrum bangense*), best hit = first line for the query. Decision code, verbatim:
`call = "form_I" if (thr or mot == "AYRGAGR(formI)") else ("form_II/other" if mot == "AYRxSFR(formII)" else
"ambiguous")` with `thr = (pid > 55.0 and qc >= 90.0)`. So **form I ⇔ (>55% identity AND ≥90% query coverage)
OR the AYRGAGR motif** — either alone suffices. The older `results/atmospheric/` method note saying "called
form I on concordance of the two" is **wrong for the current run**; M.9's "or" is correct.
Current distribution of 4,727 coxL proteins: **form_I 1,901 · ambiguous 2,815 · form_II/other 11**
(motif: AYRGAGR 1,627 · AYRxSFR 147 · none 2,953). Note the ordering artefact: 147 proteins carry the
**form-II** motif but only 11 are called form-II, because the form-I clause is evaluated first.

**Step 2c — trees.** MMseqs2 18.8cc5c `easy-cluster --min-seq-id 0.7 -c 0.8` → MAFFT 7.526 `--auto` →
FastTree 2.2.0. **No bootstrap, no model specification, no placement software (no pplacer/EPA-ng).**
Display only. M.9's "typed phylogenetically" therefore overstates: the calls are 100% DIAMOND best-hit
(+ motif for coxL). ⚠ The v2 typing job never wrote `CLASSIFY_DONE`; `classify/hydB.aln` is **0 bytes** and
`classify/hydB.nwk` is the **Aug-10 v1 tree** — the hydB tree figure input is stale. Numbers are unaffected
(typing, guild and abundance all completed before the tree step). The coxL tree is current.

**Step 3 — guild and abundance.** Guild = ≥1 high-affinity hydB and/or ≥1 form-I coxL, built with
`sort -u` / `comm -12` over the classification table. ⚠ **The 5-rep post-typing filter is in no script.**
`.prefilter` files hold 970 hydB / 1,173 coxL / 1,864 any / 279 both; the shipped 968/1,170/1,860/278 differ
by exactly the five named `_sub` reps, dropped by hand. This is the only un-scripted step in §2.2 and it
would not survive a reproducibility rerun — it also closes the README's open "small reconciliation owed".
`lifeonair_3_abundance.py` then reads the CoverM `abund.tsv` profiles, computing
`whole = Σ guild abundance` and `mapped = whole/(100 − unmapped)×100` over canonical ∩ profiled = **548**
at run time. ⚠ **That cohort is now 551**: the three missing profiles (51Sr1, F13PRr1, F5Dr1) landed
6–8 Sep 2026, so a re-verification would shift the medians slightly.

**Step 4 — energetics** (`lifeonair_4_energetics.py`; narrative = `ENERGETICS.md`; paper derivation =
Supplementary Note 1). Every constant in the script carries an inline `[lit]` / `[measured]` / `[assumed]` tag.
*Nernst:* `dG = DG0 + RT*log(...)` at T = 303.15 K, p(H₂) 0.53 µatm, p(CO) 0.10 µatm, p(CO₂) 420 µatm,
p(O₂) 0.2095 atm; ΔG°′ −237.2 / −257.2 kJ mol⁻¹ → **−198.8 / −234.2 kJ mol⁻¹**.
*Fick:* J = D₀·f_eff·C/Δx with **Millington–Quirk** `f_eff = POROSITY**(4/3)`, ε = 0.38 → D_e/D₀ = 0.275;
lower boundary **C = 0 at Δx = 0.10 m — an explicit maximum-potential supply**; D₀(H₂) 6.1×10⁻⁵,
D₀(CO) 2.0×10⁻⁵ m² s⁻¹ → J = **3.58** and **0.22** nmol m⁻² s⁻¹. ⚠ This is the *air-filled-porosity* form
with **no water-content or saturation term**: the soil is treated as fully dry.
*Power:* P = Σ J|ΔG| = **7.63×10⁻⁴ W m⁻²** over a 1 m² × 0.10 m column at dry bulk density 1.6 g cm⁻³ with
`GUILD_FR = 0.184` (the script uses 0.184; Supp Note 1 and `ENERGETICS.md` say 0.18).
*Cell density* `CELLS_G = (1e6, 1e7, 1e8)` → 2.94×10¹⁰–2.94×10¹² guild cells m⁻² → **2.59×10⁻¹⁴ to
2.59×10⁻¹⁶ W cell⁻¹**. The survival threshold ~10⁻¹⁸–10⁻¹⁹ and the active-maintenance band 10⁻¹³–10⁻¹⁵
are **hard-coded print strings, never variables** — maintenance energy is not a fitted parameter.
*Sink:* `H2_RATE = (5,25)` mg m⁻² yr⁻¹ × `AREA = 6.5e11` m² → **0.0032–0.0163 Tg yr⁻¹ = 0.004–0.027%** of
`SINK_H2 = (60,80)` Tg yr⁻¹; GWP₁₀₀ = 11 → **36–179 kt CO₂-eq yr⁻¹**; CO 0.0130–0.0650 Tg yr⁻¹ =
0.005–0.026% of `SINK_CO = 250`.
⚠ `ENERGETICS.md` records its own correction, verbatim: *"**NB Gemini attributed this to Buzzard/Conrad but
Buzzard 2022's actual semiarid flux is ~45-251 mg/m2/yr (10x higher) — 5-25 is a hyperarid downscale, state
as such.**"* and, on diffusivity, *"Gemini used D0(H2)=2.4e-5 m2/s — WRONG … H2 is light"*.
⚠ Three numbers disagree between `ENERGETICS.md` and the script (7.7e-4 vs 7.63e-4 W m⁻²; 2.7e-14 vs
2.59e-14 W cell⁻¹; H₂ sink upper bound 0.02% vs 0.027%). **The script is authoritative.**

**§2.4 antiSMASH.** v8.0.4 (`module load antismash/8.0.4`):
`antismash --minimal --genefinding-tool prodigal --cpus 8 --output-dir <mag> --output-basename <mag> <fna>`
(array 1-272%60, BATCH=20). **`--minimal` = rule-based core cluster detection only** — no ClusterBlast,
KnownClusterBlast or MIBiG comparison, which is structurally why novelty cannot be read off this run.
Contig-edge uses antiSMASH's own per-region `contig_edge` flag: **37,529 / 44,002 = 85.29%**, with
`areas_vs_gbk_mismatches = 0` (the areas-JSON and GBK region parses agree). Counting rule = per-species
presence of a product type; a multi-product region counts once per type (the product rows sum to 48,205
region-instances against 44,002 regions); terpene-precursor kept distinct from terpene.

**§2.4 AMRFinderPlus.** v4.2.7, database **2026-05-15.1**: `amrfinder -p <faa> --plus --name <mag> --threads 4`
(array 1-37, BATCH=150). `--plus` is what yields the metal/biocide **STRESS** layer; **`--organism` is NOT
passed**, so there is no organism-specific point-mutation detection and no organism blacklisting — resistance
is acquired-gene hits only. AMR and STRESS are separated on the element-type column; antibiotic hits total
**4,182 = 4,040 core + 142 plus** (verified by summing `amr_class_summary_5438.tsv`).

**§2.4 pathogen census + VFDB.** The census matches the reps' GTDB **R232** genus and species strings against
two fixed panels (see *Known defects*). Virulence factors are MMseqs2 searches of the rep proteomes against
**VFDB setB** (29,804 proteins) at **≥40% identity, ≥70% coverage, e ≤ 10⁻¹⁰**, read along identity tiers —
PI compute, reused by reference, verified by MAG naming.

## Controls

- **Clade placement, not a KO call.** The wastebasket objection is answered by typing every marker protein
  against group-labelled HydDB references. Historical-basis breakdown: 1l+1h = 6,278 of 6,602 hydB proteins,
  the only other classes being 2a (246) and 1a–1f/2b–2c (<70 combined). Current basis gives 1,124 (1l) vs
  445 (1h); the full current class table sits in `lifeonair_clean/classify/` and is not in the registry.
- **coxL form-I is a floor, not a census — confirmed on the current run.** Of 4,727 coxL proteins,
  **1,901 are form_I and 2,815 (60%) are "ambiguous"**, with 11 form-II/other. The discarded class is larger
  than the called one, on three reference proteins. The 22% is a lower bound and the paper never says so.
- **Dark-supply negative** (Box 1, "trace-gas enzymes *not* dark-supplied"): targeted HydDB DIAMOND +
  NiFe/FeFe HMMs + Cordero form-I references + AYRGAGR motif over 468,211 "reallydark" and 563,619 "isdark"
  reps returned **0 high-affinity 1h/1l and 0 form-I** — "≈ 0% of the functional dark matter is recoverable as
  a trace-gas enzyme". ⚠ computed on the June dark cohorts, **not** the current 4.29 M / 106,042 basis.
- **Aridity negatives** (registered `a1.aridity.cross` = NEGATIVE): within-biome 1l share vs precipitation
  ρ = **−0.1716, p = 0.196** (58 sites, permutation); guild size vs precipitation ρ = 0.0249, p = 0.849;
  cross-biome the 1l:1h ratio does not order by aridity (atacama ∞, EQ 5.37, bay2021 3.51, antarctic 1.22,
  mojave 0.46, arid 0.31, gtdb_soil 0.17). ⚠ **all on the 15,469 basis; never recomputed on 5,438**, yet §2.2
  cites them in the present tense. (Note: `PAPER_STRUCTURE_PROPOSED` R8 borrows ρ=−0.172/p=0.196 as an
  aridity null — it is *this* hydrogenase test, not a dark-family result.)
- **Ectoine triangulation.** antiSMASH 899 species vs KO panel 895, intersection 860 →
  Jaccard 860/(895+899−860) = 860/934 = **0.921** (arithmetic verified here). Two *independent* pipelines,
  nucleotide rule-based vs protein orthologue. The advertised **third** method — MetaCyc pathway completeness
  in pre-gap-fill draft reconstructions — is **insensitive, not concordant**: 15 ectoine calls in 5,438. It was
  dropped from M.11 along with the trehalose KO×gapsmith Jaccard 0.795 and is recoverable only from the
  workbench. So §2.3/§2.4's "two independent pipelines" is accurate; "three-method triangulation" (M.11) is not.
- **Biosafety pathogen control**: 2 named pathogen species (0.04%) + 0.9% pathogen-genus species + the VFDB
  identity-tier collapse 5,435 → 3,769 → 202 species at ≥40 → ≥70 → ≥90.
- **What has no control.** The BGC census (no negative or matched genome set, no size normalisation — 97.5% is
  a saturated statistic); the AMR census likewise; the metal gradient tests have geography partialled out but
  no null-axis control of the §2.5 vanadium kind; the energetics model has no sensitivity analysis beyond the
  three cell densities and no moisture gating; guild abundance has no biomass or sequencing-depth covariate.

## Numbers

REG = registered in `CATALOGUE_NUMBERS.md`. Sources relative to the basis folder named in *Basis*.

| § | number | value | source | mark |
|---|---|---|---|---|
| 2.2 | high-affinity hydB (HydDB 1h/1l) | **968 species (18%)** | `lifeonair_clean/README.md` | REG |
| 2.2 | form-I coxL | **1,170 (22%)** | same | REG |
| 2.2 | guild (either) / core (both) | **1,860 (34%) / 278 (5%)** | same | REG |
| 2.2 | proteins, group 1l vs 1h | **1,124 vs 445** (2.53:1) | `classify/` | REG |
| 2.2 | guild abundance | median **18.361%** whole (IQR 12.829–25.614) / **43.919%** mapped (IQR 30.764–55.813); coxL 11.186 vs hydB 9.020 whole, 25.756 vs 21.943 mapped, both 2.828/6.611 | `results/guild_abundance_summary.tsv` | REG |
| 2.2 | coxL form calls (not reported) | form_I **1,901** proteins · **ambiguous 2,815** · form_II/other 11, of 4,727 | `classify/coxL_class.tsv` | unreported |
| 2.2 | ΔG H₂ / CO at ambient | **−198.8 / −234.2 kJ mol⁻¹** | Supp Note 1 S1.1 | REG (as "≈ −200 / −230") |
| 2.2 | diffusive supply | J(H₂) 3.6, J(CO) 0.22 nmol m⁻² s⁻¹; P = 7.6×10⁻⁴ W m⁻² | Supp Note 1 S1.2–S1.3 | supp only |
| 2.2 | per-cell power | 2.6×10⁻¹⁴ … 2.6×10⁻¹⁶ W vs ~10⁻¹⁸–10⁻¹⁹ survival | Supp Table S1 | supp only |
| 2.2 | biome H₂ sink | **0.0032–0.0163 Tg yr⁻¹ = 0.004–0.027%** of the global soil sink; CO 0.0130–0.0650 = 0.005–0.026% | `lifeonair_4_energetics.py` (script is authoritative) | REG as "~0.01 Tg, <1%" |
| 2.4 | BGC census | **97.48%** of species ≥1 cluster · **44,002** regions · contig-edge **85.29%** (37,529) | `bgc_census_summary_5438.tsv` | REG |
| 2.4 | BGC products (% of species) | terpene-precursor **83.78** · terpene **62.14** · T3PKS **41.49** · RiPP-like **36.52** · NRPS-like **31.28** · NRPS **25.06** · RRE-containing **23.41** · T1PKS **22.36** · ectoine **16.53** (899) | `bgc_product_summary_5438.tsv` (82 product types) | REG (NRPS-like, RRE not in text) |
| 2.4 | BGC carriage | any-BGC median **48%** whole / **99%** mapped | `bgc_clean/results/community_summary.tsv` | REG |
| 2.4 | AMR census | any **37.4** · β-lactam **24.55** (1,335) · aminoglycoside **11.93** (649) · rifamycin **7.25** (394) · phenicol 5.02 · macrolide 4.63 · glycopeptide 4.38 | `amr_class_summary_5438.tsv` | REG |
| 2.4 | AMR carriage | median **15%** whole / **33.5%** mapped | `resistome_clean/results/community_summary.tsv` | REG |
| 2.4 | metal layer (STRESS) | arsenic **17.34** (943) · mercury **7.52** (409) · all others ≤0.2% | `stress_class_summary_5438.tsv` | REG |
| 2.4 | arsenic gradients | Ca ρ **+0.224** (q 0.011, n=248) · precipitation **−0.243** (q 2.4e-6, n=453); salinity proxy attenuates 0.24→0.12 (q 0.11) | `resistome_clean/results/gradient_tests.tsv` | REG |
| 2.4 | pathogen census | **2** named pathogen species of 5,438 (**0.04%**) · pathogen-containing genera **0.9%** | `resistome_clean/results/pathogen_census_summary.tsv` | REG — see defects |
| 2.4 | VFDB tiers | **5,435 → 3,769 → 202** species at id ≥40 → ≥70 → ≥90; 541,508 proteins with a hit | `results/vfdb/VFDB.md` | REG |
| 2.4 | ectoine triangulation | Jaccard **0.921** (860 / 895 KO / 899 antiSMASH); trehalose KO×gapsmith 0.795 | `pathway_check_clean/results/triangulation_summary.tsv` | REG (trehalose dropped from M.11) |

## Comparative context — Discussion only, never Results

| claim | number | key |
|---|---|---|
| Antarctic soil MAGs, group 1h/1l | **322/554 = 58%**; coxL **60/554 = 11%**, all 60 also carry 1h/1l; 339 (61%) any aerotrophy gene | Hutchinson 2025 *ISME J* (bioRxiv 2025.09.23.677255) — **NO BIBKEY**, PDF `trace-gas-deserts/` |
| Antarctic endolithic MAGs | **41.1%** high-affinity [NiFe]-hydrogenase | `coleine2024endolithic` ⚠ PDF on disk is the 2023 preprint, different title — re-check before citing |
| cold-desert MAGs, 1h/1l | about a third | `ray2022atmospheric` |
| 1l:1h ratio elsewhere | Mackay Glacier **2.3×** (close to our 2.5:1); Gobi likewise; Hutchinson 2025 **inverts** it | `ortiz2021antarctic`, `jordaan2020gobi` |
| our position | **18% hydrogenase / 22% CODH → CO-weighted where cold deserts are H₂-weighted** | — |
| soil H₂ sink is under-constrained, not absent | "the largest and most uncertain term in the global H₂ budget"; "paucity of observational constraints … particularly in moisture-limited, arid, and semi-arid soils"; a drier threshold shifts **4–7 pp** into arid/semi-arid | Reji et al. *Nat Commun* 10.1038/s41467-025-67208-3 — **NO BIBKEY** |
| observations "restricted to cold deserts", calls for work "extending to hot desert environments" | — | Ray et al. 2023 *MMBR* 10.1128/mmbr.00048-23 — **NO BIBKEY** |
| global soil H₂ sink denominator | **60 +30/−20 Tg yr⁻¹** soil uptake; 79 +30/−20 total sinks | `ehhalt2009h2cycle` ⚠ Supp Note 1 prints "60–80", which conflates the soil term with the total |
| pre-anthropogenic resistome | 30,000-yr permafrost resistome; pristine Antarctic soils resistant without anthropogenic selection | `dcosta2011ancient`, `vangoethem2018antarctic` |
| bioprospecting, one sentence | 97.5% carriage but 85.3% contig-edge → a **target list** for long-read/enrichment recovery, not a mineable cluster set | — |

Standing instructions from `DISCUSSION_NOTES.md`: **never write "our community is hydrogenase-rich" —
comparatively it is not**; a cross-pipeline caveat is mandatory (their MAGs vs our ANI95 reps, different
binning/annotation, 554 vs 5,438, dereplication status unverified); report 1l dominance as a description of
the pool, not an aridity signal; `bay2021trace`'s ">70% of soil taxa" is inorganic energy sources broadly,
**not** trace-gas enzymes.

## Caveats (quoted)

- §2.2: "These are genomic capacities and DNA-based abundances, evidence of potential and prevalence rather than measured flux; and the pathway is the known trace-gas machinery, not new chemistry, with the dark proteome of Section~2.5 supplying none of these enzymes."
- M.9: "the calculation is a literature-parameter model whose principal uncertainty is the assumed cell density."
- Supp S1.5: "The balance addresses energetic \emph{feasibility}, not realised flux, which is unmeasured here and is additionally moisture-gated in hyperarid soils." · "The absolute cell density is the dominant uncertainty and should be replaced by a measured Empty-Quarter value (e.g.\ 16S rRNA qPCR copies~g$^{-1}$) when available."
- Supp Note 1 opening: "the calculation is a literature-parameter model that the catalogue data scale rather than drive".
- §2.4: "85\% of the regions run off a contig edge, so cluster counts, completeness and any comparison against reference cluster families would measure assembly fragmentation more than biology, and we do not attempt them."
- M.11: "because one representative stands for each species, accessory clusters or determinants confined to non-representative conspecific genomes are not counted, and per-species prevalence is conservative with respect to pangenome content."
- Registry footnote §2.4: "BGC novelty never computed (deferred list)" · "salinity-proxy association attenuates geography-partialled (0.24→0.12) — cite Ca/Na + aridity instead."
- Box 1: "trace-gas enzymes \emph{not} dark-supplied; group-1l hydrogenases show no within-biome or cross-biome aridity gradient; VFDB ``virulence'' dual-use, collapsing at high identity".
- **Correction owed (blocker 9 / citation decision D-3):** §2.2's "effectively absent from those budgets" must become **"poorly constrained"** — `greening2022review` does not support the negative.

## Known defects to carry

1. **The pathogen census matches names without GTDB suffixes.** `code/pathogen_census.py` tests
   `if s in PATHOGEN_SPECIES` — exact set membership, although its own comment claims "substring match on
   `s__`". The genus panel *does* carry suffixed variants (`Klebsiella_A`, `Enterococcus_B`, `Bacillus_A`,
   `Burkholderia_A`); the species panel carries **none**, so `s__Escherichia_A coli`, `s__Bacillus_A
   anthracis` and `s__Enterococcus_B lactis` are silently missed. Suffix-blind the count is **5/5,438**, not 2
   (eqmags_poc/07, 2026-09-10, UNREGISTERED). **Touches**: §2.4's "two are named pathogen species … 0.04%" and
   the 0.9% genus figure. The repo copy read here is the 15,469-basis version; the PoC audited the current one
   and found the same rule.
2. **The E. coli representative is probably a reagent contaminant.** `V42PRr2__single_SemiBin_33` is detected
   in **316/550** samples; eqmags_poc/07 calls it CONTAMINATION LIKELY — 100% of trip-5 samples vs 24% of
   trip-3, negative against two independent biomass proxies in all five trips, no replicate or site structure,
   313/316 detections below 0.5% abundance. *A. baumannii*, `Enterococcus_B lactis` and the `Bacillus_A`
   cereus-group bin are UNDECIDABLE (9–28 detections, no contamination signature); the "anthracis" call is
   borderline at 95.2% ANI with no toxin genes — **never write anthrax**. **Touches**: the biosafety
   paragraph's own logic, since one of its two "named pathogen species" is likely lab contamination being read
   as evidence of pristineness. It moves **no** BGC, AMR, metal or trace-gas percentage (one rep = 0.02%).
3. **The contamination control has never been run against its own blanks.** Shotgun negatives SNTC (335,190
   reads) and SPositive exist under `trip4/controls` in `EmptyQuarter_Data` and were never profiled against the
   5,438 (blocker 3, hours of compute).
4. **`reproducibility/claims.tsv` still asserts dead 15,469-basis values** for both sections — coxL 5,645
   proteins / 3,828 MAGs, guild core 1,204, terpene 64.5, NRPS 25.3, ectoine 17.0, As 2,423, Hg 1,077,
   β-lactam 24.7, pathogens 0.05%. PI rebinds at end-game (blocker 8).
5. **§2.2 is otherwise clean.** `rep_ko_v2` was the one broken join and it was fixed, re-registered and
   re-shipped; no defect above moves a §2.2 number.

## Status

- **§2.2 DONE + SHIPPED 2026-08-12** on `rep_ko_v2` (recomputed from the 2026-08-11 5,438 run);
  **§2.4 DONE + SHIPPED 2026-08-15** on 5,438. Both registered in `CATALOGUE_NUMBERS.md` and `CORRECTIONS.md`.
- **Repo mirrors exist for only four §2.4 tables**. Community carriage, gradient tests, the pathogen census,
  the triangulation and **every current §2.2 table** live only in the workbench; `results/atmospheric/` is the
  15,469 snapshot and its README points at a `results/atmospheric_5438/` that does not exist.
- One open `\todo` in §2.2 (promote to a 5,438 guild table + a HydDB/CODH placement figure). §2.4 has none.
- Figures `fig_lifeonair_guild.pdf`, `fig_bgc_product.pdf`, `fig_amr_class.pdf` all exist, all
  `\status{candidate}`; `PAPER_STRUCTURE_PROPOSED` demotes all three to Extended Data (ED 1, ED 2).
- **Destination in the rewrite:** both sections collapse into **R7**, "The annotated fraction encodes the
  expected desert toolkit", ~200–250 words shared with §2.3's floors; biosafety is gated on blockers 1–3.

## What can and cannot be claimed

**Can.** A third of the catalogue (1,860 species) encodes high-affinity trace-gas machinery typed against
group-labelled references rather than KO calls; the carriers are a pervasive baseline, not a niche
(18.4% whole-community / 43.9% mapped, 548 samples); the pool is CO-weighted and group-1l-dominated
(1,124 vs 445 proteins); oxidation stays strongly exergonic at ambient hyperarid mixing ratios and the
diffusive supply exceeds basal maintenance by 2–4 orders of magnitude across 10⁶–10⁸ cells g⁻¹ — an
*energetic-feasibility* statement for maintenance and persistence, not growth; the biome's share of the
global soil H₂ sink is small (<1%, actually 0.004–0.03%) and the point is that hot hyperarid drylands are
poorly constrained in those budgets. Biosynthetic capacity is near-universal (97.48% of species, 44,002
regions) and broadly carried (48% whole / 99% mapped); the resistome is a host-genome trait of a biome with
minimal anthropogenic input (37.4% of species, β-lactam-led); the metal layer exceeds every antibiotic class
but β-lactam and reads geogenic (As 17.34%, Ca ρ +0.224, precipitation ρ −0.243, geography partialled); two
independent detectors agree on ectoine carriers at Jaccard 0.921.

**Cannot.** Any flux, rate or activity statement in either section. Any BGC novelty, cluster completeness or
product claim — `--minimal` mode ran no comparative module and 85.29% of regions abut a contig edge. Any
count of complete pathways. "Hydrogenase-rich" in any comparative sense. Any aridity signal from the 1l:1h
ratio (and the negatives supporting that are on the dead 15,469 basis). "Three-method triangulation" —
the third method is insensitive. A pangenome-level prevalence (one rep per species). Point-mutation
resistance of any kind (`--organism` not passed). Rifampicin resistance prevalence as a phenotype. "Effectively
absent from global budgets". And, until blockers 1–3 clear, the biosafety paragraph cannot be written as a
clean contamination control.

## Reviewer attack surface

1. **The high-affinity hydrogenase call has no identity or coverage floor.** Best hit against HydDB with
   nothing but DIAMOND's default e-value: the weakest group-1l assignments sit at **37.8% identity over 20.5%
   query coverage**, and 1l is 72% of the high-affinity pool and the whole basis of the "arid-soil 1l
   predominates" sentence. **The single sharpest attack on §2.2.** A floor (say ≥50% id / ≥70% cov) or the
   published HydDB-classifier SVM would answer it; either way the sensitivity must be reported.
2. **60% of coxL proteins are "ambiguous" and the section does not say so** (2,815 of 4,727), against three
   reference proteins. The form-I rule is also an OR, so a protein failing 55%/90% enters on motif alone —
   and 147 proteins carrying the *form-II* motif are adjudicated by a clause that never runs for them.
3. **"Typed phylogenetically" is DIAMOND best-hit.** No bootstrap, no model, no placement software; the trees
   are display-only and the hydB tree on disk is the stale v1 (0-byte alignment, v2 job never finished).
   Reword M.9 to "classified by best-hit against group-labelled references, with trees shown for display".
4. **One step of §2.2 exists in no script** — the five hand-dropped `_sub` reps between `.prefilter`
   (970/1,173) and the published 968/1,170. A data-availability rerun reproduces the wrong numbers.
5. **The energetics model is 9 literature/assumed parameters against 3 measured ones**, with a
   maximum-potential diffusive supply (C = 0 at the lower boundary) and no moisture gating — in a soil where
   Buzzard 2022, the paper's own flux anchor, reports uptake suppressed when dry. The 2–4 order margin means
   the conclusion survives, but the supply is an upper bound and should be labelled one.
6. **The Ehhalt denominator is misquoted** — `SINK_H2 = (60, 80)` is hard-coded in the script and printed in
   Table S2 as "Ehhalt & Rohrer 2009", but that paper gives soil uptake 60 +30/−20 and *total* sinks 79.
7. **No size normalisation on any census.** eqmags_poc/03: genome size alone explains R² 0.52 of KO richness.
   BGC, AMR and marker counts are all gene counts on genomes of unequal size and completeness.
8. **97.5% BGC carriage is saturated and uncontrolled** — it cannot distinguish "this soil is biosynthetically
   rich" from "antiSMASH finds a terpene cluster in almost any actinobacterial genome". There is no matched
   non-desert genome set anywhere in the section.
9. **The product ranking in the text skips two classes it outranks**: NRPS-like (31.28%, 5th) and
   RRE-containing (23.41%) are omitted while NRPS (25.06%) is reported. `fig:bgc` also says "81 product types"
   against 82 rows in the table.
10. **AMRFinderPlus without `--organism`** detects no point mutations, so rifamycin 7.25% (typically *rpoB* in
   Actinomycetota) is an acquired-gene floor, and β-lactam 24.55% is a hit class known to over-call in soil
   actinobacteria. The "intrinsic resistome" reading is right; the numbers are not phenotypes.
11. **The VFDB "secretion, nutrient uptake and adhesion" reading is selective.** At id ≥90 the top categories
    by protein count are Nutritional/Metabolic 257, Adherence 256, **Regulation 256, Motility 251**, Immune
    modulation 205, Effector delivery 202 — two categories as high as the three named are not mentioned.
12. **VFDB ran on "20.5M proteins"**, the CheckM2 gene-call denominator, where every other §2.x number uses the
    Bakta 19,186,853 — unexplained, and it is PI compute reused rather than rerun.
13. **The contamination control is argued from a census that miscounts pathogens and includes a likely
    contaminant, and has never been run against the existing blanks** (defects 1–3).
14. Minor: §2.2 cites aridity negatives computed on a superseded catalogue in the present tense; the abundance
    cohort is quoted as 548 when 551 profiles now exist; §2.3 and §2.4 both state the ectoine Jaccard.

## Open questions for Alejandra

1. **Two §2.2 methods decisions, both cheap.** (a) Put an identity/coverage floor on the HydDB best-hit, or
   report the sensitivity of 968 to one — currently the weakest 1l call is 37.8% id / 20.5% cov. (b) Say in
   §2.2 or M.9 that 2,815 of 4,727 coxL proteins are unresolved, so 1,170 is a floor. Recommendation: do both.
2. **E. coli representative** — annotate or remove from the 5,438 (blocker 1, PI call)? Removal costs nothing
   numerically (0.02%), but the biosafety paragraph needs rewriting either way.
3. **Fix `pathogen_census.py` at the root** to strip GTDB alphabetic suffixes and re-register 5/5,438, or keep
   2 and state the matching rule in M.11? Recommendation: fix at the root and report both counts.
4. **Profile SNTC/SPositive against the 5,438** before the contamination control is cited again (hours)?
5. **Supplementary Note 1** — replace "60–80 Tg yr⁻¹" with Ehhalt's 60 +30/−20 soil term, and key Reji et al.
   and Ray et al. 2023 so D-3's "poorly constrained" wording has a source. Hutchinson 2025 also needs a key.
6. **Measured cell density** — does a 16S qPCR copies g⁻¹ number exist for these samples anywhere? It is the
   model's stated dominant uncertainty and would turn Supp Note 1 from feasibility into an estimate.
7. **Size-normalise the BGC and AMR censuses** per eqmags_poc/03 before R7 is drafted?
8. **Does BGC novelty stay deferred?** The historical `RQ/F2` object exists (BiG-SLiCE v2.0.2 × MIBiG 3.1,
   119,090 regions, 42,064/42,365 = 99.29% of GCFs with no MIBiG member, 90.1% contig-edge) and is explicitly
   "not promoted"; nothing has been computed on the 44,002-region basis. Keep deferred, or state the negative?
