# Discussion material bank

Deposit box for interpretive / comparative material that must NOT go in Results. Each section
session appends its Discussion-worthy findings here as it closes; session 10 ("discussion-abstract")
generates the Discussion from this file. Keep entries short: claim → numbers → caveats → source.

## From §1.4 (deposited 2026-08-14)

- **FESNov as frame, not target.** Our purifying-selection result (median dN/dS 0.425, 96.7%<1,
  vs matched control 0.388) mirrors the global family-level surveys of unknown genes; what EQ adds
  is what a pooled catalogue cannot measure — within-biome recurrence (12 samples/9 sites),
  persistence (79.4% all five expeditions) and member-aware abundance (98/550) for the same
  families. Key: `rodriguezdelrio2024fesnov`.
- **Relic-DNA counter works only as a combination.** Persistence + selection + read recruitment
  jointly exclude relic/contaminant explanations; no single leg does (relic DNA can be abundant,
  contaminants can recur within a batch, selection alone doesn't date the DNA). Keys:
  `carini2016relicdna`, `salter2014contamination`; Lennon 2018 bib pending (citations session).
- **Matched-control lesson (methods honesty + armor).** June's larger dark-vs-control gap
  (0.391 vs 0.286) was partly a length confound; protein-level clustering preselects conserved
  sequence, so only an identically ascertained control is defensible. Designated escalation if a
  reviewer raises residual circularity: within-population pN/pS (inStrain; per-sample sorted BAMs
  exist, bounded cost, deferred to revision).
- **Unexplored biology of dominant lineages = the thesis spine.** 97.9% of attributable families
  single-phylum, led by Actinomycetota (32%) / Thermoplasmatota (19%): dark matter is the normal
  accessory biology of the community's dominant, named phyla — the Discussion's cleanest statement
  of "known phyla, unknown biology". Caveat: 74% attributable only.
- **The young class as future work.** 782,616 recurrent 1–2-variant dark clusters sit below the
  admission rule: recently emerged or rare accessory genes; deliberately unclaimed in Results,
  worth one forward-looking sentence.

## From §2.3 (deposited 2026-08-13)

- **NHEJ excess over bacteria generally.** Definition-matched: our complete Ku+LigD = **55.6%**
  (3,021/5,438 species) vs Sharda 2020's **22%** (~6,000 mostly cultured complete genomes) — 2.5×.
  Do NOT use our either-gene 74.5% next to their number. Caveats: their cultured denominator biases
  them low, our MAG incompleteness biases us low (roughly opposing); NHEJ tracks large slow-growing
  genomes (their own finding), so part of our excess is lineage composition (Actinomycetota), not
  desert selection. Full breakdown in workbench `uv_survival_clean/README.md`. Key: `sharda2020nhej`.
- **Three stressors, three response modes** as a Discussion thread: light partitions communities
  vertically (10 cm of shielding), desiccation is met by a uniform baseline, salt replaces the
  community (Halobacteriota turnover). Possible generalisation: the response mode matches the
  spatial structure of the stressor. Keep as interpretation, not a law.
- **Capacity ≠ activity** is the standing limit for all of Part 2; metatranscriptomics on the same
  samples = the natural next study. Assembly-level HMM profiling (existing 753M proteins) = the
  in-reach validation for the ~40% assembled-but-unbinned community.
- **Atacama glycogen→solutes** (Moore 2023 citing Baran 2017; secondhand, model-based, cyanobacterium):
  half-sentence hedge that synthesis capacity is not the only route to osmotic protection.

## From §2.2 comparators (deposited 2026-08-13, lit-agent ingest)

- **Hot desert is CO-weighted where cold deserts are H₂-weighted.** Genome fractions (species/MAG
  level, NOT abundance): our hydrogenase 17.8% (968/5,438) vs Hutchinson 2025 Antarctic 58%
  (322/554 MAGs; dereplication status unverified) and Coleine endolithic 41.1% (dereplicated MAGs);
  our CO 22% and CO > H₂ by abundance, consistent with 1l-over-1h. Cross-pipeline caveat mandatory;
  a same-pipeline recount would be needed for a hard claim. Never write "our community is
  hydrogenase-rich" — comparatively it is not.
- **Bay 2021 scope**: the ">70% of soil taxa" is inorganic energy sources broadly, NOT trace-gas
  enzymes specifically — phrase accordingly wherever benchmarked. Key: `bay2021trace` (consolidated
  2026-08-13; `bay2021desert` deleted).

## From §2.4 (deposited 2026-08-15, session 03-bgc-resistome)

- **Pre-anthropogenic resistome as the Discussion angle.** The intrinsic-resistome framing (β-lactam
  24.6% of species in a biome with ~zero human input) is the useful comparison seed: D'Costa 2011
  (30,000-year permafrost resistome) and Van Goethem 2018 (pristine Antarctic soils, resistance
  without anthropogenic selection) are the load-bearing prior art — bibkeys pending, \todo in §2.4.
  Claim shape: resistance genes predate and do not require the antibiotic era; the EQ adds a
  hot-hyperarid data point where the metal layer (As 17.3%) exceeds every antibiotic class except
  β-lactam and tracks mineral chemistry (Ca/Na + aridity), not human activity.
- **Biosafety/contamination control feeds C3/§3.1**: 2 named pathogen species of 5,438 (0.04%);
  VFDB collapse 5,435→202 (id≥40→90) resolving to secretion/uptake/adhesion. This licenses "the
  dark proteome is not contamination riding along with human material" — cite from §3.1/Discussion
  when the specificity claims need it, never re-argued in Results.
- **Bioprospecting line (Discussion only, one sentence):** 97.5% of species carry BGCs and the
  carriers are ~half the community by abundance, but 85.3% contig-edge fragmentation means the
  catalogue is a *target list* for recovery (long-read / enrichment), not a mineable cluster set.
  Pairs with the §1.3 fold-novelty targets in the existing "targets for discovery" Discussion
  bullet.
- **Deferred negatives worth one Discussion clause if the PI agrees** (PI list item 5): essentially
  no N-fixation capacity (nifH/D/K ~20 of 15,469 old-basis bins, vnfD = 0 — the planned V-nitrogenase
  vs XRF-vanadium survey answers "absent"); rhodopsins present but rare (341 old-basis MAGs). Both
  recompute on rep_ko_v2 in minutes if wanted; natural home = §2.2's energy story or Discussion.
- **Ectoine triangulation as a methods-confidence point**: two independent pipelines (nucleotide
  rule-based antiSMASH vs protein KO census) agree at Jaccard 0.92 on ectoine carriers; usable in
  Discussion/response-to-reviewers as evidence the marker censuses are not annotation artifacts.

## §2.5 deposits (2026-08-16, session 04-dark-function)

- **One master gradient organises known and dark biology** (the section's interpretive payoff,
  Discussion-grade): salt/water chemistry drives species turnover in the annotated community
  (§2.3, Halobacteriota) AND filters the dark-family repertoire (§2.5, ~90% negative associations
  on S/Ca/Cl/salinity). Discussion can argue the dark fraction behaves like functional biology
  under the same environmental control, completing §1.4's "real, not spurious" with "responsive,
  not inert". Do not over-read: associations, not selection or mechanism.
- **Filter vs recruitment asymmetry**: the gradient mostly EXCLUDES dark families toward the
  salty/dry end; the recruited positive minority is uncharacterized by design (no EC decode in
  this paper). The June decode suggested P-scavenging (alkaline phosphatase) among recruited
  specialists — sibling-paper territory, usable as a Discussion teaser ONLY if the annotator
  paper is cited as in-prep.
- **Trace-metal nulls at two levels**: no single dark family tracks V/Mo/Ni/Co/Mn/Zn, and a
  vanadium PERMANOVA block explains nothing (p=0.78). Pairs with the V-nitrogenase absence item
  (PI list 6, old-basis vnfD=0): vanadium chemistry is present in the sand, vanadium biology is
  not detectable. Cross-study: contrast with Atacama/Namib claims of metal-driven communities if
  the lit agent finds them.
- **Cautionary note (internal, not for the paper)**: June's Mo signal (177 families, 100%
  positive) vanished on the rebuilt universe while every other June pattern strengthened —
  small-N positive-tail associations on the old object were fragile; treat any other
  June-object-only curiosity with the same suspicion.
- **Salinity proxy behaves consistently badly across sections**: most spatially confounded in
  §2.4 (arsenic, ρ 0.24→0.12 partialled) and weakest-retained in §2.5 (41% vs 71–86% for its
  mineral components). A Methods/Discussion clause can note the composite proxy is dominated by
  between-site geography; per-element axes are the reliable objects.

## §2.1 reference bias (deposited 2026-08-16, session 5)

- **The companion amplicon paper reaches an apparently opposite conclusion, and the Discussion should
  say so explicitly rather than let a reader find it.** `empty-quarter-ecology` (`main.tex`,
  "Predicted functional potential and shotgun comparison") reports PICRUSt2-vs-shotgun agreement of
  per-sample ρ **0.743** and community-mean ρ **0.806** over **3,471 shared KOs** in **125** samples,
  concluding it "supports broad functional comparisons while leaving individual pathways as
  predictions". Read against our section title that looks like a contradiction. It is not.
  Their correlation is computed over the KOs the two methods **share** — it measures accuracy
  *conditional on a gene family already being inside the reference*. Our ceiling measures the size of
  what lies outside: a protein with no KO cannot enter that intersection at any accuracy. The two
  compose into one statement: **placement is fine, the annotated core is predicted well, and the
  majority of the gene content is not representable at all.**
- **We independently reproduce their result, which is the strongest form of this argument.** On a
  3.6× larger matched cohort (448 vs 125 samples) we get ρ **0.839** on shared KOs and NSTI **0.159**
  against their 0.165 — so we are not contesting their number, we are agreeing with it and then
  pointing at its denominator. Best rhetorical shape for the Discussion: concede the agreement first,
  in their own numbers, then give the ceiling.
- **Two distinct ceilings, worth keeping separate in prose.** (1) 55.6% of protein-coding genes carry
  no KO at all — the reference cannot describe them. (2) Of the KOs genomes *do* encode, 22%
  (2,517/11,374) lie outside the 10,543-family set PICRUSt2 can emit — a reference-vintage gap
  (eggNOG 2.1.12 vs PICRUSt2 2.4.1 defaults), not prediction error. The second is smaller but sharper,
  because it is a failure on genes we *can* name.
- **A predictor that returns a constant is worth one dry sentence.** PICRUSt2 emits all 10,543
  families in every one of the 448 samples; MAG-encoded richness varies (median 8,666, IQR
  8,169–9,003). Whatever else it estimates, it is not estimating richness. This is the honest residue
  of the withdrawn over-prediction ratio — keep the observation, never the ratio.
- **Cross-paper provenance is unusually tight and worth one line in Discussion or Data availability:**
  our PICRUSt2 input is byte-identical (sha256 `c71c2256…`) to the file the amplicon paper pins in its
  own manifest. The two papers disagree in framing while sharing the exact input bytes.
- **Standing caveat for the two-paper story:** the amplicon survey covers 1,217 non-control samples
  and the metagenome cohort 551; only 448 overlap. Neither paper's cohort is a subset of the other's,
  so per-sample claims should not be silently transferred between them.

## From §3.1 specificity (deposited 2026-08-17)

Deposited while compute was in flight; none of it depends on the final numbers. All cross-study —
**≤1 calibrating sentence each in Results**. Both key sources are preprints; don't let the framing
depend structurally on them.

- **Biome restriction of novel families is the norm, and conceding it is what makes our claim survive.**
  Aplakidou 2025 (bioRxiv 2025.12.21.695754): 385,791/608,258 = **63.4% of novel families are
  single-ecosystem**. So a large specificity fraction is unremarkable *as a fraction* — the reason §3.1
  was reframed. Discussion shape: concede the norm, then locate the contribution in *which* biome is
  missing from the references.
- **Desert convergence and our specificity concern different halves of the same proteome.** Leung 2025
  (bioRxiv 2025.11.29.691269; 25 locations, 75 metagenomes, 911 MAGs) finds deserts converge on
  composition and annotated gene-family abundance — measured on KEGG/CAZy/marker HMMs and rplB. "Dark
  matter", "UniRef", "MGnify", "no hit", "endemic" appear zero times. Deserts converge on the biology
  you can read while unreadable repertoires stay catalogue-specific. Two levers: their intro concedes
  deserts are "systematically underrepresented in global metagenome sampling" (better than asserting
  our premise ourselves), and Bay 2020 (mSystems 5:e00540-20) shows soil biogeography is strong at fine
  resolution — convergence is scale-dependent.
- **The reference-infrastructure gap is the payoff, at provenance level only.** Statable: no Arabian or
  Saharan site in any named desert comparator (Atacama, Mojave, JGI arid, Bay 2021, JGI Antarctica),
  MGnify90 dominated by aquatic/engineered biomes. Not statable: "absent from all catalogues" — MGnify90
  biome labels carry no geography. Connects §3.1 to the existing "discovery reservoir" bullet: a
  reservoir *because* the references don't contain it.
- **One method-honesty clause.** DIAMOND e-values scale with database size, so a fixed threshold is more
  permissive against smaller catalogues — the direction that would flatter a biome-matched-reference
  argument. Hence the raw cross-catalogue ordering stays descriptive and the claim rides on an
  e-value-corrected rarefaction curve. Saying so pre-empts the sharpest objection to Part 3.
- Armor already banked: the §2.4 biosafety control (2 pathogen species, 0.04%; VFDB 5,435→202) licenses
  "not contamination riding along with human material" — cite, never re-argue in Results.

## From the pipeline-guide rewrite (deposited 2026-08-17, session 7a)

**A tool default, not a chosen threshold, defines "absent" across §2.1–2.4.** Every community-level
absence statement in Part 2 rests on CoverM relative abundance > 0, and because
`--min-covered-fraction` was never passed, CoverM's genome-mode default silently set that floor at
**10% genome breadth**. The analysis is sound and the floor is a defensible one — but it was
inherited rather than selected, and it was invisible until the 2026-08-17 audit. Worth one honest
sentence in the limitations: the operational definition of detection is a mapper default, stated
explicitly (M.5) so readers can judge it, rather than a value we tuned. This is a small, concrete
instance of a general reproducibility point — that unstated defaults in widely used tools quietly
become the definitions a paper's claims depend on — and it is the kind of disclosure that costs
nothing and buys credibility with a methods-attentive referee.

**Related standing caveat, already registered:** the dark-family presence rule (`covered fraction
≥ 0.5` on any member gene) is length-biased in the opposite direction — easier to satisfy on short
genes — so any occupancy claim from `fam_presence.tsv` must be length-stratified. §1.4 survives that
stratification and the anti-relic-DNA inference is strengthened by it (see `CORRECTIONS.md`
2026-08-17). Both facts together make a fair point for the Discussion: read-mapping "presence" is a
threshold artefact as much as a biological statement, and we report which threshold and show the
claim holds without it.

### Added 2026-08-17 once §3.1's compute landed (the framing above did not depend on it; these do)

- **"Content beats depth" is now measured, and it is a statement about how reference databases should
  grow.** The five desert catalogues combined (147.5M sequences) leave less of our dark proteome
  unexplained than MGnify90 alone (717.7M): 42.60% against 52.68%. Corrected for size, the desert union
  is 24.33 points better than a size-matched random MGnify90 slice and SPIRE 21.30 points better, while
  the matched-control gaps are near 1 point. **Discussion value is generalisable, not self-congratulatory:**
  per sequence added, biome-matched sampling buys more reference coverage than generic depth. That is a
  recommendation the field can act on, and it is the strongest thing Part 3 produces.
- **The depth asymmetry sharpens the same point.** Going from 718M to 1.66G MGnify90 sequences resolves
  34.3% of the *control's* residual absences but only 5.8% of the dark set's. Deeper generic references
  keep improving coverage of proteins we already recognise, and barely touch the rest. Pairs naturally
  with the bullet above.
- **Independent concordance, and it argues for modesty.** The PI's own protein-novelty rerun searched
  299,927 UniRef50-dark family queries against the same 1.66G reference and found **48.1% still without a
  hit** (`rebuild_20260729/protein_novelty/mgnify_2026_07/`). Different object, different clustering, so
  not interchangeable with our 28.74% — but it establishes that roughly half of UniRef50-dark families
  survive the deepest reference available. Our fraction is therefore ordinary in magnitude, which is
  exactly what the reframe claims. Worth one sentence, and worth citing to the PI, whose own summary
  states our framing verbatim: "no hit means not detected in this versioned reference under these
  thresholds; it does not establish geographic endemicity."
- **A methodological caveat that is not only ours.** Calling a gene present when a fixed fraction of its
  length is read-covered is length-biased: swapping the rule for any-read presence reverses the relation
  between gene length and site occupancy in our data (median representative length 146 aa at the lowest
  occupancy against 25 aa at the highest under the covered-fraction rule, 76 aa against 175 aa under
  any-read). Any study deriving presence/absence from CoverM covered fraction over a gene catalogue
  inherits this, and short genes are exactly what dark-matter catalogues are made of. A calibrating
  sentence here is a service to readers and pre-empts the objection against us.

## Standing caveats surfaced while writing Methods (session 7b, 2026-08-17)

Deposited here because they are Discussion or limitations material, not Results. All three are things a
referee can reach on their own, so we are better off naming them first.

- **"Detected" in this paper is an inherited tool default, not a chosen threshold.** Every community
  number in §2.1–2.4 rests on CoverM relative abundance > 0, which under the never-passed
  `--min-covered-fraction` means **≥10% of the genome covered**. Three further parameters were also
  never passed (75 nt contig-end exclusion; no read-identity filter; no aligned-length filter). M.5
  now states this plainly rather than presenting 10% as a design decision. The honest framing for the
  Discussion: our absence language means "less than a tenth of the genome covered in that sample", and
  a study that chose a threshold deliberately would not necessarily choose this one. It is also worth
  noting that we *tested* the size analogue and it does not bite — median samples detected runs
  7/35/43/44/42 from <1.5 Mb to ≥5 Mb, flat above 2.5 Mb and *lower* for small genomes — so unlike the
  gene-level covered-fraction rule (see the bullet in the previous section), the genome-level rule
  carries no length artefact.

- **The paper draws taxonomy from two GTDB releases, and says so.** §1.1's rank census is R232
  (GTDB-Tk 2.7.2); the phylum labels behind §1.4's attribution, §2.3's phylum censuses and
  Halobacteriota turnover, and §2.4's phylum, genus and species matching all come from the earlier
  r226 classification, because the analysis scripts read that table. **1,898 of 5,438 representatives
  changed classification between the releases.** At phylum rank, where nearly all of the affected
  results sit, this is very likely immaterial. At **genus and species** rank it is not obviously
  immaterial, and exactly one claim depends on those ranks: the pathogen census's "2 named pathogen
  species, 0.04%". Whichever way that decision goes, the Discussion should carry one sentence
  acknowledging that a MAG catalogue's functional annotations and its headline taxonomy can drift out
  of release-sync during a long analysis, because this will be true of most large catalogues and is
  rarely disclosed.

- **The structural chain used two Foldseek builds, for a documented reason.** 3Di prediction, the
  triage search and the clustering ran on a GPU build (commit `8dc75c74`) because the packaged build
  crashes on GPU; the validation screen ran on packaged `foldseek/10-941cd33` to reproduce the control
  apparatus's calibrated invocation exactly. This is defensible and M.6 states it, but the defence
  belongs in the Discussion if a referee raises it: the two stages are separate measurements, and
  within the validation screen the 3,684 blinded known-fold controls are **co-queries under the
  identical build and command line as the targets**, so the control floor is build-matched to the
  quantity it calibrates. The triage's 68.0% no-hit rate is not, and is never reported as, the paper's
  novelty figure.

- **A small honesty point worth one clause somewhere.** Gene calling for the assembly protein catalogue
  ran over **552** assembly files against the **549** canonical assemblies — a catalogue-level basis
  difference that enters no per-sample cohort. M.4 states it. If a reviewer counts denominators, this is
  one of the places the counts will not line up, and it is cheaper to have said so.

## Deposited 2026-08-18 (session 8, citations sweep)

Cross-study comparisons surfaced while verifying claim support. **None belongs in Results.**
Evidence and sources: `resources/literature/CITATION_AUDIT.md`.

- **NHEJ prevalence needs a like-for-like caveat.** Our 74.5% sits against Sharda et al. 2020's **22%**
  across ~6,000 genomes (`survival-genome/Sharda_2020_bacterial_NHEJ_evolution_GBE.pdf`) — a 3.4-fold
  excess. Two conditions before leaning on it: (i) Sharda's 22% requires **Ku *and* LigD**; separately 2.3%
  carry Ku only and 10.3% LigD only, so a Ku-**OR**-LigD rule inflates prevalence by up to ~12 points —
  Methods must state which KOs define our 74.5%; (ii) Sharda finds NHEJ tracks **large genomes and slow
  growth**, i.e. exactly the Actinomycetota that dominate this catalogue, so the excess cannot be read as
  desiccation selection alone. A reviewer will reach for that alternative.
- **Our CO > H₂ ordering is the minority result.** Two desert MAG surveys on disk run the other way:
  Hutchinson 2025 (Antarctic) 322/554 MAGs with group 1h/1l vs **60/554** with coxL; Ray 2022 finds 1h/1l
  in about a third of cold-desert MAGs. Worth owning explicitly rather than letting a referee find it.
- **Group 1l dominance is not a universal arid rule.** It matches Ortiz 2021 (Mackay Glacier: 1l **2.3×**
  group 1h — close to our 2.5:1) and Jordaan 2020 (Gobi), but Hutchinson 2025 inverts it. Supports the
  existing decision to report the ratio as a description of the pool, not as an aridity signal.
- **The reference-coverage framing has a strong external anchor.** Aplakidou et al. 2025 (bioRxiv
  `10.64898/2025.12.21.695754`): **385,791 of 608,258 families (63.4%) are restricted to a single
  ecosystem**. Two nuances a referee will see: their "ecosystem" is GOLD's coarse 3-level/8-subtype
  scheme, and the family set was *defined* with a ≥5%-of-samples-from-one-ecosystem rule, so mild
  restriction bias is built into the denominator. Both make our concession safer, not weaker.
- **The soil H₂-sink budget gap is documentable.** Reji et al. (*Nat Commun*, 10.1038/s41467-025-67208-3):
  the soil microbial H₂ sink is "the largest and most uncertain term in the global H₂ budget", with a
  "paucity of observational constraints … particularly in moisture-limited, arid, and semi-arid soils";
  a drier threshold shifts 4–7 pp of the sink into arid/semi-arid regions. Ray et al. 2023 (*MMBR*,
  10.1128/mmbr.00048-23) says observations are "restricted to cold deserts", calling for work "extending
  to hot desert environments". This is the honest version of §2.2's current "effectively absent" claim
  (audit D-3) — **poorly constrained**, not absent.

## §1.3 structural novelty (session 6, 2026-08-18)

- **The headline for the Discussion is a methods lesson with teeth.** Our two TM normalisations give
  10.8× and 1.5× excess over their own control floors on the *same models, same references, same run*.
  A single "N novel folds" number is therefore not a property of a proteome — it is a property of the
  normalisation chosen. Any comparison we make to published fold-novelty counts has to state which
  normalisation produced them, and most do not. This is our strongest general contribution from §1.3
  and belongs in the Discussion, not in Results.
- **Novelty here is architectural, not local.** 96.7% of confident models contain an aligned region
  matching a known structure; only 63.7% are covered end-to-end. Candidates are longer (301 vs 215 aa)
  at equal confidence. The natural reading is multi-domain proteins whose parts are individually known
  but whose whole-protein architecture has no counterpart — i.e. novel *arrangements* of known
  building blocks. That is a weaker claim than "new folds" and a more interesting one for evolution:
  recombination of known domains rather than invention of new topology.
- **Cross-check worth stating once:** the PI's blinded known-fold controls, put through our own
  invocation, reproduced his false-dual-no-hit floor to the digit (101/3,013 qTM; 68/3,013 alnTM).
  Independent confirmation that the two pipelines are the same instrument, which is what licenses us
  to compare our current-basis result against his historical cohort at all.
- **Comparison to our own history (do not put in Results):** the June/August predecessor reported
  608 of 1,875 gated models on 68,527 old-basis candidates. Ours is 1,031 of 2,838 on a current-basis
  cohort with only 0.87% query overlap. The rates (32.4% vs 36.3%) are close despite almost disjoint
  inputs — worth one sentence as evidence the estimate is basis-robust, if the PI agrees.
- **Standing caveats for the limitations paragraph:** ESMFold is single-sequence (right tool for
  no-MSA dark proteins, but MSA-based refolding of the surviving candidates is the obvious revision
  request, and K is small enough to make that cheap); one representative per family speaks for its
  members; the controls were length/partiality-matched to the PI's candidate distribution, not
  re-matched to ours; the ≥3-member family admission is our choice and the pair-family sensitivity
  (+33,024 families) was not run.
- **Literature still owed** (📚 P1 in PLAN): AFESM 2025, Durairaj 2023, TED 2024. Read before drafting
  the Discussion paragraph so the comparison to published structural-dark-matter surveys is grounded
  rather than asserted.
