# Dossier 07 — darkness as a lineage property, and the Thermoplasmatota radiation

STATUS: EXPLORATORY, UNREGISTERED. Everything below comes from `/ibex/scratch/projects/c2014/alelopezv/eqmags_poc/` (01_dark_by_lineage, 03_genome_size_function, 05_rakpr_overlap), whose README says: "Numbers here are NOT registered; nothing enters `notes/CATALOGUE_NUMBERS.md` or `main.tex` from this folder until re-run and registered." No sentence of this may be drafted as it stands. Read with `notes/PAPER_STRUCTURE_PROPOSED.md` (item R4, blocker 16). Marks: VERIFIED = I checked it on disk this session; INFERENCE = mine, not in any README.

## Questions

1. Does lineage or assembly quality set a genome's dark fraction? (01)
2. Where does the dark MASS sit, once censuses are size-normalised? (01 goal C, 03)
3. Is darkness an artefact of fragmented MAGs, or does it survive in closed long-read genomes? (05 step three)

## Basis

| test | key | catalogue | n joined |
|---|---|---|---|
| 01 round one (models, residuals) | RENAME `<s>__self_SemiBin_<n>` → `selfsup__<s>__SemiBin_<n>` | per_bin_dark, 37,949 bins | **5,359 of 5,438** |
| 01 round two, all `v2_*` | NATIVE: contig → winner bin → cluster rep | Wdb, 5,438 dRep winners | **5,438 of 5,438** |
| 03 parts a–c (KO vs size) | native Wdb key (rep_ko_v2) | 5,438 | 5,438 |
| 03 part d (dark vs size, M4/M5) | RENAME, reads `01/results/reps_dark.tsv` | per_bin_dark | **5,359** |
| 05 step one (skani) | Wdb `.fa` basenames | 722 cultured × 5,438 | 203 shared pairs |
| 05 step three | catalogue ids; dark reps from `dark_reps.faa` | 651 Bakta genomes; 4,292,446 dark reps | 132 shared-species genomes |

THE CRITICAL DISTINCTION. 5,359 is a **rename, not a key**. It dropped 79 reps (44 `48PRr1__SemiBin_<n>` with no assembly token, 35 unmatched `_sub`) and — proven in round two — **for all 975 `_sub` reps it read the un-suffixed PARENT bin**: every one a strict superset, bigger by a median of 11 proteins (0.34 %), positive in 975/975. Round one's caveat had concluded the opposite ("the rename maps the same bin and not a parent bin"); round two overturned it. The *fraction* survives (rho 0.9991, median |diff| 0.0002); the *counts* do not. Round one's headline table and 03's M4/M5 sit on 5,359; every `v2_*` result sits on 5,438. INFERENCE: 03 part d is the only surviving analysis still on the rename.

## How it ran

**Join** (`v2_scan_cluster90.py`, 8 workers, 26 min): protein id `<sample>::<contig>_<gene>` → contig (id up to the last `_`) → winner bin via `contig_to_bin_refined.tsv.gz` cut to Wdb (2,485,945 contigs, none shared) → 90/80 cluster rep via `all_cluster90.tsv` (753,305,533 lines, 8 byte ranges) → darkness as a *cluster* property from `dark_rep_ids.txt.gz` (69,938,628 ids). Member logic copied from `darkclean_5_mag.sh`. `v2_boundary_check.py`: 0 lines lost at the 7 splits.

**Variance decomposition** (OLS, statsmodels 0.14.2): A = completeness + contamination + log10 N50 + log10 length; B = A + C(phylum); C = B + C(genus_grp) (genera < 5 reps pooled as `other`, blank `g__` held as one level). The lineage-first order is fitted too; unique-lineage = R²(full) − R²(A), unique-quality = R²(full) − R²(phylum+genus), shared = the remainder. An R² increment is only the extra variance a block explains *given what is already in the model*; with 51 phylum and ~215 genus levels there is no cross-validation and no permutation (adjusted R² is reported).

**Novelty gradient**: `first_unnamed_rank` = deepest Latin rank + 1 (phylum=1 … species=6; 7 = fully named). Placeholder = contains a digit, or fails `^[A-Z][a-z]+(_[A-Z]{1,3})?$`. Every placeholder call was digit-driven.

**Novelty-matched Thermoplasmatota contrast**: the stratum is just all reps sharing one `first_unnamed_rank` value. 236 of 237 Thermoplasmatota sit at `class`, so the contrast is that single stratum — 236 vs 324 others, Mann-Whitney on raw `dark_fraction`. VERIFIED: **all 324 comparators are Bacteria** (222 Actinomycetota, 28 Chloroflexota, 14 Myxococcota, …) and the stratum is **not** completeness-matched (Thermo median 85.41 vs 88.73).

**Composition enrichment, two denominators**: `enr_fam` = share of the 77,069 single-phylum dark families ÷ share of the 5,438 species reps (a **genome** denominator; numerator from §1.4's registered `fam_phylum.tsv`). `enr_prot` = share of 1,600,067 dark protein instances ÷ share of 20,379,930 protein instances (a **protein** denominator, both sides from the v2 scan). Abundance is a third axis: each sample column renormalised over the 5,438 reps, median across 549 samples (F5Dr1 all-zero, dropped) of a 550-column PI file, not the 551 canonical samples.

## Controls

REAL CONTROLS (a second measurement that would fail if the claim were wrong):
- **Matched annotated-family control**: 106,042 dark family reps + 106,042 length- and family-size-matched annotated reps from §3.1's `query_all.faa`, in the same DIAMOND run — **46.5 % vs 0.16 %, a 290-fold gap**. Proves the instrument works, so low cultured overlap is lineage restriction, not a failed search.
- **Length distributions**: cultured dark proteins are *shorter* than MAG dark reps (median 60 vs 70 aa) in genomes of median 1 contig, so shortness belongs to dark proteins, not to fragmentation.
- **Reciprocal coverage ≥ 50 %** (qcovhsp and scovhsp, applied post hoc so the relaxed count stays available).
- **Identical UniRef50 rule**: `darkclean_3b_uniref.sh` verbatim (MMseqs2 18.8cc5c, `-s 7.5 -e 1e-3 --max-seqs 10`, UniRef50 2026_01, dark = no hit), so cultured and MAG darkness are called by one rule.
- **Same-species matching** at ANI ≥ 95 and aligned fraction ≥ 30 (skani 0.3.1), so the comparison is within species.

NOT CONTROLS, ONLY COVARIATES (they adjust a model; they do not test the claim): completeness, contamination, log10 N50, log10 length enter Model A linearly and additively, and nothing tests a non-linear completeness response inside Thermoplasmatota; 03's size+quality residual column is a sensitivity. JUDGMENT: the **novelty-matched stratum** is closer to a control than a covariate, but as built it is neither completeness-, abundance-, prevalence-, size- nor domain-matched — treat it as one covariate held fixed, not as the matched null the archaeal claim needs. The boundary audit and the round-one/round-two reconciliation are diagnostics, not controls.

## Numbers (exact, from the result files; ALL UNREGISTERED)

**01 round two, n = 5,438** (`v2_model_report.txt`, `v2_scan_report.txt`): 20,379,930 proteins, 1,600,067 dark, global 0.078512; median 0.0631, IQR 0.0365–0.1038; 5,438/5,438 computed, 0 empty; 1,552,669 dark clusters touch a winner. Ladder A 0.1147 → B 0.7680 (+0.6533) → C 0.8347 (+0.0668); phylum alone 0.7394. Partition: unique lineage 0.7201, unique quality 0.0203, shared 0.0944, unexplained 0.1653. Model A slopes log10 length −0.09312, log10 N50 −0.02282, contamination +0.00178, completeness +0.00034. Novelty: N = A + first_unnamed_rank 0.4955 (+0.3808); N + phylum 0.8356 (+0.3401 — the phylum increment **survives 52.1 %**); N + phylum + genus 0.8749. Gradient medians: phylum 0.2170 (n=10), class 0.1745 (560), order 0.1084 (272), family 0.0828 (1,006), genus 0.0652 (1,613), species 0.0327 (1,855), fully named 0.0170 (122). Named genus 2,030/5,438 (37.3 %), placeholder 1,912, blank 1,496. Named vs placeholder within phylum: Actinomycetota 0.0215 vs 0.0751 (p 2.13e-138), Pseudomonadota 0.0303 vs 0.0514 (1.1e-56), Bacteroidota 0.0526 vs 0.0683 (1.24e-12); Thermoplasmatota untestable, 0 named genera.

**Thermoplasmatota**: 237 reps, median dark 0.309565, median Model-A residual +0.214914, 0 % named genus, median completeness 85.42 (64.1 % below 90), median length 3.00 Mb, median 3,302 proteins; 236/237 above the bacterial 95th percentile (0.1815). Matched stratum 236 at 0.3101 vs 324 at 0.1249, MWU p 7.17e-90. Archaea overall 591 reps (Halobacteriota 269, Thermoplasmatota 237, Thermoproteota 56, …), median dark 0.1009 vs bacteria 0.0619. Composition: species 4.36 %, abundance 0.81 %, dark families 19.84 %, protein share 3.79 %, dark proteins 15.21 %, **enr_fam 4.5519, enr_prot 4.0120**. Actinomycetota 30.03 / 52.79 / 32.76 / 29.85 / 21.46, enr 1.0908 / 0.7188. Dark families attributable 78,727, single-phylum 77,069 (97.9 %) = 72.7 % of the 106,042 universe.

**Reconciliation**: rho 0.9998 against round one; plain reps 4,384/4,384 identical; the 79 dropped reps unremarkable (44 no-token median 0.0705, 35 `_sub` median 0.0699). GFF cross-check: **exact protein-count agreement in 1 of 20 probe reps** (delta −40 to +60, id_match 97.0–99.9 %), because those proteomes were re-called with pyrodigal in Single mode.

**01 round one, n = 5,359 (SUPERSEDED, do not quote)**: A 0.1140 → B 0.7658 (+0.6518) → C 0.8325; partition 0.719 / 0.021 / 0.093 / 0.167; +phylum F 396.0 on 50 df, +genus F 9.42 on 215 df, p < 1e-235; Thermoplasmatota +0.215 on n = 230.

**03**: KO richness median 1,252 (150–2,864); corrected size median 4.03 Mb; raw 3.45 Mb. M1 R² 0.5163 (+1,370.9 KOs per 10-fold, +413 per doubling); M1raw 0.7049; M2 0.7596 (+0.2433); M3 0.9099 (+0.1503). M2 slopes size +1,347.151, completeness +14.638, contamination −4.173. Size-only phylum residuals: Bacillota +250.9, Pseudomonadota +238.3, Deinococcota +156.6; Halobacteriota −298.0, Thermoplasmatota −252.8, Planctomycetota −227.3. Size+quality: +233.0, +188.0, −214.4, −193.9; rank agreement Spearman 0.917. **M4 dark ~ size slope −0.09459, R² 0.0659, p 2.08e-81, n = 5,359**; M5 +completeness slope −0.09248, R² 0.0929; Spearman −0.1174. Bigger genomes are LESS dark.

**05 step one**: 189/722 cultured (26.2 %) share a species; 25 of 5,438 MAG species (0.46 %); 203 pairs; skani 0.3.1, job 51690633, 29 s, `search` and `dist` agree on 601 pairs. Enrichment 100/520, isolate 32/131, thermophilic soil 57/71 (80.3 %).

**05 step three**: 651 Bakta proteomes, 2,710,913 proteins, 269,871 hypotheticals. Cultured dark: all 651 46,337/2,710,913 = **1.7093 %**; 132 shared-species 5,194/525,525 = **0.9883 %**; 519 unshared 41,143/2,185,388 = 1.8826 %; hypotheticals 45,804/269,871 = 16.97 % (98.8 % of cultured dark proteins are Bakta hypotheticals). Matched MAG reps 6,267/535,613 = **1.1701 %** pooled, median 1.026 %. 1,961 of 5,194 (37.76 %) cultured dark proteins recovered as MAG-dark reps. Of 4,292,446 MAG-dark reps: 3,988 (0.0929 %) with a reciprocal-50 % cultured homolog, 25,365 (0.5909 %) relaxed, 1,275 in a shared-species genome (648 at ≥ 90 % identity), 2,713 only in unshared genomes; ≥ 100 aa 1,520/1,466,943 (0.1036 %). Controls 49,291/106,042 (46.4825 %) annotated vs 170/106,042 (0.1603 %) dark. Lengths: MAG dark median 70 aa, 65.83 % < 100 aa, 42.41 % < 60 aa, 56.40 % incomplete calls; cultured dark median 60 aa, 77.60 % < 100 aa, 49.28 % < 60 aa.

## The cultured comparison's reach

VERIFIED from `05_rakpr_overlap/results/shared_species.tsv`: the 25 shared MAG species are **11 Pseudomonadota, 11 Bacillota, 2 Actinomycetota, 1 Bacteroidota — all `d__Bacteria`, zero archaea**. Also VERIFIED from `/ibex/scratch/projects/c2014/alelopezv/eqmgc/catalogue.tsv`: **all 722 cultured genomes are `d__Bacteria`**, and the 651 with Bakta span six bacterial phyla (Pseudomonadota 292, Bacillota 288, Actinomycetota 49, Bacteroidota 15, Verrucomicrobiota 5, Spirochaetota 2). The catalogue holds 591 archaeal reps, 237 of them Thermoplasmatota. **The assembly-artefact control therefore cannot touch the archaeal claim at all.** It licenses exactly one sentence — in the species where closed long-read genomes exist, darkness is reproduced — for cultivable, well-referenced bacteria whose per-genome dark rate (~1 %) is a sixth of the catalogue median (6.3 %). Presenting it as support for the Thermoplasmatota section would be a misrepresentation a referee will catch on sight.

Also VERIFIED: under step3_04's best-ANI rule the 132 shared-species genomes resolve to only **23 distinct MAG reps** (25 across all pairs; `F27PRr2__self_SemiBin_29` and `S34PRr1__self_SemiBin_99` are never a best hit). So "0.99 % vs 1.17 %" and "cultured lower in 107 of 132 pairs" are **pseudo-replicated** — the MAG side counts 23 genomes up to 132 times. Pooled over the 23 distinct reps the MAG rate is 1.7315 % (1,513/87,381), median 1.287 %: the gap widens, but the honest n is 23.

## The unit mismatch

MAG darkness is a property of a **90/80 linclust representative**; cultured darkness is **per protein**. Five rates coexist, on four denominators:

| rate | numerator / denominator | unit | source |
|---|---|---|---|
| 15.567 % (1 % sample 15.553 %) | 69,938,628 / 449,265,306 | 90/80 cluster reps, whole catalogue | POC 02 |
| 9.885 % (sample 9.915 %) | 4,292,446 / 43,423,322 | 90/80 cluster reps, MAG-mapped (37,949 HQ+MQ bins) | POC 02 |
| 7.851 % (per-genome median 6.31 %) | 1,600,067 / 20,379,930 | protein instances on the 5,438 winners' contigs | 01 v2 |
| 1.709 % | 46,337 / 2,710,913 | Bakta proteins, 651 cultured genomes | 05 step3 |
| 0.988 % vs 1.170 % | 5,194 / 525,525 vs 6,267 / 535,613 | Bakta proteins vs matched-rep protein instances (23 genomes, 132 pairs) | 05 step3 |

Step three's denominator (4,292,446) is on the **37,949-bin** universe while the per-rep fractions it compares against are on the **5,438-winner** universe. INFERENCE: compatible only because darkness is a cluster property; the section must state which universe each number is on.

## Caveats (quoted verbatim from the workbench READMEs)

- "Two naming bases: the 5,359 rest on a rename, not a native key."
- "round one read the un-suffixed PARENT bin and not the same bin, though the fraction is unaffected"
- "Dark means no UniRef50 homolog, so the signal partly tracks how well each lineage is covered by the reference, which is itself the point to argue."
- "Novelty is near constant inside Thermoplasmatota (236 of 237 unnamed at class), so the ordinal is partly that phylum's own indicator and the shrinking phylum term is expected rather than proof of tautology; the matched-stratum test is the evidence."
- "Completeness spans 50.05 to 100 and 2,763 reps (50.8 %) sit below 90 %, so MQ genomes carry the analysis"; "27.5 % of reps have a blank GTDB `g__` … so the genus increment is a floor."
- "phyla with low median completeness (Halobacteriota 79.7, Chloroflexota 85.4, Thermoplasmatota 85.4) read as KO-poor partly for that reason … A completeness-matched repeat would settle it."; "only the extremes are safe to name."
- "KO richness is annotation-limited, so poorly referenced lineages lose twice, in KO calls and in the UniRef50 dark measure."
- "The units differ too: MAG darkness is a property of a 90 % linclust representative, cultured darkness is per protein."
- "Gene calling differs: Bakta on cultured genomes, pyrodigal on assemblies, so the per-protein sets are not identical objects."
- "the all-651 rate of 1.71 % is not comparable to the catalogue-wide MAG rate because the cultured set is well-annotated genera."
- "contamination is a HYPOTHESIS to test, not a finding, before any recruitment step." (four human-associated shared species, incl. *E. coli* in 316 of 550 samples)

## What must be rerun and registered before writing

1. **Native-key rerun of 03's dark-vs-size test (M4/M5)**: repoint `size_function.py` from `01/results/reps_dark.tsv` (5,359, parent-bin counts) to `01/results/v2_rep_dark.tsv` (5,438). ~3 lines, **< 1 min**. Until then, "darkness is not size in disguise" rests on the superseded basis.
2. **Rarity-, abundance-, prevalence-, completeness- and size-matched null for Thermoplasmatota.** All inputs are on disk (`v2_reps.tsv`, `rep_relative_abundance.tsv`), so the compute is **minutes**; `PAPER_STRUCTURE_PROPOSED.md` budgets **~4 h** for the whole item. This is the gate — without it Maeke 2025 pre-empts the claim.
3. **Domain-matched (or domain-declared) contrast**: the matched stratum's 324 comparators are all bacterial. Add Halobacteriota/Thermoproteota comparators, or state the result as archaea vs bacteria. **Minutes.**
4. **Completeness-matched repeat of 03's KO residuals**, as 03's own caveat requests. **~1 h.**
5. **De-pseudo-replicate the cultured comparison**: one row per MAG rep (23), paired test on 23 not 132, report the distinct-rep pooled rate. **Minutes.**
6. **Settle and register the protein-count basis** (contig route vs GFF route; 19/20 probes disagree). Reasoned in the README, never registered. **No compute.**
7. **Resolve the *E. coli* rep** `V42PRr2__single_SemiBin_33` (POC 07 blocker 1). If it leaves the 5,438 the whole v2 scan reruns: **26 min** plus minutes of modelling.
8. **Thermoplasmatota slice of the existing structural no-hits** (**~2 h**, no new search). It would move the claim from annotational to structural, i.e. it changes what the section *is* — run it before drafting.
9. **Register survivors** in `notes/CORRECTIONS.md` and `notes/CATALOGUE_NUMBERS.md` with basis and date. Note the dark-family share is *partly registered already*: `main.tex` §1.4 states "Actinomycetota (32 %) and Thermoplasmatota (19 %)" of 78,727 attributable families. The 19.84 % / 4.55-fold *enrichment* framing is not.

## Reviewer attack surface

- **Maeke 2025, Microbiome, 10.1186/s40168-025-02140-8** attributes high unknown-gene content in a sister Thermoplasmatota class to **MAG rarity**. It is the most dangerous paper for this claim, it is **not in `paper/references.bib`** (VERIFIED), and the null that would answer it does not exist. Do not draft before item 2.
- **Vanni 2022 eLife (AGNOSTOS)** and **Rodriguez del Rio 2024 Nature (FESNov)** already establish lineage-specific unknown space and its concentration in uncultivated taxa; **Pavlopoulos 2023 Nature** is the planetary census. All three are already cited in `main.tex`. The novelty gradient alone restates them — only the residual phylum term (+0.3401 after novelty) and the mass attribution are new.
- **Circularity**: dark = no UniRef50 homolog, so a lineage absent from UniRef50 is dark by construction. The only answers on hand are the surviving 52 % phylum term and the matched stratum, both attackable as confounded.
- **Confounds in the matched stratum**: not completeness-matched (85.41 vs 88.73), not abundance- or prevalence-matched, entirely cross-domain.
- **Phylogenetic non-independence**: 5,438 reps enter OLS as independent draws; 237 Thermoplasmatota from a few unnamed classes are not. No mixed model, no permutation, no cross-validated R².
- **The cultured control does not cover archaea**, and is pseudo-replicated (23 genomes, 132 pairs).
- **Gene-caller mismatch** (Bakta vs pyrodigal) and the contig-vs-GFF count disagreement.
- **"Radiation"** is unsupported: no phylogenomic tree of the 237 exists.

## Open questions for Alejandra

1. Is the canonical basis the 5,438 native key, and is the *E. coli* rep in or out before anything is registered?
2. Results section or Discussion paragraph for the archaeal finding, given that neither the matched null nor any archaeal cultured control exists today?
3. Write the cultured comparison at all, narrow it to "23 well-referenced bacterial species", or hold it back?
4. Does the structural slice (item 8) run before drafting, since it may change the section from annotational to structural?
5. Who runs items 2 and 3, and does the 19 % dark-family share already in `main.tex` get reused or recomputed on v2?
6. May the word "radiation" be used before a phylogenomic tree exists?
