# Dark-family pipeline — concepts in plain words (study notes)

*For defending §1.4/§2.5/§3.1. Written 2026-08-14, session 01. Numbers = current rebuild.*

## The one-paragraph big picture

The paper claims a large "dark" proteome (no database match). A skeptic's cheapest attack:
"those are not real genes." §1.4 answers with four independent measurements on 106,042
dark protein families, each killing a different artifact explanation:

| Evidence | Kills | Current number |
|---|---|---|
| Recurrence (samples/sites) | random assembly junk | median family in 12 samples / 9 sites |
| dN/dS < 1 vs control | non-coding ORFs (real DNA, fake gene) | computing |
| Abundance/prevalence (read mapping) | rare relic/dead DNA | computing |
| Persistence (5 expeditions) | transient contamination | 28.9% in all 5 trips |

---

## 1. dN/dS — "does mutation-editing show selection at work?"

**The idea.** The genetic code is redundant: some DNA mutations change the protein
(non-synonymous), others are silent (synonymous — different codon, same amino acid).
Compare two related gene copies and count both kinds, each divided by its opportunities:
dN = protein-changing rate, dS = silent rate.

- If the sequence is junk (no working protein), nothing filters mutations: both kinds
  accumulate equally → **dN/dS ≈ 1**.
- If the protein matters, individuals with protein-breaking mutations die out; silent
  changes survive fine → protein-changing differences are missing → **dN/dS < 1**
  ("purifying selection").

**How we compute it.** Within each family (≥3 diverged members): align the proteins
(MAFFT), thread the DNA back onto that alignment codon-by-codon, then for each pair count
dN and dS (NG86 counting method). Filters: ≥30 comparable codons (enough signal),
dS ≤ 2 (beyond that, silent sites are saturated and the count is unreliable), ≤80 pairs
per family (compute cap). Same code as the June run — only inputs changed.

**Why the control.** Families were *defined* by protein-level clustering, which itself
preselects sequences whose protein is conserved — pushing dN/dS down mechanically for
anything, real or not. So we run the *identical* pipeline on annotated (non-dark) proteins,
sampled to have the same length profile and matched family sizes. The claim is comparative:
dark genes behave like known-real genes *measured with the same ruler*.

**Defense line:** *"dN/dS below 1 means mutations that break the protein are being removed
from the population. Dead DNA does not get edited. Our dark families show this signature in
the same pipeline where annotated genes show theirs."*

## 2. Read mapping / recruitment — "is it there, and how much?"

**The idea.** Assembly only reconstructs a gene where coverage is high, so "assembled in
3 samples" understates presence. Mapping is the sensitive instrument: take each sample's
raw reads and align them to our 536,437 family-member genes. Two outputs per gene per
sample: **covered fraction** (are the gene's positions actually covered by reads → the
organisms carrying it live in this sample) and **TPM** (what share of the sample's DNA is
this gene, normalized by gene length and sequencing depth so samples are comparable).

**Member-aware = the honesty fix.** Short reads only align to ~90%+ identical sequence.
The old run mapped against ONE variant per family, so samples carrying a *diverged* member
looked family-negative. We map against every member and then sum members to the family.

**What it feeds.** (a) §1.4 leg 3: families are prevalent and abundant, not a rare tail;
(b) the recruitment-based persistence number (comparable to the old 70.3%); (c) the
family × sample matrix that §2.5's environment analysis runs on. This one job serves
three sections.

**Defense line:** *"Presence is measured by reads, not by assembly success, against every
variant of every family — the sensitive and unbiased version of the measurement."*

## 3. Robustness experiment (30% recluster) — "is the universe knife-edge?"

Families were built at 50% identity / 80% coverage — a standard but *chosen* threshold.
A reviewer asks: "would your 106k families evaporate at another threshold?" We recluster
the same 4.29M dark proteins at 30% identity and report the family count. Expected:
fewer, bigger families, changing smoothly — showing the universe is threshold-stable and
50/80 is the conservative choice (it *fragments* remote relatives, so counts err low).

**Defense line:** *"Family counts move smoothly with the clustering threshold; nothing in
our conclusions depends on the exact value, and our choice is the conservative one."*

## 4. Phylum attribution — "who carries the dark matter?"

Chain each family member: protein → its contig → the refined bin containing that contig →
that bin's species representative (dRep) → GTDB taxonomy. Result: **97.9% of attributable
families sit within a single phylum**, and the carriers are the dominant, ordinary lineages
(Actinomycetota 32%, Thermoplasmatota 19%). Two claims fall out: dark families are
lineage-faithful (consistent with vertical inheritance, not rampant transfer or chimerism),
and dark matter is the normal biology of the community's abundant organisms — not an exotic
fringe. Coverage caveat: 74% of families have a member inside the refined catalogue and can
be attributed; the rest are catalogue-level only.

**Defense line:** *"The novelty is carried by the same ordinary, abundant phyla that
dominate the community — dark matter is their unexplored biology, not contamination or
exotica."*

## 5. Recurrence and persistence — "independently observed, again and again"

A protein counts as recurrent when near-identical copies (≥90%) of it were independently
*assembled* in different samples — random artifacts don't replicate. Sites are the truly
independent unit (samples within a site share a community), and expeditions add the time
axis: 28.9% of families were re-assembled in **all five expeditions**, years apart — dead
DNA and transient contamination don't do that. The stricter cousin: a family's ≥3 *diverged*
members (each <90% identical to the others) cannot be produced by any assembly artifact —
only by evolution.

**Defense line:** *"Every family was observed as at least three diverged variants, across
at least three samples and two sites — an artifact would have to independently evolve."*

---

## The number funnel (memorize this spine)

```
551 samples → 753M proteins → 449M reps (90% cluster) → 69.9M dark (no UniRef50 hit)
   → 4.29M genome-resolved dark (on HQ+MQ-binned contigs)      [§1.2-§1.3's object]
   → cluster 50/80 → 3.66M clusters → ≥3 members: 108,826
   → +≥3 canonical samples & ≥2 sites: 106,042 dark families   [§1.4/§2.5/§3.1's object]
      (paralog-only dropped: 2,784 · young 1–2-variant recurrent class: 782,616)
```

One dark definition, one namespace, top to bottom — that is the sentence that makes the
whole of Part 1 defensible.
