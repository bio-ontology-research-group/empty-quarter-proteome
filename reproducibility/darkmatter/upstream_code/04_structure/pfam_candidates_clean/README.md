# pfam_candidates_clean: Pfam-A 38.2 scan of the 1,031 candidate fold families

Owner: velazqal. Turns the paper's *inference* that the novelty of the candidate architectures is
architectural rather than local into a *measurement* at the sequence level: do the K proteins carry
known Pfam domains inside longer proteins whose remainder no domain covers?

All inputs are read-only. Nothing outside this folder was written.

## STATUS (2026-09-20): COMPLETE

| step | script | how it ran | wall |
|---|---|---|---|
| build the three protein sets | `scripts/01_build_sets.py` | login node | ~4 min |
| Pfam-A 38.2 scan, `--cut_ga` | `scripts/02_hmmsearch_pfam.sh` via `scripts/02_hmmsearch_pfam.sbatch` | Slurm **52136381**, 32 CPUs, 64G, batch | 00:06:40 |
| per-protein metrics and summaries | `scripts/03_summarise.py` | login node | ~20 s |

Job 52136381 COMPLETED, exit 0, MaxRSS 1.14 GB. Marker `results/PFAM_SCAN_COMPLETE`.

Exact command lines:

```
/home/velazqal/miniconda3/envs/python-utils/bin/python scripts/01_build_sets.py

sbatch scripts/02_hmmsearch_pfam.sbatch      # thin wrapper over 02_hmmsearch_pfam.sh, which runs:
/ibex/scratch/projects/c2014/rob/miniconda3/envs/bigslice/bin/hmmsearch \
  --cut_ga --cpu 32 --noali --notextw \
  --domtblout results/pfam_domtbl.txt --tblout results/pfam_tbl.txt -o results/pfam_hmmsearch.out \
  /ibex/scratch/projects/c2014/rob/metagenomics/rebuild_20260729/references/pfam_38_2/hmm/Pfam-A.hmm \
  sequences/all_sets.faa

/home/velazqal/miniconda3/envs/python-utils/bin/python scripts/03_summarise.py
```

Versions. HMMER 3.3.2 (Nov 2020), chosen to match the version the workbench Pfam reference was built
with (`rebuild_20260729/references/pfam_38_2/hmmer_version.txt`). Pfam release 38.2, 30,134 Pfam-A
families, dated 2026-01, based on UniProtKB 2025_03; marker `PFAM_38_2_REFERENCE_COMPLETE` present.
Python `/home/velazqal/miniconda3/envs/python-utils/bin/python` (3.x, standard library only).

## Inputs, all under `W = /ibex/scratch/projects/c2014/rob/metagenomics`

| input | rows used |
|---|---|
| `04_structure/00_structure_clean/05_validate/per_rep_structures.tsv` | 17,221 admitted families; `conf_qualified_dual_novel = 1` gives 1,031 (K), `confidence_qualified = 1` gives 2,838 (Q) |
| `04_structure/00_structure_clean/05_validate/fam_reps.tsv` | 17,221 rows, `family` to `chosen_rep`; the 1,031 K families declare 4,880 members (866 FOLD, 165 REUSE) |
| `04_structure/00_structure_clean/04_cluster/clusters.tsv` | 997,646 family/member lines; the 1,031 K family keys resolve to exactly 4,880 members, matching `n_members` for every family |
| `06_darkmatter/darkfam_basis_clean/results/dark_reps.faa` | 4,292,446 sequences; supplied all 4,880 K members and all 1,807 Q-not-K representatives, 0 missing |
| `04_structure/00_structure_clean/05_validate/per_control_structures.tsv` | 3,684 blinded controls; `confidence_qualified = 1` gives 3,013 |
| `rebuild_20260729/structural_calibration/preparation/controls/controls_blinded.faa` | 3,684 sequences; supplied all 3,013 confidence-qualified controls |
| `rebuild_20260729/references/pfam_38_2/hmm/Pfam-A.hmm` | 30,134 profiles |

`06_darkmatter/darkmatter_clean/02_uncovered_reps.faa.gz` does not exist; the sequence source is
`darkfam_basis_clean/results/dark_reps.faa`.

**Membership trap, resolved.** `per_rep_structures.tsv` is keyed on the *modelled* representative,
which is `fam_reps.tsv:chosen_rep`. The cluster key in `clusters.tsv` is `fam_reps.tsv:family`, and
the two differ for the 165 REUSE families. Joining `chosen_rep` straight into `clusters.tsv` finds
only 955 of the 1,031 families and 4,366 of the 4,880 members. The script resolves
chosen_rep -> family -> members and asserts the per-family count against `n_members`.

**Length note.** The scan uses the full protein sequence. For 1,022 of the 1,031 modelled K
representatives that length equals the `length` column of `per_rep_structures.tsv` exactly; the
other 9 were truncated to 1,022 aa for folding, so their sequence is longer than the model. Mean
sequence length of the K representatives is 302.9 aa against a modelled 300.6 aa, and of the
Q-not-K representatives 216.0 aa against a modelled 214.7 aa.

## Sets scanned (disjoint; one hmmsearch over the concatenation, 9,700 targets)

| set | n | median aa | mean aa | what it is |
|---|---|---|---|---|
| `K_members` | 4,880 | 256 | 287.3 | every member of the 1,031 candidate fold families |
| of which modelled reps (`K_reps_only`) | 1,031 | 265 | 302.9 | the representatives that carry the structural claim |
| `QnotK_reps` | 1,807 | 193 | 216.0 | confident dark models that *are* covered by a known structure |
| `controls` | 3,013 | 279 | 308.3 | confidence-qualified blinded PDB chain controls |

## Metrics (per protein, `results/per_protein_pfam.tsv`)

Domains are the `--cut_ga` rows of the domtblout, so every reported domain passes its family's
gathering threshold. Coverage merges the *envelope* coordinates (`env from`/`env to`);
`covered_frac` = covered residues / protein length; `max_gap` is the longest run of residues outside
every envelope, counting the N- and C-terminal tails; `has_gap80` is `max_gap >= 80`.
`n_dom_resolved` applies greedy overlap resolution, best bit score first, dropping a domain that
overlaps a kept envelope by more than 50 % of the shorter envelope.

## Results

### Per set (`results/set_summary.tsv`)

| | K_members | K_reps_only | QnotK_reps | controls |
|---|---|---|---|---|
| n | 4,880 | 1,031 | 1,807 | 3,013 |
| with at least one Pfam domain | **427 (8.75 %)** | 93 (9.02 %) | 101 (5.59 %) | 2,732 (90.67 %) |
| Wilson 95 % | 7.99 to 9.58 | 7.42 to 10.92 | 4.62 to 6.75 | 89.58 to 91.66 |
| with a domain **and** an uncovered stretch >= 80 aa | **382 (7.83 %)** | 83 (8.05 %) | 42 (2.32 %) | 477 (15.83 %) |
| of the Pfam-positive proteins, share with a >= 80 aa gap | **89.46 %** | 89.25 % | 41.58 % | 17.46 % |
| covered fraction of Pfam-positive proteins, median (Q1 to Q3) | **0.253 (0.160 to 0.365)** | 0.223 (0.130 to 0.338) | 0.559 (0.304 to 0.818) | 0.885 (0.750 to 0.954) |
| Pfam-positive proteins with covered fraction >= 0.8 | 1.41 % | 1.08 % | 28.71 % | 67.72 % |
| Pfam-positive proteins with covered fraction < 0.5 | 90.40 % | 90.32 % | 42.57 % | 8.09 % |
| median longest uncovered stretch, Pfam-positive proteins | 159 aa | 164 aa | 65 aa | 22 aa |
| resolved domain count 0 / 1 / 2 / >= 3 | 4,453 / 388 / 36 / 3 | 938 / 84 / 7 / 2 | 1,706 / 98 / 2 / 1 | 281 / 1,510 / 740 / 482 |
| Pfam-positive proteins with >= 2 domains | 9.13 % | 9.68 % | 2.97 % | 44.73 % |
| distinct Pfam families hit | 85 | | 61 | 1,840 |

Sensitivity to the threshold convention. The table counts a protein as Pfam-positive when at least
one *domain* clears the domain gathering threshold, which is the pfam_scan convention and the basis
of every coverage number here. Counting instead at the *sequence* gathering threshold
(`results/pfam_tbl.txt`, which admits sequences whose domains sum over the sequence GA without any
single domain clearing the domain GA) raises the shares to 11.13 % for K members, 7.14 % for
Q-not-K and 93.06 % for controls. The ordering and the conclusion are unchanged.

### Families (`results/k_family_pfam.tsv`, `results/k_family_summary.json`)

| quantity | value |
|---|---|
| K families | 1,031 |
| with at least one member carrying a Pfam domain | **144 (13.97 %)**, Wilson 95 % 11.98 to 16.22 |
| with a member carrying a domain and a >= 80 aa uncovered stretch | 130 (12.61 %) |
| whose modelled representative itself carries a domain | 93 (9.02 %) |

### Hit rate by length (`results/hit_rate_by_length.tsv`)

K members are longer than the Q-not-K representatives, so the raw hit rates are not comparable.
Stratified, K is above Q-not-K in every bin from 150 to 599 aa:

| length | K_members % hit | QnotK_reps % hit | controls % hit |
|---|---|---|---|
| 100-149 | 2.66 (n=188) | 7.49 (n=427) | 84.11 (n=472) |
| 150-199 | 8.10 (n=864) | 5.60 (n=554) | 85.14 (n=498) |
| 200-299 | 9.34 (n=2,121) | 4.52 (n=598) | 91.64 (n=658) |
| 300-399 | 10.68 (n=1,133) | 5.71 (n=140) | 93.04 (n=675) |
| 400-599 | 6.37 (n=455) | 1.39 (n=72) | 95.85 (n=506) |
| 600+ | 3.39 (n=118) | 12.50 (n=16) | 95.59 (n=204) |

### Top 20 Pfam families among K members (`results/top_pfam_K_members.tsv`)

The two most frequent are themselves domains of unknown function.

| rank | Pfam | id | n proteins | % of the 427 Pfam-positive K members |
|---|---|---|---|---|
| 1 | PF05305 | DUF732 | 43 | 10.07 |
| 2 | PF14018 | DUF4234 | 38 | 8.90 |
| 3 | PF13202 | EF-hand_5 | 28 | 6.56 |
| 4 | PF13490 | zf-HC2 | 25 | 5.85 |
| 5 | PF14520 | HHH_5 | 21 | 4.92 |
| 6 | PF12728 | HTH_17 | 20 | 4.68 |
| 7 | PF05103 | DivIVA | 19 | 4.45 |
| 8 | PF20901 | Sf6_terminase | 13 | 3.04 |
| 9 | PF01391 | Collagen | 13 | 3.04 |
| 10 | PF13240 | Zn_Ribbon_1 | 12 | 2.81 |
| 11 | PF27063 | RsmA_C | 11 | 2.58 |
| 12 | PF12836 | HHH_3 | 11 | 2.58 |
| 13 | PF13559 | DUF4129 | 10 | 2.34 |
| 14 | PF21983 | NikA-like | 10 | 2.34 |
| 15 | PF06568 | YjiS-like | 7 | 1.64 |
| 16 | PF07411 | DUF1508 | 7 | 1.64 |
| 17 | PF19102 | DUF5789 | 6 | 1.41 |
| 18 | PF04659 | Arch_fla_DE | 6 | 1.41 |
| 19 | PF02517 | Rce1-like | 5 | 1.17 |
| 20 | PF12802 | MarR_2 | 5 | 1.17 |

### Unknown-function share of the hits (`results/duf_share.tsv`)

| set | Pfam-positive | all hits are DUF/UPF | % of hits | carry a named domain | % of whole set |
|---|---|---|---|---|---|
| K_members | 427 | 124 | 29.04 | 303 | 6.21 |
| QnotK_reps | 101 | 25 | 24.75 | 76 | 4.21 |
| controls | 2,732 | 10 | 0.37 | 2,722 | 90.34 |

## What this licenses about "architectural novelty"

The multi-domain reading is **not** supported as a description of the K set. 91.25 % of the 4,880
candidate-family member proteins, and 90.98 % of the 1,031 modelled representatives, carry no Pfam
domain at all at the gathering threshold (88.87 % under the looser sequence-level threshold),
against 90.67 % of the confidence-qualified controls that do carry one.
Only 8.75 % of K proteins carry at least one Pfam domain, and only 6.21 % carry a domain that is
not itself a domain of unknown function. At family level, 86.03 % of the 1,031 candidate families
contain no member with any Pfam domain. Whatever locally familiar structure the
alignment-normalised TM-score sees in 96.7 % of confident models, it is a structural resemblance and
is not, for the large majority of these proteins, a recognisable sequence domain.

Inside the Pfam-positive minority the architecture picture is exactly the predicted one, and
strongly so. Of the 427 K proteins with a domain, 89.46 % still carry an uncovered stretch of at
least 80 aa, the median covered fraction is 0.253, 90.40 % are less than half covered, and the
median longest uncovered stretch is 159 aa. The same quantities in the controls are 17.46 %, 0.885,
8.09 % and 22 aa, and in the Q-not-K comparator 41.58 %, 0.559, 42.57 % and 65 aa. So where a known
part is present in a K protein, it sits in a protein most of which no known domain covers.

The K-versus-control contrast is partly built in and cannot carry the argument on its own: the K
proteins are sequence-dark by construction, and Pfam is a sequence-profile method, so a low Pfam
rate was expected. The informative contrast is K against Q-not-K, both dark and both
confidence-qualified. There K is the more domain-bearing set, 8.75 % against 5.59 %, and that holds
after length stratification in every bin from 150 to 599 aa. The K proteins are also the longer
ones (287.3 against 216.0 aa mean). That is the direction the multi-domain reading predicts, but
the effect is small in absolute terms.

What the paper may say, on this evidence: the candidate architectures are longer proteins whose
Pfam-annotatable parts, where any exist, cover a quarter of the chain and leave a long unmatched
remainder. What the paper may not say on this evidence: that the candidates are generally
assemblies of known parts. At the sequence level the parts are not known either.

Direction if this is to be pushed further: the negative is a sensitivity limit as much as a result,
so the test that would separate the two is a profile search with sensitivity beyond Pfam GA (for
example HHsearch against Pfam, or a structure-level domain parse of the ESMFold models against
ECOD or CATH) restricted to the 382 K proteins that already show a domain plus a long remainder.

## Files

```
scripts/01_build_sets.py          set construction and FASTA extraction
scripts/02_hmmsearch_pfam.sh      plain-bash Pfam scan worker
scripts/02_hmmsearch_pfam.sbatch  thin Slurm wrapper (32 CPUs, 64G, 12h)
scripts/03_summarise.py           per-protein metrics and all summary tables
sequences/{k_members,qnotk_reps,controls,all_sets}.faa
results/ids_*.txt                 the id lists for the four sets
results/k_family_members.tsv      modelled_rep, family_key, member (4,880 rows)
results/set_membership.tsv        protein, set, is_family_rep, length (9,700 rows)
results/pfam_domtbl.txt           hmmsearch --domtblout, 7,575 domain rows
results/pfam_tbl.txt              hmmsearch --tblout, 8,100 sequence rows
results/per_protein_pfam.tsv      the per-protein metrics
results/set_summary.tsv           the per-set table above
results/covered_fraction_hist.tsv covered-fraction deciles per set
results/hit_rate_by_length.tsv    length-stratified hit rates
results/duf_share.tsv             unknown-function share of the hits
results/top_pfam_{K_members,QnotK_reps,controls}.tsv
results/k_family_pfam.tsv         per-family rollup for the 1,031 K families
results/k_family_summary.json     the family-level headline numbers
results/missing_sequences.txt     empty; every requested sequence was found
logs/                             build log, summary log, Slurm out/err
```
