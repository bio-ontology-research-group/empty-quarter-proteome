# contamination_clean: are the culture-shared MAG species contamination, and what do the blanks say?

Owner: velazqal. Promoted 2026-09-20 from the exploratory folder
`alelopezv/eqmags_poc/07_contamination_check`, then extended with the two shotgun control libraries,
which had never been profiled against the catalogue. This folder is authoritative; the exploratory
copy is superseded.

## STATUS (2026-09-20): promotion reproduces exactly, both controls profiled, all questions answered

| step | script | how it ran | wall |
|---|---|---|---|
| abundance statistics | `01_contamination_tests.py` | login | ~2 min |
| within-trip biomass | `02_biomass_within_trip.py` | login | seconds |
| census + control inventory | `03_census_and_controls.py` | login | seconds |
| CoverM, both controls vs the 5,438 | `04_coverm_controls.sbatch` | Slurm **52135671** array 1-2, 16 CPUs, 250G | 00:01:54 / 00:10:28 |
| GTDB-Tk R232 + skani + contig mapping | `05_identify_controls.sbatch` | Slurm **52135758**, 32 CPUs, 300G | 00:03:45 |
| blank genes vs UniRef50 | `06_sntc_uniref.sbatch` | Slurm **52135781**, 16 CPUs, 256G | 00:09:21 |
| tables | `07_summarise_controls.py` | login | seconds |
| blank counts, breadth gate off | `08_coverm_sntc_counts.sbatch` | Slurm **52136180**, 16 CPUs, 250G | 00:03:41 |
| check of the SPositive zero | `09_sanity_spositive_zero.sbatch` | Slurm **52136296**, 16 CPUs, 250G | 00:10:01 |

Job **52135929** is the first attempt at step 08 and failed in 5 s: CoverM refuses `-m count` while a
covered-fraction gate is set. The rerun sets `--min-covered-fraction 0`, which is why step 08 is a
separate ungated run and not a column of step 04.

**Promotion: no number moved.** All six `results/*.tsv` and `figures/abundance_by_trip.png` are
byte-identical to the exploratory outputs (`diff`, `cmp`). Only the output root changed, at the top of
each script; there is no shim. E. coli 316/550 (57.46 %), trips 72.4 / 56.8 / 23.6 / 25.9 / 100.0 %,
DNA rho -0.4900, replicate-concordance excess 0.0045, census 2 strict and 5 suffix-blind of 5,438,
22 sequenced libraries outside the canonical 551. The four-species table stands as written.

Python `/home/velazqal/miniconda3/envs/python-utils/bin/python` (pandas 1.5.0, numpy 1.26.1,
scipy 1.11.3, matplotlib 3.9.4). CoverM 0.7.0, minimap2 2.31-r1302, samtools 1.23.1 from
`rob/miniconda3/envs/eco_probes`. GTDB-Tk 2.7.2 (module `gtdb-tk/2.7.2`) against the project's own
R232 reference. skani 0.3.1. mmseqs2 and prodigal as in the scripts.

## The two controls

| | SNTC (no-template) | SPositive |
|---|---|---|
| reads before fastp | 610,596 | 125,688,774 |
| reads after fastp, paired | 335,190 | 119,683,914 |
| bases after fastp | 29,757,805 | 11,765,177,280 |
| Q30 after / duplication | 81.43 % / 2.37 % | 85.88 % / 2.74 % |
| reads assigned by CoverM | **28,694 (8.56 %)** | **0 (0.00 %)** |
| representatives with abundance > 0 | **2** | **0** |
| assembly | 11 contigs, 11,046 bp | 166 contigs, 5,194,724 bp, N50 119,198 |

CoverM was reused verbatim from the cohort run (`193_current_coverm_full_recompute.sbatch`): same
prebuilt `reps.mmi`, same `genome_def.tsv`, `-p minimap2-sr -m relative_abundance -t 16`, coupled
reads, nothing else set. "Assigned" is therefore post-gate: CoverM 0.7.0 defaults to
`--min-covered-fraction 10`, so a representative covered over less than 10 % of its length reports
zero. That default was an open caveat in the exploratory check and is now confirmed from the binary.

### SNTC: the blank contains the E. coli representative

| representative | taxon (R232) | abundance | share of assigned reads | in the 550 soil samples |
|---|---|---|---|---|
| `V42PRr2__single_SemiBin_33` | s__Escherichia coli | **4.296 %** | 50.2 % | 316/550, DNA rho **-0.490** |
| `51Dr3__single_SemiBin_390` | s__Cutibacterium acnes | **4.264 %** | 49.8 % | 290/550, DNA rho **-0.469** |

Those are the only two. C. acnes was never flagged before: it is a human skin commensal, it is not on
the pathogen list `03_census_and_controls.py` reproduces, and its soil behaviour is the E. coli
signature again (52.7 % prevalence, trips 67 / 35 / 47 / 60 / 29 %, anti-correlated with input DNA).
The blank therefore names a second contaminant that the abundance statistics alone did not reach.

With the breadth gate off (step 08) 2,968 representatives hold at least one alignment record, but the
records sum to 130.4 % of the library, so they are alignment records and not reads. They pile onto
`F4Sr2__self_SemiBin_184` (42.9 %) and `24PRr2__self_SemiBin_320` (41.1 %), neither of which reaches
10 % breadth. That is a repeat pileup, and it is the reason the gate exists; it is recorded in
`results/sntc_ungated_counts.tsv` as a diagnostic, not as blank composition.

### SNTC's assembly is laboratory plasmid, not an organism

GTDB-Tk returns `s__Neopusillimonas minor` at ANI 97.77 and AF 0.534, and that call must not be used:
the same ANI screen lists **352** further references at 98.98 to 100 % ANI with AF 0.16 to 0.36. An
11 kb assembly cannot identify an organism. UniRef50 on its 14 predicted proteins can:

| contig | best UniRef50 hit | identity |
|---|---|---|
| SNTC_72 (3,171 bp) | Aminoglycoside 3'-phosphotransferase (`KKA1_ECOLX`) | 1.000 over 271 aa |
| SNTC_72 | Beta-lactamase TEM (`BLAT_ECOLX`) | 0.938 over 103 aa |
| SNTC_72 | Regulatory protein Rop (`ROP_ECOLX`, ColE1 copy control) | 0.944 over 63 aa |
| SNTC_69 | Chloramphenicol acetyltransferase | 1.000 over 46 aa |
| SNTC_32 | Single-stranded DNA-binding protein (`SSB_BPT4`, phage T4) | 1.000 over 176 aa |

A ColE1 backbone carrying kanamycin, ampicillin and chloramphenicol markers, plus phage T4 protein.
Ten representatives carry sequence matching these contigs (`results/control_contig_hits.tsv`); the
strongest is SNTC_72 end to end at **99.94 % over 3,171 bp** against a contig of
`T5PRr1__self_SemiBin_614` (s__Mixta alhagi). Vector sequence is in the catalogue.

### SPositive is Hafnia paralvei and is not in the catalogue

| test | result |
|---|---|
| GTDB-Tk 2.7.2 / R232 ANI screen | `s__Hafnia paralvei` GCF_001655005.1, ANI **98.99**, AF **0.871** |
| next relatives | H. alvei 84.56, H. alvei_B 84.41, H. protea 84.35 |
| skani vs the 5,438, `-s 75` | **no hit** |
| CoverM reads assigned | **0 of 119,683,914** |
| best contig overlap with the catalogue | 69,247 bp of 5,194,724 (1.3 %) at 87.4 % mean identity, to `21PRr3__self_SemiBin_37` s__Enterobacter cloacae |

One organism, cleanly called, absent from the catalogue. Step 09 checks the zero rather than
reporting it: mapped ungated, 41,032 of 1,000,000 SPositive reads (4.1 %) produce a primary
alignment, and the SPositive assembly cut into 98 bp fragments gives the same rate (1,649 of 53,007,
3.1 %). So 4 % of the library does touch the catalogue, spread so thinly across Enterobacteriaceae
conserved genes that no representative reaches 10 % breadth. The zero means no catalogue organism is
present, not that no read aligns. The 47 kb of SPositive contig aligning to the E. coli
representative at 89.6 % identity is family-level conservation, not E. coli.

## Verdicts

**(a) Does SNTC contain the E. coli representative?** Yes, and it is half the blank. The E. coli
representative is one of only two catalogue genomes detectable in a no-template control, at 4.30 %
relative abundance and 50.2 % of assigned reads. This closes the exploratory check's open question in
the direction it feared: the blank confirms contamination rather than certifying cleanliness.

**(b) Is SPositive a single known organism, and is it among the 5,438?** Yes to the first, no to the
second. Hafnia paralvei at 98.99 % ANI over 87 % of the genome, next relative 84.6 %, and absent from
the catalogue by skani, by read mapping and by contig overlap.

**(c) What the blanks license the paper to say.** They license a stated, quantified contamination
control in place of an absence of one, and they change three things.
1. The E. coli representative is a reagent contaminant, now shown directly and not only inferred from
   trip skew and biomass anti-correlation. It should carry that annotation wherever the census names
   it. Whether it is removed from the 5,438 remains the PI's call, as before.
2. C. acnes (`51Dr3__single_SemiBin_390`) is a second contaminant on the same evidence and is not on
   the pathogen list, so the biosafety paragraph cannot be rewritten around E. coli alone.
3. The blank carries a ColE1 vector with three antibiotic-resistance markers, and ten
   representatives carry matching sequence, one of them over a full 3.2 kb at 99.94 %. Any resistome
   count that includes those regions is counting laboratory plasmid. This needs checking against the
   AMR census before that section is defended.
The positive control is the clean half of the result: a 120 M read library of an organism absent from
the catalogue puts zero reads on any representative, so there is no index hopping between the control
lane and the soil libraries at a level that would matter, in either direction.

## Files

`scripts/01_contamination_tests.py`, `02_biomass_within_trip.py`, `03_census_and_controls.py`
(promoted, paths repointed only), `04_coverm_controls.sbatch`, `05_identify_controls.sbatch`,
`06_sntc_uniref.sbatch`, `07_summarise_controls.py`, `08_coverm_sntc_counts.sbatch`,
`09_sanity_spositive_zero.sbatch`.

`results/` per_genome_stats.tsv, target_vs_null.tsv, by_trip.tsv, biomass_within_trip.tsv,
pathogen_census_5438.tsv, noncanonical_sequenced_samples.tsv (promoted, unchanged);
control_read_stats.tsv, control_profiles.tsv, control_contig_hits.tsv, sntc_ungated_counts.tsv;
`coverm/{SNTC,SPositive,SNTC_counts}/`; `identify/` gtdbtk/, skani_controls_vs_reps.tsv,
{SNTC,SPositive}_contigs_vs_reps.paf, sntc_uniref50_hits.tsv. `figures/abundance_by_trip.png`.

## Caveats

The 10 % covered-fraction gate is CoverM's default and applies to the cohort profiles as well as to
these controls, so "absent" throughout means "below 10 % breadth", not "no read aligned". Only one
no-template control and one positive control exist for the whole shotgun cohort, both filed under
trip 4, so the blanks cannot be resolved by batch and cannot bound contamination per trip. The four
amplicon controls under trip1 and trip2 are 16S libraries and say nothing about these MAGs. SNTC's
335,190 reads are not a deep blank: absence of a third contaminant in it is weak evidence. The
resistome implication in verdict (c) is a flag raised by contig overlap, not a recount of the AMR
census, which has not been rerun here.
