# Dossier 02 — the dark definition and the dark protein families (§1.2, §1.4; M.4, M.7)

Scope: how "dark" is called; family construction, admission, control, attribution. PERMANOVA / gene–environment (last paragraph of M.7, `darkfam_ecology_clean/`) → **separate dossier**. `M` = `/ibex/scratch/projects/c2014/rob/metagenomics`. ✅ = recomputed on disk 2026-09-17.

## Question

Is the sequence-dark proteome real coding biology or spurious ORFs from a fragmented assembly? §1.2 sizes it; §1.4 answers with four legs (recurrence · selection · abundance · persistence) + phylum attribution, all on **one** object — §1.2's own dark call carried into §1.4.

## Basis — every denominator, and the file that produces it

| Object | Count | Rule applies to | Produced by |
|---|---|---|---|
| Predicted proteins, 552 assemblies | 753,305,533 | protein (pyrodigal `-p meta --no-stop-codon`) | `M/04_protein_catalog/assembly_proteins/` |
| Cluster-90 representatives | 449,265,306 | cluster (90 % id / 80 % cov) | `M/scripts/05_pool_and_cluster.sh` → `all_linclust90/` |
| UniRef50-dark reps | 69,938,628 (15.6 %) | **cluster**, called from its *members* | `darkclean_4_combine.sh` → `04_catalog_dark/dark_rep_ids.txt.gz` |
| Genome-resolved reps | 43,423,322 | cluster, **any member** on an HQ+MQ contig | `darkclean_5_mag.sh` → `05_mag_dark/mag_dark_summary.tsv` |
| Genome-resolved dark reps | 4,292,446 (9.89 %) | cluster | `darkclean_5_mag.sh` → `mag_dark_rep_ids.txt.gz` |
| Dark families admitted | 106,042 | family (cluster of dark reps) | `fam_2_recurrence.sh` → `FAMILIES_FINAL.txt` ✅ |
| Member CDS in them | 536,437 | member gene, nucleotide | `fam_3_cds.sh` → `darkfam_members.fna`, `member_fam_final.tsv` |

§1.4 workbench: `M/06_darkmatter/darkfam_basis_clean/` (owner velazqal; README = audit + job log). Not in the paper: 3,658,797 total 50/80 clusters; 39,130,876 non-dark MAG-mapped reps (= the control pool).

## How it ran

| Step | Parameters (read from the script, not from Methods) |
|---|---|
| Catalogue clustering | `mmseqs linclust --min-seq-id 0.9 -c 0.8 --cov-mode 0`, v18.8cc5c (PI's `05_pool_and_cluster.sh`) |
| UniRef50 search | `mmseqs search -s 7.5 -e 1e-3 --max-seqs 10`, 120-task array (`darkclean_3b_uniref.sh`). Pure E-value — **no identity, no coverage filter** |
| Darkness level | **Member-level call, cluster-level verdict**: dark iff the cluster has ≥1 *searched* member and **no member returned a hit** (`darkclean_4_combine.sh`). Not a rep-only call |
| Genome-resolved gate | **Per-member**: ANY member's contig in one of the **37,949 HQ+MQ bins of the two UNREFINED binnings** (`results/catalogue_provenance.tsv` + `contig_to_bin.tsv.gz`). Membership, not read mapping |
| Family clustering | `mmseqs easy-cluster --min-seq-id 0.5 -c 0.8 --cov-mode 0` on the 4,292,446 dark reps (`fam_1_cluster.sh`); sensitivity left at the easy-cluster default **-s 4**, cluster-mode 0, 3 cascaded steps (log) |
| Admission | ≥3 members **AND** ≥3 samples **AND** ≥2 sites (`fam_2_recurrence.sh`), from a streaming pass over all 753,305,533 rows of `all_cluster90.tsv`, restricted to the **551 canonical samples**; 548 contribute instances |
| Member CDS | prodigal start/end/strand from the `.faa.gz` headers, sliced from `<sample>_contigs.fa` (`fam_3_cds.sh`) |
| Presence / abundance | CoverM `contig -p minimap2-sr -m tpm covered_fraction`, 551-task array vs **all 536,437 member CDS** (`fam_4_coverm_array.sh`); present ⇔ **any member covered_fraction ≥ 0.5**; family TPM = **sum of member TPM** (`fam_8_aggregate.sh`) |
| Phylum | member → contig → **refined** bin → dRep Cdb/Wdb species rep → GTDB-Tk **R232** (`fam_9_phylum.sh`) |
| Robustness | same reps reclustered at `--min-seq-id 0.3 -c 0.8` (`fam_7_robustness.sh`) |

Bin-set switch: admission and the control use the **37,949 unrefined** HQ+MQ bins; phylum attribution uses the **refined** catalogue (12,433 bins → 5,438 reps). Hence 74 % attribution coverage. Deliberate, undeclared in M.7.

## Controls — the length- and size-matched annotated control

Built by **Alejandra** (`fam_5_control.sh`, owner velazqal, 2026-08-13; a *rewrite* after a multi-agent review — its header names fixes F1–F5). The dN/dS **worker** is the PI's June algorithm (`94_darkfam_dnds.sh`) reused verbatim and repointed; both arms run it via `ARM=dark|ctrl`.

| Question | Answer (from the script) |
|---|---|
| Pool | All **non-dark** cluster-90 reps passing the *identical* HQ+MQ gate (log: 15,366,244 binned contigs, same as darkclean_5); darkness excluded against the **full 69.9 M** dark set, not the 4.29 M subset → 39,130,876 reps |
| Matched on (1) | **Representative length, 10-aa bins** (437 bins), sampling the pool to reproduce the dark length histogram — **distributional, not paired** |
| Matched on (2) | **Family size**, bins 3–9 exact / 10–19 / 20+; the same *number* of control families drawn per bin — **count-matched per bin, not 1:1 paired families** |
| Tolerance | 10 aa; exact for sizes 3–9, then two coarse classes |
| Clustering | identical: `easy-cluster --min-seq-id 0.5 -c 0.8 --cov-mode 0` |
| Control families | **106,042** (`CONTROL_FAMILIES.txt` ✅); size-matching shortfall **0** |
| Same admission rule? | **NO.** Control admission is `m >= 3` members only — **no ≥3-samples / ≥2-sites gate, no paralog-only exclusion** (the 2,784 drops hit the dark arm alone) |
| "Annotated"? | Means **has a UniRef50 hit**. No eggNOG/KEGG/Bakta requirement anywhere in the script |
| Shortfall | Pool ran **12.47 %** short: 3,757,129 of 4,292,446 (`ctrl_sample_n.txt`), in 20–29 aa (409,948) and 30–39 aa (125,366) |
| Seed | `random.seed(20260813)` in both sampling stages — deterministic |

dN/dS (`fam_6c_dnds.sh`, identical both arms): entry needs ≥3 usable CDS, each ≥90 nt, codon-trimmed, translatable under **table 11** with no internal stop, protein ≥25 aa. Aligner **MAFFT v7.526 `--auto`** on proteins, **back-translated** to a codon alignment. Model **NG86** (Biopython `cal_dn_ds`). Pairwise gap/N codon stripping; ≥30 clean codons; identical pairs dropped; `0 < dS ≤ 2`; ≤80 **deterministically thinned** pairs per family. Family statistic = **median** pairwise dN/dS. n = 97,200 dark · 95,647 control (91.7 % / 90.2 % of 106,042 — near-symmetric attrition).

**JUDGMENT (mine).** Fit for its stated purpose and only that. It is a **scale anchor, not a null**: it answers "what does 50 %-identity clustering plus this pipeline return on sequences we know are genes?", and it removes the length confound that inflated the June gap (0.391 vs 0.286). It does not answer "what would this pipeline return on non-coding sequence" — no shuffled, frameshifted or intergenic arm exists. Three real weaknesses: (a) the control skips the recurrence gate, so the arms are unmatched on ecological spread and the control keeps within-genome paralog sets the dark arm drops; (b) matching is on *representative* length at input, not the aligned CDS length that enters dN/dS — mitigating: the largest shortfall bin (20–29 aa ≈ 60–89 nt) lies entirely below the 90-nt entry gate, but the 30–39 aa shortfall (125,366) is inside it; (c) no taxonomic, GC or codon-usage matching. None plausibly reverses the observed direction; all are rebuttal material, not paper material.

## Numbers (exact, with source; REGISTERED = in `CATALOGUE_NUMBERS.md §1.4`)

| Quantity | Value | Source |
|---|---|---|
| Families / member CDS | **106,042 / 536,437** | `FAMILIES_FINAL.txt`, `member_fam_final.tsv` ✅ — REGISTERED |
| Excluded classes | 2,784 paralog-only · 782,616 young (1–2 variants) | `logs/fam2_50439896.out` — REGISTERED |
| Recurrence | median **12 samples / 9 sites** | admitted rows of `fam_recurrence.tsv` ✅ — REGISTERED |
| dN/dS dark | median **0.4253**, 96.69 % < 1, n = **97,200** | `results/selection/dnds.tsv` ✅ — REGISTERED as 0.425 |
| dN/dS control | median **0.3884**, 99.61 % < 1, n = **95,647** | same file ✅ — REGISTERED as 0.388; MWU p < 1e-300 (Supp. S5 only) |
| Fast tail > 1 | 3.3 % dark vs 0.4 % control | `dnds.tsv` — REGISTERED, hypothesis-generating only |
| Prevalence / abundance | median **98 / 551**, **4** single-sample, median peak **99.4 TPM** | `logs/fam8_*.out` (job 50531118) — REGISTERED |
| Persistence, recruitment | **79.4 %** all five (84,174) · **94.8 %** ≥4 (100,489) | `fam_presence.tsv` × `sample_site_trip.tsv` ✅ — REGISTERED |
| Persistence, re-assembly | **28.9 %** all five (30,609 / 106,042) | `fam_recurrence.tsv` col `n_trips` ✅ — REGISTERED |
| Length-stratified check | ≥250 aa: median 55/551, 67.1 % all five (<60 aa: 201, 89.7 %) | `CORRECTIONS.md` 2026-08-17 |
| Phylum | 78,727 (74 %) attributable, **97.9 %** single-phylum; Actinomycetota 32 %, Thermoplasmatota 19 % | `fam_phylum.tsv` — REGISTERED |
| Threshold robustness | ≥3-member clusters 108,826 (50 %) → **130,843** (30 %) | `logs/fam7_50480205.out` ✅ — REGISTERED |

**How presence per sample was decided.** Not by assembly. Each of the 551 canonical samples' reads were mapped (CoverM/minimap2-sr) against all 536,437 member CDS; a family is **present in a sample iff at least one member gene reaches covered_fraction ≥ 0.5** there (`fam_8_aggregate.sh`). Expedition persistence = number of the five trips with ≥1 present sample. The 28.9 % is a **different rule**: trips in which an instance was independently **assembled** (`n_trips`, instance pass) — it never touches read mapping.

## Caveats (quoted)

- M.7: "Clustering at 50\% identity makes any set of families look conserved, so the selection test below and the catalogue search of M.12 both run against a matched annotated control that inherits the same artefact." · "The sampled pool ran 12.5\% short of the target, concentrated in the shortest bins." · "The estimate is family-scoped and supports no catalogue-wide statement about selection regimes."
- M.4: "Genome-resolved here means protein-to-contig membership, not read mapping."
- Abstract: "…is the purifying signature of coding genes measured with one instrument, **not stronger constraint**."
- §1.4: "A tail of 3.3\% of dark families (0.4\% of the control) lies above one; we treat these as hypothesis-generating candidates for diversifying selection and draw no conclusion from them here."
- `CORRECTIONS.md` 2026-08-17: "**Any distribution claim from this matrix must be length-stratified.**" · "§1.4 is the only affected claim and it survives."
- Fig. caption (main.tex:1288): "The line at dN/dS${}=1$ is a reference value, not a pipeline-calibrated null."

## Status — current vs dead bases

| | Current | Dead |
|---|---|---|
| Dark call | darkclean: `-s 7.5 -e 1e-3 --max-seqs 10`, dark = **no hit** | June `76_candidate_isdark.sh`: `-s 5.7 --max-seqs 5`, dark = no hit **OR** pident < 30 |
| Family object | **106,042** / 536,437 CDS | **23,382** June orthogroups — 20.1 % dark only via the `<30 %` clause; 67.7 % of eligible clusters never screened |
| dN/dS | 0.425 vs matched control 0.388 | 0.391 vs 0.286 (unmatched on length) |
| Genome basis | 5,438 dRep species reps (R232); 37,949 unrefined HQ+MQ bins for the protein gate | **15,469** single-sample MAG snapshot / 51.9 M Bakta proteins / 2.84 M dark |
| Manuscript | `tab:ladder` at main.tex:1033; Supp. Table S5 `tab:darkfam-evidence` | R226 `tab:ladder` (main.tex:1065) and the Part B scaffold (main.tex:398–417, 23,382 / 0.391) sit inside `\iffalse` — fenced, **not printed**, must not be deleted |

## What can and cannot be claimed

Dark 0.425 vs control 0.388 = dark families are measurably **less** constrained than the matched annotated control.

- **CAN**: the dark families carry a purifying-selection signature (median far below 1; 96.7 % < 1) of the same kind and order as known genes measured with the identical instrument — i.e. translated coding sequence, not non-coding ORFs. The control fixes the scale on which 0.425 is read.
- **CAN**: the earlier dark-vs-control gap was substantially a length artefact; matching removes most of it and the paper reports the residual honestly.
- **CANNOT**: "as constrained as annotated genes" (they are not; the gap is highly significant), and emphatically not "more constrained" or "under stronger selection".
- **CANNOT**: "the below-1 value is not an artefact of 50 %-identity clustering" — the control *shares* that artefact by design; it bounds it at 0.388, it does not remove it. With **no ~1.0 null arm**, the absolute distance from 1 is uncalibrated; only the distance from the annotated arm is measured.
- **CANNOT**: anything catalogue-wide about selection regimes, or anything from the 3.3 % fast tail.
- The abstract already says "not stronger constraint"; **§1.4 Results never states the direction** — it prints both medians side by side. INFERENCE (mine): that is the most exposed sentence in §1.4.

## Reviewer attack surface

1. **No non-coding null** — "your control is another set of genes; where is the arm expected at ≈1?" Strongest attack; answerable today only by the pre-argued M.7 sentence.
2. **Direction** — "dark is *less* constrained than your own control, consistent with real genes diluted by spurious ORFs." No quantitative decomposition exists.
3. **Control admission asymmetry** — no ≥3-samples/≥2-sites gate and no paralog-only exclusion on the control; M.7 is silent and `fam_5_control.sh` shows it in four lines.
4. **"Annotated" = UniRef50 hit**, not functionally annotated — the word overstates.
5. **Presence rule is length-biased** (internal, 2026-08-17); §1.4 pre-empts with the ≥250 aa stratum, but 98/551 is composition-weighted and spans 201→55.
6. **Bin-set switch** (37,949 unrefined for admission vs refined for phylum) and the 74 % coverage are unexplained in M.7.
7. **Two persistence rules in one sentence** (79.4 % recruitment vs 28.9 % re-assembly); a fast reader takes 28.9 % as a failure of 79.4 %.
8. **Site count** — the gate's labels give **61 distinct sites** across the 551 canonical samples (site52 has one sample); §1.1 and M.1 report **60**.
9. **easy-cluster `-s 4` default** unstated; the 30 % recluster tests threshold stability, not sensitivity stability.

## Open questions for Alejandra

1. Add one explicit direction sentence to §1.4 ("slightly lower than, not higher than, the control"), matching the abstract? Cheap; closes attack 2.
2. Disclose in M.7 that the control was admitted on **member count only** (and why that is conservative), or rebuild it under the full admission rule (re-run fam_5 + 6a/6b/6c)?
3. Is a shuffled/frameshifted null arm worth one job, or is the pre-argued M.7 sentence the final answer?
4. Replace "annotated families" with "UniRef50-matched families" in §1.4 and M.7?
5. 61 vs 60 sites — which is canonical for the ≥2-sites gate, and does M.7 need a footnote?
6. Should M.7 state the unrefined-vs-refined bin-set switch and the `-s 4` easy-cluster default?
