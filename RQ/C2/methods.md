# RQ-C2 — Historical exploratory pangenome workflow

**Status.** This is the exact archived workflow behind the exploratory RQ-C2 output. Its raw-union dRep inputs
and intermediate family membership have not been promoted into the active checksum-closed chain, so its results
are not used as a mechanism for RQ-C1. The tested hypothesis was whether a historical core/accessory partition
was associated with UniRef50 no-hit status; association would not establish why taxonomy-based prediction errs.

Framework: Tettelin et al. 2005 (PNAS; "pan-genome", open vs closed) [`tettelin2005pangenome`]
and McInerney, McNally & O'Connell 2017 (Nat Microbiol; open-pangenome interpretation)
[`mcinerney2017openpangenome`].

## Historical inputs on IBEX

| Input | Path |
|---|---|
| Species clusters (dRep secondary clusters ≈ 95% ANI) | `$M/04_drep_self_vs_single/out/data_tables/Cdb.csv` (col `secondary_cluster`) |
| Member MAG bins (pre-dRep, 37,085) | `$M/04_drep_self_vs_single/inputs/*.fa.gz` (contig FASTA per bin) |
| Per-assembly proteins (Prodigal) | `$M/04_protein_catalog/assembly_proteins/<SAMPLE>/<SAMPLE>.faa.gz` |
| Dark proteome / UniRef50 status | `$M/06_darkmatter/isdark/is_dark.tsv` (`candidate_protein`,`best_uniref_pident`,`is_dark`) |

`$M = /ibex/scratch/projects/c2014/rob/metagenomics`

**ID join keys (verified):**
- Bin contig header `>F6Dr2_2241028` ↔ protein header `>F6Dr2::F6Dr2_2241028_1`.
  Protein→contig: take the string after `::`, drop the trailing `_<geneidx>`.
- Dark table id `F6Dr2__single_SemiBin_54::F6Dr2_1498071_39`; join on the string after `::`
  (`contig_gene`), which is globally unique because contig names are sample-prefixed.

## Pipeline (`140_pangenome_openness.py`, driven by `140_pangenome_openness.sbatch`)

Real compute runs as a **batch** sbatch job (`--account=pi-hohndor --partition=batch`,
32 cpus, 256G, 24h). MMseqs2 16.747c6 from the `eggnog` conda env. **No clustering on login nodes.**

1. **Species & non-redundant members.** Parse `Cdb.csv` → species = `secondary_cluster`.
   The catalogue contains *two* binning passes of each assembly (`ss__` single-sample,
   `selfsup__` self-supervised), which can place a near-duplicate of the same genome into one
   species and inflate the core. To prevent this we keep **one genome per (species, sample)**
   (prefer `ss__`; lexicographic tie-break) → distinct field samples = true pangenome members.
   Keep species with **≥ 10** members.
   - Raw: **600** species ≥10 / **17,215** member bins.
   - Non-redundant (one per sample): **345** species / **8,635** genomes.
2. **Per-genome proteomes.** Build a contig→genome map from member-bin headers; stream each
   needed sample `.faa.gz` once (parallel) and route every protein to its genome's FASTA,
   re-headed `>{genome}|{origid}`.
3. **Gene families.** Per species, pool member proteins and cluster with
   `mmseqs easy-cluster --min-seq-id 0.5 -c 0.8 --cov-mode 1 -e 1e-3`. Each cluster = a gene
   family; presence/absence = which genomes contribute ≥1 member.
   Classify: **core** ≥95% of genomes, **soft-core** ≥90%, **shell** 15–90%, **cloud/accessory** <15%.
4. **Openness.** Over 50 random genome orderings, accumulate pan-genome size and new genes.
   Fit Heaps' law new-genes ∝ n^(−α) (Tettelin): **α < 1 ⇒ OPEN**. Also report the pan-size
   exponent γ and the core/pan ratio.
5. **Functional split / dark overlap.** Tag each family against `is_dark.tsv`
   (dark = no UniRef50 hit) and by mean UniRef50 best-pident ("knownness"; high pident ≈
   conserved housekeeping). Pool a **(core/soft-core vs shell/cloud) × (dark vs known)** 2×2
   across species → odds ratio (+ Fisher p if scipy present), and test whether accessory
   families have lower UniRef pident than core (Mann–Whitney). Report the dark overlap of the
   open accessory fraction vs the core.
   *Note:* `is_dark.tsv` covers the dRep representative proteomes (one per species), so the
   dark/known enrichment is evaluated on annotated family members; the Heaps/openness result
   uses the full non-redundant member set.

## Historical outputs (land in `$M/06_pangenome/out/`)

- `species_openness.tsv` — per species: n_genomes, n_families, n_core, core_pan_ratio,
  heaps_alpha, pan_exponent_gamma, heaps_R2, open/closed.
- `family_table.tsv` — per family: class, dark_frac, mean_pident.
- `species_members.tsv` — the non-redundant member set per species.
- `enrichment_summary.json` — openness counts, dark odds ratio (accessory vs core),
  accessory/core dark overlap, mean pident core vs accessory, Fisher/Mann–Whitney p.

## Parameters (exact)

`MIN_GENOMES=10  CORE_FRAC=0.95  SOFT_FRAC=0.90  SHELL_FRAC=0.15  PERMS=50`
`mmseqs easy-cluster --min-seq-id 0.5 -c 0.8 --cov-mode 1 -e 1e-3 --threads 32`.

Script (archived IBEX path): `$M/scripts/ibex/140_pangenome_openness.{py,sbatch}` (copies in this folder).
Re-running it reproduces that historical universe; it does not update the analysis to the current catalogue.
